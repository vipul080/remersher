#include "quad_cleanup.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <unordered_map>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

Vec3 polyNormal(const std::vector<Vec3>& P, const std::vector<int>& f) {
    Vec3 n{0, 0, 0};
    for (size_t i = 0; i < f.size(); ++i) n = add(n, cross(P[f[i]], P[f[(i + 1) % f.size()]]));
    return n;
}

// A quad is acceptable if it is non-degenerate and every corner turns the same way as the face.
bool goodQuad(const std::vector<Vec3>& P, const std::vector<int>& q, const Vec3& refNormal) {
    Vec3 n = polyNormal(P, q);
    double ln = norm(n);
    if (ln <= 0 || dot(n, refNormal) <= 0.3 * ln * norm(refNormal)) return false;
    for (int i = 0; i < 4; ++i) {
        Vec3 c = cross(sub(P[q[(i + 1) % 4]], P[q[i]]), sub(P[q[(i + 3) % 4]], P[q[i]]));
        if (dot(c, n) <= 0) return false;
    }
    return true;
}

// Mean |corner angle - 90| (degrees) of a quad, or 0 for other polygons.
double quadAngleError(const std::vector<Vec3>& P, const std::vector<int>& q) {
    if (q.size() != 4) return 0;
    double sum = 0;
    for (int i = 0; i < 4; ++i) {
        Vec3 u = sub(P[q[(i + 3) % 4]], P[q[i]]), w = sub(P[q[(i + 1) % 4]], P[q[i]]);
        double d = norm(u) * norm(w);
        double c = d > 0 ? std::clamp(dot(u, w) / d, -1.0, 1.0) : 1.0;
        sum += std::abs(std::acos(c) * 180.0 / 3.14159265358979323846 - 90.0);
    }
    return sum / 4;
}

struct Op {
    int gain;
    int kind;  // 0 = diagonal collapse, 1 = edge rotation
    int f1, f2;
    int a, c;                      // collapse: merged corners (c into a)
    std::array<int, 4> q1, q2;     // rotation: replacement quads
};

}  // namespace

PoleCleanupStats cancelPoles(Mesh& mesh, const std::vector<char>& locked, int rounds,
                             std::vector<char>* changed) {
    PoleCleanupStats st;
    if (changed) changed->assign(mesh.vertices.size(), 0);
    auto& P = mesh.vertices;
    auto& F = mesh.faces;
    std::vector<char> alive(F.size(), 1);

    for (int round = 0; round < rounds; ++round) {
        const size_t nv = P.size();
        std::vector<std::vector<int>> vfaces(nv), nbrs(nv);
        std::unordered_map<uint64_t, std::vector<int>> edgeFaces;
        for (int fi = 0; fi < (int)F.size(); ++fi) {
            if (!alive[fi]) continue;
            const auto& f = F[fi];
            for (size_t i = 0; i < f.size(); ++i) {
                int a = f[i], b = f[(i + 1) % f.size()];
                vfaces[a].push_back(fi);
                nbrs[a].push_back(b), nbrs[b].push_back(a);
                edgeFaces[edgeKey(a, b)].push_back(fi);
            }
        }
        std::vector<char> border(nv, 0);
        for (const auto& [k, fs] : edgeFaces)
            if (fs.size() != 2) border[k >> 32] = border[k & 0xffffffff] = 1;
        std::vector<int> val(nv, 0);
        for (size_t v = 0; v < nv; ++v) {
            auto& n = nbrs[v];
            std::sort(n.begin(), n.end());
            n.erase(std::unique(n.begin(), n.end()), n.end());
            val[v] = (int)n.size();
        }
        auto adjacent = [&](int a, int b) { return std::binary_search(nbrs[a].begin(), nbrs[a].end(), b); };
        auto free = [&](int v) { return !border[v] && !(v < (int)locked.size() && locked[v]); };
        // Cost of a valence: mostly whether the vertex is a pole at all, then how irregular it is,
        // so an operation never trades one pole for two milder ones.
        auto dev = [](int valence) { return valence == 4 ? 0 : 10 + std::abs(valence - 4); };

        std::vector<Op> ops;
        // Diagonal collapses.
        for (int fi = 0; fi < (int)F.size(); ++fi) {
            const auto& f = F[fi];
            if (!alive[fi] || f.size() != 4) continue;
            for (int i = 0; i < 2; ++i) {
                const int a = f[i], b = f[i + 1], c = f[i + 2], d = f[(i + 3) % 4];
                if (!free(a) || !free(c) || border[b] || border[d] || adjacent(a, c)) continue;
                const int before = dev(val[a]) + dev(val[c]) + dev(val[b]) + dev(val[d]);
                const int after = dev(val[a] + val[c] - 2) + dev(val[b] - 1) + dev(val[d] - 1);
                if (after >= before || val[b] <= 3 || val[d] <= 3) continue;
                ops.push_back({before - after, 0, fi, -1, a, c, {}, {}});
            }
        }
        // Edge rotations.
        for (const auto& [k, fs] : edgeFaces) {
            if (fs.size() != 2 || F[fs[0]].size() != 4 || F[fs[1]].size() != 4) continue;
            int u = (int)(k >> 32), v = (int)(k & 0xffffffff);
            if (!free(u) || !free(v)) continue;
            // Name the quads [u, v, x1, x2] and [v, u, y1, y2] in winding order.
            auto rotateTo = [](std::vector<int> q, int first, int second, std::array<int, 4>& out) {
                for (int r = 0; r < 4; ++r) {
                    if (q[0] == first && q[1] == second) {
                        std::copy(q.begin(), q.end(), out.begin());
                        return true;
                    }
                    std::rotate(q.begin(), q.begin() + 1, q.end());
                }
                return false;
            };
            std::array<int, 4> A, B;
            if (!rotateTo(F[fs[0]], u, v, A)) std::swap(u, v);
            if (!rotateTo(F[fs[0]], u, v, A) || !rotateTo(F[fs[1]], v, u, B)) continue;
            const int x1 = A[2], x2 = A[3], y1 = B[2], y2 = B[3];
            if (border[x1] || border[x2] || border[y1] || border[y2]) continue;
            if (val[u] <= 3 || val[v] <= 3) continue;
            struct Option { int p, q; std::array<int, 4> q1, q2; };
            const Option options[2] = {{x1, y1, {x1, x2, u, y1}, {y1, y2, v, x1}},
                                       {x2, y2, {x2, u, y1, y2}, {y2, v, x1, x2}}};
            for (const auto& o : options) {
                if (o.p == o.q || adjacent(o.p, o.q)) continue;
                const int before = dev(val[u]) + dev(val[v]) + dev(val[o.p]) + dev(val[o.q]);
                const int after = dev(val[u] - 1) + dev(val[v] - 1) + dev(val[o.p] + 1) + dev(val[o.q] + 1);
                if (after >= before) continue;
                ops.push_back({before - after, 1, fs[0], fs[1], -1, -1, o.q1, o.q2});
            }
        }
        if (ops.empty()) break;
        std::stable_sort(ops.begin(), ops.end(), [](const Op& x, const Op& y) { return x.gain > y.gain; });

        std::vector<char> touched(nv, 0);
        int applied = 0;
        for (const Op& op : ops) {
            if (op.kind == 0) {
                const int a = op.a, c = op.c;
                if (!alive[op.f1]) continue;
                bool busy = false;
                for (int fi : vfaces[a]) for (int w : F[fi]) busy |= touched[w] != 0;
                for (int fi : vfaces[c]) for (int w : F[fi]) busy |= touched[w] != 0;
                if (busy) continue;
                // Link condition: a and c may only share the other two corners of the quad.
                std::vector<int> common;
                std::set_intersection(nbrs[a].begin(), nbrs[a].end(), nbrs[c].begin(), nbrs[c].end(),
                                      std::back_inserter(common));
                if (common.size() != 2) continue;
                const Vec3 merged = {0.5 * (P[a][0] + P[c][0]), 0.5 * (P[a][1] + P[c][1]), 0.5 * (P[a][2] + P[c][2])};
                std::vector<Vec3> Q = P;
                Q[a] = Q[c] = merged;
                bool ok = true;
                double errBefore = quadAngleError(P, F[op.f1]), errAfter = 0;
                int counted = 0;
                for (int v : {a, c})
                    for (int fi : vfaces[v]) {
                        if (fi == op.f1) continue;
                        Vec3 before = polyNormal(P, F[fi]), after = polyNormal(Q, F[fi]);
                        double lb = norm(before), la = norm(after);
                        if (la < 1e-3 * lb || dot(after, before) <= 0.3 * la * lb) ok = false;
                        errBefore += quadAngleError(P, F[fi]);
                        errAfter += quadAngleError(Q, F[fi]);
                        ++counted;
                    }
                // Removing two poles is not worth visibly distorting the quads around them.
                if (!ok || errAfter / std::max(1, counted) > errBefore / (counted + 1) + 3.0) continue;
                for (int v : {a, c})
                    for (int fi : vfaces[v])
                        for (int w : F[fi]) {
                            touched[w] = 1;
                            if (changed) (*changed)[w] = 1;
                        }
                P[a] = merged;
                alive[op.f1] = 0;
                for (int fi : vfaces[c])
                    if (fi != op.f1)
                        for (int& w : F[fi])
                            if (w == c) w = a;
                ++st.diagonalCollapses;
                ++applied;
            } else {
                if (!alive[op.f1] || !alive[op.f2]) continue;
                bool busy = false;
                for (int fi : {op.f1, op.f2}) for (int w : F[fi]) busy |= touched[w] != 0;
                if (busy) continue;
                const Vec3 ref = add(polyNormal(P, F[op.f1]), polyNormal(P, F[op.f2]));
                std::vector<int> q1(op.q1.begin(), op.q1.end()), q2(op.q2.begin(), op.q2.end());
                if (!goodQuad(P, q1, ref) || !goodQuad(P, q2, ref)) continue;
                for (int fi : {op.f1, op.f2})
                    for (int w : F[fi]) {
                        touched[w] = 1;
                        if (changed) (*changed)[w] = 1;
                    }
                F[op.f1] = q1;
                F[op.f2] = q2;
                ++st.edgeRotations;
                ++applied;
            }
        }
        if (!applied) break;
    }

    std::vector<std::vector<int>> kept;
    kept.reserve(F.size());
    for (size_t fi = 0; fi < F.size(); ++fi)
        if (alive[fi]) kept.push_back(std::move(F[fi]));
    F = std::move(kept);
    return st;
}

}  // namespace remersher
