#include "crease.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdint>
#include <unordered_map>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

}  // namespace

std::vector<std::array<Vec3, 2>> detectCreases(const TriangleMesh& mesh, double angleDegrees,
                                                double spacing) {
    std::vector<std::array<Vec3, 2>> out;
    if (angleDegrees <= 0 || angleDegrees >= 180) return out;
    const auto& V = mesh.vertices;
    const auto& T = mesh.triangles;

    std::vector<Vec3> normal(T.size());
    for (size_t f = 0; f < T.size(); ++f) {
        const Vec3 &a = V[T[f][0]], &b = V[T[f][1]], &c = V[T[f][2]];
        const Vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        Vec3 n{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
        const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        normal[f] = l > 0 ? Vec3{n[0] / l, n[1] / l, n[2] / l} : Vec3{0, 0, 0};
    }
    std::unordered_map<uint64_t, std::vector<int>> edgeFaces;
    for (int f = 0; f < (int)T.size(); ++f)
        for (int i = 0; i < 3; ++i) edgeFaces[edgeKey(T[f][i], T[f][(i + 1) % 3])].push_back(f);

    // Dihedral angle (degrees) of every manifold edge; borders and non-manifold edges get none.
    std::unordered_map<uint64_t, double> dihedral;
    dihedral.reserve(edgeFaces.size());
    for (const auto& [k, fs] : edgeFaces) {
        if (fs.size() != 2) continue;
        const Vec3 &n0 = normal[fs[0]], &n1 = normal[fs[1]];
        const double c = std::max(-1.0, std::min(1.0, n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]));
        dihedral[k] = std::acos(c) * 180.0 / 3.14159265358979323846;
    }

    const bool bevels = spacing > 0;
    const double lowAngle = bevels ? std::min(angleDegrees, 20.0) : angleDegrees;
    std::vector<uint64_t> candidates;
    for (const auto& [k, theta] : dihedral)
        if (theta > lowAngle) candidates.push_back(k);
    std::sort(candidates.begin(), candidates.end());  // deterministic grouping
    std::unordered_map<int, std::vector<int>> incident;  // vertex -> candidate indices
    for (int i = 0; i < (int)candidates.size(); ++i) {
        incident[(int)(candidates[i] >> 32)].push_back(i);
        incident[(int)(candidates[i] & 0xffffffff)].push_back(i);
    }

    // Line-like connected pieces.
    struct Piece {
        std::vector<int> edges;
        double meanAngle = 0;
    };
    std::vector<Piece> pieces;
    std::vector<int> visited(candidates.size(), 0);
    for (int seed = 0; seed < (int)candidates.size(); ++seed) {
        if (visited[seed]) continue;
        Piece p;
        std::vector<int> stack{seed};
        visited[seed] = 1;
        while (!stack.empty()) {
            const int e = stack.back();
            stack.pop_back();
            p.edges.push_back(e);
            for (int v : {(int)(candidates[e] >> 32), (int)(candidates[e] & 0xffffffff)})
                for (int other : incident[v])
                    if (!visited[other]) visited[other] = 1, stack.push_back(other);
        }
        std::unordered_map<int, int> degree;
        for (int e : p.edges) {
            ++degree[(int)(candidates[e] >> 32)], ++degree[(int)(candidates[e] & 0xffffffff)];
            p.meanAngle += dihedral[candidates[e]];
        }
        p.meanAngle /= p.edges.size();
        int loose = 0;
        for (const auto& [v, d] : degree) loose += d == 1;
        const double meanDegree = 2.0 * p.edges.size() / degree.size();
        if (p.edges.size() < 4 || loose > 0.25 * degree.size() || meanDegree > 3.2) continue;
        pieces.push_back(std::move(p));
    }

    auto emit = [&](const Piece& p) {
        for (int e : p.edges) out.push_back({V[candidates[e] >> 32], V[candidates[e] & 0xffffffff]});
    };
    std::vector<int> soft;  // pieces below the angle: possible bevel segments
    for (int i = 0; i < (int)pieces.size(); ++i) {
        if (pieces[i].meanAngle > angleDegrees) emit(pieces[i]);
        else if (bevels) soft.push_back(i);
    }
    if (soft.empty()) return out;

    // Group soft pieces that run parallel within `spacing` (segments of one bevel).
    struct Sample { Vec3 mid, dir; int piece; };
    std::vector<Sample> samples;
    for (int i : soft)
        for (int e : pieces[i].edges) {
            const Vec3 &a = V[candidates[e] >> 32], &b = V[candidates[e] & 0xffffffff];
            Vec3 d{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (l <= 0) continue;
            samples.push_back({{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])},
                               {d[0] / l, d[1] / l, d[2] / l}, i});
        }
    auto cellOf = [&](const Vec3& p) {
        return std::array<int64_t, 3>{(int64_t)std::floor(p[0] / spacing), (int64_t)std::floor(p[1] / spacing),
                                      (int64_t)std::floor(p[2] / spacing)};
    };
    auto hash = [](const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) | (uint64_t)(c[2] & 0x1FFFFF);
    };
    std::unordered_map<uint64_t, std::vector<int>> grid;
    for (int i = 0; i < (int)samples.size(); ++i) grid[hash(cellOf(samples[i].mid))].push_back(i);
    std::unordered_map<uint64_t, int> links;  // (piece a, piece b) -> parallel close pairs
    for (int i = 0; i < (int)samples.size(); ++i) {
        const auto c = cellOf(samples[i].mid);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = grid.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        const Sample &x = samples[i], &y = samples[j];
                        if (x.piece >= y.piece) continue;
                        double d2 = 0, dot = 0;
                        for (int k = 0; k < 3; ++k)
                            d2 += (x.mid[k] - y.mid[k]) * (x.mid[k] - y.mid[k]), dot += x.dir[k] * y.dir[k];
                        if (d2 <= spacing * spacing && std::abs(dot) > 0.9)
                            ++links[((uint64_t)(uint32_t)x.piece << 32) | (uint32_t)y.piece];
                    }
                }
    }
    std::unordered_map<int, int> parent;
    for (int i : soft) parent[i] = i;
    std::function<int(int)> find = [&](int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); };
    for (const auto& [k, n] : links) {
        const int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
        // Enough close parallel pairs that the pieces really run side by side.
        const size_t shorter = std::min(pieces[a].edges.size(), pieces[b].edges.size());
        if (n >= 3 && n >= 0.3 * shorter) parent[find(a)] = find(b);
    }
    std::unordered_map<int, std::vector<int>> clusters;
    for (int i : soft) clusters[find(i)].push_back(i);

    for (auto& [root, members] : clusters) {
        if (members.size() < 2) continue;  // a lone soft piece is just a gentle bend
        double total = 0;
        for (int i : members) total += pieces[i].meanAngle;
        if (total <= angleDegrees) continue;
        // Middle piece: the one with the smallest summed distance to the others' samples.
        int best = members[0];
        double bestScore = INFINITY;
        for (int i : members) {
            double score = 0;
            for (int e : pieces[i].edges) {
                const Vec3 &a = V[candidates[e] >> 32], &b = V[candidates[e] & 0xffffffff];
                const Vec3 m{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
                for (int j : members) {
                    if (j == i) continue;
                    double nearest = INFINITY;
                    for (int f : pieces[j].edges) {
                        const Vec3 &c = V[candidates[f] >> 32], &d = V[candidates[f] & 0xffffffff];
                        double d2 = 0;
                        for (int k = 0; k < 3; ++k) d2 += std::pow(m[k] - 0.5 * (c[k] + d[k]), 2);
                        nearest = std::min(nearest, d2);
                    }
                    score += std::sqrt(nearest);
                }
            }
            score /= pieces[i].edges.size();
            if (score < bestScore) bestScore = score, best = i;
        }
        emit(pieces[best]);
    }
    return out;
}

}  // namespace remersher
