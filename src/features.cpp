#include "features.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

struct Segment {
    Vec3 a, b;
    Vec3 dir;  // unit
};

Vec3 closestOnSegment(const Segment& s, const Vec3& p) {
    Vec3 ab = sub(s.b, s.a);
    double len2 = dot(ab, ab);
    double t = len2 > 0 ? std::clamp(dot(sub(p, s.a), ab) / len2, 0.0, 1.0) : 0.0;
    return add(s.a, scale(ab, t));
}

// Uniform grid over segments for nearest-point queries.
class SegmentGrid {
   public:
    SegmentGrid(std::vector<Segment> segs, double cell) : segs_(std::move(segs)), cell_(cell) {
        for (int i = 0; i < (int)segs_.size(); ++i) {
            auto c0 = cellOf(segs_[i].a), c1 = cellOf(segs_[i].b);
            for (int j = 0; j < 3; ++j)
                if (c0[j] > c1[j]) std::swap(c0[j], c1[j]);
            double count = (double)(c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
            if (count > 4096) { large_.push_back(i); continue; }
            for (int64_t x = c0[0]; x <= c1[0]; ++x)
                for (int64_t y = c0[1]; y <= c1[1]; ++y)
                    for (int64_t z = c0[2]; z <= c1[2]; ++z) cells_[hash({x, y, z})].push_back(i);
        }
    }

    bool empty() const { return segs_.empty(); }

    // Nearest point on any segment within maxDist; returns the segment index or -1.
    int nearest(const Vec3& p, double maxDist, Vec3& out) const {
        double best = maxDist * maxDist;
        int bestSeg = -1;
        auto consider = [&](int i) {
            Vec3 q = closestOnSegment(segs_[i], p);
            Vec3 d = sub(q, p);
            if (dot(d, d) <= best) best = dot(d, d), bestSeg = i, out = q;
        };
        const auto c = cellOf(p);
        const int r = std::max(1, (int)std::ceil(maxDist / cell_));
        for (int dx = -r; dx <= r; ++dx)
            for (int dy = -r; dy <= r; ++dy)
                for (int dz = -r; dz <= r; ++dz) {
                    auto it = cells_.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                    if (it != cells_.end())
                        for (int i : it->second) consider(i);
                }
        for (int i : large_) consider(i);
        return bestSeg;
    }

    const Segment& operator[](int i) const { return segs_[i]; }

   private:
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {(int64_t)std::floor(p[0] / cell_), (int64_t)std::floor(p[1] / cell_),
                (int64_t)std::floor(p[2] / cell_)};
    }
    static uint64_t hash(const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) |
               (uint64_t)(c[2] & 0x1FFFFF);
    }
    std::vector<Segment> segs_;
    double cell_;
    std::unordered_map<uint64_t, std::vector<int>> cells_;
    std::vector<int> large_;
};

Vec3 newellNormal(const Mesh& m, const std::vector<int>& f) {
    Vec3 n{0, 0, 0};
    for (size_t i = 0; i < f.size(); ++i) n = add(n, cross(m.vertices[f[i]], m.vertices[f[(i + 1) % f.size()]]));
    return n;
}

}  // namespace

SnapStats snapToFeatures(Mesh& output, const TriangleMesh& input, double hardEdgeAngle) {
    SnapStats st;
    if (output.faces.empty() || input.triangles.empty()) return st;

    // --- input features ---------------------------------------------------------------------
    std::unordered_map<uint64_t, std::vector<int>> edgeTris;
    for (int t = 0; t < (int)input.triangles.size(); ++t)
        for (int i = 0; i < 3; ++i)
            edgeTris[edgeKey(input.triangles[t][i], input.triangles[t][(i + 1) % 3])].push_back(t);
    auto triNormal = [&](int t) {
        const auto& T = input.triangles[t];
        Vec3 n = cross(sub(input.vertices[T[1]], input.vertices[T[0]]), sub(input.vertices[T[2]], input.vertices[T[0]]));
        double l = norm(n);
        return l > 0 ? scale(n, 1 / l) : n;
    };
    const double cosHard = hardEdgeAngle > 0 ? std::cos(hardEdgeAngle * 3.14159265358979323846 / 180.0) : -2;

    std::vector<Segment> creases, borders;
    std::vector<std::vector<Vec3>> featureDirs(input.vertices.size());
    for (const auto& [k, ts] : edgeTris) {
        const int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
        const bool border = ts.size() == 1;
        const bool crease = ts.size() == 2 && dot(triNormal(ts[0]), triNormal(ts[1])) < cosHard;
        if (!border && !crease) continue;
        Vec3 d = sub(input.vertices[b], input.vertices[a]);
        double l = norm(d);
        if (l <= 0) continue;
        Segment s{input.vertices[a], input.vertices[b], scale(d, 1 / l)};
        (border ? borders : creases).push_back(s);
        featureDirs[a].push_back(s.dir);
        featureDirs[b].push_back(scale(s.dir, -1));
    }
    if (creases.empty() && borders.empty()) return st;

    // Feature corners: where feature lines meet, end, or turn sharply.
    std::vector<Vec3> corners;
    for (size_t v = 0; v < featureDirs.size(); ++v) {
        const auto& d = featureDirs[v];
        if (d.empty()) continue;
        if (d.size() != 2 || dot(d[0], d[1]) > -std::cos(30.0 * 3.14159265358979323846 / 180.0))
            corners.push_back(input.vertices[v]);
    }

    // --- output adjacency -------------------------------------------------------------------
    const size_t nv = output.vertices.size();
    std::vector<std::vector<int>> nbrs(nv), vfaces(nv);
    std::unordered_map<uint64_t, int> edgeCount;
    for (int fi = 0; fi < (int)output.faces.size(); ++fi) {
        const auto& f = output.faces[fi];
        for (size_t i = 0; i < f.size(); ++i) {
            int a = f[i], b = f[(i + 1) % f.size()];
            nbrs[a].push_back(b), nbrs[b].push_back(a);
            vfaces[a].push_back(fi);
            ++edgeCount[edgeKey(a, b)];
        }
    }
    std::vector<char> onBorder(nv, 0);
    for (const auto& [k, c] : edgeCount)
        if (c == 1) onBorder[k >> 32] = onBorder[k & 0xffffffff] = 1;
    std::vector<double> h(nv, 0);
    double meanEdge = 0;
    size_t edges = 0;
    for (size_t v = 0; v < nv; ++v) {
        auto& n = nbrs[v];
        std::sort(n.begin(), n.end());
        n.erase(std::unique(n.begin(), n.end()), n.end());
        for (int u : n) h[v] += norm(sub(output.vertices[u], output.vertices[v]));
        if (!n.empty()) h[v] /= n.size(), meanEdge += h[v], ++edges;
    }
    if (!edges) return st;
    meanEdge /= edges;

    const std::vector<Vec3> original = output.vertices;
    std::vector<char> moved(nv, 0);

    // Corners: each corner pulls in its nearest output vertex (one vertex per corner).
    {
        std::vector<Segment> pts;  // degenerate segments = points, to reuse the grid
        for (size_t v = 0; v < nv; ++v) pts.push_back({original[v], original[v], {0, 0, 0}});
        SegmentGrid outGrid(pts, meanEdge);
        std::vector<double> claim(nv, INFINITY);
        std::vector<Vec3> target(nv);
        for (const auto& c : corners) {
            Vec3 q;
            int v = outGrid.nearest(c, 0.8 * meanEdge, q);
            if (v < 0 || nbrs[v].empty()) continue;
            double d = norm(sub(original[v], c));
            if (d > 0.75 * h[v] || d >= claim[v]) continue;
            claim[v] = d;
            target[v] = c;
        }
        for (size_t v = 0; v < nv; ++v)
            if (claim[v] < INFINITY) output.vertices[v] = target[v], moved[v] = 1, ++st.corners;
    }

    // Creases and borders.
    const SegmentGrid creaseGrid(creases, meanEdge), borderGrid(borders, meanEdge);
    const double cos25 = std::cos(25.0 * 3.14159265358979323846 / 180.0);
    for (size_t v = 0; v < nv; ++v) {
        if (moved[v] || nbrs[v].empty()) continue;
        Vec3 q;
        if (onBorder[v]) {
            if (!borderGrid.empty() && borderGrid.nearest(original[v], 0.5 * h[v], q) >= 0) {
                output.vertices[v] = q, moved[v] = 1, ++st.borders;
            }
            continue;
        }
        if (creaseGrid.empty()) continue;
        int s = creaseGrid.nearest(original[v], 0.4 * h[v], q);
        if (s < 0) continue;
        // Only vertices on an edge loop running along the crease: two incident edges parallel to it.
        int parallel = 0;
        for (int u : nbrs[v]) {
            Vec3 e = sub(original[u], original[v]);
            double l = norm(e);
            if (l > 0 && std::abs(dot(e, creaseGrid[s].dir)) > cos25 * l) ++parallel;
        }
        if (parallel >= 2) output.vertices[v] = q, moved[v] = 1, ++st.creases;
    }

    // Grow snapped chains along creases: a neighbour of a snapped vertex joins the crease if it is
    // close to it and the edge between them runs along the crease. This pulls in the vertices
    // next to corners, whose edges were not parallel to the crease before the corner snapped.
    if (!creaseGrid.empty()) {
        std::vector<int> queue;
        for (size_t v = 0; v < nv; ++v)
            if (moved[v] && !onBorder[v]) queue.push_back((int)v);
        while (!queue.empty()) {
            const int v = queue.back();
            queue.pop_back();
            for (int u : nbrs[v]) {
                if (moved[u] || onBorder[u]) continue;
                Vec3 q;
                int s = creaseGrid.nearest(original[u], 0.5 * h[u], q);
                if (s < 0) continue;
                Vec3 e = sub(q, output.vertices[v]);
                double l = norm(e);
                if (l < 0.3 * h[u] || std::abs(dot(e, creaseGrid[s].dir)) < cos25 * l) continue;
                output.vertices[u] = q, moved[u] = 1, ++st.creases;
                queue.push_back(u);
            }
        }
    }

    // Undo snaps that fold or collapse an incident face.
    for (size_t v = 0; v < nv; ++v) {
        if (!moved[v]) continue;
        bool bad = false;
        for (int fi : vfaces[v]) {
            const auto& f = output.faces[fi];
            Vec3 after = newellNormal(output, f);
            Vec3 before{0, 0, 0};
            for (size_t i = 0; i < f.size(); ++i) {
                const Vec3& a = f[i] == (int)v ? original[v] : output.vertices[f[i]];
                const Vec3& b = f[(i + 1) % f.size()] == (int)v ? original[v] : output.vertices[f[(i + 1) % f.size()]];
                before = add(before, cross(a, b));
            }
            double la = norm(after), lb = norm(before);
            if (la < 1e-3 * lb || dot(after, before) <= 0.2 * la * lb) { bad = true; break; }
        }
        if (bad) {
            output.vertices[v] = original[v];
            ++st.reverted;
        }
    }
    return st;
}

}  // namespace remersher
