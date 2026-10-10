#pragma once

// Point-to-surface projection onto a triangle mesh, shared by resampling and relaxation.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace remersher {
namespace geom {

using Vec3 = std::array<double, 3>;
using Tri = std::array<int, 3>;

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// Closest point on triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5).
inline Vec3 closestOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return a;
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return b;
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return c;
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

// Projects points onto the original surface. Triangles are binned into a uniform grid; a query
// looks at the point's cell and its 26 neighbours, which is enough because smoothing moves a
// vertex by less than one target edge length (the cell size is two of them).
class SurfaceProjector {
   public:
    SurfaceProjector(const std::vector<Vec3>& V, const std::vector<Tri>& F, double cell)
        : V_(V), F_(F), cell_(cell) {
        for (int f = 0; f < (int)F.size(); ++f) {
            Vec3 lo = V[F[f][0]], hi = lo;
            for (int k = 1; k < 3; ++k)
                for (int j = 0; j < 3; ++j)
                    lo[j] = std::min(lo[j], V[F[f][k]][j]), hi[j] = std::max(hi[j], V[F[f][k]][j]);
            auto c0 = cellOf(lo), c1 = cellOf(hi);
            double count = (double)(c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
            if (count > 4096) { large_.push_back(f); continue; }  // checked on every query
            for (int64_t x = c0[0]; x <= c1[0]; ++x)
                for (int64_t y = c0[1]; y <= c1[1]; ++y)
                    for (int64_t z = c0[2]; z <= c1[2]; ++z) cells_[hash({x, y, z})].push_back(f);
        }
    }

    bool project(const Vec3& p, Vec3& out) const {
        double best = INFINITY;
        auto consider = [&](int f) {
            Vec3 q = closestOnTriangle(p, V_[F_[f][0]], V_[F_[f][1]], V_[F_[f][2]]);
            Vec3 d = q - p;
            double d2 = dot(d, d);
            if (d2 < best) best = d2, out = q;
        };
        auto c = cellOf(p);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = cells_.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                    if (it != cells_.end())
                        for (int f : it->second) consider(f);
                }
        for (int f : large_) consider(f);
        return best < INFINITY;
    }

   private:
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {(int64_t)std::floor(p[0] / cell_), (int64_t)std::floor(p[1] / cell_),
                (int64_t)std::floor(p[2] / cell_)};
    }
    static uint64_t hash(const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) |
               (uint64_t)(c[2] & 0x1FFFFF);
    }
    const std::vector<Vec3>& V_;
    const std::vector<Tri>& F_;
    double cell_;
    std::unordered_map<uint64_t, std::vector<int>> cells_;
    std::vector<int> large_;
};

}  // namespace geom
}  // namespace remersher
