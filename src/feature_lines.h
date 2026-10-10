#pragma once

// User-specified feature lines (material borders, sharp/split-normal edges) kept as positions, so
// they can be recognised again after clean-up, resampling and the solver's own subdivision have
// renumbered and split the mesh.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace remersher {

class FeatureLines {
   public:
    using Vec3 = std::array<double, 3>;
    struct Segment {
        Vec3 a, b, dir;  // dir is unit length
    };

    FeatureLines() = default;

    // tolerance: how far (absolute) a point may be from a line and still count as on it.
    FeatureLines(const std::vector<std::array<Vec3, 2>>& segments, double tolerance) : tol_(tolerance) {
        double total = 0;
        for (const auto& s : segments) {
            Vec3 d{s[1][0] - s[0][0], s[1][1] - s[0][1], s[1][2] - s[0][2]};
            const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (l <= 0) continue;
            segs_.push_back({s[0], s[1], {d[0] / l, d[1] / l, d[2] / l}});
            total += l;
        }
        if (segs_.empty()) return;
        cell_ = std::max(total / segs_.size(), 4 * tol_);
        for (int i = 0; i < (int)segs_.size(); ++i) {
            auto c0 = cellOf(segs_[i].a), c1 = cellOf(segs_[i].b);
            for (int j = 0; j < 3; ++j)
                if (c0[j] > c1[j]) std::swap(c0[j], c1[j]);
            for (int64_t x = c0[0] - 1; x <= c1[0] + 1; ++x)
                for (int64_t y = c0[1] - 1; y <= c1[1] + 1; ++y)
                    for (int64_t z = c0[2] - 1; z <= c1[2] + 1; ++z) cells_[hash({x, y, z})].push_back(i);
        }
    }

    bool empty() const { return segs_.empty(); }
    const std::vector<Segment>& segments() const { return segs_; }

    // True if the edge a-b runs along a feature line: both endpoints and the midpoint lie on the
    // lines and the edge is parallel to the line under its midpoint.
    bool covers(const Vec3& a, const Vec3& b) const {
        if (segs_.empty()) return false;
        const Vec3 m{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
        const int s = nearest(m);
        if (s < 0 || nearest(a) < 0 || nearest(b) < 0) return false;
        Vec3 d{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (l <= 0) return false;
        const double c = (d[0] * segs_[s].dir[0] + d[1] * segs_[s].dir[1] + d[2] * segs_[s].dir[2]) / l;
        return std::abs(c) > 0.95;
    }

   private:
    // Index of a segment within tol_ of p, or -1.
    int nearest(const Vec3& p) const {
        auto it = cells_.find(hash(cellOf(p)));
        if (it == cells_.end()) return -1;
        int best = -1;
        double bestD = tol_ * tol_;
        for (int i : it->second) {
            const auto& s = segs_[i];
            Vec3 ab{s.b[0] - s.a[0], s.b[1] - s.a[1], s.b[2] - s.a[2]};
            Vec3 ap{p[0] - s.a[0], p[1] - s.a[1], p[2] - s.a[2]};
            const double len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
            const double t = std::clamp((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / len2, 0.0, 1.0);
            double d = 0;
            for (int j = 0; j < 3; ++j) d += (ap[j] - t * ab[j]) * (ap[j] - t * ab[j]);
            if (d <= bestD) bestD = d, best = i;
        }
        return best;
    }
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {(int64_t)std::floor(p[0] / cell_), (int64_t)std::floor(p[1] / cell_),
                (int64_t)std::floor(p[2] / cell_)};
    }
    static uint64_t hash(const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) |
               (uint64_t)(c[2] & 0x1FFFFF);
    }

    std::vector<Segment> segs_;
    double tol_ = 0;
    double cell_ = 1;
    std::unordered_map<uint64_t, std::vector<int>> cells_;
};

}  // namespace remersher
