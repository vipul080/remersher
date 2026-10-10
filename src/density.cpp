#include "density.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "curvature.h"

namespace remersher {

DensitySource::DensitySource(std::vector<std::array<double, 3>> points, std::vector<float> values)
    : points_(std::move(points)), values_(std::move(values)) {
    if (points_.empty()) return;
    std::array<double, 3> hi = points_[0];
    lo_ = points_[0];
    for (const auto& p : points_)
        for (int j = 0; j < 3; ++j) lo_[j] = std::min(lo_[j], p[j]), hi[j] = std::max(hi[j], p[j]);
    // About two points per cell on average.
    double volume = 1;
    for (int j = 0; j < 3; ++j) volume *= std::max(hi[j] - lo_[j], 1e-9);
    cell_ = std::cbrt(volume * 2.0 / points_.size());
    for (int j = 0; j < 3; ++j)
        dims_[j] = std::clamp((int)std::ceil((hi[j] - lo_[j]) / cell_) + 1, 1, 512);
    for (int j = 0; j < 3; ++j) cell_ = std::max(cell_, (hi[j] - lo_[j]) / (dims_[j] - 0.5));
    buckets_.assign((size_t)dims_[0] * dims_[1] * dims_[2], {});
    for (int i = 0; i < (int)points_.size(); ++i) {
        int c[3];
        for (int j = 0; j < 3; ++j) c[j] = std::clamp((int)((points_[i][j] - lo_[j]) / cell_), 0, dims_[j] - 1);
        buckets_[bucketOf(c[0], c[1], c[2])].push_back(i);
    }
}

float DensitySource::sample(const std::array<double, 3>& p) const {
    int c[3];
    for (int j = 0; j < 3; ++j) c[j] = std::clamp((int)((p[j] - lo_[j]) / cell_), 0, dims_[j] - 1);
    double best = std::numeric_limits<double>::infinity();
    int bestIdx = -1;
    const int maxRing = std::max({dims_[0], dims_[1], dims_[2]});
    for (int r = 0; r <= maxRing; ++r) {
        if (bestIdx >= 0 && (r - 1) * cell_ > std::sqrt(best)) break;
        for (int x = c[0] - r; x <= c[0] + r; ++x)
            for (int y = c[1] - r; y <= c[1] + r; ++y)
                for (int z = c[2] - r; z <= c[2] + r; ++z) {
                    if (std::max({std::abs(x - c[0]), std::abs(y - c[1]), std::abs(z - c[2])}) != r) continue;
                    if (x < 0 || y < 0 || z < 0 || x >= dims_[0] || y >= dims_[1] || z >= dims_[2]) continue;
                    for (int i : buckets_[bucketOf(x, y, z)]) {
                        double d = 0;
                        for (int j = 0; j < 3; ++j) d += (points_[i][j] - p[j]) * (points_[i][j] - p[j]);
                        if (d < best) best = d, bestIdx = i;
                    }
                }
    }
    return bestIdx >= 0 ? values_[bestIdx] : 0.5f;
}

void applyDensityPaint(qflow::Hierarchy& h, const DensitySource& src, double normScale,
                       const qflow::Vector3d& normOffset) {
    if (src.empty() || h.mV.empty() || h.mS.empty()) return;
    const auto& V = h.mV[0];
    const int n = (int)V.cols();
    std::vector<double> logSize(n);
    for (int v = 0; v < n; ++v) {
        const qflow::Vector3d p = V.col(v) * normScale + normOffset;
        const double paint = std::clamp((double)src.sample({p[0], p[1], p[2]}), 0.0, 1.0);
        logSize[v] = (1.0 - 2.0 * paint) * std::log(2.0);  // edge length 2^(1 - 2 * paint)
    }
    // A little smoothing so density changes over a few quads instead of in one step.
    std::vector<std::vector<int>> nbrs(n);
    for (int f = 0; f < h.mF.cols(); ++f)
        for (int i = 0; i < 3; ++i) {
            const int a = h.mF(i, f), b = h.mF((i + 1) % 3, f);
            nbrs[a].push_back(b), nbrs[b].push_back(a);
        }
    for (int round = 0; round < 3; ++round) {
        std::vector<double> next = logSize;
        for (int v = 0; v < n; ++v) {
            double sum = logSize[v];
            for (int j : nbrs[v]) sum += logSize[j];
            next[v] = sum / (1 + nbrs[v].size());
        }
        logSize.swap(next);
    }
    double mean = 0;
    for (double x : logSize) mean += std::exp(x);
    mean /= std::max(1, n);
    auto& S = h.mS[0];
    for (int v = 0; v < n; ++v) {
        const double k = std::exp(logSize[v]) / mean;
        S(0, v) *= k;
        S(1, v) *= k;
    }
    propagateSizing(h);
}

}  // namespace remersher
