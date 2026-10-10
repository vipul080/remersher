#pragma once

#include <array>
#include <vector>

#include "parametrizer.hpp"

namespace remersher {

// Density paint sampled from the original input, so it survives clean-up, resampling and the
// solver's own subdivision: every query takes the value of the nearest painted input vertex.
class DensitySource {
   public:
    DensitySource(std::vector<std::array<double, 3>> points, std::vector<float> values);
    float sample(const std::array<double, 3>& p) const;
    bool empty() const { return points_.empty(); }

   private:
    std::vector<std::array<double, 3>> points_;
    std::vector<float> values_;
    double cell_ = 1;
    std::vector<std::vector<int>> buckets_;
    std::array<double, 3> lo_{};
    std::array<int, 3> dims_{};
    int bucketOf(int x, int y, int z) const { return (z * dims_[1] + y) * dims_[0] + x; }
};

// Scales the solver's sizing field by painted density: paint 0.5 is neutral, 1 makes quads half
// as long (4x the density), 0 twice as long. The field keeps its mean, so the overall quad count
// is unchanged and density only moves between regions. normScale / normOffset map the solver's
// normalized coordinates back to input space.
void applyDensityPaint(qflow::Hierarchy& h, const DensitySource& src, double normScale,
                       const qflow::Vector3d& normOffset);

}  // namespace remersher
