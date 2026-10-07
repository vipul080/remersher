#pragma once

#include <string>

#include "mesh_cleanup.h"
#include "remersher/mesh.h"

namespace remersher {

// Cheap sanity metrics of a remeshed result against its (cleaned) input. Distances are fractions
// of the input's bounding-box diagonal and use nearest-vertex queries, so they are approximate
// to within the input's edge length; good enough to catch collapsed or folded solves.
struct QualityCheck {
    double farOutputFraction = 0;  // output vertices far from the input surface
    double uncoveredFraction = 0;  // input vertices far from any output vertex (holes)
    double flippedFraction = 0;    // output faces facing against the input normal
    double angleDeviation = 0;     // mean |corner angle - 90°| over quads, in degrees

    bool acceptable(std::string* why = nullptr) const;
};

QualityCheck checkQuality(const TriangleMesh& input, const Mesh& output);

}  // namespace remersher
