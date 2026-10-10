#pragma once

#include <array>
#include <vector>

#include "mesh_cleanup.h"

namespace remersher {

// Hard edges of a triangle mesh, as segments. Edges whose dihedral angle exceeds angleDegrees are
// grouped into connected pieces, and a piece is kept only if it looks like a crease: at least 4
// edges and loose ends at no more than a quarter of its vertices. Real creases (cube edges, the
// rim of a boolean-cut dent) form lines and closed loops; the strong facet angles of a coarse
// smooth surface form scattered short pieces, which would otherwise all become hard edges at a
// low threshold.
std::vector<std::array<std::array<double, 3>, 2>> detectCreases(const TriangleMesh& mesh,
                                                                 double angleDegrees);

}  // namespace remersher
