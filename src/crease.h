#pragma once

#include <array>
#include <vector>

#include "mesh_cleanup.h"

namespace remersher {

// Hard edges of a triangle mesh, as segments, judged at the scale of the output quads.
//
// Candidate edges are grouped into connected pieces, and a piece must look like a line: at least
// 4 edges, loose ends at no more than a quarter of its vertices and a low mean vertex degree.
// Real creases form lines and closed loops (cube edges, the rim of a boolean-cut dent); the facet
// angles of a coarse smooth surface form scattered bits or a dense net.
//
// A piece sharper than angleDegrees is a hard edge. With spacing > 0, pieces down to a lower
// angle are also considered, and pieces running parallel within `spacing` of each other are
// treated as the segments of one small bevel or fillet: their angles add up, and if the total
// exceeds angleDegrees only the middle piece is returned, so the bevel gets one edge loop rather
// than a band of thin quads.
std::vector<std::array<std::array<double, 3>, 2>> detectCreases(const TriangleMesh& mesh,
                                                                 double angleDegrees,
                                                                 double spacing = 0);

}  // namespace remersher
