#pragma once

#include "mesh_cleanup.h"

namespace remersher {

struct ResampleStats {
    int iterations = 0;
    int splits = 0, collapses = 0, flips = 0;
};

// Share of triangles that are slivers (smallest angle below 15 degrees). The field solver
// becomes unstable on inputs with many slivers, e.g. the pole regions of dense UV spheres.
double sliverFraction(const TriangleMesh& mesh);

// Isotropic remeshing (Botsch & Kobbelt 2004): repeated edge splits, collapses, valence-improving
// flips and tangential smoothing projected back onto the input surface, towards edge length
// targetEdge. Boundary edges and edges with a dihedral angle above hardEdgeAngle (degrees;
// <= 0 disables) are kept as features: they are split but never flipped, and their vertices are
// never smoothed or removed. Non-manifold edges and their vertices are left untouched.
TriangleMesh resampleIsotropic(const TriangleMesh& mesh, double targetEdge, double hardEdgeAngle,
                               int iterations = 5, ResampleStats* stats = nullptr);

}  // namespace remersher
