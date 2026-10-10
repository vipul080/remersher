#pragma once

#include "mesh_cleanup.h"
#include "remersher/mesh.h"

namespace remersher {

// Keeps the part of the mesh with coordinate[axis] >= 0, cutting triangles that cross the plane.
// New vertices on the cut lie exactly on the plane and are shared between neighbouring triangles.
TriangleMesh clipToHalfSpace(const TriangleMesh& mesh, int axis);

// Mirrors a mesh built on the positive side of the plane coordinate[axis] = 0 and joins the two
// halves. Border vertices within seamTolerance of the plane are moved onto it and shared with
// their mirror image; the mirrored faces get reversed winding so normals stay consistent.
Mesh mirrorAcross(const Mesh& half, int axis, double seamTolerance);

}  // namespace remersher
