#pragma once

#include <array>
#include <vector>

#include "remersher/mesh.h"

namespace remersher {

struct TriangleMesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
};

struct CleanupStats {
    int weldedVertices = 0;     // input vertices merged into another vertex
    int droppedTriangles = 0;   // degenerate or duplicate triangles removed
    int flippedTriangles = 0;   // triangles re-wound to match their neighbours
    int components = 0;         // connected pieces of the cleaned mesh
};

// Prepares an arbitrary polygon mesh for field-aligned remeshing:
//  - fan-triangulates polygons,
//  - welds vertices closer than weldTolerance * bounding-box diagonal (open seams, split normals),
//  - drops degenerate and duplicate triangles,
//  - makes winding consistent within each connected piece; closed pieces end up facing outwards,
//    open pieces keep the orientation most of their input faces had.
// Inconsistent winding matters: the solver pairs half-edges by direction, so a flipped face is
// seen as a hole border and drags the whole field towards it.
TriangleMesh cleanupForRemeshing(const Mesh& mesh, double weldTolerance = 1e-6,
                                 CleanupStats* stats = nullptr);

// Repairs solver output: welds vertices closer than weldTolerance * bounding-box diagonal (the
// solver duplicates vertices in place when it untangles non-manifold spots), collapses repeated
// consecutive corners, and removes folded faces (a corner repeated non-consecutively), duplicate
// faces and, on edges still shared by more than two faces, the smallest of those faces.
// Returns the number of faces removed.
int repairPolygonMesh(Mesh& mesh, double weldTolerance = 1e-7);

}  // namespace remersher
