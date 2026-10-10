#pragma once

#include "mesh_cleanup.h"
#include "remersher/mesh.h"

namespace remersher {

struct SnapStats {
    int corners = 0;   // output vertices moved onto input feature corners
    int creases = 0;   // output vertices moved onto input hard edges
    int borders = 0;   // output border vertices moved onto input borders
    int reverted = 0;  // snaps undone because they would fold a face
    int diagonalCollapses = 0, edgeRotations = 0;  // pole clean-up operations
};

// Moves output vertices that lie near input features exactly onto them: feature corners first,
// then hard edges (only for vertices whose edge loop runs along the crease) and open borders.
// The solver aligns edge loops with features but leaves them up to about half a quad away,
// which rounds off corners and creases. hardEdgeAngle <= 0 snaps borders only.
// Then (optionally) cancels pole pairs with local quad operations (see cancelPoles) and relaxes
// the quad mesh for relaxIterations rounds (vertices kept on the surface and on their
// features), which evens out quad shapes left distorted by the solver.
SnapStats snapToFeatures(Mesh& output, const TriangleMesh& input, double hardEdgeAngle,
                         int relaxIterations = 0, bool cleanupPoles = false);

}  // namespace remersher
