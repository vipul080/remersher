#pragma once

#include "mesh_cleanup.h"
#include "remersher/mesh.h"

namespace remersher {

struct SnapStats {
    int corners = 0;   // output vertices moved onto input feature corners
    int creases = 0;   // output vertices moved onto input hard edges
    int borders = 0;   // output border vertices moved onto input borders
    int reverted = 0;  // snaps undone because they would fold a face
};

// Moves output vertices that lie near input features exactly onto them: feature corners first,
// then hard edges (only for vertices whose edge loop runs along the crease) and open borders.
// The solver aligns edge loops with features but leaves them up to about half a quad away,
// which rounds off corners and creases. hardEdgeAngle <= 0 snaps borders only.
SnapStats snapToFeatures(Mesh& output, const TriangleMesh& input, double hardEdgeAngle);

}  // namespace remersher
