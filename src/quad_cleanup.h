#pragma once

#include <vector>

#include "remersher/mesh.h"

namespace remersher {

struct PoleCleanupStats {
    int diagonalCollapses = 0;
    int edgeRotations = 0;
};

// Greedy local clean-up of irregular vertices in a quad(-dominant) mesh. Applies diagonal
// collapses (merging opposite corners of a quad, e.g. two valence-3 poles into one regular
// vertex) and edge rotations (re-pairing two adjacent quads) whenever they lower the total
// pole count (then sum |valence - 4|) of the vertices involved and keep the faces unfolded.
// Vertices flagged in `locked` (features, borders) are never moved, merged, or have their edges
// rotated. Merged-away vertices stay in the vertex list unreferenced.
// If `changed` is given, it is sized to the vertex count and flags every vertex whose faces were
// modified, so callers can re-smooth just those areas.
PoleCleanupStats cancelPoles(Mesh& mesh, const std::vector<char>& locked, int rounds = 10,
                             std::vector<char>* changed = nullptr);

}  // namespace remersher
