#pragma once

#include "parametrizer.hpp"

namespace remersher {

// Adds soft orientation constraints to the finest level of the solver hierarchy that pull the
// 4-RoSy field towards principal curvature directions. QuadriFlow's field only minimises
// variation, so on surfaces like a torus its edge loops drift diagonally and need extra poles;
// curvature-aligned loops follow the shape the way an artist would lay them out.
//
// Each vertex's weight grows with how anisotropic its curvature is (|k1 - k2| relative to
// |k1| + |k2|) and how strong it is at the target quad size, so flat and umbilic (sphere-like)
// regions stay unconstrained. Existing constraints (open borders) are kept. `strength` scales
// all weights (0 disables). Returns the number of constrained vertices.
int addCurvatureConstraints(qflow::Hierarchy& h, double strength);

// Shrinks the solver's sizing field where the surface is strongly curved, so quads get smaller
// there and larger in flat regions. adaptivity is the 0..100 slider value; up to 50 the field is
// left untouched, above 50 curvature sizing is blended in. The field keeps a mean of 1, so the face budget still sets the overall density.
void applyCurvatureSizing(qflow::Hierarchy& h, double adaptivity);

// Copies the finest level's sizing field (mS[0]) down the hierarchy to the coarser levels.
void propagateSizing(qflow::Hierarchy& h);

}  // namespace remersher
