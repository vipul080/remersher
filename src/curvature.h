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

}  // namespace remersher
