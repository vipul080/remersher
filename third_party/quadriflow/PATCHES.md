# Local changes to QuadriFlow

QuadriFlow is vendored (see UPSTREAM.txt) because remersher modifies the solver. Changes:

- `src/main.cpp`, `src/Optimizer.cu`: removed; remersher has its own CLI and no CUDA path.
- `src/parametrizer-mesh.cpp` (`FixValence`): the vertex fan walk stops at already-visited
  edges, fixing an infinite loop on non-manifold intermediate connectivity.
- `src/parametrizer-flip.cpp` (`FixHoles`): the boundary fan walk is capped; a corrupted loop is
  skipped instead of spinning forever.
- `src/optimizer.cpp` (`optimize_scale`): the sizing-field system (SPD, well conditioned) is
  solved with conjugate gradients instead of sparse LU, which was ~65% of solve time.
- `src/hierarchy.cpp` (`propagateConstraints`): soft constraint weights are averaged onto
  coarser levels instead of being forced to 1; hard weights (borders) behave as before.
- `src/parametrizer.hpp`, `src/parametrizer-mesh.cpp`: the hard-edge angle (was fixed at 60 degrees)
  is a member, `sharp_angle_degrees`.
- `3rd/MapleCOMSPS_LRB`: not vendored (only used by an optional SAT path that is not built).
