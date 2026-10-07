# Local changes to QuadriFlow

QuadriFlow is vendored (see UPSTREAM.txt) because remersher modifies the solver. Changes:

- `src/main.cpp`, `src/Optimizer.cu`: removed; remersher has its own CLI and no CUDA path.
- `src/parametrizer-mesh.cpp` (`FixValence`): the vertex fan walk stops at already-visited
  edges, fixing an infinite loop on non-manifold intermediate connectivity.
- `src/parametrizer-flip.cpp` (`FixHoles`): the boundary fan walk is capped; a corrupted loop is
  skipped instead of spinning forever.
- `3rd/MapleCOMSPS_LRB`: not vendored (only used by an optional SAT path that is not built).
