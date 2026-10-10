# Remersher

An open-source automatic quad remesher. Give it a triangle (or mixed) mesh and a target quad
count; it returns a clean, all-quad mesh whose edge loops follow the shape's curvature and hard
edges.

> Status: early (v0.1). The engine works on the bundled test meshes; the Blender add-on,
> vertex-color density, symmetry and material/normal splitting are on the roadmap below.

## Build

Requirements: CMake ≥ 3.16, a C++17 compiler, Eigen 3 and Boost headers.

```bash
# macOS
brew install cmake eigen boost
# Debian / Ubuntu
sudo apt-get install cmake libeigen3-dev libboost-dev
```

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Usage

```bash
./build/remersher -i meshes/bumpy.obj -o bumpy_quads.obj --target 2000
```

| option | default | meaning |
|---|---|---|
| `-t, --target <n>` | 5000 | target quad count |
| `--passes <n>` | 6 | extra solves used to land close to the target count |
| `-j, --jobs <n>` | 3 | calibration solves run in parallel (macOS / Linux) |
| `--tolerance <f>` | 0.03 | stop calibrating once within this fraction of the target |
| `--no-adaptive` | off | uniform quad size instead of curvature-adaptive |
| `--no-hard-edges` | off | do not align edge loops to sharp creases |
| `--no-boundary` | off | do not constrain open borders |
| `--resample <mode>` | auto | `auto`, `always` or `never`: isotropic resampling of the input before solving; `auto` only resamples sliver-heavy input |
| `--seed <n>` | 0 | solver seed; results are deterministic per seed |
| `--timeout <s>` | 120 | per-solve time limit (0 = none) |

On macOS and Linux each solve runs in an isolated child process: if a solve hangs, crashes or
returns invalid geometry it is killed and retried with another seed, so the caller (for example
Blender) never freezes.

## Blender add-on

```bash
./scripts/package_addon.sh        # -> dist/remersher_blender-<platform>.zip (CLI bundled inside)
```

In Blender: *Edit › Preferences › Add-ons › Install from Disk* and pick the zip (Blender 3.6+).
The **Remersher** tab in the 3D view sidebar has Target Quad Count, Adaptive Size, Detect Hard
Edges, Preserve Borders and Seed. *Remesh* runs in the background (Esc cancels), remeshes the
active object with its modifiers applied, and adds the result as a new object with the same
transform and materials. CI builds the zip for Linux and macOS on every push.

## How it works

1. **Input clean-up**: concave-safe ngon triangulation (ear clipping), tolerance welding of
   seams, removal of degenerate/duplicate faces, and consistent winding per connected piece.
2. **Isotropic resampling** (Botsch & Kobbelt) when the input has many slivers: edge splits,
   collapses, valence flips and smoothing projected onto the input; hard edges and borders are
   kept. This is what makes dense UV-sphere-like inputs stable.
3. **Field-aligned parametrization** (QuadriFlow, Huang et al. 2018): a 4-RoSy orientation field
   aligned to principal curvature and sharp edges, an optional curvature-adaptive scale field,
   and a position field, solved on a multi-resolution hierarchy.
4. **Quad extraction** with network-flow based singularity and flip removal.
5. **Feature snapping**: output vertices near input corners, hard edges and borders are moved
   exactly onto them, growing chains along creases from the corners.
6. **Validation**: every result is checked against the input surface (stray vertices, holes,
   folded faces); bad solves are rejected and retried with another seed; if every attempt fails, the
   original triangulation is tried and finally the least-bad result is returned with a warning.
7. **Count calibration**: the solve is repeated with a corrected face budget until the result is
   within `--tolerance` of the target.

## Benchmark

`bench/` scores remersher's output on objective quality metrics, and can score any other
remesher's output on the same cases for a head-to-head comparison:

| metric | meaning |
|---|---|
| `target_err_pct` | distance of the face count from the target |
| `nonquad_pct` | share of faces that are not quads |
| `irregular_pct` | share of vertices that are poles (valence ≠ 4 inside, ≠ 3 on borders) |
| `angle_dev_mean`, `angle_dev_p95` | how far quad corners are from 90° |
| `dev_mean_pct`, `hausdorff_pct` | mean / max distance to the input surface (% of bbox diagonal) |
| `sharp_dev_pct` | how far the input's hard edges are from the nearest output edge |
| `folded_pct`, `nonmanifold_edges` | broken output: edges where the mesh folds over itself, edges with more than two faces |

```bash
python3 -m venv .venv && .venv/bin/pip install numpy scipy
.venv/bin/python bench/run_bench.py            # writes bench/results/<timestamp>/report.md
```

To compare against another remesher, run it on `meshes/` with the settings in
`bench/cases.py` and put its outputs in `bench/reference/<mesh>__<paramset>.obj`
(`bench/import_reference.py` copies them from a `<mesh>__<paramset>/out.obj` folder layout).
Reference outputs are git-ignored and never committed.

## Roadmap

- [x] Engine + CLI, count calibration, crash/hang isolation, input clean-up, result validation
- [x] Benchmark harness
- [x] Stable solves on sliver-heavy input (isotropic resampling)
- [x] Remove non-manifold edges and folded faces from the output
- [x] Snap edge loops onto hard edges, corners and borders
- [ ] Proper corner singularities (valence-3 poles on cube-like corners)
- [x] Speed: CG sizing solve, faster resampling, parallel calibration (median case 3 s)
- [x] Blender add-on (target count, adaptive size, hard edges, one-click remesh)
- [ ] Continuous adaptivity (0–100) instead of on/off
- [ ] Vertex-color density painting
- [ ] Symmetry (X / Y / Z)
- [ ] Split by materials / normals
- [ ] Windows build and prebuilt release binaries

## License

MIT, see [LICENSE](LICENSE). Bundled third-party code is listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
