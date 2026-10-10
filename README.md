# Remersher

An open-source automatic quad remesher. Give it a triangle (or mixed) mesh and a target quad
count; it returns a clean, all-quad mesh whose edge loops follow the shape's curvature and hard
edges.

> Status: early (v0.1). CLI and Blender add-on work and are tested in CI (including inside
> Blender); quality is being compared case by case against a commercial reference.

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
| `-j, --jobs <n>` | 3 | calibration solves run in parallel (macOS, Linux) |
| `--tolerance <f>` | 0.03 | stop calibrating once within this fraction of the target |
| `--curvature <f>` | 0.5 | how strongly edge loops follow principal curvature (0 = off, 1 = strongest) |
| `--relax <n>` | 1 | tangential relaxation rounds on the final mesh (more = more even quads, less adaptive sizing) |
| `--no-pole-cleanup` | off | skip the local pole-pair cancellation on the final mesh |
| `--symmetry <axes>` | – | mirror symmetry across the object-space X / Y / Z planes, e.g. `x` or `xz` |
| `--adaptivity <0-100>` | 50 | quad size variation with curvature: 0 = uniform, 50 = mild default, 51–100 = increasingly denser in curved areas |
| `--hard-angle <deg>` | 50 | dihedral angle above which an edge can be hard; edges are kept as hard only where they form crease lines or loops, not where a coarse smooth surface is just faceted |
| `--no-hard-edges` | off | do not align edge loops to sharp creases |
| `--vertex-color` | off | use vertex colors as density paint: white = 4× denser, black = 4× sparser, mid-grey unchanged (OBJ `v x y z r g b`) |
| `--materials` | off | keep borders between OBJ materials (`usemtl`) as edge loops; edges given as OBJ `l` elements are always kept |
| `--no-boundary` | off | do not constrain open borders |
| `--resample <mode>` | auto | `auto`, `always` or `never`: isotropic resampling of the input before solving; `auto` only resamples sliver-heavy input |
| `--seed <n>` | 0 | solver seed; results are deterministic per seed |
| `--timeout <s>` | 120 | per-solve time limit (0 = none) |

On macOS and Linux each solve runs in an isolated child process: if a solve hangs, crashes or
returns invalid geometry it is killed and retried with another seed, so the caller (for example
Blender) never freezes.

## Blender add-on

Download `remersher_blender-<platform>.zip` from the
[releases page](https://github.com/vipul080/remersher/releases), or build it yourself:

```bash
python3 scripts/package_addon.py  # -> dist/remersher_blender-<platform>.zip (CLI bundled inside)
```

In Blender: *Edit › Preferences › Add-ons › Install from Disk* and pick the zip (Blender 3.6+).
The **Remersher** tab in the 3D view sidebar mirrors the familiar retopology workflow: Quad Count,
Adaptive Size, Detect Hard Edges by angle, Preserve Borders, Use Vertex Color (density paint with
the active color attribute), Use Materials, Use Normals Splitting (edges marked sharp), Symmetry
X/Y/Z and Seed, with a **Remesh It** button. It runs in the background (Esc cancels), remeshes
the active object with its modifiers applied, and adds the result as `Retopo_<name>` with the same
transform and materials (each new face gets the material of the face under it).

CI builds the add-on for macOS, Linux and Windows on every push and runs it inside Blender
(`tests/blender_smoke.py`).

## How it works

1. **Input clean-up**: concave-safe ngon triangulation (ear clipping), tolerance welding of
   seams, removal of degenerate/duplicate faces, and consistent winding per connected piece.
   Hard edges are detected here once: edges sharper than `--hard-angle` that form crease lines
   or loops (scattered facet edges of coarse smooth surfaces are ignored). Together with material
   borders and user edges they travel as feature lines through every later stage.
2. **Isotropic resampling** (Botsch & Kobbelt) when the input has many slivers: edge splits,
   collapses, valence flips and smoothing projected onto the input; hard edges and borders are
   kept. This is what makes dense UV-sphere-like inputs stable.
3. **Field-aligned parametrization** (QuadriFlow, Huang et al. 2018): a 4-RoSy orientation field
   aligned to sharp edges and, where the curvature direction is strong and coherent over a few
   quads (tubes, tori, limbs), softly to principal curvature; a sizing field (curvature,
   density paint); and a position field, solved on a multi-resolution hierarchy.
4. **Quad extraction** with network-flow based singularity and flip removal.
5. **Feature snapping, pole clean-up and relaxation**: output vertices near input corners, hard
   edges and borders are moved exactly onto them (growing chains along creases from the corners);
   pole pairs are cancelled with diagonal collapses and edge rotations where that does not distort
   the surrounding quads; then the quads are relaxed tangentially, with vertices kept on the
   surface and on their features.
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
- [x] Continuous adaptivity (0–100)
- [x] Vertex-color density painting
- [x] Symmetry (X / Y / Z)
- [x] Keep material borders and sharp (split-normal) edges as edge loops
- [x] Windows build (CI) and a release workflow publishing add-on zips for macOS, Linux and Windows

## License

MIT, see [LICENSE](LICENSE). Bundled third-party code is listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
