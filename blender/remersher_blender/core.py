"""Blender-independent helpers for the Remersher add-on: locating the CLI, building its command
line, and moving meshes in and out of it as OBJ files. Kept free of `bpy` so it can be tested
with plain Python."""
from __future__ import annotations

import os
import shutil
import sys
from dataclasses import dataclass

BINARY_NAME = "remersher.exe" if sys.platform == "win32" else "remersher"


@dataclass
class RemeshSettings:
    target_quad_count: int = 5000
    adaptivity: int = 50  # 0 = uniform quad size, 100 = strongest curvature adaptivity
    detect_hard_edges: bool = True
    hard_edge_angle: float = 45.0
    preserve_boundary: bool = True
    seed: int = 0
    symmetry: str = ""  # any of "x", "y", "z"
    use_vertex_color: bool = False  # vertex colors as density paint
    timeout: float = 120.0


def find_binary(addon_dir: str, preferred: str = "") -> str | None:
    """Return the remersher executable to use: the path set in preferences, then a copy bundled
    in the add-on's bin/ folder, then one on PATH."""
    candidates = []
    if preferred:
        candidates.append(os.path.expanduser(preferred))
    candidates.append(os.path.join(addon_dir, "bin", BINARY_NAME))
    for path in candidates:
        if not os.path.isfile(path):
            continue
        if not os.access(path, os.X_OK):
            # Blender installs add-ons with zipfile, which drops the executable bit.
            try:
                os.chmod(path, os.stat(path).st_mode | 0o111)
            except OSError:
                continue
        if os.access(path, os.X_OK):
            return path
    return shutil.which("remersher")


def build_command(binary: str, in_path: str, out_path: str, s: RemeshSettings) -> list[str]:
    cmd = [binary, "-i", in_path, "-o", out_path, "--target", str(int(s.target_quad_count)),
           "--seed", str(int(s.seed)), "--timeout", str(float(s.timeout))]
    axes = "".join(a for a in "xyz" if a in s.symmetry.lower())
    if axes:
        cmd += ["--symmetry", axes]
    cmd += ["--adaptivity", str(int(s.adaptivity))]
    if not s.detect_hard_edges:
        cmd.append("--no-hard-edges")
    else:
        cmd += ["--hard-angle", "%g" % s.hard_edge_angle]
    if not s.preserve_boundary:
        cmd.append("--no-boundary")
    if s.use_vertex_color:
        cmd.append("--vertex-color")
    return cmd


def write_obj(path: str, vertices, faces, colors=None, feature_edges=None) -> None:
    """vertices: iterable of (x, y, z); faces: iterable of vertex-index sequences (0-based);
    colors: optional per-vertex (r, g, b) in 0..1, written as `v x y z r g b`;
    feature_edges: optional (a, b) vertex pairs that must become edge loops, written as `l a b`."""
    with open(path, "w") as fh:
        fh.write("# remersher blender export\n")
        if colors is None:
            fh.writelines("v %.9g %.9g %.9g\n" % (v[0], v[1], v[2]) for v in vertices)
        else:
            fh.writelines("v %.9g %.9g %.9g %.4f %.4f %.4f\n" % (v[0], v[1], v[2], c[0], c[1], c[2])
                          for v, c in zip(vertices, colors))
        fh.writelines("f " + " ".join(str(i + 1) for i in f) + "\n" for f in faces)
        if feature_edges is not None:
            fh.writelines("l %d %d\n" % (a + 1, b + 1) for a, b in feature_edges)


def read_obj(path: str) -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    vertices, faces = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                x, y, z = line.split()[1:4]
                vertices.append((float(x), float(y), float(z)))
            elif line.startswith("f "):
                face = []
                for tok in line.split()[1:]:
                    i = int(tok.split("/")[0])
                    face.append(i - 1 if i > 0 else len(vertices) + i)
                faces.append(face)
    return vertices, faces


def summarize_stderr(text: str, limit: int = 200) -> str:
    """Last meaningful line of the CLI's stderr, for error reports in the UI."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return (lines[-1] if lines else "unknown error")[:limit]
