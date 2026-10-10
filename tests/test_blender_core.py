"""Tests for the Blender add-on's bpy-free helpers, including a round trip through the real CLI.

usage: python3 tests/test_blender_core.py [path/to/remersher]
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Load core.py directly: importing the package would run __init__.py, which needs bpy.
_spec = importlib.util.spec_from_file_location(
    "remersher_core_helpers", os.path.join(ROOT, "blender", "remersher_blender", "core.py"))
core = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = core  # dataclasses look the module up here
_spec.loader.exec_module(core)


def test_obj_round_trip(tmp):
    verts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)]
    faces = [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4]]
    path = os.path.join(tmp, "rt.obj")
    core.write_obj(path, verts, faces)
    v, f = core.read_obj(path)
    assert v == [tuple(map(float, p)) for p in verts], v
    assert f == faces, f


def test_build_command():
    s = core.RemeshSettings(target_quad_count=1234, adaptivity=0, detect_hard_edges=True,
                            preserve_boundary=False, seed=7, timeout=30)
    cmd = core.build_command("/bin/remersher", "a.obj", "b.obj", s)
    assert cmd[:5] == ["/bin/remersher", "-i", "a.obj", "-o", "b.obj"], cmd
    assert cmd[cmd.index("--target") + 1] == "1234"
    assert cmd[cmd.index("--seed") + 1] == "7"
    assert cmd[cmd.index("--adaptivity") + 1] == "0"
    assert "--no-boundary" in cmd
    assert cmd[cmd.index("--hard-angle") + 1] == "45" and "--no-hard-edges" not in cmd
    assert "--symmetry" not in cmd
    cmd = core.build_command("/bin/remersher", "a.obj", "b.obj", core.RemeshSettings(symmetry="Xz"))
    assert cmd[cmd.index("--symmetry") + 1] == "xz"


def test_find_binary(tmp, binary):
    assert core.find_binary(tmp, binary) == binary
    # A bundled copy that lost its executable bit (zip install) is repaired and used.
    bundled = os.path.join(tmp, "bin", core.BINARY_NAME)
    os.makedirs(os.path.dirname(bundled), exist_ok=True)
    with open(binary, "rb") as src, open(bundled, "wb") as dst:
        dst.write(src.read())
    os.chmod(bundled, 0o644)
    assert core.find_binary(tmp) == bundled and os.access(bundled, os.X_OK)
    os.remove(bundled)
    assert core.find_binary(tmp, os.path.join(tmp, "missing")) in (None, core.shutil.which("remersher"))


def test_cli_round_trip(tmp, binary):
    src = os.path.join(ROOT, "meshes", "torus.obj")
    verts, faces = core.read_obj(src)
    in_path, out_path = os.path.join(tmp, "in.obj"), os.path.join(tmp, "out.obj")
    core.write_obj(in_path, verts, faces)
    cmd = core.build_command(binary, in_path, out_path, core.RemeshSettings(target_quad_count=800))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, core.summarize_stderr(proc.stderr)
    v, f = core.read_obj(out_path)
    assert len(f) > 400 and all(len(x) == 4 for x in f), (len(f), {len(x) for x in f})
    assert all(0 <= i < len(v) for x in f for i in x)


def main():
    binary = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "remersher")
    if not os.path.isfile(binary):
        sys.exit(f"remersher binary not found at {binary}; build it first")
    with tempfile.TemporaryDirectory() as tmp:
        test_obj_round_trip(tmp)
        test_build_command()
        test_find_binary(tmp, binary)
        test_cli_round_trip(tmp, binary)
    print("blender core tests passed")


if __name__ == "__main__":
    main()
