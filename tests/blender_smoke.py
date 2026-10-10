"""End-to-end check of the Blender add-on inside Blender (background mode).

usage: blender -b --factory-startup --python tests/blender_smoke.py -- <add-on zip>
Exits non-zero on failure.
"""
import sys
import traceback

import bpy


def main():
    zip_path = sys.argv[sys.argv.index("--") + 1]
    bpy.ops.preferences.addon_install(filepath=zip_path, overwrite=True)
    bpy.ops.preferences.addon_enable(module="remersher_blender")

    # Suzanne with a subdivision modifier: the add-on must remesh the evaluated (modified) mesh.
    bpy.ops.mesh.primitive_monkey_add()
    src = bpy.context.active_object
    mod = src.modifiers.new("subd", "SUBSURF")
    mod.levels = 2
    src.location = (1.0, 2.0, 3.0)

    s = bpy.context.scene.remersher
    s.target_quad_count = 1500
    s.symmetry_x = True
    result = bpy.ops.remersher.remesh()
    assert result == {"FINISHED"}, result

    out = bpy.context.active_object
    assert out is not src and out.name == "Retopo_" + src.name, out.name
    faces = out.data.polygons
    quads = sum(1 for p in faces if p.loop_total == 4)
    assert len(faces) > 700 and quads / len(faces) > 0.98, (len(faces), quads)
    assert tuple(out.location) == tuple(src.location), (out.location, src.location)
    assert src.hide_get(), "source should be hidden"
    xs = [v.co.x for v in out.data.vertices]
    assert abs(max(xs) + min(xs)) < 1e-3, "symmetric result expected"
    print(f"BLENDER SMOKE OK: {len(faces)} faces, {quads} quads")


try:
    main()
except Exception:
    traceback.print_exc()
    sys.exit(1)
