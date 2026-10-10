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
    # Density paint: a sphere painted white on top and black below gets more quads on top.
    bpy.ops.mesh.primitive_uv_sphere_add(segments=48, ring_count=24, location=(5, 0, 0))
    ball = bpy.context.active_object
    col = ball.data.color_attributes.new("paint", "FLOAT_COLOR", "POINT")
    for i, v in enumerate(ball.data.vertices):
        c = 1.0 if v.co.z > 0 else 0.0
        col.data[i].color_srgb = (c, c, c, 1.0)
    ball.data.color_attributes.active_color = col
    s.symmetry_x = False
    s.use_vertex_color = True
    s.target_quad_count = 1200
    assert bpy.ops.remersher.remesh() == {"FINISHED"}
    painted = bpy.context.active_object
    top = sum(1 for p in painted.data.polygons if p.center.z > 0)
    bottom = len(painted.data.polygons) - top
    assert top > 2 * bottom, (top, bottom)
    # Materials: a sphere whose top cap has a second material keeps the border as an edge loop.
    bpy.ops.mesh.primitive_uv_sphere_add(segments=48, ring_count=24, location=(10, 0, 0))
    mat_ball = bpy.context.active_object
    for name in ("A", "B"):
        mat_ball.data.materials.append(bpy.data.materials.new(name))
    for p in mat_ball.data.polygons:
        p.material_index = 1 if p.center.z > 0.45 else 0
    owners = {}
    for p in mat_ball.data.polygons:
        for v in p.vertices:
            owners.setdefault(v, set()).add(p.material_index)
    border_z = sum(mat_ball.data.vertices[v].co.z for v, m in owners.items() if len(m) == 2)
    border_z /= sum(1 for m in owners.values() if len(m) == 2)
    s.use_vertex_color = False
    s.use_materials = True
    s.target_quad_count = 1200
    assert bpy.ops.remersher.remesh() == {"FINISHED"}
    split = bpy.context.active_object
    assert len(split.data.materials) == 2
    used = {p.material_index for p in split.data.polygons}
    assert used == {0, 1}, used
    # Faces with the cap material sit above the border, the others below it.
    wrong = sum(1 for p in split.data.polygons if (p.material_index == 1) != (p.center.z > border_z))
    assert wrong < 0.03 * len(split.data.polygons), wrong
    print(f"BLENDER SMOKE OK: {len(faces)} faces, {quads} quads; density paint {top} top / {bottom} bottom; "
          f"materials {len(split.data.polygons)} faces")


try:
    main()
except Exception:
    traceback.print_exc()
    sys.exit(1)
