"""Builds a die (bevelled cube with spherical pips cut by booleans) and exports it as OBJ.

usage: blender -b --factory-startup --python scripts/make_dice.py -- <out.obj>
"""
import sys

import bpy

out = sys.argv[sys.argv.index("--") + 1]
bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete()

bpy.ops.mesh.primitive_cube_add(size=2)
die = bpy.context.active_object
bev = die.modifiers.new("bevel", "BEVEL")
bev.width, bev.segments = 0.25, 6
bpy.ops.object.modifier_apply(modifier="bevel")

# Pip layout per face (in face-local u, v), like a real die.
layouts = {1: [(0, 0)], 2: [(-1, -1), (1, 1)], 3: [(-1, -1), (0, 0), (1, 1)],
           4: [(-1, -1), (-1, 1), (1, -1), (1, 1)], 5: [(-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1)],
           6: [(-1, -1), (-1, 0), (-1, 1), (1, -1), (1, 0), (1, 1)]}
faces = {1: ((0, 0, 1), (1, 0, 0), (0, 1, 0)), 6: ((0, 0, -1), (1, 0, 0), (0, 1, 0)),
         2: ((1, 0, 0), (0, 1, 0), (0, 0, 1)), 5: ((-1, 0, 0), (0, 1, 0), (0, 0, 1)),
         3: ((0, 1, 0), (1, 0, 0), (0, 0, 1)), 4: ((0, -1, 0), (1, 0, 0), (0, 0, 1))}
for count, (n, u, v) in faces.items():
    for a, b in layouts[count]:
        p = [n[i] * 1.12 + (u[i] * a + v[i] * b) * 0.5 for i in range(3)]
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.2, location=p, segments=24, ring_count=12)
        pip = bpy.context.active_object
        boolean = die.modifiers.new("pip", "BOOLEAN")
        boolean.operation, boolean.object = "DIFFERENCE", pip
        bpy.context.view_layer.objects.active = die
        bpy.ops.object.modifier_apply(modifier="pip")
        bpy.data.objects.remove(pip)

me = die.data
with open(out, "w") as fh:
    for v in me.vertices:
        fh.write("v %.6f %.6f %.6f\n" % tuple(v.co))
    for p in me.polygons:
        fh.write("f " + " ".join(str(i + 1) for i in p.vertices) + "\n")
print("dice:", len(me.vertices), "verts", len(me.polygons), "faces")
