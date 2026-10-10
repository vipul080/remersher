"""Builds extra, more production-like test meshes in Blender and exports them as OBJ:
suzanne (organic, subdivided) and mech_part (bevelled disc with through-holes and a slot).

usage: blender -b --factory-startup --python scripts/make_test_models.py -- <out_dir>
"""
import math
import os
import sys

import bpy

out_dir = sys.argv[sys.argv.index("--") + 1]


def export(obj, name):
    bpy.context.view_layer.objects.active = obj
    me = obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).to_mesh()
    with open(os.path.join(out_dir, name + ".obj"), "w") as fh:
        for v in me.vertices:
            fh.write("v %.6f %.6f %.6f\n" % tuple(v.co))
        for p in me.polygons:
            fh.write("f " + " ".join(str(i + 1) for i in p.vertices) + "\n")
    print(name, len(me.vertices), "verts", len(me.polygons), "faces")


def clear():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()


def cut(target, cutter):
    mod = target.modifiers.new("cut", "BOOLEAN")
    mod.operation, mod.object = "DIFFERENCE", cutter
    bpy.context.view_layer.objects.active = target
    bpy.ops.object.modifier_apply(modifier="cut")
    bpy.data.objects.remove(cutter)


clear()
bpy.ops.mesh.primitive_monkey_add()
monkey = bpy.context.active_object
monkey.modifiers.new("subd", "SUBSURF").levels = 2
export(monkey, "suzanne")

clear()
bpy.ops.mesh.primitive_cylinder_add(vertices=64, radius=1.0, depth=0.5)
disc = bpy.context.active_object
bpy.ops.mesh.primitive_cylinder_add(vertices=48, radius=0.35, depth=1.0)
cut(disc, bpy.context.active_object)
for k in range(4):
    a = k * math.pi / 2 + math.pi / 4
    bpy.ops.mesh.primitive_cylinder_add(vertices=24, radius=0.1, depth=1.0, location=(0.7 * math.cos(a), 0.7 * math.sin(a), 0))
    cut(disc, bpy.context.active_object)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.75, 0.25))
slot = bpy.context.active_object
slot.scale = (0.6, 0.12, 0.2)
cut(disc, slot)
bev = disc.modifiers.new("bevel", "BEVEL")
bev.width, bev.segments, bev.limit_method = 0.03, 3, "ANGLE"
bpy.context.view_layer.objects.active = disc
bpy.ops.object.modifier_apply(modifier="bevel")
export(disc, "mech_part")
