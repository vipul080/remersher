bl_info = {
    "name": "Remersher",
    "author": "The Remersher contributors",
    "version": (0, 1, 0),
    "blender": (3, 6, 0),
    "location": "View3D > Sidebar > Remersher",
    "description": "Automatic quad remeshing with the open-source remersher engine",
    "category": "Mesh",
}

import os
import shutil
import subprocess
import tempfile

import bpy
import numpy as np

from . import core

ADDON_DIR = os.path.dirname(os.path.abspath(__file__))


class RemersherPreferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    binary_path: bpy.props.StringProperty(
        name="remersher executable",
        description="Path to the remersher CLI. Leave empty to use a bundled copy or one on PATH",
        subtype="FILE_PATH",
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "binary_path")
        found = core.find_binary(ADDON_DIR, self.binary_path)
        layout.label(text=f"Using: {found}" if found else "remersher executable not found",
                     icon="CHECKMARK" if found else "ERROR")


class RemersherSettings(bpy.types.PropertyGroup):
    target_quad_count: bpy.props.IntProperty(
        name="Quad Count", default=5000, min=10, soft_max=100000,
        description="Approximate number of quads in the result")
    adaptivity: bpy.props.IntProperty(
        name="Adaptive Size", default=50, min=0, max=100, subtype="PERCENTAGE",
        description="How much quads shrink in curved regions and grow in flat ones (0 = uniform)")
    detect_hard_edges: bpy.props.BoolProperty(
        name="Detect Hard Edges by angle", default=True,
        description="Align edge loops to creases sharper than the angle below")
    hard_edge_angle: bpy.props.FloatProperty(
        name="Angle", default=45.0, min=1.0, max=179.0, subtype="NONE",
        description="Dihedral angle (degrees) above which an edge counts as hard")
    preserve_boundary: bpy.props.BoolProperty(
        name="Preserve Borders", default=True,
        description="Keep open borders of the mesh")
    symmetry_x: bpy.props.BoolProperty(name="X", description="Mirror symmetry across the object's X plane")
    symmetry_y: bpy.props.BoolProperty(name="Y", description="Mirror symmetry across the object's Y plane")
    symmetry_z: bpy.props.BoolProperty(name="Z", description="Mirror symmetry across the object's Z plane")
    seed: bpy.props.IntProperty(
        name="Seed", default=0, min=0,
        description="Try another seed for a different layout")
    use_vertex_color: bpy.props.BoolProperty(
        name="Use Vertex Color", default=False,
        description="Paint density with the active color attribute: white = denser, black = sparser, "
                    "mid-grey = unchanged")
    hide_original: bpy.props.BoolProperty(
        name="Hide Original", default=True,
        description="Hide the source object after remeshing")


def _vertex_colors(mesh, loops):
    """Per-vertex sRGB colors of the active color attribute (corner colors are averaged), or None."""
    attr = getattr(mesh.color_attributes, "active_color", None) if hasattr(mesh, "color_attributes") else None
    if attr is None or len(attr.data) == 0:
        return None
    rgba = np.empty(len(attr.data) * 4, dtype=np.float64)
    prop = "color_srgb" if hasattr(attr.data[0], "color_srgb") else "color"
    attr.data.foreach_get(prop, rgba)
    rgba = rgba.reshape(-1, 4)[:, :3]
    if attr.domain == "POINT":
        return rgba
    # CORNER domain: average the corners of each vertex.
    sums = np.zeros((len(mesh.vertices), 3))
    counts = np.zeros(len(mesh.vertices))
    np.add.at(sums, loops, rgba)
    np.add.at(counts, loops, 1)
    return sums / np.maximum(counts, 1)[:, None]


def _mesh_arrays(obj, depsgraph):
    """World-independent (object space) vertices and faces of obj with modifiers applied."""
    evaluated = obj.evaluated_get(depsgraph)
    mesh = evaluated.to_mesh()
    try:
        verts = np.empty(len(mesh.vertices) * 3, dtype=np.float64)
        mesh.vertices.foreach_get("co", verts)
        loops = np.empty(len(mesh.loops), dtype=np.int64)
        mesh.loops.foreach_get("vertex_index", loops)
        starts = np.empty(len(mesh.polygons), dtype=np.int64)
        totals = np.empty(len(mesh.polygons), dtype=np.int64)
        mesh.polygons.foreach_get("loop_start", starts)
        mesh.polygons.foreach_get("loop_total", totals)
        faces = [loops[s:s + t].tolist() for s, t in zip(starts, totals)]
        return verts.reshape(-1, 3), faces, _vertex_colors(mesh, loops)
    finally:
        evaluated.to_mesh_clear()


class REMERSHER_OT_remesh(bpy.types.Operator):
    bl_idname = "remersher.remesh"
    bl_label = "Remesh"
    bl_description = "Quad-remesh the active mesh object into a new object"
    bl_options = {"REGISTER", "UNDO"}

    _proc = None
    _timer = None
    _log = None
    _workdir = None
    _source_name = ""

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == "MESH" and context.mode == "OBJECT"

    def _start(self, context):
        """Exports the mesh and launches the CLI. Returns an error message, or None."""
        prefs = context.preferences.addons[__package__].preferences
        binary = core.find_binary(ADDON_DIR, prefs.binary_path)
        if not binary:
            return "remersher executable not found; set it in the add-on preferences"
        obj = context.active_object
        verts, faces, colors = _mesh_arrays(obj, context.evaluated_depsgraph_get())
        if not faces:
            return "Mesh has no faces"

        s = context.scene.remersher
        settings = core.RemeshSettings(
            target_quad_count=s.target_quad_count, adaptivity=s.adaptivity,
            detect_hard_edges=s.detect_hard_edges, hard_edge_angle=s.hard_edge_angle,
            preserve_boundary=s.preserve_boundary,
            seed=s.seed, use_vertex_color=s.use_vertex_color,
            symmetry=("x" if s.symmetry_x else "") + ("y" if s.symmetry_y else "") + ("z" if s.symmetry_z else ""))

        self._workdir = tempfile.mkdtemp(prefix="remersher_")
        in_path = os.path.join(self._workdir, "in.obj")
        self._out_path = os.path.join(self._workdir, "out.obj")
        self._log_path = os.path.join(self._workdir, "log.txt")
        s_use_colors = context.scene.remersher.use_vertex_color and colors is not None
        core.write_obj(in_path, verts, faces, colors if s_use_colors else None)

        self._log = open(self._log_path, "w")
        self._proc = subprocess.Popen(core.build_command(binary, in_path, self._out_path, settings),
                                      stdout=subprocess.DEVNULL, stderr=self._log)
        self._source_name = obj.name
        return None

    def _complete(self, context):
        """Reads the CLI's result once the process has exited and builds the new object."""
        code = self._proc.returncode
        self._log.close()
        with open(self._log_path) as fh:
            log = fh.read()
        if code != 0 or not os.path.exists(self._out_path):
            self._finish(context)
            self.report({"ERROR"}, "Remesh failed: " + core.summarize_stderr(log))
            return {"CANCELLED"}

        verts, faces = core.read_obj(self._out_path)
        self._finish(context)
        source = bpy.data.objects.get(self._source_name)
        if source is None:
            self.report({"ERROR"}, "Source object was removed during remeshing")
            return {"CANCELLED"}
        result = self._create_object(context, source, verts, faces)
        quads = sum(1 for f in faces if len(f) == 4)
        self.report({"INFO"}, f"Remersher: {len(faces)} faces ({quads} quads) → {result.name}")
        return {"FINISHED"}

    def execute(self, context):
        # Blocking path, used from scripts and in background mode where modal timers do not run.
        error = self._start(context)
        if error:
            self.report({"ERROR"}, error)
            return {"CANCELLED"}
        self._proc.wait()
        return self._complete(context)

    def invoke(self, context, event):
        error = self._start(context)
        if error:
            self.report({"ERROR"}, error)
            return {"CANCELLED"}
        self._timer = context.window_manager.event_timer_add(0.2, window=context.window)
        context.window_manager.modal_handler_add(self)
        context.workspace.status_text_set("Remersher: remeshing… (Esc to cancel)")
        return {"RUNNING_MODAL"}

    def modal(self, context, event):
        if event.type == "ESC":
            self._proc.kill()
            self._proc.wait()
            self._finish(context)
            self.report({"WARNING"}, "Remesh cancelled")
            return {"CANCELLED"}
        if event.type != "TIMER" or self._proc.poll() is None:
            return {"PASS_THROUGH"}
        return self._complete(context)

    def _create_object(self, context, source, verts, faces):
        mesh = bpy.data.meshes.new("Retopo_" + source.name)
        mesh.from_pydata(verts, [], faces)
        mesh.validate()
        mesh.update()
        for mat in source.data.materials:
            mesh.materials.append(mat)
        result = bpy.data.objects.new(mesh.name, mesh)
        result.matrix_world = source.matrix_world.copy()
        for coll in source.users_collection:
            coll.objects.link(result)
        if not source.users_collection:
            context.scene.collection.objects.link(result)

        if context.scene.remersher.hide_original:
            source.hide_set(True)
        for o in context.selected_objects:
            o.select_set(False)
        result.select_set(True)
        context.view_layer.objects.active = result
        return result

    def _finish(self, context):
        if self._timer:
            context.window_manager.event_timer_remove(self._timer)
            self._timer = None
        if context.workspace:
            context.workspace.status_text_set(None)
        if self._log and not self._log.closed:
            self._log.close()
        if self._workdir:
            shutil.rmtree(self._workdir, ignore_errors=True)
            self._workdir = None

    def cancel(self, context):
        if self._proc and self._proc.poll() is None:
            self._proc.kill()
        self._finish(context)


class REMERSHER_PT_panel(bpy.types.Panel):
    bl_label = "Remersher"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Remersher"

    def draw(self, context):
        layout = self.layout
        s = context.scene.remersher
        layout.operator(REMERSHER_OT_remesh.bl_idname, text="Remesh It", icon="MOD_REMESH")
        layout.prop(s, "target_quad_count")
        layout.prop(s, "adaptivity", text="Adaptive Size", slider=True)

        box = layout.box()
        row = box.row()
        row.prop(s, "detect_hard_edges")
        sub = row.row()
        sub.enabled = s.detect_hard_edges
        sub.prop(s, "hard_edge_angle", text="")
        box.prop(s, "preserve_boundary")
        layout.prop(s, "use_vertex_color")

        row = layout.row(align=True)
        row.label(text="Symmetry")
        row.prop(s, "symmetry_x", toggle=True)
        row.prop(s, "symmetry_y", toggle=True)
        row.prop(s, "symmetry_z", toggle=True)

        layout.prop(s, "hide_original")
        layout.prop(s, "seed")


classes = (RemersherPreferences, RemersherSettings, REMERSHER_OT_remesh, REMERSHER_PT_panel)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.remersher = bpy.props.PointerProperty(type=RemersherSettings)


def unregister():
    del bpy.types.Scene.remersher
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
