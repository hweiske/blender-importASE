"""Render collections separately -- one still each, or one whole animation each.

Part of the ASE Importer add-on (was a separate render_vpts.py add-on
that had to be installed alongside it).

Two operators, both built on the same idea: show exactly one collection at a
time and render it, so a scene holding many structures yields one image set per
structure with nothing else in frame.

``render.render_vpts`` -- "Render structure vpts"
    One still per collection per camera, named ``<collection>_<camera>.png``
    directly in the chosen folder.

``render.render_animations`` -- "Render multiple animations"
    The whole animation of every collection, into one subfolder per collection::

        <folder>/<collection>/<camera>_0000.png
        <folder>/<collection>/<camera>_0001.png

    Written for trajectories imported by this add-on: each imported trajectory
    is its own collection, and its length is read from its own shape-key
    ``eval_time`` f-curve rather than from the scene frame range, because the
    importer rewrites ``scene.frame_end`` on every import -- in a file holding
    several trajectories the scene range belongs to whichever one was imported
    last, and it is one frame too long (the last shape key sits at
    ``n_images - 1``).

    It runs as a modal timer operator, one image per tick, so Blender keeps
    redrawing and Esc cancels a long job. For headless use
    (``blender -b --python ...``), where modal operators do not run, call
    :func:`render_animations` directly.
"""
import time
from os.path import join
from pathlib import Path

import bpy
from bpy_extras.io_utils import ExportHelper


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------
def visibility_snapshot(scene):
    """Baseline ``hide_render`` of every top-level collection and its objects.

    Restoring this baseline is not the same as clearing ``hide_render``: the
    importer deliberately leaves its ``*_pair_table`` / ``*_element_table``
    helper meshes hidden, and force-showing everything in a collection drags
    those into the image as stray geometry.
    """
    return {
        coll.name: {
            'collection': coll.hide_render,
            'objects': {ob.name: ob.hide_render for ob in coll.objects},
        }
        for coll in scene.collection.children
    }


def set_collection_visible(coll, base, visible):
    """``visible=True`` restores the baseline visibility, ``False`` hides all."""
    if visible:
        coll.hide_render = base['collection']
        for ob in coll.objects:
            ob.hide_render = base['objects'].get(ob.name, ob.hide_render)
    else:
        coll.hide_render = True
        for ob in coll.objects:
            ob.hide_render = True


def iter_fcurves(action):
    """F-curves of an action, for the slotted (4.4+) and the legacy layout."""
    if getattr(action, 'layers', None):
        for layer in action.layers:
            for strip in layer.strips:
                for bag in getattr(strip, 'channelbags', []):
                    for fcurve in bag.fcurves:
                        yield fcurve
    else:
        for fcurve in action.fcurves:
            yield fcurve


def collection_frame_range(coll):
    """(first, last) animated frame of a collection, or None if it is static.

    Read from the shape-key ``eval_time`` f-curves the trajectory importer
    writes -- see the module docstring for why not ``scene.frame_start/end``.
    """
    lo = hi = None
    for ob in coll.objects:
        if ob.type != 'MESH' or ob.data.shape_keys is None:
            continue
        anim = ob.data.shape_keys.animation_data
        if anim is None or anim.action is None:
            continue
        for fcurve in iter_fcurves(anim.action):
            if not fcurve.keyframe_points:
                continue
            first = fcurve.keyframe_points[0].co[0]
            last = fcurve.keyframe_points[-1].co[0]
            lo = first if lo is None else min(lo, first)
            hi = last if hi is None else max(hi, last)
    if lo is None:
        return None
    return int(round(lo)), int(round(hi))


def scene_cameras(scene, names=''):
    """Cameras to render with; ``names`` is an optional comma-separated filter."""
    cameras = [ob for ob in scene.objects if ob.type == 'CAMERA']
    wanted = [n.strip() for n in names.split(',') if n.strip()]
    if wanted:
        by_name = {c.name: c for c in cameras}
        missing = [n for n in wanted if n not in by_name]
        if missing:
            raise ValueError('no such camera(s): %s; scene has %s'
                             % (', '.join(missing), ', '.join(by_name) or 'none'))
        cameras = [by_name[n] for n in wanted]
    return cameras


def output_state(scene):
    """Save the render output settings this module overwrites, so they can go back.

    ``use_file_extension`` has to be off: these operators write one explicit
    filename per frame, already carrying the extension, and Blender would append
    a second one. ``filepath`` is a user setting and should not be left pointing
    at the last frame of a batch job.
    """
    return (scene.render.filepath, scene.render.use_file_extension)


def restore_output_state(scene, state):
    scene.render.filepath, scene.render.use_file_extension = state


def build_jobs(scene, directory, cameras, stride=1, use_scene_range=False,
               start=None, end=None, collections=None, skip_existing=False):
    """Flat list of (collection, camera, frame, filepath) for the whole run."""
    stride = max(1, int(stride))
    ext = scene.render.file_extension
    root = Path(bpy.path.abspath(directory))
    prefix_camera = len(cameras) > 1
    jobs = []
    for coll in (collections if collections is not None
                 else list(scene.collection.children)):
        rng = collection_frame_range(coll)
        if rng is None or use_scene_range:
            first, last = scene.frame_start, scene.frame_end
        else:
            first, last = rng
        if start is not None:
            first = start
        if end is not None:
            last = end
        folder = root / coll.name
        for camera in cameras:
            tag = '%s_' % camera.name if prefix_camera else ''
            for frame in range(first, last + 1, stride):
                path = folder / ('%s%04d%s' % (tag, frame, ext))
                if skip_existing and path.exists():
                    continue
                jobs.append((coll, camera, frame, path))
    return jobs


def render_animations(scene=None, directory='//', cameras='', stride=1,
                      use_scene_range=False, start=None, end=None,
                      skip_existing=False, collections=None, report=print):
    """Blocking version, for ``blender -b --python`` where modal cannot run.

    Renders every collection's animation into ``<directory>/<collection>/``
    and returns the number of images written. ``collections`` restricts the run
    to a subset; the default is every top-level collection in the scene.
    """
    scene = scene or bpy.context.scene
    cams = scene_cameras(scene, cameras)
    if not cams:
        raise ValueError('the scene has no camera')

    jobs = build_jobs(scene, directory, cams, stride, use_scene_range,
                      start, end, collections=collections,
                      skip_existing=skip_existing)
    base = visibility_snapshot(scene)
    saved_output = output_state(scene)
    scene.render.use_file_extension = False
    for coll in scene.collection.children:
        set_collection_visible(coll, base[coll.name], False)

    shown = None
    written = 0
    started = time.time()
    try:
        for n, (coll, camera, frame, path) in enumerate(jobs, 1):
            if shown is not coll:
                if shown is not None:
                    set_collection_visible(shown, base[shown.name], False)
                set_collection_visible(coll, base[coll.name], True)
                shown = coll
            path.parent.mkdir(parents=True, exist_ok=True)
            scene.camera = camera
            scene.frame_set(frame)
            scene.render.filepath = str(path)
            bpy.ops.render.render(write_still=True)
            written += 1
            if n % 10 == 0 or n == len(jobs):
                rate = n / max(time.time() - started, 1e-9)
                report('[render] %d/%d  %s %s f%d  %.2f img/s  ETA %.1f min'
                       % (n, len(jobs), coll.name, camera.name, frame, rate,
                          (len(jobs) - n) / rate / 60))
    finally:
        for coll in scene.collection.children:
            set_collection_visible(coll, base[coll.name], True)
        restore_output_state(scene, saved_output)
    return written


# --------------------------------------------------------------------------
# operator 1: one still per collection per camera (unchanged behaviour)
# --------------------------------------------------------------------------
class RenderImageOperator(bpy.types.Operator, ExportHelper):
    bl_idname = "render.render_vpts"
    bl_label = "Render structure vpts"
    bl_description = "Render every collection separately for every camera"

    # ExportHelper expects this; the operator builds its own per-image
    # filenames from the collection and camera names
    filename_ext = ".png"

    directory: bpy.props.StringProperty(
        name='folder',
        description='where to put the images',
        subtype='DIR_PATH'
    )
    imagepath: bpy.props.StringProperty(
        name='imagepath',
        description="path for the image",
        default='',
        subtype='FILE_PATH'
    )

    def execute(self, context):
        # Set up scene
        scene = bpy.context.scene
        cameras=self.get_camera_list()
        collections=scene.collection.children
        #print(collections)
        for camera in cameras:
            self.toggle(camera,SET=True)
        for collection in collections:
            self.toggle_collection(collection,SET=True)
        for collection in collections:
            self.toggle_collection(collection,SET=False)
            #print(collection.name)
            for camera in cameras:
                self.toggle(camera,SET=False)
                bpy.context.scene.camera = camera
                self.RENDER(FILEPATH=str(Path(f'{join(self.directory,collection.name)}_{camera.name}.png')))
                #print(camera.name)
                self.toggle(camera,SET=True)
            self.toggle_collection(collection,SET=True)
        for camera in cameras:
            self.toggle(camera,SET=False)
        for collection in collections:
            self.toggle_collection(collection,SET=False)
        return {'FINISHED'}
    #def draw_func(self, context):
    #    layout = self.layout
    #    layout.operator("render.render_image", text="Render Image")
    def toggle(self,obj,SET=True):
        obj.hide_render = SET
        obj.hide_viewport = SET  # Optional: hide in the viewport as well
        for child in obj.children:
                child.hide_render = SET
                child.hide_viewport = SET  # Optional: hide in the viewport as well
        return(None)
    def toggle_collection(self,coll_obj,SET=True):
        coll_obj.hide_render = SET
        #coll_obj.hide_viewport = SET
        #bpy.context.view_layer.active_layer_collection =
        for obj in coll_obj.objects:
            obj.hide_render = SET
        #     obj.hide_viewport = SET  # Optional: hide in the viewport as well
    def RENDER(self,FILEPATH='./'):
        render_settings = bpy.context.scene.render
        render_settings.filepath = FILEPATH  # Replace with desired output file path
#        render_settings.engine = 'CYCLES'  # Replace with desired render engine
        #render_settings.resolution_x = 2000  # Replace with desired resolution
        #render_settings.resolution_y = 2000
        bpy.ops.render.render(write_still=True)
    def get_camera_list(self):
        camera_list = []
        scene = bpy.context.scene
        for obj in scene.objects:
            if obj.type == 'CAMERA':
                camera_list.append(obj)
        return camera_list
    @classmethod
    def poll(cls, context):
        return context.scene is not None
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "imagepath")

    def check(self, context):
        return self.filepath != ""


# --------------------------------------------------------------------------
# operator 2: the whole animation of every collection
# --------------------------------------------------------------------------
class RenderAnimationsOperator(bpy.types.Operator, ExportHelper):
    bl_idname = "render.render_animations"
    bl_label = "Render multiple animations"
    bl_description = ("Render the full animation of every collection separately, "
                      "one subfolder per collection. Esc cancels")

    filename_ext = ".png"

    directory: bpy.props.StringProperty(
        name='folder',
        description='where to put the per-collection subfolders',
        subtype='DIR_PATH'
    )
    cameras: bpy.props.StringProperty(
        name='cameras',
        description="comma-separated camera names; empty means every camera in "
                    "the scene. With more than one camera the file names are "
                    "prefixed <camera>_",
        default='',
    )
    stride: bpy.props.IntProperty(
        name='every nth frame',
        description='render only every nth frame of each animation',
        default=1, min=1,
    )
    use_scene_range: bpy.props.BoolProperty(
        name='use scene frame range',
        description="ignore each collection's own animated range and use the "
                    "scene's frame_start/frame_end for all of them",
        default=False,
    )
    skip_existing: bpy.props.BoolProperty(
        name='skip existing',
        description='leave frames that are already on disk alone, so an '
                    'interrupted run can be restarted where it stopped',
        default=True,
    )

    _timer = None
    _jobs = []
    _index = 0
    _shown = None
    _base = None
    _output = None
    _started = 0.0

    @classmethod
    def poll(cls, context):
        return context.scene is not None

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'cameras')
        layout.prop(self, 'stride')
        layout.prop(self, 'use_scene_range')
        layout.prop(self, 'skip_existing')

    def check(self, context):
        return True

    def execute(self, context):
        scene = context.scene
        try:
            cameras = scene_cameras(scene, self.cameras)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        if not cameras:
            self.report({'ERROR'}, 'the scene has no camera')
            return {'CANCELLED'}
        if not scene.collection.children:
            self.report({'ERROR'}, 'the scene has no collections to render')
            return {'CANCELLED'}

        self._jobs = build_jobs(scene, self.directory, cameras, self.stride,
                                self.use_scene_range,
                                skip_existing=self.skip_existing)
        if not self._jobs:
            self.report({'INFO'}, 'nothing to render (all frames already exist?)')
            return {'FINISHED'}

        self._base = visibility_snapshot(scene)
        self._output = output_state(scene)
        scene.render.use_file_extension = False
        for coll in scene.collection.children:
            set_collection_visible(coll, self._base[coll.name], False)
        self._index = 0
        self._shown = None
        self._started = time.time()

        self.report({'INFO'}, 'rendering %d images into %s - Esc to cancel'
                    % (len(self._jobs), self.directory))
        wm = context.window_manager
        wm.progress_begin(0, len(self._jobs))
        # one image per timer tick: the render call itself blocks, but Blender
        # gets to redraw and handle Esc between frames
        self._timer = wm.event_timer_add(0.01, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        if event.type in {'ESC'}:
            self.report({'WARNING'}, 'cancelled after %d/%d images'
                        % (self._index, len(self._jobs)))
            return self._finish(context, {'CANCELLED'})
        if event.type != 'TIMER':
            return {'PASS_THROUGH'}

        scene = context.scene
        coll, camera, frame, path = self._jobs[self._index]
        if self._shown is not coll:
            if self._shown is not None:
                set_collection_visible(self._shown,
                                       self._base[self._shown.name], False)
            set_collection_visible(coll, self._base[coll.name], True)
            self._shown = coll
        path.parent.mkdir(parents=True, exist_ok=True)
        scene.camera = camera
        scene.frame_set(frame)
        scene.render.filepath = str(path)
        bpy.ops.render.render(write_still=True)

        self._index += 1
        context.window_manager.progress_update(self._index)
        if self._index >= len(self._jobs):
            elapsed = (time.time() - self._started) / 60.0
            self.report({'INFO'}, 'rendered %d images in %.1f min'
                        % (self._index, elapsed))
            return self._finish(context, {'FINISHED'})
        return {'RUNNING_MODAL'}

    def _finish(self, context, status):
        wm = context.window_manager
        if self._timer is not None:
            wm.event_timer_remove(self._timer)
            self._timer = None
        wm.progress_end()
        if self._base is not None:
            for coll in context.scene.collection.children:
                if coll.name in self._base:
                    set_collection_visible(coll, self._base[coll.name], True)
        if self._output is not None:
            restore_output_state(context.scene, self._output)
        self._jobs = []
        self._shown = None
        self._base = None
        self._output = None
        return status


_CLASSES = (RenderImageOperator, RenderAnimationsOperator)


def menu_func(self, context):
    self.layout.operator(RenderImageOperator.bl_idname)
    self.layout.operator(RenderAnimationsOperator.bl_idname)


def register():
    for cls in _CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.TOPBAR_MT_render.append(menu_func)


def unregister():
    bpy.types.TOPBAR_MT_render.remove(menu_func)
    for cls in reversed(_CLASSES):
        bpy.utils.unregister_class(cls)
