"""Render the README trajectory-panel GIF from Cu-block-TMA.traj.

    blender -b --factory-startup --python docs/render_trajectory_gif.py

Recipe: nodes representation, Cu-Cu bond hidden, Cu shown as vdW
spheres, 1x3x1 supercell, bond distance 0.67, orthographic view
straight down -X centered on the supercell. Rendered from
blender_startup.blend. Frames are assembled into docs/images/
trajectory.gif when imageio is available (otherwise left as PNGs).
"""
import math
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
import bpy
import numpy as np
from mathutils import Euler, Vector

import blender_importASE
try:
    blender_importASE.register()
except Exception:
    pass
from blender_importASE.ui import import_ase_molecule
from blender_importASE import controls
from blender_importASE.controls import pair_id
from blender_importASE.node_networks.compat import set_mod_input, get_mod_input

bpy.ops.wm.open_mainfile(filepath=f'{REPO}/blender_startup.blend')
import_ase_molecule(f'{REPO}/test/fixtures/Cu-block-TMA.traj', 'Cu-block-TMA.traj',
                    representation='nodes', animate=True, read_density=False,
                    outline=False, add_supercell=True, imageslice=6)
struct = next(o for o in bpy.data.objects if o.type == 'MESH' and 'element' in o.data.attributes)
bpy.context.view_layer.objects.active = struct

set_mod_input(struct.modifiers['GeometryNodes.001'], 'Socket_3', 3)     # 1x3x1 supercell (repeat_y)
set_mod_input(struct.modifiers['GeometryNodes.002'], 'Socket_2', 0.67)  # bond distance
mod, idents = controls.find_ase_modifier(struct)
pair_table = get_mod_input(mod, idents['pair_table'])
element_table = get_mod_input(mod, idents['element_table'])
pair_table.data.attributes['cut'].data[pair_id(29, 29)].value = True   # hide Cu-Cu
element_table.data.attributes['radius_mode'].data[29].value = 1        # Cu -> vdW
pair_table.data.update()
element_table.data.update()
struct.update_tag()

frames = bpy.data.scenes['Scene'].frame_end
bpy.context.scene.frame_set(1)
deps = bpy.context.evaluated_depsgraph_get()
ev = struct.evaluated_get(deps)
me = ev.to_mesh()
co = np.array([v.co[:] for v in me.vertices])
co = co @ np.array(struct.matrix_world.to_3x3()).T + np.array(struct.matrix_world.translation)
lo, hi = co.min(0), co.max(0)
center = (lo + hi) / 2
ev.to_mesh_clear()

cam = bpy.data.objects['Camera']
cam.data.type = 'ORTHO'
cam.data.ortho_scale = 22.0
cam.data.clip_start = 0.01
cam.data.clip_end = 100000
euler = Euler((math.radians(90), 0.0, math.radians(45)), 'XYZ')
cam.rotation_euler = euler
view_dir = (euler.to_matrix() @ Vector((0, 0, -1))).normalized()
cam.location = Vector(center) - view_dir * (float(np.linalg.norm(hi - lo)) + 50)
bpy.context.scene.camera = cam

sc = bpy.context.scene
sc.cycles.samples = 48
sc.cycles.use_denoising = True
sc.render.resolution_x = 480
sc.render.resolution_y = 480
sc.render.image_settings.file_format = 'PNG'
frame_dir = f'{REPO}/docs/_traj_frames'
os.makedirs(frame_dir, exist_ok=True)
for fr in range(1, frames + 1):
    sc.frame_set(fr)
    sc.render.filepath = f'{frame_dir}/{fr:03d}.png'
    bpy.ops.render.render(write_still=True)

try:
    import imageio.v2 as imageio
    imgs = [imageio.imread(f'{frame_dir}/{fr:03d}.png')[:, :, :3] for fr in range(1, frames + 1)]
    imageio.mimsave(f'{REPO}/docs/images/trajectory.gif', imgs, duration=0.09, loop=0)
    print('wrote docs/images/trajectory.gif')
except ImportError:
    print(f'frames in {frame_dir}; assemble with any GIF tool')
