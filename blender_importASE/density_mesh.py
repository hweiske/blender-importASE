"""Import a density isosurface as a real mesh via marching cubes.

Unlike the volume-based density import (import_cubefiles), this computes
the +/- isosurfaces in Python with skimage's marching cubes and creates
an ordinary mesh object. A second, optional density file can be used to
color the surface: its values are sampled at every vertex and stored in
a 'density_color' color attribute (used by the generated material).

Ported from import_cubefiles_marchingcube.ipynb.
"""
import os
import bpy
import numpy as np
from ase.io.cube import read_cube
from ase.calculators.vasp import VaspChargeDensity

from .import_cubefiles import is_vasp_density
from .node_networks.electron_density_nodes import newMaterial

DENSITY_MESH_MATERIAL = 'density_mesh material'
COLOR_ATTRIBUTE = 'density_color'


def _ensure_skimage():
    """Import skimage, installing scikit-image on demand like
    check_dependency() does for ase."""
    try:
        from skimage.measure import marching_cubes
        return marching_cubes
    except ImportError:
        pass
    import sys
    import subprocess
    import importlib
    print("scikit-image not present in Blender python. Attempting install...")
    install_path = os.path.join(bpy.utils.script_path_user(), "modules")
    subprocess.check_call([sys.executable, "-m", "pip", "install",
                           "--target", install_path, "scikit-image"])
    if install_path not in sys.path:
        sys.path.append(install_path)
    importlib.invalidate_caches()
    from skimage.measure import marching_cubes
    return marching_cubes


def read_density_grid(filepath):
    """Read a volumetric file (.cube or VASP CHGCAR-like) and return
    (volume, spacing, origin) with spacing as a 3x3 matrix of grid step
    vectors."""
    if is_vasp_density(os.path.basename(filepath)):
        density = VaspChargeDensity(filepath)
        volume = density.chg[-1]
        cell = np.array(density.atoms[-1].get_cell())
        spacing = cell / np.array(volume.shape)[:, None]
        origin = np.zeros(3)
    else:
        with open(filepath, 'r') as f:
            data = read_cube(f, read_data=True)
        volume = data['data']
        spacing = np.array(data['spacing'])
        origin = np.array(data['origin'])
    return volume, spacing, origin


def density_to_mesh_data(filepath, color_filepath=None, iso_value=0.03):
    """Run marching cubes on the +/- isosurfaces of a density file.

    Returns (vertices, faces, colors): cartesian vertex positions, face
    index triples, and one RGB color per vertex. Without a color file the
    positive surface is white and the negative one black; with a color
    file its values are sampled at each vertex (nearest voxel) and
    normalized to a black-to-white gradient.
    """
    marching_cubes = _ensure_skimage()
    volume, spacing, origin = read_density_grid(filepath)

    surfaces = []  # (index-space verts, faces, constant color)
    for level, const_color in ((abs(iso_value), 1.0), (-abs(iso_value), 0.0)):
        if not (volume.min() < level < volume.max()):
            continue  # e.g. total densities have no negative lobe
        verts, faces, _, _ = marching_cubes(volume, level=level)
        surfaces.append((verts, faces, const_color))
    if not surfaces:
        raise ValueError(
            f'isovalue {iso_value} is outside the data range '
            f'[{volume.min():.3g}, {volume.max():.3g}] of {filepath}')

    offset = 0
    all_verts, all_faces, const_colors = [], [], []
    for verts, faces, const_color in surfaces:
        all_verts.append(verts)
        all_faces.append(faces + offset)
        const_colors.append(np.full(len(verts), const_color))
        offset += len(verts)
    verts_index = np.vstack(all_verts)
    faces = np.vstack(all_faces)

    # marching_cubes returns vertices in grid-index space; grid point i
    # sits at origin + i @ spacing
    verts_cart = origin + verts_index @ spacing

    if color_filepath:
        color_volume, _, _ = read_density_grid(color_filepath)
        idx = np.round(verts_index).astype(int)
        idx = np.clip(idx, 0, np.array(color_volume.shape) - 1)
        vals = color_volume[idx[:, 0], idx[:, 1], idx[:, 2]]
        vals = (vals - vals.min()) / (vals.max() - vals.min() + 1e-12)
    else:
        vals = np.concatenate(const_colors)
    colors = np.repeat(vals[:, None], 3, axis=1)
    return verts_cart, faces, colors


def _density_mesh_material():
    mat = newMaterial(DENSITY_MESH_MATERIAL)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    principled = nodes.get('Principled BSDF')
    color_attr = nodes.get('Color Attribute')
    if color_attr is None:
        color_attr = nodes.new('ShaderNodeVertexColor')
        color_attr.name = 'Color Attribute'
        color_attr.location = (-300, 200)
    color_attr.layer_name = COLOR_ATTRIBUTE
    links.new(color_attr.outputs['Color'], principled.inputs['Base Color'])
    return mat


def import_density_mesh(filepath, filename, color_filepath=None,
                        iso_value=0.03, shade_smooth=True, **kwargs):
    verts, faces, colors = density_to_mesh_data(
        filepath, color_filepath=color_filepath, iso_value=iso_value)
    print(f'density mesh: {len(verts)} verts, {len(faces)} faces')

    name = filename.split('.')[0] + '_isomesh'
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts.tolist(), [], faces.tolist())

    attr = mesh.color_attributes.new(name=COLOR_ATTRIBUTE,
                                     type='FLOAT_COLOR', domain='POINT')
    rgba = np.concatenate([colors, np.ones((len(colors), 1))], axis=1)
    attr.data.foreach_set('color', rgba.ravel())

    if shade_smooth:
        mesh.polygons.foreach_set('use_smooth', [True] * len(mesh.polygons))
    mesh.update()

    obj = bpy.data.objects.new(name, mesh)
    # context.collection can be None (e.g. headless after scene cleanup)
    collection = bpy.context.collection or bpy.context.scene.collection
    collection.objects.link(obj)
    obj.data.materials.append(_density_mesh_material())
    return obj
