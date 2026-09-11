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

# shader presets: material name + color-ramp stops (position, rgba).
# The ramp maps the normalized color-density values (0..1).
SHADER_PRESETS = {
    'DEFAULT': (DENSITY_MESH_MATERIAL,
                [(0.0, (0.85, 0.10, 0.10, 1)),
                 (0.5, (0.95, 0.95, 0.95, 1)),
                 (1.0, (0.05, 0.15, 0.85, 1))]),
    # electrostatic potential: pure blue -> white -> red
    'ELSTAT': ('elstat_potential material',
               [(0.0, (0.0, 0.0, 1.0, 1)),
                (0.5, (1.0, 1.0, 1.0, 1)),
                (1.0, (1.0, 0.0, 0.0, 1))]),
    # LED analysis: red/green/blue concentrated at the top of the range
    'LED': ('LED material',
            [(0.8, (1.0, 0.0, 0.0, 1)),
             (0.9, (0.0, 1.0, 0.0, 1)),
             (1.0, (0.0, 0.0, 1.0, 1))]),
    # isovalue shells: cold outside (weak) to hot inside (strong), and for a
    # signed density the negative side mirrored through the middle
    'SHELLS': ('density_shells material',
               [(0.0, (0.10, 0.20, 0.90, 1)),
                (0.25, (0.10, 0.75, 0.90, 1)),
                (0.5, (0.95, 0.95, 0.95, 1)),
                (0.75, (0.98, 0.70, 0.10, 1)),
                (1.0, (0.90, 0.10, 0.05, 1))]),
}

# how transparent the shells get: the outermost (weakest) shell at the low
# end, the innermost at the high end. Without this a nest of closed
# surfaces would show nothing but its outside.
SHELL_ALPHA = (0.10, 0.75)


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


def shell_levels(iso_value, shells, shell_max=None, spacing='LOG', limit=None):
    """The isovalues of a nest of shells, outermost (weakest) first.

    Densities fall off exponentially, so the levels are spaced
    geometrically by default: linear spacing bunches every shell against
    the outer surface. `shell_max` is the innermost level; left out, it is
    half the largest value in the data (`limit`), which is deep enough to
    sit inside the outer surface without collapsing to a point.
    """
    if shells <= 1:
        return [abs(iso_value)]
    low = abs(iso_value)
    high = abs(shell_max) if shell_max else (abs(limit) * 0.5 if limit else low * 10)
    if high <= low:
        high = low * 10
    if spacing == 'LINEAR':
        return list(np.linspace(low, high, shells))
    return list(np.geomspace(low, high, shells))


def density_to_mesh_data(filepath, color_filepath=None, iso_value=0.03,
                         color_min=None, color_max=None, sample_interior=False,
                         shells=1, shell_max=None, shell_spacing='LOG'):
    """Run marching cubes on the +/- isosurfaces of a density file.

    Returns (vertices, faces, colors): cartesian vertex positions, face
    index triples, and one RGBA color per vertex - the alpha channel
    carries how strong that vertex's shell is (0 outermost, 1 innermost),
    which the 'SHELLS' material turns into transparency. Without a color file the
    positive surface is white and the negative one black; with a color
    file its values are sampled at each vertex (nearest voxel) and
    normalized to a black-to-white gradient.

    color_min/color_max: normalize the color-density values between these
    two values (clamped) instead of the sampled min/max - lets colors stay
    comparable between imports. Left equal/unset, the sampled range is
    used.

    sample_interior: instead of the color value directly on the surface,
    use the strongest (largest magnitude) value found anywhere along the
    surface normal through the volume - projects features buried inside
    the isosurface (e.g. LED energies) onto it.
    """
    marching_cubes = _ensure_skimage()
    volume, spacing, origin = read_density_grid(filepath)

    # one surface per level per sign. With shells=1 this is the familiar
    # pair at +/- iso_value; beyond that the levels nest inwards, and each
    # shell is colored by the level it stands for (see shell_levels).
    levels = shell_levels(iso_value, shells, shell_max, shell_spacing,
                          limit=np.abs(volume).max())
    signs = [sign for sign in (1.0, -1.0)
             if volume.min() < sign * abs(iso_value) < volume.max()]
    both_signs = len(signs) == 2

    surfaces = []  # (index-space verts, faces, normals, color, strength)
    for index, level in enumerate(levels):
        strength = index / (len(levels) - 1) if len(levels) > 1 else 1.0
        for sign in signs:
            if not (volume.min() < sign * level < volume.max()):
                continue  # this shell is deeper than the data goes
            verts, faces, normals, _ = marching_cubes(volume, level=sign * level)
            if both_signs:
                # 0.5 is the weakest level, the two signs run out to the
                # ends of the ramp from there
                color = 0.5 + sign * 0.5 * strength
            else:
                color = strength if shells > 1 else (1.0 if sign > 0 else 0.0)
            surfaces.append((verts, faces, normals, color, strength))
    if not surfaces:
        raise ValueError(
            f'isovalue {iso_value} is outside the data range '
            f'[{volume.min():.3g}, {volume.max():.3g}] of {filepath}')

    offset = 0
    all_verts, all_faces, all_normals, const_colors, strengths = [], [], [], [], []
    for verts, faces, normals, const_color, strength in surfaces:
        all_verts.append(verts)
        all_faces.append(faces + offset)
        all_normals.append(normals)
        const_colors.append(np.full(len(verts), const_color))
        strengths.append(np.full(len(verts), strength))
        offset += len(verts)
    verts_index = np.vstack(all_verts)
    faces = np.vstack(all_faces)
    normals_index = np.vstack(all_normals)

    # marching_cubes returns vertices in grid-index space; grid point i
    # sits at origin + i @ spacing
    verts_cart = origin + verts_index @ spacing

    if color_filepath:
        color_volume, _, _ = read_density_grid(color_filepath)
        shape_max = np.array(color_volume.shape) - 1

        def sample_at(points):
            idx = np.round(points).astype(int)
            inside = np.all((idx >= 0) & (idx <= shape_max), axis=1)
            idx = np.clip(idx, 0, shape_max)
            values = color_volume[idx[:, 0], idx[:, 1], idx[:, 2]]
            values[~inside] = 0.0  # never beats an in-volume maximum
            return values

        vals = sample_at(verts_index)
        if sample_interior:
            # strongest value anywhere along the +/- normal ray through the
            # whole volume
            lengths = np.linalg.norm(normals_index, axis=1, keepdims=True)
            directions = normals_index / np.maximum(lengths, 1e-12)
            best_abs = np.abs(vals)
            max_steps = int(np.ceil(np.linalg.norm(color_volume.shape)))
            for step in range(1, max_steps + 1):
                for sign in (1.0, -1.0):
                    probed = sample_at(verts_index + directions * (sign * step))
                    stronger = np.abs(probed) > best_abs
                    vals[stronger] = probed[stronger]
                    best_abs[stronger] = np.abs(probed[stronger])
        if color_min is not None and color_max is not None and color_min != color_max:
            vals = np.clip((vals - color_min) / (color_max - color_min), 0.0, 1.0)
        else:
            vals = (vals - vals.min()) / (vals.max() - vals.min() + 1e-12)
    else:
        vals = np.concatenate(const_colors)
    colors = np.repeat(vals[:, None], 3, axis=1)
    # alpha carries the shell strength; a single surface stays opaque, so
    # every path that is not a shell nest renders exactly as before
    alpha = np.concatenate(strengths)[:, None] if shells > 1 else np.ones((len(colors), 1))
    return verts_cart, faces, np.concatenate([colors, alpha], axis=1)


def _density_mesh_material(preset='DEFAULT'):
    """Material mapping the density_color attribute through a color ramp.
    One material per preset; the ramp is only initialized on creation, so
    user edits survive re-imports."""
    name, stops = SHADER_PRESETS[preset]
    mat = newMaterial(name)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    principled = nodes.get('Principled BSDF')
    color_attr = nodes.get('Color Attribute')
    if color_attr is None:
        color_attr = nodes.new('ShaderNodeVertexColor')
        color_attr.name = 'Color Attribute'
        color_attr.location = (-500, 200)
    color_attr.layer_name = COLOR_ATTRIBUTE
    ramp = nodes.get('Color Ramp')
    if ramp is None:
        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.name = 'Color Ramp'
        ramp.location = (-300, 200)
        elements = ramp.color_ramp.elements
        elements[0].position, elements[0].color = stops[0]
        elements[1].position, elements[1].color = stops[-1]
        for position, color in stops[1:-1]:
            el = elements.new(position)
            el.color = color
    links.new(color_attr.outputs['Color'], ramp.inputs['Fac'])
    links.new(ramp.outputs['Color'], principled.inputs['Base Color'])
    if preset == 'SHELLS' and not principled.inputs['Alpha'].is_linked:
        # the attribute's alpha channel holds how strong each shell is;
        # mapped to transparency so the outer ones let you see the inner
        alpha_range = nodes.new('ShaderNodeMapRange')
        alpha_range.name = 'Shell Alpha'
        alpha_range.location = (-300, -100)
        alpha_range.inputs['To Min'].default_value = SHELL_ALPHA[0]
        alpha_range.inputs['To Max'].default_value = SHELL_ALPHA[1]
        links.new(color_attr.outputs['Alpha'], alpha_range.inputs['Value'])
        links.new(alpha_range.outputs['Result'], principled.inputs['Alpha'])
        # EEVEE needs to be told to blend; the property moved in 4.2
        if hasattr(mat, 'surface_render_method'):
            mat.surface_render_method = 'BLENDED'
        elif hasattr(mat, 'blend_method'):
            mat.blend_method = 'BLEND'
    return mat


def import_density_mesh(filepath, filename, color_filepath=None,
                        iso_value=0.03, shade_smooth=True, preset='DEFAULT',
                        import_atoms=True, color_min=None, color_max=None,
                        sample_interior=False, outline=True, shells=1,
                        shell_max=None, shell_spacing='LOG', **kwargs):
    """Import a density isosurface as a mesh.

    shells > 1 imports a nest of isosurfaces instead of one, each colored
    by the isovalue it stands for and made more transparent the further
    out it is - a density colored by its own value. The levels run from
    iso_value inwards to shell_max (see shell_levels), and a nest defaults
    to the 'SHELLS' preset since the other ramps are opaque.
    """
    if shells > 1 and preset == 'DEFAULT':
        preset = 'SHELLS'
    if import_atoms:
        # the structure from the same file, as the nodes representation;
        # this also creates the collection the isomesh is linked into
        from .ui import import_ase_molecule
        import_ase_molecule(filepath, filename, representation='nodes',
                            read_density=False, animate=False, outline=outline,
                            add_supercell=False)

    verts, faces, colors = density_to_mesh_data(
        filepath, color_filepath=color_filepath, iso_value=iso_value,
        color_min=color_min, color_max=color_max,
        sample_interior=sample_interior, shells=shells, shell_max=shell_max,
        shell_spacing=shell_spacing)
    print(f'density mesh: {len(verts)} verts, {len(faces)} faces'
          + (f', {shells} shells' if shells > 1 else ''))

    name = filename.split('.')[0] + '_isomesh'
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts.tolist(), [], faces.tolist())

    attr = mesh.color_attributes.new(name=COLOR_ATTRIBUTE,
                                     type='FLOAT_COLOR', domain='POINT')
    attr.data.foreach_set('color', np.ascontiguousarray(colors).ravel())

    if shade_smooth:
        mesh.polygons.foreach_set('use_smooth', [True] * len(mesh.polygons))
    mesh.update()

    obj = bpy.data.objects.new(name, mesh)
    # context.collection can be None (e.g. headless after scene cleanup)
    collection = bpy.context.collection or bpy.context.scene.collection
    collection.objects.link(obj)
    obj.data.materials.append(_density_mesh_material(preset))
    return obj
