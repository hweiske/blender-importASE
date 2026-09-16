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

from .import_cubefiles import is_vasp_density, read_vasp_density
from .node_networks.compat import (compositor_tree, compositor_output,
                                   alpha_over_sockets)
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
    # isovalue shells: a hue sweep from blue at the outermost (weakest)
    # level to red at the innermost, which is what an HSV ramp gives with
    # just these two stops - blue, cyan, green, yellow, red, matplotlib's
    # jet without having to spell every stop out. Positions start at 0.1
    # so the outermost shell is solidly blue rather than half-faded.
    'SHELLS': ('density_jet material',
               [(0.1, (0.0, 0.0, 1.0, 1)),
                (1.0, (1.0, 0.0, 0.0, 1))]),
}

# The shells ramp is an HSV sweep rather than a list of RGB stops, so the
# hue has to walk the long way round the circle - blue through cyan, green
# and yellow to red. Blender calls that 'FAR' (NEAR takes the short arc,
# blue straight through magenta to red).
SHELL_HUE_PATH = 'FAR'


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
        # the add-on's own reader, not ase's VaspChargeDensity directly:
        # that one silently returns no grids at all for perfectly good
        # files (gzipped ones, or a POTCAR-style species line like 'Fe/'),
        # which is exactly what read_vasp_density works around
        density = read_vasp_density(filepath)
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
                         shells=1, shell_max=None, shell_spacing='LOG',
                         per_shell=False, repeat=(1, 1, 1)):
    """Run marching cubes on the +/- isosurfaces of a density file.

    Returns (vertices, faces, colors): cartesian vertex positions, face
    index triples, and one RGBA color per vertex - the alpha channel
    carries how strong that vertex's shell is (0 outermost, 1 innermost),
    which the 'SHELLS' material turns into transparency.

    per_shell returns a list of (vertices, faces, colors, level) instead,
    one per level (both signs of a level together, since they are equally
    strong), weakest first - what the compositor router needs to put each
    shell on its own layer, and what lets a shell be named after the
    isovalue it stands for. Without a color file the
    positive surface is white and the negative one black; with a color
    file its values are sampled at each vertex (nearest voxel) and
    normalized to a black-to-white gradient.

    color_min/color_max: normalize the color-density values between these
    two values (clamped) instead of the sampled min/max - lets colors stay
    comparable between imports. Left equal/unset, the sampled range is
    used.

    repeat: tile the grid first, so the isosurfaces span a supercell.

    sample_interior: instead of the color value directly on the surface,
    use the strongest (largest magnitude) value found anywhere along the
    surface normal through the volume - projects features buried inside
    the isosurface (e.g. LED energies) onto it.
    """
    marching_cubes = _ensure_skimage()
    volume, spacing, origin = read_density_grid(filepath)
    if tuple(repeat) != (1, 1, 1):
        # the grid of a periodic calculation is periodic itself, so tiling
        # it is the supercell's density - and the marching cubes then runs
        # across the interior boundaries instead of closing every copy off
        # at them, which is what tiling the finished mesh would do
        volume = np.tile(volume, tuple(int(n) for n in repeat))

    # one surface per level per sign. With shells=1 this is the familiar
    # pair at +/- iso_value; beyond that the levels nest inwards, and each
    # shell is colored by the level it stands for (see shell_levels).
    levels = shell_levels(iso_value, shells, shell_max, shell_spacing,
                          limit=np.abs(volume).max())
    signs = [sign for sign in (1.0, -1.0)
             if volume.min() < sign * abs(iso_value) < volume.max()]

    surfaces = []  # (index-space verts, faces, normals, color, strength, level)
    for index, level in enumerate(levels):
        strength = index / (len(levels) - 1) if len(levels) > 1 else 1.0
        for sign in signs:
            if not (volume.min() < sign * level < volume.max()):
                continue  # this shell is deeper than the data goes
            verts, faces, normals, _ = marching_cubes(volume, level=sign * level)
            # marching cubes winds its triangles by the gradient, which
            # points the other way for the negative lobe. Make every shell
            # face outwards - measured from the winding itself, the way
            # Blender reads it, not from the gradient normals skimage
            # returns (those point into a positive lobe, i.e. the other
            # way again). Without this the 'see inside' material culls the
            # far wall of half the shells instead of the near one, and
            # that half renders as a solid blob.
            triangles = verts[faces]
            winding = np.cross(triangles[:, 1] - triangles[:, 0],
                               triangles[:, 2] - triangles[:, 0])
            outward = np.einsum('ij,ij->i', winding,
                                triangles.mean(axis=1) - verts.mean(axis=0)).mean()
            if outward < 0:
                faces = faces[:, ::-1]
            if shells > 1:
                # the shell's own strength, so the ramp is read as a color
                # map: 0 the outermost (weakest) level, 1 the innermost.
                # Both signs use the same scale - a lobe's sign is not in
                # the color any more, it is in where the lobe is.
                color = strength
            else:
                color = 1.0 if sign > 0 else 0.0
            surfaces.append((verts, faces, normals, color, strength, index))
    if not surfaces:
        raise ValueError(
            f'isovalue {iso_value} is outside the data range '
            f'[{volume.min():.3g}, {volume.max():.3g}] of {filepath}')

    if per_shell:
        # one triple per level, so each can become its own object and its
        # own render layer; the colors are constant within a shell, which
        # is why they can skip the color-file path below
        grouped = []
        for index in range(len(levels)):
            same = [entry for entry in surfaces if entry[5] == index]
            if not same:
                continue
            shell_offset = 0
            verts_list, faces_list, colors_list = [], [], []
            for verts, faces, _normals, const_color, strength, _ in same:
                verts_list.append(origin + verts @ spacing)
                faces_list.append(faces + shell_offset)
                colors_list.append(np.tile([const_color] * 3 + [strength],
                                           (len(verts), 1)))
                shell_offset += len(verts)
            grouped.append((np.vstack(verts_list), np.vstack(faces_list),
                            np.vstack(colors_list), levels[index]))
        return grouped

    offset = 0
    all_verts, all_faces, all_normals, const_colors, strengths = [], [], [], [], []
    for verts, faces, normals, const_color, strength, _level_index in surfaces:
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
    user edits survive re-imports.

    The shells preset is the whole shader in four nodes - attribute, ramp,
    **Emission**, output. A contour map is a color map, not a lit surface,
    so emission keeps every band the color the map says it is from any
    angle and under any lighting, and nothing else is needed: with every
    shell on its own render pass (see shell_compositor) no shell can
    occlude another, so there is nothing to see through and no
    transparency to fade.
    """
    name, stops = SHADER_PRESETS[preset]
    mat = newMaterial(name)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
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
        if preset == 'SHELLS':
            ramp.color_ramp.color_mode = 'HSV'
            ramp.color_ramp.hue_interpolation = SHELL_HUE_PATH
        elements = ramp.color_ramp.elements
        elements[0].position, elements[0].color = stops[0]
        elements[1].position, elements[1].color = stops[-1]
        for position, color in stops[1:-1]:
            el = elements.new(position)
            el.color = color
    links.new(color_attr.outputs['Color'], ramp.inputs['Fac'])

    if preset != 'SHELLS':
        principled = next((n for n in nodes
                           if n.bl_idname == 'ShaderNodeBsdfPrincipled'), None)
        if principled is None:
            # use_nodes=True does not reliably leave a Principled behind -
            # newShader creates one explicitly for the same reason
            principled = nodes.new('ShaderNodeBsdfPrincipled')
            principled.location = (-60, 200)
            output = next(n for n in nodes
                          if n.bl_idname == 'ShaderNodeOutputMaterial')
            links.new(principled.outputs['BSDF'], output.inputs['Surface'])
        links.new(ramp.outputs['Color'], principled.inputs['Base Color'])
        return mat

    emission = nodes.get('Shell Emission')
    if emission is None:
        for node in list(nodes):
            if node.bl_idname == 'ShaderNodeBsdfPrincipled':
                nodes.remove(node)          # the default shader, unused here
        emission = nodes.new('ShaderNodeEmission')
        emission.name = 'Shell Emission'
        emission.location = (-60, 200)
    links.new(ramp.outputs['Color'], emission.inputs['Color'])
    output = next(n for n in nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
    links.new(emission.outputs['Emission'], output.inputs['Surface'])
    return mat


def _isomesh_object(name, verts, faces, colors, shade_smooth, preset,
                    source=None):
    """One isosurface mesh object, colored by the density_color attribute.

    `source` records what it was built from - the density file and the
    level - so a supercell can rebuild it from a tiled grid later.
    """
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts.tolist(), [], faces.tolist())
    attr = mesh.color_attributes.new(name=COLOR_ATTRIBUTE,
                                     type='FLOAT_COLOR', domain='POINT')
    attr.data.foreach_set('color', np.ascontiguousarray(colors).ravel())
    if shade_smooth:
        mesh.polygons.foreach_set('use_smooth', [True] * len(mesh.polygons))
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    collection = bpy.context.collection or bpy.context.scene.collection
    collection.objects.link(obj)
    obj.data.materials.append(_density_mesh_material(preset))
    for key, value in (source or {}).items():
        obj[key] = value
    return obj


def density_mesh_supercell(mesh_obj, repeat=(1, 1, 1)):
    """Rebuild an isosurface mesh with its grid tiled into a supercell.

    The surface has to be recomputed rather than the mesh repeated: a copy
    of the mesh would still be closed off at the cell face it was
    generated in, while marching cubes over the tiled grid runs straight
    through. The object, its material, its collection and its render layer
    all stay as they are - only the mesh data is replaced - so a shell
    keeps its place in the compositor.
    """
    path = mesh_obj.get('ase_density_file')
    level = mesh_obj.get('ase_shell_level')
    if path is None or level is None:
        raise ValueError(f'{mesh_obj.name} was not imported as an ASE density mesh')
    path = bpy.path.abspath(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f'the density file {path} is gone - re-import it')

    groups = density_to_mesh_data(path, iso_value=level, shells=1,
                                  per_shell=True, repeat=repeat)
    verts, faces, colors, _level = groups[0]
    if 'ase_shell_color' in mesh_obj:
        # one shell of a nest: it is one flat color, the level it stands for
        colors = np.tile([mesh_obj['ase_shell_color']] * 3
                         + [mesh_obj['ase_shell_alpha']], (len(verts), 1))

    old_mesh = mesh_obj.data
    smooth = bool(old_mesh.polygons and old_mesh.polygons[0].use_smooth)
    mesh = bpy.data.meshes.new(old_mesh.name)
    mesh.from_pydata(verts.tolist(), [], faces.tolist())
    attr = mesh.color_attributes.new(name=COLOR_ATTRIBUTE,
                                     type='FLOAT_COLOR', domain='POINT')
    attr.data.foreach_set('color', np.ascontiguousarray(colors).ravel())
    if smooth:
        mesh.polygons.foreach_set('use_smooth', [True] * len(mesh.polygons))
    mesh.update()
    for material in old_mesh.materials:
        mesh.materials.append(material)
    mesh_obj.data = mesh
    bpy.data.meshes.remove(old_mesh)
    mesh_obj['ase_density_repeat'] = [int(n) for n in repeat]
    return mesh_obj


def shell_compositor(shell_objects, scene=None, structure_on_top=True):
    """Put every shell on its own render pass and composite them by value.

    Within one render pass two surfaces are ordered by where they are in
    space. For a nest that *is* the value order, but two separate lobes
    can overlap the other way round - a weak band of one in front of a
    strong band of another - and no shader can fix that, because the
    ordering happens per ray. Compositing can: each shell renders on its
    own view layer and they are alpha-overed weakest first, so a stronger
    value always lands on top of a weaker one no matter where either sits
    in space.

    The structure gets a pass of its own, composited last, so the atoms
    and bonds sit over the contour bands the way such a map is normally
    drawn - and it is excluded from the shell passes so it cannot occlude
    the bands it is drawn over. Whatever else is in the scene keeps the
    original view layer and goes at the bottom, as the backdrop.

    So the stack, bottom to top, is

        rest of the scene -> weakest shell .. strongest shell -> structure

    `shell_objects` must be ordered weakest first. Returns the view layer
    names in that order.
    """
    scene = scene or bpy.context.scene
    master = scene.collection

    # the shells have to be alone in their collections to be separable
    origins = []
    shell_collections = []
    for obj in shell_objects:
        origins.extend(obj.users_collection)
        collection = bpy.data.collections.new(f'{obj.name}_layer')
        master.children.link(collection)
        for previous in list(obj.users_collection):
            previous.objects.unlink(obj)
        collection.objects.link(obj)
        shell_collections.append(collection)
    # where the shells came from is where the structure is - unless they
    # were linked straight into the scene's master collection, which is
    # the view layer's root and cannot be excluded from anything
    structure_collections = [c for c in dict.fromkeys(origins)
                             if c not in shell_collections and c is not master]
    # the shell collections sit at the scene root, not inside the
    # structure's collection, because a nested one cannot be excluded from
    # a view layer on its own. This says which structure they belong to,
    # so a supercell can find them again.
    for collection in shell_collections:
        collection['ase_structure'] = (structure_collections[0].name
                                       if structure_collections else '')

    def isolate(view_layer, keep):
        """Leave only the collections in `keep` of the ones this router
        owns; everything it does not own is left alone."""
        for child in view_layer.layer_collection.children:
            if child.collection in shell_collections or \
                    child.collection in structure_collections:
                child.exclude = child.collection not in keep

    def layer_for(name, keep):
        view_layer = scene.view_layers.get(name) or scene.view_layers.new(name)
        view_layer.use = True
        isolate(view_layer, keep)
        return view_layer.name

    # the scene's first view layer is what the viewport shows and what you
    # work in, so it keeps showing everything - isolating passes there
    # would empty the viewport of the very structure being imported. It is
    # taken out of the render instead, or its contents would render twice.
    base = scene.view_layers[0]
    for child in base.layer_collection.children:
        child.exclude = False
    base.use = False

    order = [layer_for(collection.name, [collection])
             for collection in shell_collections]
    if structure_collections:
        structure_layer = layer_for(f'{structure_collections[0].name}_structure',
                                    structure_collections)
        order.insert(len(order) if structure_on_top else 0, structure_layer)
    # a backdrop pass, but only when there is anything else to render
    owned = shell_collections + structure_collections
    rest = [ob for ob in scene.objects
            if ob.type in {'MESH', 'VOLUME', 'CURVE', 'SURFACE', 'META', 'FONT'}
            and not any(coll in owned for coll in ob.users_collection)]
    if rest:
        order.insert(0, layer_for(f'{scene.name}_backdrop', []))

    # alpha over needs something to composite onto
    scene.render.film_transparent = True
    tree = compositor_tree(scene)
    for node in list(tree.nodes):
        tree.nodes.remove(node)

    stack = None
    for height, layer_name in enumerate(order):
        render_layer = tree.nodes.new('CompositorNodeRLayers')
        render_layer.scene = scene
        render_layer.layer = layer_name
        render_layer.location = (-400, -220 * height)
        if stack is None:
            stack = render_layer.outputs['Image']
            continue
        over = tree.nodes.new('CompositorNodeAlphaOver')
        over.location = (-100, -220 * height + 100)
        background, foreground, factor = alpha_over_sockets(over)
        factor.default_value = 1.0
        tree.links.new(stack, background)
        tree.links.new(render_layer.outputs['Image'], foreground)
        stack = over.outputs['Image']
    output, image_socket = compositor_output(tree)
    output.location = (200, 0)
    tree.links.new(stack, image_socket)
    return order


def import_density_mesh(filepath, filename, color_filepath=None,
                        iso_value=0.03, shade_smooth=True, preset='DEFAULT',
                        import_atoms=True, color_min=None, color_max=None,
                        sample_interior=False, outline=True, shells=1,
                        shell_max=None, shell_spacing='LOG', layered=True,
                        **kwargs):
    """Import a density isosurface as a mesh.

    shells > 1 imports a nest of isosurfaces instead of one, each colored
    by the isovalue it stands for - a density colored by its own value.
    Each shell is flat-shaded by its own color from the map, and the
    ordering - strongest on top - comes from putting every shell on its
    own render pass, which is what layered does and why it defaults to on
    for a nest. The levels run from iso_value inwards to shell_max (see
    shell_levels), and a nest defaults to the 'SHELLS' preset, the only
    one built as a color map.

    layered (on by default for a nest) puts every shell on its own view
    layer and composites them by value (shell_compositor), so a stronger
    value lands on top of a weaker one however the two sit in space. It
    returns the list of shell objects, weakest first, instead of a single
    object - and costs one render pass per shell plus one for the
    structure. Turned off, the shells are opaque surfaces in one pass, so
    the outermost hides the rest.
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

    name = filename.split('.')[0] + '_isomesh'
    if layered and shells > 1:
        groups = density_to_mesh_data(
            filepath, iso_value=iso_value, shells=shells, shell_max=shell_max,
            shell_spacing=shell_spacing, per_shell=True)
        # named after the isovalue each one stands for, so the outliner,
        # the collections and the render layers all say which level they are
        objects = [_isomesh_object(f'{name}_shell_{level:.4g}', verts, faces,
                                   colors, shade_smooth, preset,
                                   source={'ase_density_file': filepath,
                                           'ase_shell_level': float(level),
                                           'ase_shell_color': float(colors[0][0]),
                                           'ase_shell_alpha': float(colors[0][3])})
                   for verts, faces, colors, level in groups]
        print(f'density mesh: {shells} shells on their own render layers, '
              f'{sum(len(o.data.vertices) for o in objects)} verts')
        shell_compositor(objects)
        return objects

    verts, faces, colors = density_to_mesh_data(
        filepath, color_filepath=color_filepath, iso_value=iso_value,
        color_min=color_min, color_max=color_max,
        sample_interior=sample_interior, shells=shells, shell_max=shell_max,
        shell_spacing=shell_spacing)
    print(f'density mesh: {len(verts)} verts, {len(faces)} faces'
          + (f', {shells} shells' if shells > 1 else ''))
    return _isomesh_object(name, verts, faces, colors, shade_smooth, preset,
                           source={'ase_density_file': filepath,
                                   'ase_shell_level': float(iso_value)})
