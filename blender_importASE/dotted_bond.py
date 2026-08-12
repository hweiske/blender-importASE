"""Custom bonds between two atoms of a structure.

Three styles are available for interactions the distance-based bond
search does not draw (a partial bond in a transition state, a hydrogen
bond, ...):

- dotted: a row of spheres along the line between the two atoms
- scaled: a solid bond whose thickness shrinks as the bond gets longer,
  capped at the chosen radius
- dashed: alternating cylinder segments along the line

Each style is its own geometry node group, all sharing the same
front-end: the two atom positions are sampled straight from the
structure object, so the bond follows the structure (including
trajectory animation) instead of being baked in. They are colored like a
real bond - the group writes the same 'COLOR_CURVE' attribute the bond
material reads, blended between the two atoms' colors - and set the
material index of the structure's bond material slot.

Dashes are built by resampling the line, dropping every other segment
and running the result through Curve to Mesh with a circle profile: the
segments come out oriented along the bond for free, with no rotation
maths, and the color gradient survives the conversion.

'Replace' additionally hides the solid bond between the two atoms. The
structure mesh carries a per-atom 'dotted_partner' int attribute holding
*the other atom's index + 1* (0, and therefore a missing attribute,
means "none"), which atoms_and_bonds ORs into its bond delete selection.
The +1 offset is what makes an absent attribute safe: it reads as 0,
which never matches a real index.

Because that is one partner per atom, an atom can only have one of its
solid bonds replaced at a time: replacing 0-1 and then 0-2 leaves the
0-1 bond visible again. Adding the custom bond itself is unaffected -
only the hiding of the underlying solid bond is one-per-atom.
"""
import contextlib

import bpy

from . import __version__
from .node_networks.compat import cin, set_mod_input
from .node_networks.outline import outline_objects

DOTTED_PARTNER_ATTRIBUTE = 'dotted_partner'
OUTLINE_MATERIAL = 'outline_color'

# style -> node group name
BOND_GROUPS = {
    'DOTTED': 'dotted_bond',
    'SCALED': 'scaled_bond',
    'DASHED': 'dashed_bond',
}
DOTTED_BOND_GROUP = BOND_GROUPS['DOTTED']  # backwards compatible alias

BOND_STYLE_ITEMS = (
    ('DOTTED', 'Dotted', 'A row of spheres between the two atoms'),
    ('SCALED', 'Scaled',
     'A solid bond that gets thinner the longer it is, up to the chosen radius'),
    ('DASHED', 'Dashed', 'Alternating cylinder segments between the two atoms'),
)


def _bond_front_end(group):
    """Nodes every style shares: sample both atom positions and colors from
    the structure and build the line between them.

    Returns (group_input, group_output, curve_line, color_a, color_b).
    """
    nodes, links = group.nodes, group.links

    group_input = nodes.new('NodeGroupInput')
    group_input.location = (-900, 0)
    group_output = nodes.new('NodeGroupOutput')
    group_output.location = (900, 0)

    object_info = nodes.new('GeometryNodeObjectInfo')
    object_info.name = 'Structure Info'
    object_info.location = (-700, -200)
    object_info.inputs['As Instance'].default_value = False
    links.new(group_input.outputs['structure'], object_info.inputs['Object'])

    position = nodes.new('GeometryNodeInputPosition')
    position.location = (-700, -400)

    sample_a = nodes.new('GeometryNodeSampleIndex')
    sample_a.name = 'Sample Position A'
    sample_a.data_type = 'FLOAT_VECTOR'
    sample_a.domain = 'POINT'
    sample_a.location = (-500, -100)
    sample_b = nodes.new('GeometryNodeSampleIndex')
    sample_b.name = 'Sample Position B'
    sample_b.data_type = 'FLOAT_VECTOR'
    sample_b.domain = 'POINT'
    sample_b.location = (-500, -300)
    for node, socket in ((sample_a, 'atom A'), (sample_b, 'atom B')):
        links.new(object_info.outputs['Geometry'], node.inputs['Geometry'])
        links.new(position.outputs['Position'], node.inputs['Value'])
        links.new(group_input.outputs[socket], node.inputs['Index'])

    curve_line = nodes.new('GeometryNodeCurvePrimitiveLine')
    curve_line.name = 'Bond Line'
    curve_line.mode = 'POINTS'
    curve_line.location = (-300, -100)
    links.new(sample_a.outputs['Value'], curve_line.inputs['Start'])
    links.new(sample_b.outputs['Value'], curve_line.inputs['End'])

    atom_color = nodes.new('GeometryNodeInputNamedAttribute')
    atom_color.name = 'atom_color'
    atom_color.data_type = 'FLOAT_COLOR'
    atom_color.inputs[0].default_value = 'atom_color'
    atom_color.location = (-700, -600)

    color_a = nodes.new('GeometryNodeSampleIndex')
    color_a.name = 'Sample Color A'
    color_a.data_type = 'FLOAT_COLOR'
    color_a.domain = 'POINT'
    color_a.location = (-500, -520)
    color_b = nodes.new('GeometryNodeSampleIndex')
    color_b.name = 'Sample Color B'
    color_b.data_type = 'FLOAT_COLOR'
    color_b.domain = 'POINT'
    color_b.location = (-500, -700)
    for node, socket in ((color_a, 'atom A'), (color_b, 'atom B')):
        links.new(object_info.outputs['Geometry'], node.inputs['Geometry'])
        links.new(atom_color.outputs[0], node.inputs['Value'])
        links.new(group_input.outputs[socket], node.inputs['Index'])

    return group_input, group_output, curve_line, color_a, color_b


def _store_bond_color(group, curve_socket, color_a, color_b, location):
    """Blend the two atom colors along the curve into 'COLOR_CURVE', the
    attribute the bond material reads. Returns the store node."""
    nodes, links = group.nodes, group.links

    spline_parameter = nodes.new('GeometryNodeSplineParameter')
    spline_parameter.location = (location[0] - 200, location[1] - 460)

    mix_color = nodes.new('ShaderNodeMix')
    mix_color.name = 'Bond Color'
    mix_color.data_type = 'RGBA'
    mix_color.blend_type = 'MIX'
    mix_color.location = (location[0] - 40, location[1] - 460)
    links.new(spline_parameter.outputs['Factor'], mix_color.inputs['Factor'])
    links.new(color_a.outputs['Value'], mix_color.inputs[6])
    links.new(color_b.outputs['Value'], mix_color.inputs[7])

    store_color = nodes.new('GeometryNodeStoreNamedAttribute')
    store_color.name = 'Store COLOR_CURVE'
    store_color.data_type = 'FLOAT_COLOR'
    store_color.domain = 'POINT'
    store_color.inputs['Name'].default_value = 'COLOR_CURVE'
    store_color.location = location
    links.new(curve_socket, store_color.inputs['Geometry'])
    links.new(mix_color.outputs[2], store_color.inputs['Value'])
    return store_color


def _bond_tail(group, group_input, group_output, geometry_socket):
    """Join the host mesh so the object's material slot list travels with
    the geometry (geometry built from scratch inside a tree has an empty
    material list, and Set Material Index would then be out of range and
    clamp to slot 0 as soon as something - e.g. the outline modifier -
    joins it with other geometry), then set the bond material slot."""
    nodes, links = group.nodes, group.links

    join = nodes.new('GeometryNodeJoinGeometry')
    join.name = 'Join Host Mesh'
    join.location = (700, -100)
    links.new(geometry_socket, join.inputs['Geometry'])
    links.new(group_input.outputs['Geometry'], join.inputs['Geometry'])

    set_material_index = nodes.new('GeometryNodeSetMaterialIndex')
    set_material_index.name = 'Bond Material Slot'
    set_material_index.location = (800, -100)
    links.new(join.outputs['Geometry'], set_material_index.inputs['Geometry'])
    links.new(group_input.outputs['material_index'],
              set_material_index.inputs['Material Index'])
    links.new(set_material_index.outputs['Geometry'], group_output.inputs['Geometry'])


def _new_group(style):
    """Create the node group for one bond style, replacing any group of that
    name built by an older add-on version (or by hand - these started life as
    manual node setups, so a name clash is likely)."""
    name = BOND_GROUPS[style]
    existing = bpy.data.node_groups.get(name)
    if existing is not None:
        if existing.description == __version__:
            return existing, False
        existing.name = f'{name}_old'

    group = bpy.data.node_groups.new(type='GeometryNodeTree', name=name)
    group.description = __version__
    group.is_modifier = True

    group.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
    group.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
    group.interface.new_socket('structure', in_out='INPUT', socket_type='NodeSocketObject')
    group.interface.new_socket('atom A', in_out='INPUT', socket_type='NodeSocketInt')
    group.interface.new_socket('atom B', in_out='INPUT', socket_type='NodeSocketInt')
    return group, True


def _finish_group(group):
    mat_index = group.interface.new_socket('material_index', in_out='INPUT',
                                           socket_type='NodeSocketInt')
    mat_index.default_value = 0
    mat_index.min_value = 0


def dotted_bond_node_group():
    """Row of spheres along the bond."""
    group, is_new = _new_group('DOTTED')
    if not is_new:
        return group

    count = group.interface.new_socket('dots', in_out='INPUT', socket_type='NodeSocketInt')
    count.default_value = 10
    count.min_value = 2
    radius = group.interface.new_socket('dot radius', in_out='INPUT',
                                        socket_type='NodeSocketFloat')
    radius.default_value = 0.08
    radius.min_value = 0.0
    resolution = group.interface.new_socket('resolution', in_out='INPUT',
                                            socket_type='NodeSocketInt')
    resolution.default_value = 2
    resolution.min_value = 1
    _finish_group(group)

    nodes, links = group.nodes, group.links
    group_input, group_output, curve_line, color_a, color_b = _bond_front_end(group)

    resample = nodes.new('GeometryNodeResampleCurve')
    resample.name = 'Dot Positions'
    resample.location = (-120, -100)
    links.new(curve_line.outputs['Curve'], resample.inputs['Curve'])
    links.new(group_input.outputs['dots'], resample.inputs['Count'])

    store_color = _store_bond_color(group, resample.outputs['Curve'],
                                    color_a, color_b, (80, -100))

    ico = nodes.new('GeometryNodeMeshIcoSphere')
    ico.name = 'Dot'
    ico.location = (80, -350)
    links.new(group_input.outputs['dot radius'], ico.inputs['Radius'])
    links.new(group_input.outputs['resolution'], ico.inputs['Subdivisions'])

    instance = nodes.new('GeometryNodeInstanceOnPoints')
    instance.name = 'Dots'
    instance.location = (400, -100)
    links.new(store_color.outputs['Geometry'], instance.inputs['Points'])
    links.new(ico.outputs['Mesh'], instance.inputs['Instance'])

    realize = nodes.new('GeometryNodeRealizeInstances')
    realize.location = (560, -100)
    links.new(instance.outputs['Instances'], realize.inputs['Geometry'])

    _bond_tail(group, group_input, group_output, realize.outputs['Geometry'])
    return group


def scaled_bond_node_group():
    """Solid bond whose radius shrinks with increasing bond length."""
    group, is_new = _new_group('SCALED')
    if not is_new:
        return group

    radius = group.interface.new_socket('bond radius', in_out='INPUT',
                                        socket_type='NodeSocketFloat')
    radius.default_value = 0.1
    radius.min_value = 0.0
    reference = group.interface.new_socket('reference length', in_out='INPUT',
                                           socket_type='NodeSocketFloat')
    reference.default_value = 1.5
    reference.min_value = 0.0001
    resolution = group.interface.new_socket('resolution', in_out='INPUT',
                                            socket_type='NodeSocketInt')
    resolution.default_value = 16
    resolution.min_value = 3
    _finish_group(group)

    nodes, links = group.nodes, group.links
    group_input, group_output, curve_line, color_a, color_b = _bond_front_end(group)

    # radius = min(bond radius, bond radius * reference length / bond length):
    # at or below the reference length the bond has its full radius, beyond it
    # it thins out in proportion
    endpoints = nodes.new('ShaderNodeVectorMath')
    endpoints.name = 'Bond Length'
    endpoints.operation = 'DISTANCE'
    endpoints.location = (-300, -700)
    links.new(curve_line.inputs['Start'].links[0].from_socket, endpoints.inputs[0])
    links.new(curve_line.inputs['End'].links[0].from_socket, endpoints.inputs[1])

    scale_by_reference = nodes.new('ShaderNodeMath')
    scale_by_reference.name = 'Radius x reference'
    scale_by_reference.operation = 'MULTIPLY'
    scale_by_reference.location = (-120, -700)
    links.new(group_input.outputs['bond radius'], scale_by_reference.inputs[0])
    links.new(group_input.outputs['reference length'], scale_by_reference.inputs[1])

    divide_by_length = nodes.new('ShaderNodeMath')
    divide_by_length.name = 'Scaled radius'
    divide_by_length.operation = 'DIVIDE'
    divide_by_length.location = (40, -700)
    links.new(scale_by_reference.outputs[0], divide_by_length.inputs[0])
    links.new(endpoints.outputs['Value'], divide_by_length.inputs[1])

    cap_radius = nodes.new('ShaderNodeMath')
    cap_radius.name = 'Cap at bond radius'
    cap_radius.operation = 'MINIMUM'
    cap_radius.location = (200, -700)
    links.new(divide_by_length.outputs[0], cap_radius.inputs[0])
    links.new(group_input.outputs['bond radius'], cap_radius.inputs[1])

    store_color = _store_bond_color(group, curve_line.outputs['Curve'],
                                    color_a, color_b, (80, -100))

    profile = nodes.new('GeometryNodeCurvePrimitiveCircle')
    profile.name = 'Bond Profile'
    profile.location = (400, -400)
    links.new(group_input.outputs['resolution'], profile.inputs['Resolution'])
    links.new(cap_radius.outputs[0], profile.inputs['Radius'])

    to_mesh = nodes.new('GeometryNodeCurveToMesh')
    to_mesh.name = 'Bond Mesh'
    to_mesh.location = (560, -100)
    to_mesh.inputs['Fill Caps'].default_value = True
    links.new(store_color.outputs['Geometry'], to_mesh.inputs['Curve'])
    links.new(profile.outputs['Curve'], to_mesh.inputs['Profile Curve'])

    _bond_tail(group, group_input, group_output, to_mesh.outputs['Mesh'])
    return group


def dashed_bond_node_group():
    """Alternating cylinder segments along the bond.

    Built by resampling the line, converting it to a wire mesh, deleting
    every other edge and running the remainder through Curve to Mesh: the
    dashes come out aligned with the bond without any rotation maths.
    """
    group, is_new = _new_group('DASHED')
    if not is_new:
        return group

    dashes = group.interface.new_socket('dashes', in_out='INPUT',
                                        socket_type='NodeSocketInt')
    dashes.default_value = 6
    dashes.min_value = 1
    radius = group.interface.new_socket('bond radius', in_out='INPUT',
                                        socket_type='NodeSocketFloat')
    radius.default_value = 0.08
    radius.min_value = 0.0
    resolution = group.interface.new_socket('resolution', in_out='INPUT',
                                            socket_type='NodeSocketInt')
    resolution.default_value = 12
    resolution.min_value = 3
    _finish_group(group)

    nodes, links = group.nodes, group.links
    group_input, group_output, curve_line, color_a, color_b = _bond_front_end(group)

    # 2 points per dash: one segment drawn, one left as the gap
    segments = nodes.new('ShaderNodeMath')
    segments.name = 'Dash segments'
    segments.operation = 'MULTIPLY'
    segments.inputs[1].default_value = 2.0
    segments.location = (-300, -450)
    links.new(group_input.outputs['dashes'], segments.inputs[0])

    resample = nodes.new('GeometryNodeResampleCurve')
    resample.name = 'Dash Positions'
    resample.location = (-120, -100)
    links.new(curve_line.outputs['Curve'], resample.inputs['Curve'])
    links.new(segments.outputs[0], resample.inputs['Count'])

    store_color = _store_bond_color(group, resample.outputs['Curve'],
                                    color_a, color_b, (80, -100))

    # curve -> wire mesh, so the segments become edges we can drop
    wire = nodes.new('GeometryNodeCurveToMesh')
    wire.name = 'Wire'
    wire.location = (260, -100)
    links.new(store_color.outputs['Geometry'], wire.inputs['Curve'])

    edge_index = nodes.new('GeometryNodeInputIndex')
    edge_index.location = (200, -420)
    parity = nodes.new('ShaderNodeMath')
    parity.name = 'Edge parity'
    parity.operation = 'MODULO'
    parity.inputs[1].default_value = 2.0
    parity.location = (360, -420)
    links.new(edge_index.outputs[0], parity.inputs[0])

    is_gap = nodes.new('FunctionNodeCompare')
    is_gap.name = 'Is gap'
    is_gap.data_type = 'FLOAT'
    is_gap.operation = 'EQUAL'
    is_gap.location = (520, -420)
    links.new(parity.outputs[0], cin(is_gap, 0))
    cin(is_gap, 1).default_value = 1.0

    drop_gaps = nodes.new('GeometryNodeDeleteGeometry')
    drop_gaps.name = 'Drop gaps'
    drop_gaps.domain = 'EDGE'
    drop_gaps.mode = 'ALL'
    drop_gaps.location = (420, -100)
    links.new(wire.outputs['Mesh'], drop_gaps.inputs['Geometry'])
    links.new(is_gap.outputs[0], drop_gaps.inputs['Selection'])

    back_to_curve = nodes.new('GeometryNodeMeshToCurve')
    back_to_curve.name = 'Dashes'
    back_to_curve.location = (560, -100)
    links.new(drop_gaps.outputs['Geometry'], back_to_curve.inputs['Mesh'])

    profile = nodes.new('GeometryNodeCurvePrimitiveCircle')
    profile.name = 'Dash Profile'
    profile.location = (560, -400)
    links.new(group_input.outputs['resolution'], profile.inputs['Resolution'])
    links.new(group_input.outputs['bond radius'], profile.inputs['Radius'])

    to_mesh = nodes.new('GeometryNodeCurveToMesh')
    to_mesh.name = 'Dash Mesh'
    to_mesh.location = (640, -100)
    to_mesh.inputs['Fill Caps'].default_value = True
    links.new(back_to_curve.outputs['Curve'], to_mesh.inputs['Curve'])
    links.new(profile.outputs['Curve'], to_mesh.inputs['Profile Curve'])

    _bond_tail(group, group_input, group_output, to_mesh.outputs['Mesh'])
    return group


_GROUP_BUILDERS = {
    'DOTTED': dotted_bond_node_group,
    'SCALED': scaled_bond_node_group,
    'DASHED': dashed_bond_node_group,
}


def bond_node_group(style):
    return _GROUP_BUILDERS[style]()


def _socket_identifiers(group):
    return {item.name: item.identifier
            for item in group.interface.items_tree
            if getattr(item, 'in_out', None) == 'INPUT'}


def set_dotted_partner(structure_obj, index_a, index_b, connected=True):
    """Mark (or unmark) the solid bond between two atoms as replaced.

    Stores the partner index + 1 on both atoms, so atoms_and_bonds can
    delete that bond. 0 means "no dotted partner", which is also what a
    missing attribute reads as.
    """
    mesh = structure_obj.data
    # attribute writes do not stick while the mesh is open in the edit cage
    with _object_mode():
        if DOTTED_PARTNER_ATTRIBUTE not in mesh.attributes:
            if not connected:
                return
            mesh.attributes.new(name=DOTTED_PARTNER_ATTRIBUTE, type='INT',
                                domain='POINT')
        data = mesh.attributes[DOTTED_PARTNER_ATTRIBUTE].data
        if max(index_a, index_b) >= len(data):
            raise ValueError('atom index outside the structure')
        data[index_a].value = index_b + 1 if connected else 0
        data[index_b].value = index_a + 1 if connected else 0
        mesh.update()
        structure_obj.update_tag()


def custom_bond_objects(structure_obj):
    """The custom-bond objects belonging to a structure."""
    found = []
    for collection in structure_obj.users_collection:
        for ob in collection.objects:
            if ob.get('ase_bond_structure') == structure_obj.name:
                found.append(ob)
    return found


def reset_custom_bonds(structure_obj, remove_objects=True):
    """Restore every solid bond that was replaced, and (by default) delete
    the custom-bond objects of this structure. Returns (restored, removed)."""
    mesh = structure_obj.data
    restored = 0
    # same as set_dotted_partner: clearing the attribute needs object mode
    with _object_mode():
        if DOTTED_PARTNER_ATTRIBUTE in mesh.attributes:
            data = mesh.attributes[DOTTED_PARTNER_ATTRIBUTE].data
            for entry in data:
                if entry.value:
                    restored += 1
                    entry.value = 0
            mesh.update()
            structure_obj.update_tag()

        removed = 0
        if remove_objects:
            for ob in custom_bond_objects(structure_obj):
                bpy.data.objects.remove(ob, do_unlink=True)
                removed += 1
    return restored, removed


def add_bond(structure_obj, index_a, index_b, style='DOTTED', segments=None,
             radius=None, resolution=None, reference_length=1.5,
             replace=False, outline=True):
    """Create a custom-bond object between two atoms of structure_obj.

    style: 'DOTTED', 'SCALED' or 'DASHED'. segments is the number of dots
    (dotted) or dashes (dashed) and is ignored for the scaled style;
    passing None for segments/radius/resolution keeps the style's default.
    """
    if style not in BOND_GROUPS:
        raise ValueError(f'unknown bond style {style!r}')
    if index_a == index_b:
        raise ValueError('pick two different atoms')
    n_atoms = len(structure_obj.data.vertices)
    if not (0 <= index_a < n_atoms and 0 <= index_b < n_atoms):
        raise ValueError(f'atom index outside the structure (0-{n_atoms - 1})')

    # The natural way to pick the two atoms is in edit mode, but object-level
    # operators (outline_objects) refuse to run there, and mesh attribute
    # writes (the 'replace' bookkeeping) do not stick while the mesh is open
    # in the edit cage. Drop to object mode for the duration and go back
    # afterwards, so the user stays where they were.
    with _object_mode():
        return _add_bond(structure_obj, index_a, index_b, style, segments,
                         radius, resolution, reference_length, replace, outline)


@contextlib.contextmanager
def _object_mode():
    obj = bpy.context.view_layer.objects.active if bpy.context.view_layer else None
    previous = obj.mode if obj is not None else 'OBJECT'
    if previous != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    try:
        yield
    finally:
        if previous != 'OBJECT':
            restore = bpy.context.view_layer.objects.active
            if restore is not None and restore.mode == 'OBJECT':
                try:
                    bpy.ops.object.mode_set(mode=previous)
                except RuntimeError:
                    pass


def _add_bond(structure_obj, index_a, index_b, style, segments, radius,
              resolution, reference_length, replace, outline):
    group = bond_node_group(style)

    mesh = bpy.data.meshes.new(f'{structure_obj.name}_{style.lower()}_bond')
    # a single loose vertex: the geometry comes from the node group, but the
    # object still needs mesh data to carry the material slots
    mesh.from_pydata([(0.0, 0.0, 0.0)], [], [])
    obj = bpy.data.objects.new(mesh.name, mesh)
    for collection in structure_obj.users_collection:
        collection.objects.link(obj)
    if not obj.users_collection:
        bpy.context.scene.collection.objects.link(obj)
    obj.matrix_world = structure_obj.matrix_world.copy()

    # same material slots as the structure, so the bond material keeps its
    # index and can be swapped in the object's material tab like any other.
    # The outline material is skipped: outline_objects() appends its own, and
    # a duplicate slot would make Blender resolve its Set Material node to the
    # wrong index.
    bond_slot = 0
    for slot in structure_obj.data.materials:
        if slot is not None and slot.name.startswith(OUTLINE_MATERIAL):
            continue
        if slot is not None and slot.name.startswith('BOND'):
            bond_slot = len(mesh.materials)
        mesh.materials.append(slot)

    modifier = obj.modifiers.new(name='GeometryNodes', type='NODES')
    modifier.node_group = group
    idents = _socket_identifiers(group)
    set_mod_input(modifier, idents['structure'], structure_obj)
    set_mod_input(modifier, idents['atom A'], index_a)
    set_mod_input(modifier, idents['atom B'], index_b)
    set_mod_input(modifier, idents['material_index'], bond_slot)

    # style specific inputs, each optional
    segment_socket = {'DOTTED': 'dots', 'DASHED': 'dashes'}.get(style)
    if segment_socket and segments is not None:
        set_mod_input(modifier, idents[segment_socket], segments)
    radius_socket = 'dot radius' if style == 'DOTTED' else 'bond radius'
    if radius is not None:
        set_mod_input(modifier, idents[radius_socket], radius)
    if resolution is not None:
        set_mod_input(modifier, idents['resolution'], resolution)
    if style == 'SCALED':
        set_mod_input(modifier, idents['reference length'], reference_length)

    if outline:
        # outline_objects() makes its target the active object; put the
        # structure back so the ASE panel stays on it and another bond
        # can be added straight away
        view_layer = bpy.context.view_layer
        previous_active = view_layer.objects.active
        previous_selection = [ob for ob in view_layer.objects if ob.select_get()]
        outline_objects([obj], modifier='GeometryNodes.001')
        for ob in view_layer.objects:
            ob.select_set(ob in previous_selection)
        if previous_active is not None:
            view_layer.objects.active = previous_active

    if replace:
        set_dotted_partner(structure_obj, index_a, index_b, connected=True)

    obj['ase_dotted_bond'] = [index_a, index_b]
    obj['ase_bond_structure'] = structure_obj.name
    obj['ase_bond_style'] = style
    return obj


def add_dotted_bond(structure_obj, index_a, index_b, dots=10, radius=0.08,
                    resolution=2, replace=False, outline=True):
    """Backwards-compatible wrapper for the dotted style."""
    return add_bond(structure_obj, index_a, index_b, style='DOTTED',
                    segments=dots, radius=radius, resolution=resolution,
                    replace=replace, outline=outline)
