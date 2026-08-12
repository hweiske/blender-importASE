"""Dotted bonds between two atoms of a structure.

A dotted bond is a row of small spheres along the line between two atoms
of an imported structure - useful for partial bonds in a transition
state, hydrogen bonds, or any interaction the distance-based bond search
does not draw.

The dots live in their own object whose geometry-node group samples the
two atom positions straight from the structure object, so they follow
the structure (including trajectory animation) without being baked in.
They are colored like a real bond: the group writes the same
'COLOR_CURVE' attribute the bond material reads, blended from the two
atoms' colors, and sets the material index of the structure's bond
material slot.

'Replace' additionally hides the solid bond between the two atoms. The
structure mesh carries a per-atom 'dotted_partner' int attribute holding
*the other atom's index + 1* (0, and therefore a missing attribute,
means "none"), which atoms_and_bonds ORs into its bond delete selection.
The +1 offset is what makes an absent attribute safe: it reads as 0,
which never matches a real index.
"""
import bpy

from . import __version__
from .node_networks.compat import set_mod_input
from .node_networks.outline import outline_objects

DOTTED_BOND_GROUP = 'dotted_bond'
DOTTED_PARTNER_ATTRIBUTE = 'dotted_partner'
OUTLINE_MATERIAL = 'outline_color'


def dotted_bond_node_group():
    """Build (or return) the dotted-bond geometry node group.

    Like the outline and bond groups, the tree is stamped with the add-on
    version: a group of this name from an older add-on - or a hand-built
    one that happens to share the name, which is likely since this started
    life as a manual node setup - is renamed out of the way and rebuilt,
    instead of being reused with an interface we do not recognise.
    """
    existing = bpy.data.node_groups.get(DOTTED_BOND_GROUP)
    if existing is not None:
        if existing.description == __version__:
            return existing
        existing.name = f'{DOTTED_BOND_GROUP}_old'

    group = bpy.data.node_groups.new(type='GeometryNodeTree', name=DOTTED_BOND_GROUP)
    group.description = __version__
    group.is_modifier = True

    group.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
    group.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')

    group.interface.new_socket('structure', in_out='INPUT',
                               socket_type='NodeSocketObject')
    group.interface.new_socket('atom A', in_out='INPUT', socket_type='NodeSocketInt')
    group.interface.new_socket('atom B', in_out='INPUT', socket_type='NodeSocketInt')
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
    mat_index = group.interface.new_socket('material_index', in_out='INPUT',
                                           socket_type='NodeSocketInt')
    mat_index.default_value = 0
    mat_index.min_value = 0

    nodes = group.nodes
    links = group.links

    group_input = nodes.new('NodeGroupInput')
    group_input.location = (-900, 0)
    group_output = nodes.new('NodeGroupOutput')
    group_output.location = (700, 0)

    object_info = nodes.new('GeometryNodeObjectInfo')
    object_info.name = 'Structure Info'
    object_info.location = (-700, -200)
    object_info.inputs['As Instance'].default_value = False
    links.new(group_input.outputs['structure'], object_info.inputs['Object'])

    position = nodes.new('GeometryNodeInputPosition')
    position.location = (-700, -400)

    # positions of the two atoms
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

    # the line between them, resampled into the dot positions
    curve_line = nodes.new('GeometryNodeCurvePrimitiveLine')
    curve_line.name = 'Bond Line'
    curve_line.mode = 'POINTS'
    curve_line.location = (-300, -100)
    links.new(sample_a.outputs['Value'], curve_line.inputs['Start'])
    links.new(sample_b.outputs['Value'], curve_line.inputs['End'])

    resample = nodes.new('GeometryNodeResampleCurve')
    resample.name = 'Dot Positions'
    resample.location = (-120, -100)
    links.new(curve_line.outputs['Curve'], resample.inputs['Curve'])
    links.new(group_input.outputs['dots'], resample.inputs['Count'])

    # color the dots like a real bond: blend the two atom colors along the
    # line into the 'COLOR_CURVE' attribute the bond material reads
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

    spline_parameter = nodes.new('GeometryNodeSplineParameter')
    spline_parameter.location = (-300, -560)

    mix_color = nodes.new('ShaderNodeMix')
    mix_color.name = 'Bond Color'
    mix_color.data_type = 'RGBA'
    mix_color.blend_type = 'MIX'
    mix_color.location = (-120, -560)
    links.new(spline_parameter.outputs['Factor'], mix_color.inputs['Factor'])
    links.new(color_a.outputs['Value'], mix_color.inputs[6])
    links.new(color_b.outputs['Value'], mix_color.inputs[7])

    store_color = nodes.new('GeometryNodeStoreNamedAttribute')
    store_color.name = 'Store COLOR_CURVE'
    store_color.data_type = 'FLOAT_COLOR'
    store_color.domain = 'POINT'
    store_color.inputs['Name'].default_value = 'COLOR_CURVE'
    store_color.location = (80, -100)
    links.new(resample.outputs['Curve'], store_color.inputs['Geometry'])
    links.new(mix_color.outputs[2], store_color.inputs['Value'])

    # a sphere per dot
    ico = nodes.new('GeometryNodeMeshIcoSphere')
    ico.name = 'Dot'
    ico.location = (80, -350)
    links.new(group_input.outputs['dot radius'], ico.inputs['Radius'])
    links.new(group_input.outputs['resolution'], ico.inputs['Subdivisions'])

    instance = nodes.new('GeometryNodeInstanceOnPoints')
    instance.name = 'Dots'
    instance.location = (280, -100)
    links.new(store_color.outputs['Geometry'], instance.inputs['Points'])
    links.new(ico.outputs['Mesh'], instance.inputs['Instance'])

    realize = nodes.new('GeometryNodeRealizeInstances')
    realize.location = (440, -100)
    links.new(instance.outputs['Instances'], realize.inputs['Geometry'])

    # Join the host mesh (a single loose vertex, renders nothing) so the
    # object's material slot list travels with the geometry - geometry built
    # from scratch inside the tree has an empty material list, and Set
    # Material Index would then be out of range and clamp to slot 0 as soon
    # as something (e.g. the outline modifier) joins it with other geometry.
    join = nodes.new('GeometryNodeJoinGeometry')
    join.name = 'Join Host Mesh'
    join.location = (500, -100)
    links.new(realize.outputs['Geometry'], join.inputs['Geometry'])
    links.new(group_input.outputs['Geometry'], join.inputs['Geometry'])

    # use the structure's bond material by slot index
    set_material_index = nodes.new('GeometryNodeSetMaterialIndex')
    set_material_index.name = 'Bond Material Slot'
    set_material_index.location = (560, -100)
    links.new(join.outputs['Geometry'], set_material_index.inputs['Geometry'])
    links.new(group_input.outputs['material_index'],
              set_material_index.inputs['Material Index'])
    links.new(set_material_index.outputs['Geometry'], group_output.inputs['Geometry'])

    return group


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
    if DOTTED_PARTNER_ATTRIBUTE not in mesh.attributes:
        if not connected:
            return
        mesh.attributes.new(name=DOTTED_PARTNER_ATTRIBUTE, type='INT', domain='POINT')
    data = mesh.attributes[DOTTED_PARTNER_ATTRIBUTE].data
    if max(index_a, index_b) >= len(data):
        raise ValueError('atom index outside the structure')
    data[index_a].value = index_b + 1 if connected else 0
    data[index_b].value = index_a + 1 if connected else 0
    mesh.update()
    structure_obj.update_tag()


def add_dotted_bond(structure_obj, index_a, index_b, dots=10, radius=0.08,
                    resolution=2, replace=False, outline=True):
    """Create a dotted-bond object between two atoms of structure_obj."""
    if index_a == index_b:
        raise ValueError('pick two different atoms')
    n_atoms = len(structure_obj.data.vertices)
    if not (0 <= index_a < n_atoms and 0 <= index_b < n_atoms):
        raise ValueError(f'atom index outside the structure (0-{n_atoms - 1})')

    group = dotted_bond_node_group()

    mesh = bpy.data.meshes.new(f'{structure_obj.name}_dotted_bond')
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
    set_mod_input(modifier, idents['dots'], dots)
    set_mod_input(modifier, idents['dot radius'], radius)
    set_mod_input(modifier, idents['resolution'], resolution)
    set_mod_input(modifier, idents['material_index'], bond_slot)

    if outline:
        # outline_objects() makes its target the active object; put the
        # structure back so the ASE panel stays on it and another dotted bond
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
    return obj
