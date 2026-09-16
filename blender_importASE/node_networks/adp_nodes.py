"""Thermal ellipsoids with principal-axis rings in the atoms_and_bonds tree.

Reads the per-atom attributes written by the polyhedra importer (see
..adp): 'adp_rotation' (quaternion onto the principal axes), 'adp_axes'
(RMS displacements along them) and 'adp_valid'. Adds three modifier inputs:

- adps: switch between the ellipsoids and the normal spheres
- hydrogen_adps: ellipsoids for hydrogens too (off: the usual small spheres)
- adp_scale: RMS displacement -> ellipsoid semi-axis factor, set from the
  probability level at import (1.538 for 50 %)

The rings leave the tree as curves; adp_rings_node_group is the modifier
after the outline that turns them into tubes.
"""
import math

import bpy

import numpy as np

from ..adp import ADP_RINGS_MATERIAL, ellipsoids, probability_scale
from .compat import cin, setup_curve_to_mesh, set_mod_input
from .electron_density_nodes import newMaterial

# the set_atoms sphere has radius 0.5 before its per-atom scale
SPHERE_RADIUS = 0.5
RING_PROFILE_RESOLUTION = 8


def _input(node, name):
    """Input socket by name, the enabled one where a node has several
    (lookup by name alone can land on a disabled socket in Blender 4.x)."""
    sockets = [socket for socket in node.inputs if socket.name == name]
    if not sockets:
        raise KeyError(f"node {node.name!r} has no input {name!r}")
    return next((socket for socket in sockets if socket.enabled), sockets[0])


def _named(tree, name, data_type, location):
    node = tree.nodes.new("GeometryNodeInputNamedAttribute")
    node.name = f"Named Attribute.{name}"
    node.data_type = data_type
    node.inputs[0].default_value = name
    node.location = location
    return node


def add_adp_nodes(tree, atoms_geometry, points_geometry):
    """Add the ADP sockets and nodes to an atoms_and_bonds tree.

    atoms_geometry: output socket carrying the atom sphere instances
    points_geometry: output socket carrying the structure's points

    Returns (atoms, rings): the atom instances, turned into ellipsoids while
    the adps switch is on, and the rings as curves (empty while it is off),
    swept into tubes only after the outline by adp_rings_node_group.
    """
    nodes, links = tree.nodes, tree.links

    adps_socket = tree.interface.new_socket(
        name="adps", in_out='INPUT', socket_type='NodeSocketBool')
    adps_socket.default_value = True
    hydrogen_socket = tree.interface.new_socket(
        name="hydrogen_adps", in_out='INPUT', socket_type='NodeSocketBool')
    hydrogen_socket.description = ("draw hydrogens as ellipsoids too; off keeps "
                                   "them as the usual small spheres without rings")
    hydrogen_socket.default_value = False
    scale_socket = tree.interface.new_socket(
        name="adp_scale", in_out='INPUT', socket_type='NodeSocketFloat')
    scale_socket.description = ("RMS displacement to ellipsoid semi-axis: "
                                "1.538 draws 50 % probability ellipsoids")
    scale_socket.default_value = 1.5382
    scale_socket.min_value = 0.0
    scale_socket.max_value = 10.0

    x0, y0 = 300.0, -2600.0
    group_input = nodes.new("NodeGroupInput")
    group_input.name = "Group Input ADP"
    group_input.location = (x0 - 900, y0)

    rotation = _named(tree, "adp_rotation", 'QUATERNION', (x0 - 900, y0 - 250))
    axes = _named(tree, "adp_axes", 'FLOAT_VECTOR', (x0 - 900, y0 - 400))
    has_tensor = _named(tree, "adp_valid", 'BOOLEAN', (x0 - 900, y0 - 550))

    # an atom gets an ellipsoid when it has a usable tensor, and, if it is a
    # hydrogen, only while hydrogen_adps is on
    element = _named(tree, "element", 'INT', (x0 - 900, y0 - 700))
    is_not_hydrogen = nodes.new("FunctionNodeCompare")
    is_not_hydrogen.name = "ADP Not Hydrogen"
    is_not_hydrogen.data_type = 'INT'
    is_not_hydrogen.operation = 'NOT_EQUAL'
    is_not_hydrogen.location = (x0 - 700, y0 - 700)
    links.new(element.outputs[0], cin(is_not_hydrogen, 2))
    cin(is_not_hydrogen, 3).default_value = 1
    hydrogen_allowed = nodes.new("FunctionNodeBooleanMath")
    hydrogen_allowed.name = "ADP Hydrogen Allowed"
    hydrogen_allowed.operation = 'OR'
    hydrogen_allowed.location = (x0 - 500, y0 - 700)
    links.new(is_not_hydrogen.outputs[0], hydrogen_allowed.inputs[0])
    links.new(group_input.outputs['hydrogen_adps'], hydrogen_allowed.inputs[1])
    valid = nodes.new("FunctionNodeBooleanMath")
    valid.name = "ADP Selection"
    valid.operation = 'AND'
    valid.location = (x0 - 300, y0 - 600)
    links.new(has_tensor.outputs[0], valid.inputs[0])
    links.new(hydrogen_allowed.outputs[0], valid.inputs[1])

    # semi-axes = adp_scale * RMS displacement, relative to the sphere radius
    to_sphere = nodes.new("ShaderNodeMath")
    to_sphere.name = "ADP Scale To Sphere"
    to_sphere.operation = 'DIVIDE'
    to_sphere.inputs[1].default_value = SPHERE_RADIUS
    to_sphere.location = (x0 - 650, y0 - 100)
    links.new(group_input.outputs['adp_scale'], to_sphere.inputs[0])

    semi_axes = nodes.new("ShaderNodeVectorMath")
    semi_axes.name = "ADP Semi Axes"
    semi_axes.operation = 'SCALE'
    semi_axes.location = (x0 - 450, y0 - 400)
    links.new(axes.outputs[0], semi_axes.inputs[0])
    links.new(to_sphere.outputs[0], _input(semi_axes, 'Scale'))

    # --- ellipsoids: replace each sphere's uniform scale by the tensor's
    # rotation and semi-axes, keeping its position
    position = nodes.new("GeometryNodeInputPosition")
    position.location = (x0 - 450, y0 + 150)
    transform = nodes.new("FunctionNodeCombineTransform")
    transform.name = "ADP Transform"
    transform.location = (x0 - 200, y0)
    links.new(position.outputs[0], _input(transform, 'Translation'))
    links.new(rotation.outputs[0], _input(transform, 'Rotation'))
    links.new(semi_axes.outputs[0], _input(transform, 'Scale'))

    set_transform = nodes.new("GeometryNodeSetInstanceTransform")
    set_transform.name = "ADP Ellipsoids"
    set_transform.location = (x0, y0 + 100)
    links.new(atoms_geometry, _input(set_transform, 'Instances'))
    # atoms without a usable tensor stay spheres
    links.new(valid.outputs[0], _input(set_transform, 'Selection'))
    links.new(transform.outputs[0], _input(set_transform, 'Transform'))

    switch_atoms = nodes.new("GeometryNodeSwitch")
    switch_atoms.name = "ADP Switch Atoms"
    switch_atoms.input_type = 'GEOMETRY'
    switch_atoms.location = (x0 + 250, y0 + 100)
    links.new(group_input.outputs['adps'], switch_atoms.inputs[0])
    links.new(atoms_geometry, switch_atoms.inputs[1])
    links.new(set_transform.outputs[0], switch_atoms.inputs[2])

    # --- principal-axis rings: the three principal sections of the
    # ellipsoid, as unit-sphere great circles in the xy, xz and yz planes
    # that get the same rotation and scale. They stay curves, swept after
    # the outline (adp_rings_node_group), so the tube keeps a constant
    # thickness on a stretched ellipse and gets no outline shell
    ring_resolution = nodes.new("ShaderNodeMath")
    ring_resolution.name = "ADP Ring Resolution"
    ring_resolution.operation = 'MULTIPLY'
    ring_resolution.inputs[1].default_value = 2.0
    ring_resolution.location = (x0 - 650, y0 - 800)
    links.new(group_input.outputs['RESOLUTION'], ring_resolution.inputs[0])

    circle = nodes.new("GeometryNodeCurvePrimitiveCircle")
    circle.name = "ADP Ring"
    circle.mode = 'RADIUS'
    _input(circle, 'Radius').default_value = SPHERE_RADIUS
    circle.location = (x0 - 450, y0 - 800)
    links.new(ring_resolution.outputs[0], _input(circle, 'Resolution'))

    join_rings = nodes.new("GeometryNodeJoinGeometry")
    join_rings.name = "ADP Join Rings"
    join_rings.location = (x0 - 50, y0 - 800)
    links.new(circle.outputs[0], join_rings.inputs[0])
    for n, rotation_euler in enumerate(((math.pi / 2, 0.0, 0.0),
                                        (0.0, math.pi / 2, 0.0))):
        turn = nodes.new("GeometryNodeTransform")
        turn.name = f"ADP Ring Plane {n}"
        _input(turn, 'Rotation').default_value = rotation_euler
        turn.location = (x0 - 250, y0 - 900 - 150 * n)
        links.new(circle.outputs[0], _input(turn, 'Geometry'))
        links.new(turn.outputs[0], join_rings.inputs[0])

    ring_instances = nodes.new("GeometryNodeInstanceOnPoints")
    ring_instances.name = "ADP Ring Instances"
    ring_instances.location = (x0 + 150, y0 - 500)
    links.new(points_geometry, _input(ring_instances, 'Points'))
    links.new(valid.outputs[0], _input(ring_instances, 'Selection'))
    links.new(join_rings.outputs[0], _input(ring_instances, 'Instance'))
    links.new(rotation.outputs[0], _input(ring_instances, 'Rotation'))
    links.new(semi_axes.outputs[0], _input(ring_instances, 'Scale'))

    realize = nodes.new("GeometryNodeRealizeInstances")
    realize.name = "ADP Realize Rings"
    realize.location = (x0 + 350, y0 - 500)
    links.new(ring_instances.outputs[0], realize.inputs[0])

    switch_rings = nodes.new("GeometryNodeSwitch")
    switch_rings.name = "ADP Switch Rings"
    switch_rings.input_type = 'GEOMETRY'
    switch_rings.location = (x0 + 1150, y0 - 500)
    links.new(group_input.outputs['adps'], switch_rings.inputs[0])
    links.new(realize.outputs[0], switch_rings.inputs[2])

    return switch_atoms.outputs[0], switch_rings.outputs[0]


def adp_rings_node_group(material):
    """Modifier after the outline that sweeps the ring curves into tubes.

    The rings leave atoms_and_bonds as bare curves, so the outline (which
    only shells the mesh) passes them through untouched; here they become
    tubes of a constant `ring_radius` with the `material` set on them. Every
    other component of the geometry passes through as it is.
    """
    tree = bpy.data.node_groups.new(type='GeometryNodeTree', name="adp_rings")
    tree.is_modifier = True
    nodes, links = tree.nodes, tree.links

    tree.interface.new_socket(name="Geometry", in_out='OUTPUT',
                              socket_type='NodeSocketGeometry')
    tree.interface.new_socket(name="Geometry", in_out='INPUT',
                              socket_type='NodeSocketGeometry')
    radius_socket = tree.interface.new_socket(
        name="ring_radius", in_out='INPUT', socket_type='NodeSocketFloat')
    radius_socket.default_value = 0.012
    radius_socket.min_value = 0.0
    radius_socket.max_value = 1.0
    radius_socket.subtype = 'DISTANCE'
    resolution_socket = tree.interface.new_socket(
        name="profile_resolution", in_out='INPUT', socket_type='NodeSocketInt')
    resolution_socket.default_value = RING_PROFILE_RESOLUTION
    resolution_socket.min_value = 3
    resolution_socket.max_value = 64
    material_socket = tree.interface.new_socket(
        name="ring_material", in_out='INPUT', socket_type='NodeSocketMaterial')
    material_socket.default_value = material

    group_input = nodes.new("NodeGroupInput")
    group_input.location = (-600, 0)
    group_output = nodes.new("NodeGroupOutput")
    group_output.is_active_output = True
    group_output.location = (800, 0)

    separate = nodes.new("GeometryNodeSeparateComponents")
    separate.name = "Separate Components"
    separate.location = (-350, 0)
    links.new(group_input.outputs['Geometry'], separate.inputs[0])

    profile = nodes.new("GeometryNodeCurvePrimitiveCircle")
    profile.name = "Ring Profile"
    profile.mode = 'RADIUS'
    profile.location = (-350, -300)
    links.new(group_input.outputs['ring_radius'], _input(profile, 'Radius'))
    links.new(group_input.outputs['profile_resolution'], _input(profile, 'Resolution'))

    sweep = nodes.new("GeometryNodeCurveToMesh")
    sweep.name = "Ring Tubes"
    sweep.location = (-100, -200)
    setup_curve_to_mesh(tree, sweep, fill_caps=False, use_radius=False)
    links.new(separate.outputs['Curve'], _input(sweep, 'Curve'))
    links.new(profile.outputs[0], _input(sweep, 'Profile Curve'))

    set_material = nodes.new("GeometryNodeSetMaterial")
    set_material.name = "Ring Material"
    set_material.location = (150, -200)
    links.new(sweep.outputs[0], set_material.inputs['Geometry'])
    links.new(group_input.outputs['ring_material'], _input(set_material, 'Material'))

    smooth = nodes.new("GeometryNodeSetShadeSmooth")
    smooth.name = "Ring Smooth"
    smooth.location = (350, -200)
    links.new(set_material.outputs[0], smooth.inputs[0])

    # Set Material (unlike a raw material index) survives the join: the
    # material lists of the two meshes are merged and the indices remapped
    join = nodes.new("GeometryNodeJoinGeometry")
    join.name = "Join Rings"
    join.location = (600, 0)
    links.new(smooth.outputs[0], join.inputs[0])
    for socket in separate.outputs:
        if socket.name != 'Curve':
            links.new(socket, join.inputs[0])
    links.new(join.outputs[0], group_output.inputs[0])
    return tree


def store_adps(mesh, U):
    """Write the thermal ellipsoids of the per-vertex displacement tensors
    U as the adp_rotation / adp_axes / adp_valid point attributes the
    atoms_and_bonds tree reads."""
    rotations, axes, valid = ellipsoids(U)
    missing = int((~np.isfinite(U).all(axis=(1, 2))).sum())
    npd = int(len(U) - missing - valid.sum())
    if missing or npd:
        print(f'adp: {missing} atoms without displacement parameters, '
              f'{npd} with a non-positive-definite tensor - drawn as spheres')
    for attribute, kind, prop, values in (
            ('adp_rotation', 'QUATERNION', 'value', rotations),
            ('adp_axes', 'FLOAT_VECTOR', 'vector', axes),
            ('adp_valid', 'BOOLEAN', 'value', valid)):
        if attribute in mesh.attributes:
            mesh.attributes.remove(mesh.attributes[attribute])
        mesh.attributes.new(name=attribute, type=kind, domain='POINT')
        mesh.attributes[attribute].data.foreach_set(prop, values.ravel())
    mesh.update()


def adp_rings_material():
    """Flat black material of the principal-axis rings: a black Emission
    shader, so the rings read as ink lines under any lighting."""
    mat = newMaterial(ADP_RINGS_MATERIAL)
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    if any(node.bl_idname == 'ShaderNodeEmission' for node in nodes):
        return mat
    for node in list(nodes):
        nodes.remove(node)
    emission = nodes.new('ShaderNodeEmission')
    emission.inputs['Color'].default_value = (0.0, 0.0, 0.0, 1.0)
    emission.inputs['Strength'].default_value = 1.0
    output = nodes.new('ShaderNodeOutputMaterial')
    output.location = (250, 0)
    links.new(emission.outputs['Emission'], output.inputs['Surface'])
    return mat


def setup_adp_inputs(atoms_mod, probability=0.5, hydrogen_adps=False):
    """Set the ellipsoid size for a probability level and the hydrogen
    switch on an atoms_and_bonds modifier built with_adps."""
    tree = atoms_mod.node_group
    set_mod_input(atoms_mod, tree.interface.items_tree['adp_scale'].identifier,
                  probability_scale(probability))
    set_mod_input(atoms_mod, tree.interface.items_tree['hydrogen_adps'].identifier,
                  hydrogen_adps)


def add_adp_rings(obj):
    """Append the modifier sweeping the ring curves into black tubes - after
    the outline, which lets curves through unshelled."""
    material = adp_rings_material()
    obj.data.materials.append(material)
    rings = obj.modifiers.new(name='adp_rings', type='NODES')
    rings.node_group = adp_rings_node_group(material)
    set_mod_input(rings, rings.node_group.interface.items_tree['ring_material'].identifier,
                  material)
    return rings
