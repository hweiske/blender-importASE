"""Helpers for API differences between Blender versions."""
import bpy


def cin(node, old_index):
    """Operand socket of a FunctionNodeCompare, addressed by its Blender <=5.1
    index (INT operands were A=2, B=3).

    Before Blender 5.2 the Compare node exposed an A/B pair for every data type
    and only enabled the active pair; 5.2 collapses it to a single active A/B
    pair (A=0, B=1), so the old fixed indices raise IndexError. Selecting the
    enabled 'A'/'B' socket by name works on every version.
    """
    name = {0: 'A', 1: 'B', 2: 'A', 3: 'B'}[old_index]
    for socket in node.inputs:
        if socket.enabled and socket.name == name:
            return socket
    return node.inputs[old_index]


def pin(node, old_index):
    """Principled BSDF input socket addressed by its Blender <=5.1 index.

    Blender 5.2 inserted a 'Thin Wall' socket at index 5, shifting every socket
    after it by one (0-4 are unchanged). Add that offset on 5.2+.
    """
    if old_index >= 5 and bpy.app.version >= (5, 2, 0):
        return node.inputs[old_index + 1]
    return node.inputs[old_index]


def _input(node, name):
    # collection lookup by name skips disabled sockets in Blender 4.x
    # (e.g. 'Voxel Size' on a Volume to Mesh node in GRID mode), so iterate
    for socket in node.inputs:
        if socket.name == name:
            return socket
    raise KeyError(f"node {node.name!r} has no input {name!r}")


def setup_merge_by_distance(node, mode='ALL', selection=True, distance=0.001):
    """Configure a Merge by Distance node on Blender 4.x and 5.x.

    In Blender 5.0 the node's ``mode`` property became a 'Mode' menu input
    socket, which also shifted the socket indices.
    """
    if bpy.app.version >= (5, 0, 0):
        # the menu socket uses item names, not the old enum identifiers
        node.inputs['Mode'].default_value = {'ALL': 'All', 'CONNECTED': 'Connected'}[mode]
    else:
        node.mode = mode
    _input(node, 'Selection').default_value = selection
    _input(node, 'Distance').default_value = distance


def setup_curve_to_mesh(tree, node, fill_caps=False, use_radius=False):
    """Configure a Curve to Mesh node on Blender 4.4 and >= 4.5.

    Blender 4.5 added a 'Scale' input (shifting 'Fill Caps' from index 2 to
    3) and stopped multiplying the profile by the curve's Radius attribute
    implicitly. Pass use_radius=True when a Set Curve Radius node feeds this
    node: it wires the Radius attribute into 'Scale' to keep the old
    behavior. Leave it False otherwise - the Radius input node returns 0 for
    curves without a stored radius attribute, which would collapse the mesh.
    """
    _input(node, 'Fill Caps').default_value = fill_caps
    if use_radius and bpy.app.version >= (4, 5, 0):
        radius = tree.nodes.new('GeometryNodeInputRadius')
        radius.label = f"radius for {node.name}"
        tree.links.new(radius.outputs[0], _input(node, 'Scale'))
        return radius
    return None


def setup_volume_to_mesh(node, resolution_mode='GRID', voxel_size=0.3,
                         voxel_amount=64.0, adaptivity=0.0):
    """Configure a Volume to Mesh node on Blender 4.x and 5.x.

    In Blender 5.0 the node's ``resolution_mode`` property became a
    'Resolution Mode' menu input socket, which also shifted the socket
    indices.
    """
    if bpy.app.version >= (5, 0, 0):
        names = {'GRID': 'Grid', 'VOXEL_AMOUNT': 'Amount', 'VOXEL_SIZE': 'Size'}
        node.inputs['Resolution Mode'].default_value = names[resolution_mode]
    else:
        node.resolution_mode = resolution_mode
    _input(node, 'Voxel Size').default_value = voxel_size
    _input(node, 'Voxel Amount').default_value = voxel_amount
    _input(node, 'Adaptivity').default_value = adaptivity
