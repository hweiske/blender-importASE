"""Helpers for API differences between Blender versions."""
import bpy


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
    node.inputs['Selection'].default_value = selection
    node.inputs['Distance'].default_value = distance
