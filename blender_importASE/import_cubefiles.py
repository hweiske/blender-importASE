import bpy
from ase.io.cube import read_cube
if bpy.app.version[1] < 4:
    import pyopenvdb as vdb
else:
    import openvdb as vdb
import os
from .node_networks.electron_density_nodes import visualize_edensity_node_group
from .utils import toggle
import os.path


def cube2vol(filename, filepath=os.environ.get('HOME'),modifier='GeometryNodes'):
    with open(filename, 'r') as f:
        atoms = read_cube(f, read_data=True, verbose=True)
        ORIGIN = atoms['origin']
    VOLUME = atoms['data']
    GRID = vdb.FloatGrid()
    GRID.copyFromArray(VOLUME.astype(float))
    # SPACING=np.array(atoms['spacing']).T
    SPACING = atoms['spacing']
    SX = list(SPACING[0]) + [0.]
    SY = list(SPACING[1]) + [0.]
    SZ = list(SPACING[2]) + [0.]
    #    GRID.transform=vdb.createLinearTransform(np.array([SX,SY,SZ,[0,0,0,1]]).reshape((4,4)))
    GRID.transform = vdb.createLinearTransform(
        [[SX[0], SX[1], SX[2], SX[3]], [SY[0], SY[1], SY[2], SY[3]], [SZ[0], SZ[1], SZ[2], SZ[3]], [0, 0, 0, 1]])

    GRID.gridClass = vdb.GridClass.FOG_VOLUME
    GRID.name = 'density'
    TMPFILE = filename.split('.')[-2] + '_density.vdb'
    if not os.path.isfile(TMPFILE):
        vdb.write(TMPFILE, GRID)
    _ = bpy.ops.object.volume_import(filepath=TMPFILE, location=ORIGIN)
    #    os.remove(TMPFILE)
    # bpy.data.objects[TMPFILE.split('.')[-2].split('/')[-1]].select_set(True)
    #    for n,color in enumerate([[1,0,0,1],[0,0,1,1]]):
    #        if n == 0:
    #            name='+ material'
    #        else:
    #            name='- material'
    #        mat=bpy.data.materials.new(name = name)
    #        mat.use_nodes=True
    #        tree=mat.node_tree
    #        shader=tree.nodes['Principled BSDF']
    #        shader.inputs[0].default_value=color
    density_obj = bpy.context.active_object
    visualize_edensity_node_group()
    bpy.ops.object.modifier_add(type='NODES')
    node = bpy.data.node_groups["visualize_edensity"]
    bpy.context.object.modifiers[modifier].node_group = node
    bpy.context.object.modifiers[modifier]["Socket_9"] = bpy.data.materials["+ material"]
    bpy.context.object.modifiers[modifier]["Socket_10"] = bpy.data.materials["- material"]
    toggle(bpy.context.object, SET=False)
    return density_obj
