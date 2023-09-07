import bpy
import numpy as np
from ase.io.cube import read_cube
import pyopenvdb as vdb
import os
from os import path
from .setup_nodetree import visualize_edensity_node_group,newShader
from .utils import toggle
def cube2vol(filename,filepath=os.environ.get('HOME')):
    with open(filename, 'r') as f:
        atoms=read_cube(f,read_data=True, verbose=True)
    SX=atoms['spacing'][0][0]
    SY=atoms['spacing'][1][1]
    SZ=atoms['spacing'][2][2]
    ORIGIN=atoms['origin']
    VOLUME= atoms['data']  # Replace this with your own 3D NumPy array
    GRID=vdb.FloatGrid()
    GRID.copyFromArray(VOLUME.astype(float))
    GRID.transform=vdb.createLinearTransform([[SX,0,0,0],[0,SY,0,0],[0,0,SZ,0],[0,0,0,1]])
    GRID.gridClass = vdb.GridClass.FOG_VOLUME
    GRID.name='density'
    TMPFILE=filename.split('.')[-2]+'_density.vdb'
    vdb.write(TMPFILE,GRID)
    VOL=bpy.ops.object.volume_import(filepath=TMPFILE,location=ORIGIN)
#    os.remove(TMPFILE)
    #bpy.data.objects[TMPFILE.split('.')[-2].split('/')[-1]].select_set(True)
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
    visualize_edensity_node_group()
    bpy.ops.object.modifier_add(type='NODES')
    node = bpy.data.node_groups["visualize_edensity"]
    bpy.context.object.modifiers['GeometryNodes'].node_group=node
    mat = newShader("+ material", "diffuse", 1, 0, 0)
    bpy.context.active_object.data.materials.append(mat)
    mat = newShader("- material", "diffuse", 0, 0, 1)
    bpy.context.active_object.data.materials.append(mat)
    bpy.context.object.modifiers["GeometryNodes"]["Input_9"] = bpy.data.materials["+ material"]
    bpy.context.object.modifiers["GeometryNodes"]["Input_10"] = bpy.data.materials["- material"]
    toggle(bpy.context.object,SET=False)
    return(atoms['atoms'])

