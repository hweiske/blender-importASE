import bpy
from bpy_extras.io_utils import ImportHelper
from ase import io
import ase
from ase.data import covalent_radii,colors
from ase.build import make_supercell
import numpy as np
from ase import Atoms
from os.path import join
import os
import ase.neighborlist

def draw_atoms(atoms,scale=1,representation="Balls'n'Sticks"):
    cnt = 0
    #bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0),segments = 16 ,ring_count = 16)
    bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    for n,atom in enumerate(atoms):
        ob = sphere.copy()
        ob.data = sphere.data.copy()
        ob.location = atom.position
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob) 
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].name = atom.symbol
        if representation == "Balls'n'Sticks":
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [covalent_radii[atom.number]*0.5*scale,]*3
        elif representation == 'Licorice':
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [0.1]*3
        else:
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [covalent_radii[atom.number]*scale,]*3
        print(bpy.data.node_groups)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials[atom.symbol])
        cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_sphere'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None

def draw_bonds(atoms):
    nl = ase.neighborlist.NeighborList([covalent_radii[atomic_number]*0.9  for atomic_number in atoms.numbers], self_interaction = False,bothways=True)
    nl.update(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='bonds')
    except:
        None
#    bpy.ops.surface.primitive_nurbs_surface_cylinder_add(radius=1.0, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_cylinder_add(vertices=16)
    bpy.ops.object.shade_smooth()
    bond = bpy.context.object
    bond.name = 'ref_bond'
    cnt = 0
    for atom in atoms:
        if nl.get_neighbors(atom.index)[0].size > 0:
            neighbors,offsets=nl.get_neighbors(atom.index)
            for neighbor,offset in zip(neighbors, offsets):
                displacements =[ 0.5*(atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)), 0.5*( atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor] )]
                for n,displacement in enumerate(displacements):
                    if n == 0:
                        location= atom.position+ (displacement/2)
    #                else:
     #                   location=atoms[neighbor].position + (displacement/2)
                    distance = atoms.get_distance(atom.index,neighbor,mic=True)/2
                    ob = bond.copy()
                    ob.data = bond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob) 
                    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}'
                    bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials[f'{atom.symbol}-bond'])
                    ob.location = location
                    ob.scale = (1,1,1)
                    #if n == 0:
                    ob.dimensions = (0.2,0.2,distance)
                    phi = np.arctan2(displacement[1], displacement[0]) 
                    theta = np.arccos(displacement[2] / distance) 
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    break
                cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_bond'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None
def draw_unit_cell(atoms):
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='cell')
    except:
        None
    
    #SETUP MATERIAL
    matu=bpy.data.materials.new(name = 'unit_cell')
    matu.use_nodes=True
    tu=matu.node_tree
    su=tu.nodes['Principled BSDF']
    COL=(0.1,0.1,0.1,1)
    su.inputs[0].default_value=COL
    bpy.ops.mesh.primitive_cylinder_add(vertices=16)
    bpy.ops.object.shade_smooth()
    cell = bpy.context.object
    cell.name = 'ref_cell'
    cnt = 0
    X=[0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0]
    Y=[0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,0]
    Z=[0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1]
    for n in range(1,len(X)):
        pos1=np.array([X[n-1],Y[n-1],Z[n-1]])
        pos2=np.array([X[n],Y[n],Z[n]])
        location1=np.dot(pos1,atoms.cell)
        location2=np.dot(pos2,atoms.cell)
        print(n,location1,location2,pos1,pos2)
        displacement=location2-location1
        distance=np.linalg.norm(location1-location2)
        ob=cell.copy()
        ob.data=cell.data.copy()
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
        ob.name = f'unitcell-cylinder'
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials['unit_cell'])
        ob.location=location1+(displacement/2)
        ob.scale=(1,1,1)
        ob.dimensions=(0.1,0.1,distance)
        phi=np.arctan2(displacement[1],displacement[0])
        print(displacement,distance)
        theta=np.arccos(displacement[2]/distance)
        ob.rotation_euler[1]=theta
        ob.rotation_euler[2]=phi
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_cell'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None
