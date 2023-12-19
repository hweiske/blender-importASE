import bpy
from bpy_extras.io_utils import ImportHelper
from ase import io
import ase
from ase.data import covalent_radii, colors
from ase.build import make_supercell
import numpy as np
from ase import Atoms
from os.path import join
import os
import ase.neighborlist


def move_atoms(trajectory,list_of_atoms,imageslice):
#    view_layer=bpy.context.view_layer
    trajectory=trajectory[::imageslice]
    for ni,image in enumerate(trajectory):
        print(ni*imageslice)
        for n, atom_ob in enumerate(list_of_atoms):
            #print(ni,atom_ob)
            bpy.ops.object.select_all(action='DESELECT')
 #           bpy.context.active_object = atom_ob
            atom_ob.select_set(True)
            bpy.context.scene.frame_set(ni+1)
  #          bpy.context.active_object.location = image[n].position
            atom_ob.location=image[n].position
            if n == 1:
                print(ni,n,image[n].position,atom_ob.location)
            atom_ob.keyframe_insert(data_path='location')
            atom_ob.select_set(False)
    return None
def move_bonds(trajectory,list_of_bonds,NEIGHBORLIST,imageslice):
    trajectory=trajectory[::imageslice]
    for ni,image in enumerate(trajectory):
        cnt=0
        print(ni*imageslice)
        for na,atom in enumerate(image):
            neighbors, offsets = NEIGHBORLIST.get_neighbors(atom.index)
            for neighbor, offset in zip(neighbors, offsets):
                    displacements = [0.5 * (image.positions[neighbor] - atom.position + np.dot(offset, image.cell)),
                                     0.5 * (atom.position - np.dot(offset, image.cell) - image.positions[neighbor])]
                    for n, displacement in enumerate(displacements):
                        if n == 0:
                            location = atom.position + (displacement / 2)
                        distance = image.get_distance(atom.index, neighbor, mic=True) / 2
###ANIM###                        
                        ob=list_of_bonds[cnt]
                        ob.select_set(True)
                        bpy.context.scene.frame_set(ni+1)
                        ob.location = location
                        ob.dimensions = (0.2, 0.2, distance)
                        phi = np.arctan2(displacement[1], displacement[0])
                        theta = np.arccos(displacement[2] / distance)
                        ob.rotation_euler[1] = theta
                        ob.rotation_euler[2] = phi
                        cnt+=1
                        bpy.ops.anim.keyframe_insert(type='LocRotScale')
#                        ob.keyframe_insert(data_path='LocRotScale')
                        break
    print(f'plotted {cnt*len(trajectory)} bonds')
    return None
def move_longbonds(trajectory,list_of_bonds,NEIGHBORLIST,bondlengths,imageslice):
    print("Using Longbond mech")
    trajectory = trajectory[::imageslice]
    for ni,image in enumerate(trajectory):
        cnt=0
        print(ni*imageslice,cnt)
        print(cnt,len(list_of_bonds))
        for na,atom in enumerate(image):
            neighbors, offsets = NEIGHBORLIST.get_neighbors(atom.index)
            for neighbor, offset in zip(neighbors, offsets):
                neighbor_pos = image.positions[neighbor]
                neighbor_position = image.positions[neighbor] + offset.dot(image.cell)
                pos = atom.position
                ob=list_of_bonds[cnt]
                ob.select_set(True)
                bpy.context.scene.frame_set(ni+1)
                if bondlengths[cnt] == 'long': 
                    displacements = [0.5 * (image.positions[neighbor] - atom.position + np.dot(offset, image.cell)),
                                     0.5 * (atom.position - np.dot(offset, image.cell) - image.positions[neighbor])]
                    location = atom.position + (displacements[0] / 2)
                    distance = image.get_distance(atom.index, neighbor, mic=True)
                    ob.location = location
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / (distance / 2))
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    bpy.ops.anim.keyframe_insert(type='LocRotScale')
                    cnt+=1
                else:
                    ob=list_of_bonds[cnt]
                    ob.select_set(True)
                    bpy.context.scene.frame_set(ni+1)
                    displacements = [0.5 * (image.positions[neighbor] - atom.position + np.dot(offset, image.cell)),
                                     0.5 * (atom.position - np.dot(offset, image.cell) - image.positions[neighbor])]
                    location = atom.position + (displacements[0] / 2)
                    distance = image.get_distance(atom.index, neighbor, mic=True) / 2
                    ob.location = location
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / distance)
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    bpy.ops.anim.keyframe_insert(type='LocRotScale')
                    cnt += 1
                    # neighbor to atom
                    ob = list_of_bonds[cnt]
                    ob.select_set(True)
                    bpy.context.scene.frame_set(ni+1)
                    displacements = [0.5 * (atom.position - image.positions[neighbor] - np.dot(offset, image.cell)),
                                     0.5 * (image.positions[neighbor] + np.dot(offset, image.cell) + atom.position)]
                    location = image.positions[neighbor] + (displacements[0] / 2)
                    distance = image.get_distance(atom.index, neighbor, mic=True) / 2
                    ob.location = location
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / distance)
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    bpy.ops.anim.keyframe_insert(type='LocRotScale')
                    cnt += 1
    print(f'plotted {cnt*len(trajectory)} bonds')
    return None
