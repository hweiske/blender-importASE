import bpy
import ase
from ase import Atoms
from .import_cubefiles import cube2vol
from .utils import setup_materials, group_atoms
from .drawobjects import draw_atoms, draw_bonds, draw_unit_cell, draw_bonds_new
from .trajectory import move_atoms, move_bonds,move_longbonds


def import_ase_molecule(filepath, filename, matrix, resolution=16, colorbonds=False, fix_bonds=False, color=0.2, scale=1,
                        unit_cell=False,
                        representation="Balls'n'Sticks", separate_collections=False,
                        read_density=True, SUPERCELL=True, shift_cell=False, 
                        imageslice=1, animate = True, **kwargs):
    
    atoms = ase.io.read(filepath,index = ':')
    if isinstance(atoms[0],Atoms) and len(atoms) > 1:
        trajectory = True
        TRAJECTORY=atoms.copy()[1:]
        atoms=atoms[0]
        if animate is False:
            atoms = TRAJECTORY[-1]
    elif len(atoms) == 1:
        trajectory = False
        atoms=atoms[0]
    # When importing molecules from AMS, the resulting atoms do not lie in the unit cell since AMS uses unit cells centered around 0
    cell = atoms.cell
    shift_vector = 0.5 * cell[0] + 0.5 * cell[1] + 0.5 * cell[2]
    if shift_cell:
        atoms.positions += shift_vector
    #    if SUPERCELL == True:
    #        atoms=make_supercell(atoms,matrix)
    setup_materials(atoms, colorbonds=colorbonds, color=color)
    if separate_collections:
        my_coll = bpy.data.collections.new(name=atoms.get_chemical_formula() + '_' + filename.split('.')[0] + '_atoms')
    else:
        my_coll = bpy.data.collections.new(name=atoms.get_chemical_formula() + '_' + filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection
    group_atoms(atoms)
    list_of_atoms=draw_atoms(atoms, scale=scale,resolution=resolution ,representation=representation)
    if representation != 'VDW':
        if separate_collections:
            my_coll = bpy.data.collections.new(
                name=atoms.get_chemical_formula() + '_' + filename.split('.')[0] + '_bonds')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        if fix_bonds:
            list_of_bonds,nl,bondlengths=draw_bonds_new(atoms,resolution=resolution)
        else:
            list_of_bonds,nl=draw_bonds(atoms,resolution=resolution)
    if unit_cell is True and atoms.pbc.all() is not False:
        if separate_collections:
            my_coll = bpy.data.collections.new(
                name=atoms.get_chemical_formula() + '_' + filename.split('.')[0] + '_cell')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        draw_unit_cell(atoms)
    if read_density:
        if 'cube' in filename:
            density_obj = cube2vol(filepath)
            #print(density_obj)
            if shift_cell is True:
                density_obj.location.x += shift_vector[0]
                density_obj.location.y += shift_vector[1]
                density_obj.location.z += shift_vector[2]
            # name = cube2vol(filepath)
            # name = name.split('.')[0].split("/")[-1]
            # bpy.data.objects[name].location.x += shift_vector[0]
            # bpy.data.objects[name].location.y += shift_vector[1]
            # bpy.data.objects[name].location.z += shift_vector[2]
    if trajectory is True and animate is True:

        move_atoms(TRAJECTORY,list_of_atoms,imageslice)
        if representation != 'VDW':
            if fix_bonds is True:
                move_longbonds(TRAJECTORY,list_of_bonds,nl,bondlengths,imageslice)
            else:
                move_bonds(TRAJECTORY,list_of_bonds,nl,imageslice)

            
