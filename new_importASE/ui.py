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
from .import_cubefiles import cube2vol
from .utils import setup_materials, group_atoms
from .drawobjects import draw_atoms, draw_bonds, draw_unit_cell, draw_bonds_new


def import_ase_molecule(filepath, filename, matrix, colorbonds=False, fixbonds=False, color=0.2, scale=1,
                        unit_cell=False,
                        representation="Balls'n'Sticks", separate_collections=False,
                        read_density=True, SUPERCELL=True, shift_cell=False
                        ):
    atoms = ase.io.read(filepath)
    # When importing molecules from AMS, the resulting atoms do not lie in the unit cell since AMS uses unit cells centered around 0
    shift_cell = True
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
    draw_atoms(atoms, scale=scale, representation=representation)
    if representation != 'VDW':
        if separate_collections:
            my_coll = bpy.data.collections.new(
                name=atoms.get_chemical_formula() + '_' + filename.split('.')[0] + '_bonds')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        if fixbonds:
            draw_bonds_new(atoms)
        else:
            draw_bonds(atoms)
    if unit_cell == True and atoms.pbc.all() != False:
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
            print(density_obj)
            density_obj.location.x += shift_vector[0]
            density_obj.location.y += shift_vector[1]
            density_obj.location.z += shift_vector[2]
            # name = cube2vol(filepath)
            # name = name.split('.')[0].split("/")[-1]
            # bpy.data.objects[name].location.x += shift_vector[0]
            # bpy.data.objects[name].location.y += shift_vector[1]
            # bpy.data.objects[name].location.z += shift_vector[2]


