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
from .import_cubefiles import cube2vol
from .utils import setup_materials,group_atoms
from .drawobjects import draw_atoms,draw_bonds,draw_unit_cell

def import_ase_molecule(filepath,filename,matrix,colorbonds=False,color=0.2,scale=1,unit_cell=False, representation="Balls'n'Sticks",separate_collections=False,
        read_density=True,SUPERCELL=True
        ):
    atoms = ase.io.read(filepath)
#    if SUPERCELL == True:
#        atoms=make_supercell(atoms,matrix)
    setup_materials(atoms,colorbonds=colorbonds,color=color)
    if separate_collections == True:
        my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_atoms')
    else:
        my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection
    group_atoms(atoms)
    draw_atoms(atoms,scale=scale,representation=representation)
    if representation != 'VDW':
        if separate_collections == True:
            my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_bonds')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        draw_bonds(atoms)
    if unit_cell == True and atoms.pbc.all() != False:
        if separate_collections == True:
            my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_cell')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        draw_unit_cell(atoms)
    if read_density == True:
        if 'cube' in filename:
            atoms=cube2vol(filepath)
