import bpy
import ase
from ase import Atoms
from .import_cubefiles import cube2vol
from .utils import atomcolors, group_atoms
from .drawobjects import draw_atoms, draw_bonds, draw_unit_cell, draw_bonds_new
from .trajectory import move_atoms, move_bonds,move_longbonds
from .node_networks.nodes_atoms_and_bonds import set_atoms_node_group, atoms_and_bonds, read_structure
from .node_networks.supercell import make_supercell
from .node_networks.bond_mat import create_bondmat
from .node_networks.outline import outline_objects
from .node_networks.bond_node import make_bonds
from .node_networks.hide_atoms import hide_atoms
import time
def import_ase_molecule(filepath, filename, overwrite=True, add_supercell=True, resolution=16, colorbonds=False, fix_bonds=False, color=0.2, scale=1,
                        unit_cell=False,
                        representation="Balls'n'Sticks",
                        read_density=True, shift_cell=False, 
                        imageslice=1, animate = True, outline = True, **kwargs):
    start=time.time()
    modifier_counter = 0
    modifier_chosen=''
    atoms = ase.io.read(filepath,index = ':')
    end_read=time.time()
    print('Time to read file: ',end_read-start)
    trajectory = False
    if isinstance(atoms[0],Atoms) and len(atoms) > 1:
        trajectory = True
        TRAJECTORY=atoms.copy()[1:]
        atoms=atoms[0]
        if animate is False:
            atoms = TRAJECTORY[-1]
        else:
            bpy.data.scenes['Scene'].frame_end = len(TRAJECTORY[::imageslice])
            bpy.data.scenes['Scene'].frame_start = 0
            bpy.data.scenes['Scene'].frame_current = 0
    elif len(atoms) == 1:
        atoms=atoms[0]
    if overwrite and representation != 'nodes' and animate and trajectory:
        self.representation = 'nodes'
    # When importing molecules from AMS, the resulting atoms do not lie in the unit cell since AMS uses unit cells centered around 0
    cell = atoms.cell
    shift_vector = 0.5 * cell[0] + 0.5 * cell[1] + 0.5 * cell[2]
    if shift_cell:
        atoms.positions += shift_vector
    #    if SUPERCELL == True:
    #        atoms=make_supercell(atoms,matrix)
    atomcolor=atomcolors()

    atomcolor.setup_materials(atoms, colorbonds=colorbonds, color=color)        
    my_coll = bpy.data.collections.new(name=atoms.get_chemical_formula() + '_' + filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection
    if representation != 'bonds_fromnodes' and representation != 'nodes':
        group_atoms(atoms)
        list_of_atoms=draw_atoms(atoms, scale=scale,resolution=resolution ,representation=representation)
        if representation != 'VDW':
            if fix_bonds:
                list_of_bonds,nl,bondlengths=draw_bonds_new(atoms,resolution=resolution)
            else:
                list_of_bonds,nl=draw_bonds(atoms,resolution=resolution)
        group_atoms(atoms)
        list_of_atoms=draw_atoms(atoms, scale=scale,resolution=resolution ,representation=representation)
        if representation != 'VDW':
            if fix_bonds:
                list_of_bonds,nl,bondlengths=draw_bonds_new(atoms,resolution=resolution)
            else:
                list_of_bonds,nl=draw_bonds(atoms,resolution=resolution)
    if representation == 'nodes':
        
        if animate and trajectory:
            obj,mesh=read_structure(TRAJECTORY[::imageslice],atoms.get_chemical_formula() + '_' + filename.split('.')[0],animate=True)

            

        else:
            obj,mesh=read_structure(atoms,atoms.get_chemical_formula() + '_' + filename.split('.')[0],animate=False)
        print(f'add hide modifier to GeometryNodes{modifier_chosen}')
        hide_atoms(obj,atoms,modifier='GeometryNodes'+modifier_chosen)
        modifier_counter += 1
        modifier_chosen=f'.00{modifier_counter}'
        if add_supercell:
            print(f'add supercell modifier to GeometryNodes{modifier_chosen}')
            added=make_supercell([obj], atoms, 'GeometryNodes'+modifier_chosen,representation=representation)
            if added:
                modifier_counter += 1
                modifier_chosen=f'.00{modifier_counter}'
        set_atoms_node_group()
        create_bondmat()
        print(f'add atoms_and_bonds modifier to GeometryNodes{modifier_chosen}')
        atoms_from_verts = atoms_and_bonds(obj,atoms,'GeometryNodes'+modifier_chosen)
       
        bpy.context.object.modifiers['GeometryNodes'+modifier_chosen].node_group = atoms_from_verts
        bpy.context.object.modifiers['GeometryNodes'+modifier_chosen]["Socket_2"] = 0.66
        bpy.context.object.modifiers['GeometryNodes'+modifier_chosen]["Socket_3"] = 0.1
        bpy.context.object.modifiers['GeometryNodes'+modifier_chosen]["Socket_4"] = 16
        modifier_counter += 1
        modifier_chosen=f'.00{modifier_counter}'
        

        #bond_nodes = bond_nodes_node_group(atoms, atoms_from_verts)
    if representation == 'bonds_fromnodes':
       
        sec_coll = bpy.data.collections.new(name='atoms')
        my_coll.children.link(sec_coll)
        seclayer_collection = layer_collection.children[sec_coll.name]
        bpy.context.view_layer.active_layer_collection = seclayer_collection
        group_atoms(atoms)
        list_of_atoms=draw_atoms(atoms, scale=scale,resolution=resolution ,representation=representation)
       
        

        bpy.context.view_layer.active_layer_collection = layer_collection
        bonds_obj = make_bonds(modifier='GeometryNodes')
        bpy.context.object.modifiers['GeometryNodes']["Socket_1"] = 0.66
        bpy.context.object.modifiers['GeometryNodes']["Socket_2"] = 0.1
        bpy.context.object.modifiers['GeometryNodes']["Socket_3"] = sec_coll
        


    if unit_cell is True and atoms.pbc.all() is not False:
        draw_unit_cell(atoms)
    if read_density:
        if 'cube' in filename:
            density_obj = cube2vol(filepath,modifier='GeometryNodes'+modifier_chosen)
            modifier_counter += 1
            modifier_chosen=f'.00{modifier_counter}'
            #print(density_obj)
            if shift_cell is True:
                density_obj.location.x += shift_vector[0]
                density_obj.location.y += shift_vector[1]
                density_obj.location.z += shift_vector[2]
    if trajectory is True and animate is True:
        if representation != 'nodes' and representation != 'bonds_fromnodes':
            
            move_atoms(TRAJECTORY,list_of_atoms,imageslice)
            if representation != 'VDW':
                if fix_bonds is True:
                    move_longbonds(TRAJECTORY,list_of_bonds,nl,bondlengths,imageslice)
                else:
                    move_bonds(TRAJECTORY,list_of_bonds,nl,imageslice)  
    end=time.time()
    print('Time to import atoms_object: ',end-end_read)

    if outline:
        print(f'add outline modifier to GeometryNodes{modifier_chosen}')
        if representation != 'nodes' and representation != 'bonds_fromnodes' and representation != 'VDW':
            outline_objects(list_of_atoms + list_of_bonds,modifier='GeometryNodes'+modifier_chosen)
            modifier_counter += 1
            modifier_chosen=f'.00{modifier_counter}'
             
        if representation == 'bonds_fromnodes':
            outline_objects([bonds_obj],modifier='GeometryNodes.001')
            outline_objects(list_of_atoms,modifier='GeometryNodes')
            modifier_counter += 1
            modifier_chosen=f'.00{modifier_counter}'
            
        if representation == 'nodes':
            outline_objects([obj],modifier='GeometryNodes'+modifier_chosen)
        if representation == 'VDW':
            outline_objects(list_of_atoms,modifier='GeometryNodes'+modifier_chosen)
            modifier_counter += 1
            modifier_chosen=f'.00{modifier_counter}'

    if add_supercell:
        if representation == 'bonds_fromnodes':
            added=make_supercell(list_of_atoms, atoms, 'GeometryNodes'+modifier_chosen,representation=representation)
            if added:
                print(f'added supercell to GeometryNodes{modifier_chosen}')
                modifier_counter += 1
                modifier_chosen=f'.00{modifier_counter}'
        
        if representation == 'VDW':
            added=make_supercell(list_of_atoms, atoms, 'GeometryNodes'+modifier_chosen,representation=representation)       
        
        if representation == 'licorice' or representation == "Balls'n'Sticks":
            added=make_supercell(list_of_atoms, atoms, 'GeometryNodes'+modifier_chosen,representation=representation) 
    
    
    END_END=time.time()
    print('modifier time: ',END_END-end)

            
