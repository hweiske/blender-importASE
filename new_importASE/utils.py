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

def group_atoms(atoms):
    atom_types = set(atoms.get_chemical_symbols())
    bpy.ops.object.select_all(action='DESELECT')
    for atom_type in atom_types:
        try:
            bpy.ops.group.create(name=atom_type)
        except:
            None
    return None    
def setup_materials(atoms,colorbonds=False,color=0.2):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0),segments = 16 ,ring_count = 16)
    bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    bpy.data.objects['ref_sphere'].select_set(True)

    color_dict={  'H'     :(      1.00, 1.00, 1.00        ),
        'C'     :(      0.05, 0.05, 0.05        ),
        'Si'    :(      0.001, 0.093, 0.314        ),
        'Ge'    :(      0.024,0.212,0.212  ),
        'Ga'    :(      0.33, 0.71, 0.09        ),
        'In'    :(     0,0,0            ),
        'N'     :(      0.00, 0.00, 1.00        ),
        'P'     :(      1.00, 0.50, 0.00        ),
        'As'    :(      0.75, 0.54, 0.00        ),
        'Sb'    :(      0.74, 0.46, 0.17        ),
        'Bi'    :(      0.82, 0.71, 0.55        ),
        'O'     :(      1.00, 0.00, 0.00        ),
        'S'     :(      1.00, 1.00, 0.00        ),
        'F'     :(      0.00, 1.00, 0.00        ),
        'Cl'    :(      0.50, 1.00, 0.00        ),
        'Br'    :(      0.39, 0.15, 0.03        ),
        'I'     :(      1.00, 0.00, 1.00        ),
        'Ti'    :(      0.25, 1.75, 0.75        )}
    atom_types = set(atoms.get_chemical_symbols())
    atom_n=list(set(atoms.numbers))
    for n,atom_type in enumerate(atom_types):
            matat=bpy.data.materials.new(name = str(atom_type))
            matb=bpy.data.materials.new(name = f'{atom_type}-bond')
            matb.use_nodes=True
            matat.use_nodes=True
            ta=matat.node_tree
            tb=matb.node_tree
            sa=ta.nodes['Principled BSDF']
            sb=tb.nodes['Principled BSDF']
            if atom_type in color_dict:
                COL=list(color_dict[atom_type]) + [1]
            else:
                COL=list(colors.jmol_colors[atom_n[n]]) + [1]
#            bpy.data.materials[str(atom_type)].diffuse_color = COL
            sa.inputs[0].default_value=COL
            if colorbonds == False:
 #               bpy.data.materials[f'{atom_type}-bond'].diffuse_color = [color,color,color,1]
                sb.inputs[0].default_value=[color,color,color,1]
            else:
                sb.inputs[0].default_value=COL
  #              bpy.data.materials[f'{atom_type}-bond'].diffuse_color = COL
    #        bpy.data.materials[str(atom_type)].metallic = 0.2
   #         bpy.data.materials[f'{atom_type}-bond'].metallic = 0.2

    bpy.data.objects['ref_sphere'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT') 
    return None 
def toggle(obj,SET=True):
        obj.hide_render = SET
        obj.hide_viewport = SET  # Optional: hide in the viewport as well
        for child in obj.children:
                child.hide_render = SET
                child.hide_viewport = SET  # Optional: hide in the viewport as well
        return(None)

