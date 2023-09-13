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
def setup_materials(atoms,colorbonds=False,color=0.6):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0),segments = 16 ,ring_count = 16)
    bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    bpy.data.objects['ref_sphere'].select_set(True)

    color_dict={  'H'     :(      1.00, 1.00, 1.00        ),
        'C'     :(      0.05, 0.05, 0.05        ),
        'Si'    :(      0.001, 0.093, 0.314        ),
        'Ge'    :(      0.024,0.22,0.22  ),
        'Ga'    :(      0.33, 0.71, 0.09        ),
        'In'    :(     0,0,0.01            ),
        'N'     :(      0.00, 0.00, 1.00        ),
        'P'     :(      0.413, 0.013, 0.0004        ),
        'As'    :(      0.482, 0.378, 0.00        ),
        'Sb'    :(      0.74, 0.46, 0.17        ),
        'Bi'    :(      0.82, 0.71, 0.55        ),
        'O'     :(      1.00, 0.00, 0.00        ),
        'S'     :(      1.00, 1.00, 0.00        ),
        'F'     :(      0.00, 1.00, 0.00        ),
        'Cl'    :(      0.50, 1.00, 0.00        ),
        'Br'    :(      0.39, 0.15, 0.03        ),
        'I'     :(      0.257, 0.00, 0.257        ),
        'Ti'    :(      0.25, 1.75, 0.75        ),
        'Au'    :(0.518,0.312,0.006),
        'Cu'    :(0.594,0.1,0.048),
        }
    roughness_dict={  'H'     : 0.7,
        'C'     :0.7,
        'Si'    :0.7,
        'Ge'    :0.7,
        'Ga'    :0.6,
        'In'    :0.6,
        'N'     :0.7,
        'P'     :0.7,
        'As'    :0.7,
        'Sb'    :0.6,
        'Bi'    :0.5,
        'O'     :0.7,
        'S'     :0.7,
        'F'     :0.7,
        'Cl'    :0.7,
        'Br'    :0.7,
        'I'     :0.7,
        'Ti'    :0.5,
        'Cu'    :0.5,
        'Au'    :0.5,
        'Fe'    :0.5,
        'Ag'    :0.5,
        }
    metal_dict={  'H'     :0.1,
        'C'     :0.5,
        'Si'    :0.6,
        'Ge'    :0.6,
        'Ga'    :0.6,
        'In'    :0.6,
        'N'     :0.1,
        'P'     :0.7,
        'As'    :0.8,
        'Sb'    :0.7,
        'Bi'    :0.9,
        'O'     :0.1,
        'S'     :0.1,
        'F'     :0.1,
        'Cl'    :0.1,
        'Br'    :0.1,
        'I'     :0.1,
        'Ti'    :1,
        'Cu'    : 1,
        'Au'    : 1,
        'Fe'    : 1,
        'Ag'    : 1,
        }
    default_metal=0.2
    default_roughness=0.7
    atom_types = set(atoms.get_chemical_symbols())
    atom_n=list(set(atoms.numbers))
    #base color 0
    #subsurface 1
    #subsurface radius 2
    #subsurface color 3
    #subsurface IOR 4
    #subsurface Anisortropy 5
    #Metallic 6 
    #specular 7
    #specular tint 8
    #roughness 9
    #anisotropic 10
    #anisotropic rotation 11
    #Sheen 12
    #sheen tint 13
    #clearcoat 14
    #clearcoat roughness 15
    #IOR 16
    #Transmission 17
    #transmission roughness 18
    #emission 19
    #emission strength 20
    #alpha 21
    #normal 22
    #Clearcoat Normal 23
    #Tangent 24
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
            if atom_type in metal_dict:
                metal=metal_dict[atom_type]
            else:
                metal=default_metal
            if atom_type in roughness_dict:
                rough=roughness_dict[atom_type]
            else:
                rough=default_roughness
#            bpy.data.materials[str(atom_type)].diffuse_color = COL
            sa.inputs[0].default_value=COL
            sa.inputs[9].default_value=rough
            sa.inputs[6].default_value=metal
            if colorbonds == False:
 #               bpy.data.materials[f'{atom_type}-bond'].diffuse_color = [color,color,color,1]
                sb.inputs[0].default_value=[color,color,color,1]
                sb.inputs[9].default_value=rough
                sb.inputs[6].default_value=metal
            else:
                sb.inputs[0].default_value=COL
                sb.inputs[9].default_value=0.7
                sb.inputs[6].default_value=0.2
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

