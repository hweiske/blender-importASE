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
    bpy.ops.object.shade_smooth(use_auto_smooth=True)
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    bpy.data.objects['ref_sphere'].select_set(True)
    color_dict={
    'H'     :(      1,1,1                   ),
    'C'     :(      0.05, 0.05, 0.05        ),
    'Al'    :(      0.6,0.42,0.42           ),#
    'B'     :(      0.8,0.33,0.43           ),#
    'Si'    :(      0.001, 0.093, 0.314     ),
    'Ge'    :(      0.024,0.22,0.22         ),
    'Ga'    :(      0.33, 0.71, 0.09        ),
    'In'    :(      0,0,0.01                 ),
    'N'     :(      0.00, 0.00, 1.00        ),
    'P'     :(      0.413, 0.013, 0.0004    ),
    'As'    :(      0.482, 0.378, 0.00      ),
    'Sb'    :(      0.74, 0.46, 0.17        ),
    'Bi'    :(      0.82, 0.71, 0.55        ),
    'O'     :(      1.00, 0.00, 0.00        ),
    'S'     :(      1.00, 1.00, 0.00        ),
    'Se'    :(      1,0.3,0                 ), #
    'Te'    :(      0.2,0.102,0.017         ),
    'F'     :(      0.00, 1.00, 0.00        ),
    'Cl'    :(      0.50, 1.00, 0.00        ),
    'Br'    :(      0.39, 0.15, 0.03        ),
    'I'     :(      0.257, 0.00, 0.257      ),
    'Ti'    :(      0.25, 1.75, 0.75        ),
    'Au'    :(      0.518,0.312,0.006       ),
    'Cu'    :(      0.594,0.1,0.048         ),
    'Ag'    :(        0.59,0.59,0.59        ),
    'Hf'    :(        0.365,0.509,0.920     ),
    }
    color_dict_VMD_CHEMCRAFT={
    'H'     :(      1, 1, 1                   ),
    'C'     :(      0.19, 0.19, 0.19        ),
    'B'     :(      0.8, 0.33, 0.43         ),
    'Al'    :(      0.72, 0.62, 0.62        ),
    'Si'    :(      0.02, 0.38,	0.67        ),
    'Ge'    :(      0.05, 0.45, 0.45        ),
    'Ga'    :(      0.33, 0.71, 0.09        ),
    'In'    :(      0, 0, 0.01              ),
    'N'     :(      0.00, 0.00, 1.00        ),
    'P'     :(      0.413, 0.013, 0.0004    ),
    'As'    :(      0.482, 0.378, 0.00      ),
    'Sb'    :(      0.74, 0.46, 0.17        ),
    'Bi'    :(      0.82, 0.71, 0.55        ),
    'O'     :(      1.00, 0.00, 0.00        ),
    'S'     :(      1.00, 1.00, 0.00        ),
    'F'     :(      0.00, 1.00, 0.00        ),
    'Cl'    :(      0.50, 1.00, 0.00        ),
    'Br'    :(      0.39, 0.15, 0.03        ),
    'I'     :(      0.257, 0.00, 0.257      ),
    'Ti'    :(      0.25, 1.75, 0.75        ),
    'Au'    :(      0.518,0.312,0.006       ),
    'Cu'    :(      0.594,0.1,0.048         ),
    'Ag'    :(        0.59,0.59,0.59        ),
    'Hf'    :(        0.365,0.509,0.920     ),
        }
    roughness_dict={  
    'H'     :   0.5,
    'C'     :   0.5,
    'B'     :   0.5,
    'Si'    :   0.7,
    'Ge'    :   0.7,
    'Ga'    :   0.6,
    'In'    :   0.4,
    'N'     :   0.5,
    'P'     :   0.5,
    'As'    :   0.7,
    'Sb'    :   0.6,
    'Bi'    :   0.5,
    'O'     :   0.5,
    'S'     :   0.5,
    'F'     :   0.5,
    'Cl'    :   0.5,
    'Br'    :   0.5,
    'I'     :   0.5,
    'Ti'    :   0.3,
    'Cu'    :   0.4,
    'Au'    :   0.4,
    'Fe'    :   0.5,
    'Ag'    :   0.4,
    'Al'    :   0.4,
    'Se'    :   0.5, #
    'Te'    :   0.4,
    'Hf'    :   0.3,
    }
    metal_dict={    
    'H'     :   0,
    'C'     :   0.5,
    'B'     :   0,
    'Si'    :   0.6,
    'Ge'    :   0.6,
    'Ga'    :   0.6,
    'In'    :   0.6,
    'N'     :   0,
    'P'     :   0.2,
    'As'    :   0.8,
    'Sb'    :   0.7,
    'Bi'    :   0.9,
    'O'     :   0,
    'S'     :   0,
    'F'     :   0,
    'Cl'    :   0,
    'Br'    :   0,
    'I'     :   0,
    'Ti'    :   1,
    'Cu'    :   1,
    'Au'    :   1,
    'Fe'    :   1,
    'Ag'    :   1,
    'Al'    :   1,
    'Se'    :   0.5,
    'Te'    :   1,
    'Hf'    :   1,
    }
    specular_dict={'H'     :0.2,
    'C'     :0,
    'B'     :0.5,
    'Si'    :0.5,
    'Ge'    :0.5,
    'Ga'    :0.5,
    'In'    :0.5,
    'N'     :0.1,
    'P'     :0.7,
    'As'    :0.8,
    'Sb'    :0.7,
    'Bi'    :0.9,
    'O'     :0.2,
    'S'     :0.2,
    'F'     :0.2,
    'Cl'    :0.2,
    'Br'    :0.2,
    'I'     :0.2,
    'Ti'    :0.5,
    'Cu'    : 0.5,
    'Au'    : 0.5,
    'Fe'    : 0.5,
    'Ag'    : 0.5,
    'Al'    :   1,
    'Se'    :   0.5,
    'Te'    :   1,
}

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
    METALS=['Li','Be','Na','Mg','K','Ca','Rb','Sr','Cs','Ba','Fr','Ra',
                        'Sc','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
                         'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                        'La','Ta', 'W','Re','Os','Ir','Pt','Au','Hg',
                        'Ac','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn',
                        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
                        'Th','Pa', 'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'
                        'Al','Sn','Tl','Pb']
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
                if atom_type in METALS:    
                    metal  = 1
                else:
                    metal=0
            if atom_type in roughness_dict:
                rough=roughness_dict[atom_type]
            else:
                rough=default_roughness
            if atom_type in specular_dict:
                specular = specular_dict[atom_type]
            else:
                specular = 0.5
            sa.inputs[0].default_value=COL
            sa.inputs[12].default_value=0
            sa.inputs[3].default_value=1.45
            sa.inputs[2].default_value=rough
            sa.inputs[1].default_value=metal
            if colorbonds == False:
                sb.inputs[0].default_value=[color,color,color,1]
                sb.inputs[12].default_value=0
                sb.inputs[3].default_value=1.45
                sb.inputs[2].default_value=0.5
                sb.inputs[1].default_value=0
            else:
                sb.inputs[0].default_value=COL
                sb.inputs[12].default_value=0
                sb.inputs[3].default_value=1.45
                sb.inputs[2].default_value=0.7
                sb.inputs[1].default_value=0.2
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

