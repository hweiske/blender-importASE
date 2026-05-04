import bpy
from ase.data import colors, atomic_numbers

class atomcolors():
    def __init__(self):
        self.default_roughness=0.7
        self.color_dict={
    'H'     :(      1,1,1                   ),
    'Na'    :(      0.67, 0.36, 0.95        ),
    'C'     :(      0.05, 0.05, 0.05        ),
    'Al'    :(      0.6,0.42,0.42           ),#
    'B'     :(      0.8,0.33,0.43           ),#
    'Si'    :(      0.001, 0.093, 0.314     ),
    'Ge'    :(      0.024,0.22,0.22         ),
    'Ga'    :(      0.33, 0.71, 0.09        ),
    'In'    :(0.016710, 0.008966, 0.049616  ),
    'N'     :(      0.00, 0.00, 1.00        ),
    'P'     :( 0.099349, 0.004972, 0.000206 ),
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
    'Ta'    :(0.214031, 0.008155, 0.000296  ),
    'Gray'  :(      0.6,0.6,0.6             ),
    }
        self.roughness_dict={  
    'H'     :   0.5,
    'C'     :   0.5,
    'B'     :   0.5,
    'Si'    :   0.5,
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
    'Ta'    :   0.3,
    'Gray'  :   0.7,
    }
        self.metal_dict={    
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
    'Ta'    :   1,
    'Gray'  :   0.5,
    }


    
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
        self.METALS=['Li','Be','Na','Mg','K','Ca','Rb','Sr','Cs','Ba','Fr','Ra',
                        'Sc','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
                         'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                        'La','Ta', 'W','Re','Os','Ir','Pt','Au','Hg',
                        'Ac','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn',
                        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
                        'Th','Pa', 'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'
                        'Al','Sn','Tl','Pb']
    
    def setup_materials(self,atoms,colorbonds=False):
        atom_types = set(atoms.get_chemical_symbols())
        atom_n=list(set(atoms.numbers))
        atom_types.add('Gray')
        bpy.ops.object.select_all(action='DESELECT')
        bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0),segments = 16 ,ring_count = 16)
        if bpy.app.version[1] == 0: #use_auto_smoot dropped after 4.0
            bpy.ops.object.shade_smooth(use_auto_smooth=True)
        else:
            bpy.ops.object.shade_smooth()
        sphere = bpy.context.object
        sphere.name = 'ref_sphere'
        bpy.data.objects['ref_sphere'].select_set(True)
        
        for n,atom_type in enumerate(atom_types):
                if atom_type not in bpy.data.materials:
                    matat=bpy.data.materials.new(name = str(atom_type))
                else:
                    matat=bpy.data.materials[atom_type]
                if atom_type+'-bond' not in bpy.data.materials:
                    matb=bpy.data.materials.new(name = str(atom_type)+'-bond')
                else:
                    matb=bpy.data.materials[atom_type+'-bond']
                matb.use_nodes=True
                matat.use_nodes=True
                ta=matat.node_tree
                tb=matb.node_tree
                sa=ta.nodes['Principled BSDF']
                sb=tb.nodes['Principled BSDF']
                if atom_type in self.color_dict:
                    COL=list(self.color_dict[atom_type]) + [1]
                else:
                    COL=list(colors.jmol_colors[atomic_numbers[atom_type]]) + [1]
                if atom_type in self.metal_dict:
                    metal=self.metal_dict[atom_type]
                else:
                    if atom_type in self.METALS:    
                        metal  = 1
                    else:
                        metal=0
                if atom_type in self.roughness_dict:
                    rough=self.roughness_dict[atom_type]
                else:
                    rough=self.default_roughness
                # if atom_type in specular_dict:
                #     specular = specular_dict[atom_type]
                # else:
                #     specular = 0.5
                sa.inputs[0].default_value=COL
                sa.inputs[12].default_value=0
                sa.inputs[3].default_value=1.45
                sa.inputs[2].default_value=rough
                sa.inputs[1].default_value=metal
                sb.inputs[0].default_value=COL
                sb.inputs[12].default_value=0
                sb.inputs[3].default_value=1.45
                sb.inputs[2].default_value=rough
                sb.inputs[1].default_value=metal
        bpy.data.objects['ref_sphere'].select_set(True)
        bpy.ops.object.delete()
        bpy.ops.object.select_all(action='DESELECT') 
        return None 
    
    def create_bondmat(self, atom_1, atom_2, smooth, name):
        if atom_1 == atom_2:
            return bpy.data.materials[atom_1+'-bond']
        if smooth:
            name = f"{atom_1}-{atom_2}-bond-smooth"
        else:
            name = f"{atom_1}-{atom_2}-bond"
        if name in bpy.data.materials:
            mat = bpy.data.materials[name]
            return(mat)
        mat = bpy.data.materials.new(name = str(name))
        mat.use_nodes = True
        """Initialize Material.001 node group"""
        material_nodes = mat.node_tree

        # Start with a clean node tree
        for node in material_nodes.nodes:
            material_nodes.nodes.remove(node)
        material_nodes.color_tag = 'NONE'
        material_nodes.description = ""
        material_nodes.default_group_node_width = 140
        # material_001 interface

        # Initialize first atom material
        atom_type = atom_1
        if atom_type in self.color_dict:
            COL=list(self.color_dict[atom_type]) + [1]
        else:
            COL=list(colors.jmol_colors[atomic_numbers[atom_type]]) + [1]
        if atom_type in self.metal_dict:
            metal=self.metal_dict[atom_type]
        else:
            if atom_type in self.METALS:    
                metal  = 1
            else:
                metal=0
        if atom_type in self.roughness_dict:
            rough=self.roughness_dict[atom_type]
        else:
            rough=self.default_roughness

        # Node Principled BSDF
        principled_bsdf = material_nodes.nodes.new("ShaderNodeBsdfPrincipled")
        principled_bsdf.name = "Atom 1 Principled BSDF"
        principled_bsdf.distribution = 'MULTI_GGX'
        principled_bsdf.subsurface_method = 'RANDOM_WALK'
        # Base Color
        principled_bsdf.inputs[0].default_value = COL
        # Metallic
        principled_bsdf.inputs[1].default_value = metal
        # Roughness
        principled_bsdf.inputs[2].default_value = rough
        # IOR
        principled_bsdf.inputs[3].default_value = 1.45
        # Subsurface Anisotropy
        principled_bsdf.inputs[12].default_value = 0.0

        # Initialize second atom material
        atom_type = atom_2
        if atom_type in self.color_dict:
            COL=list(self.color_dict[atom_type]) + [1]
        else:
            COL=list(colors.jmol_colors[atomic_numbers[atom_type]]) + [1]
        if atom_type in self.metal_dict:
            metal=self.metal_dict[atom_type]
        else:
            if atom_type in self.METALS:    
                metal  = 1
            else:
                metal=0
        if atom_type in self.roughness_dict:
            rough=self.roughness_dict[atom_type]
        else:
            rough=self.default_roughness

        # Node Principled BSDF.001
        principled_bsdf_001 = material_nodes.nodes.new("ShaderNodeBsdfPrincipled")
        principled_bsdf_001.name = "Atom 2 Principled BSDF"
        principled_bsdf_001.distribution = 'MULTI_GGX'
        principled_bsdf_001.subsurface_method = 'RANDOM_WALK'
        # Base Color
        principled_bsdf_001.inputs[0].default_value = COL
        # Metallic
        principled_bsdf_001.inputs[1].default_value = metal
        # Roughness
        principled_bsdf_001.inputs[2].default_value = rough
        # IOR
        principled_bsdf_001.inputs[3].default_value = 1.45
        # Subsurface Anisotropy
        principled_bsdf_001.inputs[12].default_value = 0.0
        
        # Node Material Output
        material_output = material_nodes.nodes.new("ShaderNodeOutputMaterial")
        material_output.name = "Material Output"
        material_output.is_active_output = True
        material_output.target = 'ALL'
        # Displacement
        material_output.inputs[2].default_value = (0.0, 0.0, 0.0)
        # Thickness
        material_output.inputs[3].default_value = 0.0

        # Node Mix Shader
        mix_shader = material_nodes.nodes.new("ShaderNodeMixShader")
        mix_shader.name = "Mix Shader"

        # Node Color Ramp
        if smooth:
            limit = [0.25, 0.75]
        else:
            limit = [0.5, 0.5]

        color_ramp = material_nodes.nodes.new("ShaderNodeValToRGB")
        color_ramp.name = "Color Ramp"
        color_ramp.color_ramp.color_mode = 'RGB'
        color_ramp.color_ramp.hue_interpolation = 'NEAR'
        color_ramp.color_ramp.interpolation = 'LINEAR'

        # Initialize color ramp elements
        color_ramp.color_ramp.elements.remove(color_ramp.color_ramp.elements[0])
        color_ramp_cre_0 = color_ramp.color_ramp.elements[0]
        color_ramp_cre_0.position = limit[0]
        color_ramp_cre_0.alpha = 1.0
        color_ramp_cre_0.color = (0.0, 0.0, 0.0, 1.0)

        color_ramp_cre_1 = color_ramp.color_ramp.elements.new(limit[1])
        color_ramp_cre_1.alpha = 1.0
        color_ramp_cre_1.color = (1.0, 1.0, 1.0, 1.0)


        # Node Texture Coordinate
        texture_coordinate = material_nodes.nodes.new("ShaderNodeTexCoord")
        texture_coordinate.name = "Texture Coordinate"
        texture_coordinate.from_instancer = False

        # Node Separate XYZ
        separate_xyz = material_nodes.nodes.new("ShaderNodeSeparateXYZ")
        separate_xyz.name = "Separate XYZ"

        # Set locations
        principled_bsdf.location = (-40.0, 380.0)
        material_output.location = (920.0, 280.0)
        mix_shader.location = (660.0, 260.0)
        principled_bsdf_001.location = (-40.0, 20.0)
        color_ramp.location = (340.0, 660.0)
        texture_coordinate.location = (-320.0, 680.0)
        separate_xyz.location = (-40.0, 700.0)

        # Set dimensions
        principled_bsdf.width, principled_bsdf.height = 240.0, 100.0
        material_output.width, material_output.height = 140.0, 100.0
        mix_shader.width, mix_shader.height = 140.0, 100.0
        principled_bsdf_001.width, principled_bsdf_001.height = 240.0, 100.0
        color_ramp.width, color_ramp.height = 240.0, 100.0
        texture_coordinate.width, texture_coordinate.height = 140.0, 100.0
        separate_xyz.width, separate_xyz.height = 140.0, 100.0

        # Initialize material_001 links

        # mix_shader.Shader -> material_output.Surface
        material_nodes.links.new(mix_shader.outputs[0], material_output.inputs[0])
        # principled_bsdf.BSDF -> mix_shader.Shader
        material_nodes.links.new(principled_bsdf.outputs[0], mix_shader.inputs[1])
        # principled_bsdf_001.BSDF -> mix_shader.Shader
        material_nodes.links.new(principled_bsdf_001.outputs[0], mix_shader.inputs[2])
        # color_ramp.Color -> mix_shader.Fac
        material_nodes.links.new(color_ramp.outputs[0], mix_shader.inputs[0])
        # separate_xyz.Z -> color_ramp.Fac
        material_nodes.links.new(separate_xyz.outputs[2], color_ramp.inputs[0])
        # texture_coordinate.Generated -> separate_xyz.Vector
        material_nodes.links.new(texture_coordinate.outputs[0], separate_xyz.inputs[0])

        return mat



def group_atoms(atoms):
    atom_types = set(atoms.get_chemical_symbols())
    bpy.ops.object.select_all(action='DESELECT')
    for atom_type in atom_types:
        try:
            bpy.ops.group.create(name=atom_type)
        except Exception:
            pass
    return None    

def toggle(obj,SET=True):
        obj.hide_render = SET
        obj.hide_viewport = SET  # Optional: hide in the viewport as well
        for child in obj.children:
                child.hide_render = SET
                child.hide_viewport = SET  # Optional: hide in the viewport as well
        return(None)

