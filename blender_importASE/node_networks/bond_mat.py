import bpy
from .compat import pin

def bond_nodes_node_group(mat,colorbonds=False):

    bond_nodes = mat.node_tree
    #start with a clean node tree
    for node in bond_nodes.nodes:
        bond_nodes.nodes.remove(node)
    bond_nodes.color_tag = 'NONE'
    bond_nodes.description = ""
    bond_nodes.default_group_node_width = 140
    

    #bond_nodes interface

    #initialize bond_nodes nodes
    #node Material Output
    material_output = bond_nodes.nodes.new("ShaderNodeOutputMaterial")
    material_output.name = "Material Output"
    material_output.is_active_output = True
    material_output.target = 'ALL'
    #Displacement
    material_output.inputs[2].default_value = (0.0, 0.0, 0.0)
    #Thickness
    material_output.inputs[3].default_value = 0.0

    #node Attribute.001
    attribute_001 = bond_nodes.nodes.new("ShaderNodeAttribute")
    attribute_001.name = "Attribute.001"
    attribute_001.attribute_name = "EMISSION_CURVE"
    attribute_001.attribute_type = 'GEOMETRY'

    #node Math
    math = bond_nodes.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = 0.10000000149011612

    #node Principled BSDF
    principled_bsdf = bond_nodes.nodes.new("ShaderNodeBsdfPrincipled")
    principled_bsdf.name = "Principled BSDF"
    principled_bsdf.distribution = 'GGX'
    principled_bsdf.subsurface_method = 'RANDOM_WALK_SKIN'
    #Metallic
    principled_bsdf.inputs[1].default_value = 0.0
    #Roughness
    principled_bsdf.inputs[2].default_value = 0.5
    #IOR
    principled_bsdf.inputs[3].default_value = 1.4500000476837158
    #Alpha
    principled_bsdf.inputs[4].default_value = 1.0
    #Normal
    pin(principled_bsdf, 5).default_value = (0.0, 0.0, 0.0)
    #Diffuse Roughness
    pin(principled_bsdf, 7).default_value = 0.0
    #Subsurface Weight
    pin(principled_bsdf, 8).default_value = 0.0
    #Subsurface Radius
    pin(principled_bsdf, 9).default_value = (1.0, 0.20000000298023224, 0.10000000149011612)
    #Subsurface Scale
    pin(principled_bsdf, 10).default_value = 0.05000000074505806
    #Subsurface IOR
    pin(principled_bsdf, 11).default_value = 1.399999976158142
    #Subsurface Anisotropy
    pin(principled_bsdf, 12).default_value = 0.0
    #Specular IOR Level
    pin(principled_bsdf, 13).default_value = 0.5
    #Specular Tint
    pin(principled_bsdf, 14).default_value = (1.0, 1.0, 1.0, 1.0)
    #Anisotropic
    pin(principled_bsdf, 15).default_value = 0.0
    #Anisotropic Rotation
    pin(principled_bsdf, 16).default_value = 0.0
    #Tangent
    pin(principled_bsdf, 17).default_value = (0.0, 0.0, 0.0)
    #Transmission Weight
    pin(principled_bsdf, 18).default_value = 0.05740181356668472
    #Coat Weight
    pin(principled_bsdf, 19).default_value = 0.0
    #Coat Roughness
    pin(principled_bsdf, 20).default_value = 0.029999999329447746
    #Coat IOR
    pin(principled_bsdf, 21).default_value = 1.5
    #Coat Tint
    pin(principled_bsdf, 22).default_value = (1.0, 1.0, 1.0, 1.0)
    #Coat Normal
    pin(principled_bsdf, 23).default_value = (0.0, 0.0, 0.0)
    #Sheen Weight
    pin(principled_bsdf, 24).default_value = 0.0
    #Sheen Roughness
    pin(principled_bsdf, 25).default_value = 0.5
    #Sheen Tint
    pin(principled_bsdf, 26).default_value = (1.0, 1.0, 1.0, 1.0)
    #Emission Color
    pin(principled_bsdf, 27).default_value = (0.0, 0.0, 0.0, 1.0)
    #Emission Strength
    pin(principled_bsdf, 28).default_value = 0.009999999776482582
    #Thin Film Thickness
    pin(principled_bsdf, 29).default_value = 0.0
    #Thin Film IOR
    pin(principled_bsdf, 30).default_value = 1.3300000429153442

    #node Attribute
    attribute = bond_nodes.nodes.new("ShaderNodeAttribute")
    attribute.name = "Attribute"
    attribute.attribute_name = "COLOR_CURVE"
    attribute.attribute_type = 'GEOMETRY'

    #node RGB
    rgb = bond_nodes.nodes.new("ShaderNodeRGB")
    rgb.name = "RGB"

    rgb.outputs[0].default_value = (0.6000000238418579, 0.6000000238418579, 0.6000000238418579, 1.0)
    #node Mix
    mix = bond_nodes.nodes.new("ShaderNodeMix")
    mix.label = "mix_atom_color"
    mix.name = "Mix"
    mix.blend_type = 'MIX'
    mix.clamp_factor = True
    mix.clamp_result = False
    mix.data_type = 'RGBA'
    mix.factor_mode = 'UNIFORM'
    #Factor_Float
    if colorbonds:
        mix.inputs[0].default_value = 1
    else:
        mix.inputs[0].default_value = 0


    #Set locations
    material_output.location = (344.0, 300.0)
    attribute_001.location = (-591.1474609375, -175.78738403320312)
    math.location = (-292.2492370605469, -216.01373291015625)
    principled_bsdf.location = (-37.92226028442383, 256.66119384765625)
    attribute.location = (-572.2391967773438, 203.6114959716797)
    rgb.location = (-565.62109375, 389.9770812988281)
    mix.location = (-361.1031799316406, 313.6980285644531)

    #Set dimensions
    material_output.width, material_output.height = 140.0, 100.0
    attribute_001.width, attribute_001.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0
    principled_bsdf.width, principled_bsdf.height = 240.0, 100.0
    attribute.width, attribute.height = 140.0, 100.0
    rgb.width, rgb.height = 140.0, 100.0
    mix.width, mix.height = 140.0, 100.0

    #initialize bond_nodes links
    #principled_bsdf.BSDF -> material_output.Surface
    bond_nodes.links.new(principled_bsdf.outputs[0], material_output.inputs[0])
    #attribute_001.Fac -> math.Value
    bond_nodes.links.new(attribute_001.outputs[2], math.inputs[0])
    #mix.Result -> principled_bsdf.Base Color
    bond_nodes.links.new(mix.outputs[2], principled_bsdf.inputs[0])
    #attribute.Color -> mix.B
    bond_nodes.links.new(attribute.outputs[0], mix.inputs[7])
    #rgb.Color -> mix.A
    bond_nodes.links.new(rgb.outputs[0], mix.inputs[6])
    return bond_nodes


def create_bondmat(colorbonds=False,name="atoms"):
    
    if "BOND_nodes" not in bpy.data.materials:
        mat = bpy.data.materials.new(name = "BOND_nodes_"+name)
    else:
        return
    mat.use_nodes = True
    #initialize BOND_nodes node group
    bond_nodes_node_group(mat,colorbonds=colorbonds)
    return(mat)

