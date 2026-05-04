import bpy
from .. import __version__
def outline_node_group(mat=None):
    """Initialize outline node group"""
    outline = bpy.data.node_groups.new(type='GeometryNodeTree', name="outline")

    outline.color_tag = 'NONE'
    outline.description = __version__
    outline.default_group_node_width = 140
    outline.is_modifier = True

    # outline interface

    # Socket Geometry
    geometry_socket = outline.interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    # Socket Geometry
    geometry_socket_1 = outline.interface.new_socket(name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry')
    geometry_socket_1.attribute_domain = 'POINT'

    # Socket Value
    value_socket = outline.interface.new_socket(name="Value", in_out='INPUT', socket_type='NodeSocketFloat')
    value_socket.default_value = 0.5
    value_socket.min_value = -10000.0
    value_socket.max_value = 10000.0
    value_socket.subtype = 'NONE'
    value_socket.attribute_domain = 'POINT'

    # Socket Outline-Mat
    outline_mat_socket = outline.interface.new_socket(name="Outline-Mat", in_out='INPUT', socket_type='NodeSocketMaterial')
    outline_mat_socket.default_value = mat
    outline_mat_socket.attribute_domain = 'POINT'
    outline_mat_socket.description = "Outline Material"

    # Initialize outline nodes

    # Node Frame.001
    frame_001 = outline.nodes.new("NodeFrame")
    frame_001.name = "Frame.001"
    frame_001.label_size = 20
    frame_001.shrink = True

    # Node Group Output.001
    group_output_001 = outline.nodes.new("NodeGroupOutput")
    group_output_001.name = "Group Output.001"
    group_output_001.is_active_output = True

    # Node Group Input.001
    group_input_001 = outline.nodes.new("NodeGroupInput")
    group_input_001.name = "Group Input.001"

    # Node Join Geometry.001
    join_geometry_001 = outline.nodes.new("GeometryNodeJoinGeometry")
    join_geometry_001.name = "Join Geometry.001"

    # Node Set Shade Smooth
    set_shade_smooth = outline.nodes.new("GeometryNodeSetShadeSmooth")
    set_shade_smooth.name = "Set Shade Smooth"
    set_shade_smooth.domain = 'FACE'
    # Selection
    set_shade_smooth.inputs[1].default_value = True
    # Shade Smooth
    set_shade_smooth.inputs[2].default_value = True

    # Node Set Material.001
    set_material_001 = outline.nodes.new("GeometryNodeSetMaterial")
    set_material_001.name = "Set Material.001"
    # Selection
    set_material_001.inputs[1].default_value = True

    # Node Vector Math
    vector_math = outline.nodes.new("ShaderNodeVectorMath")
    vector_math.name = "Vector Math"
    vector_math.operation = 'SCALE'

    # Node Set Position
    set_position = outline.nodes.new("GeometryNodeSetPosition")
    set_position.name = "Set Position"
    # Selection
    set_position.inputs[1].default_value = True
    # Position
    set_position.inputs[2].default_value = (0.0, 0.0, 0.0)

    # Node Normal
    normal = outline.nodes.new("GeometryNodeInputNormal")
    normal.name = "Normal"
    if hasattr(normal, 'legacy_corner_normals'):
        normal.legacy_corner_normals = True

    # Node Group Input
    group_input = outline.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    # Node Math
    math = outline.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    # Value_001
    math.inputs[1].default_value = 0.05000000074505806

    # Set parents
    group_input_001.parent = frame_001
    set_shade_smooth.parent = frame_001

    # Set locations
    frame_001.location = (-325.0, -1055.0)
    group_output_001.location = (1067.1932373046875, -1736.042724609375)
    group_input_001.location = (30.27032470703125, -30.378662109375)
    join_geometry_001.location = (880.1171264648438, -1655.3466796875)
    set_shade_smooth.location = (632.2693481445312, -153.72265625)
    set_material_001.location = (554.8684692382812, -1089.7171630859375)
    vector_math.location = (96.314208984375, -812.8528442382812)
    set_position.location = (278.2383117675781, -792.7901611328125)
    normal.location = (-317.62713623046875, -734.376220703125)
    group_input.location = (-532.8179321289062, -956.6204223632812)
    math.location = (-299.7628479003906, -840.5046997070312)

    # Set dimensions
    frame_001.width, frame_001.height = 802.0, 327.0
    group_output_001.width, group_output_001.height = 140.0, 100.0
    group_input_001.width, group_input_001.height = 140.0, 100.0
    join_geometry_001.width, join_geometry_001.height = 140.0, 100.0
    set_shade_smooth.width, set_shade_smooth.height = 140.0, 100.0
    set_material_001.width, set_material_001.height = 140.0, 100.0
    vector_math.width, vector_math.height = 140.0, 100.0
    set_position.width, set_position.height = 140.0, 100.0
    normal.width, normal.height = 140.0, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0

    # Initialize outline links

    # group_input_001.Geometry -> set_shade_smooth.Geometry
    outline.links.new(group_input_001.outputs[0], set_shade_smooth.inputs[0])
    # vector_math.Vector -> set_position.Offset
    outline.links.new(vector_math.outputs[0], set_position.inputs[3])
    # set_shade_smooth.Geometry -> set_position.Geometry
    outline.links.new(set_shade_smooth.outputs[0], set_position.inputs[0])
    # normal.Normal -> vector_math.Vector
    outline.links.new(normal.outputs[0], vector_math.inputs[0])
    # join_geometry_001.Geometry -> group_output_001.Geometry
    outline.links.new(join_geometry_001.outputs[0], group_output_001.inputs[0])
    # set_material_001.Geometry -> join_geometry_001.Geometry
    outline.links.new(set_material_001.outputs[0], join_geometry_001.inputs[0])
    # set_position.Geometry -> set_material_001.Geometry
    outline.links.new(set_position.outputs[0], set_material_001.inputs[0])
    # math.Value -> vector_math.Scale
    outline.links.new(math.outputs[0], vector_math.inputs[3])
    # group_input.Value -> math.Value
    outline.links.new(group_input.outputs[1], math.inputs[0])
    # group_input.Outline-Mat -> set_material_001.Material
    outline.links.new(group_input.outputs[2], set_material_001.inputs[2])
    # set_shade_smooth.Geometry -> join_geometry_001.Geometry
    outline.links.new(set_shade_smooth.outputs[0], join_geometry_001.inputs[0])

    return outline



def outline_color_node_group(mat):

    outline_color = mat.node_tree
    print(outline_color)
    #start with a clean node tree
    for node in outline_color.nodes:
        outline_color.nodes.remove(node)
    outline_color.color_tag = 'NONE'
    outline_color.description = ""
    outline_color.default_group_node_width = 140
    

    #outline_color interface

    #initialize outline_color nodes
    #node Math
    math = outline_color.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False

    #node Mix Shader
    mix_shader = outline_color.nodes.new("ShaderNodeMixShader")
    mix_shader.name = "Mix Shader"

    #node Emission
    emission = outline_color.nodes.new("ShaderNodeEmission")
    emission.name = "Emission"
    #Color
    emission.inputs[0].default_value = (0.0, 0.0, 0.0, 1.0)
    #Strength
    emission.inputs[1].default_value = 1.0

    #node Transparent BSDF
    transparent_bsdf = outline_color.nodes.new("ShaderNodeBsdfTransparent")
    transparent_bsdf.name = "Transparent BSDF"
    #Color
    transparent_bsdf.inputs[0].default_value = (1.0, 1.0, 1.0, 1.0)

    #node Geometry
    geometry = outline_color.nodes.new("ShaderNodeNewGeometry")
    geometry.name = "Geometry"

    #node Light Path
    light_path = outline_color.nodes.new("ShaderNodeLightPath")
    light_path.name = "Light Path"

    #node Material Output
    material_output = outline_color.nodes.new("ShaderNodeOutputMaterial")
    material_output.name = "Material Output"
    material_output.is_active_output = True
    material_output.target = 'ALL'
    #Displacement
    material_output.inputs[2].default_value = (0.0, 0.0, 0.0)

    #node Value
    value = outline_color.nodes.new("ShaderNodeValue")
    value.name = "Value"

    value.outputs[0].default_value = 19.999998092651367

    #Set locations
    math.location = (-116.0, 283.03466796875)
    mix_shader.location = (104.21751403808594, 118.79178619384766)
    emission.location = (-407.0325622558594, 73.0439453125)
    transparent_bsdf.location = (-356.71563720703125, 186.7467498779297)
    geometry.location = (-586.4742431640625, 464.6316223144531)
    light_path.location = (-408.0494384765625, 627.9410400390625)
    material_output.location = (812.8287353515625, 152.38636779785156)
    value.location = (812.8287353515625, -7.6136322021484375)

    #Set dimensions
    math.width, math.height = 140.0, 100.0
    mix_shader.width, mix_shader.height = 140.0, 100.0
    emission.width, emission.height = 140.0, 100.0
    transparent_bsdf.width, transparent_bsdf.height = 140.0, 100.0
    geometry.width, geometry.height = 140.0, 100.0
    light_path.width, light_path.height = 140.0, 100.0
    material_output.width, material_output.height = 140.0, 100.0
    value.width, value.height = 140.0, 100.0

    #initialize outline_color links
    #math.Value -> mix_shader.Fac
    outline_color.links.new(math.outputs[0], mix_shader.inputs[0])
    #geometry.Backfacing -> math.Value
    outline_color.links.new(geometry.outputs[6], math.inputs[0])
    #light_path.Is Camera Ray -> math.Value
    outline_color.links.new(light_path.outputs[0], math.inputs[1])
    #emission.Emission -> mix_shader.Shader
    outline_color.links.new(emission.outputs[0], mix_shader.inputs[2])
    #transparent_bsdf.BSDF -> mix_shader.Shader
    outline_color.links.new(transparent_bsdf.outputs[0], mix_shader.inputs[1])
    #value.Value -> material_output.Thickness
    outline_color.links.new(value.outputs[0], material_output.inputs[3])
    #mix_shader.Shader -> material_output.Surface
    outline_color.links.new(mix_shader.outputs[0], material_output.inputs[0])
    return outline_color



def outline_objects(list_of_objects,modifier='GeometryNodes'):
    if 'outline_color' not in bpy.data.materials:
        mat = bpy.data.materials.new(name = "outline_color")
    else:
        mat = bpy.data.materials["outline_color"]
    mat.use_nodes = True
    #outline_color =   # ruff unhappy
    outline_color_node_group(mat)
    
    if 'outline' in bpy.data.node_groups:
        node = bpy.data.node_groups["outline"]
        # check if it is the current version otherwise rename and create new
        desc = node.description
        if desc != __version__:
            node.name = "outline_old"
            node = outline_node_group(mat=mat)
    else:
        node = outline_node_group(mat=mat)
    # Always ensure the material socket is set correctly
    node.interface.items_tree["Outline-Mat"].default_value = mat

    node["Socket_1"] = 1
    node["Socket_2"] = False
    bpy.ops.object.select_all(action='DESELECT')
    
    for obj in list_of_objects:
        if obj == list_of_objects[0]:
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.modifier_add(type='NODES')
            #node = bpy.data.node_groups["outline"]
            bpy.context.object.modifiers[modifier].node_group = node

        obj.select_set(True)
        if not obj.data.materials:
            obj.data.materials[0] = mat
        else:
            obj.data.materials.append(mat)
        
    if len(list_of_objects) > 1:
        bpy.ops.object.make_links_data(type='MODIFIERS')

        
