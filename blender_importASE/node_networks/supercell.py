import bpy
from ase.data import chemical_symbols

#initialize cutoff_group node group
def cutoff_group_node_group():

    cutoff_group = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "cutoff_group")

    cutoff_group.color_tag = 'NONE'
    cutoff_group.description = ""
    cutoff_group.default_group_node_width = 140
    


    #cutoff_group interface
    #Socket Geometry
    geometry_socket = cutoff_group.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_1 = cutoff_group.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_1.attribute_domain = 'POINT'

    #Socket +x
    _x_socket = cutoff_group.interface.new_socket(name = "+x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket.default_value = 16.0
    _x_socket.min_value = -10000.0
    _x_socket.max_value = 10000.0
    _x_socket.subtype = 'NONE'
    _x_socket.attribute_domain = 'POINT'

    #Socket +y
    _y_socket = cutoff_group.interface.new_socket(name = "+y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket.default_value = 16.0
    _y_socket.min_value = -10000.0
    _y_socket.max_value = 10000.0
    _y_socket.subtype = 'NONE'
    _y_socket.attribute_domain = 'POINT'

    #Socket +z
    _z_socket = cutoff_group.interface.new_socket(name = "+z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket.default_value = 16.0
    _z_socket.min_value = -10000.0
    _z_socket.max_value = 10000.0
    _z_socket.subtype = 'NONE'
    _z_socket.attribute_domain = 'POINT'

    #Socket -x
    _x_socket_1 = cutoff_group.interface.new_socket(name = "-x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket_1.default_value = -16.0
    _x_socket_1.min_value = -10000.0
    _x_socket_1.max_value = 10000.0
    _x_socket_1.subtype = 'NONE'
    _x_socket_1.attribute_domain = 'POINT'

    #Socket -y
    _y_socket_1 = cutoff_group.interface.new_socket(name = "-y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket_1.default_value = -16.0
    _y_socket_1.min_value = -10000.0
    _y_socket_1.max_value = 10000.0
    _y_socket_1.subtype = 'NONE'
    _y_socket_1.attribute_domain = 'POINT'

    #Socket -z
    _z_socket_1 = cutoff_group.interface.new_socket(name = "-z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket_1.default_value = -16.0
    _z_socket_1.min_value = -10000.0
    _z_socket_1.max_value = 10000.0
    _z_socket_1.subtype = 'NONE'
    _z_socket_1.attribute_domain = 'POINT'


    #initialize cutoff_group nodes
    #node Group Output
    group_output = cutoff_group.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"
    group_output.is_active_output = True

    #node Compare.002
    compare_002 = cutoff_group.nodes.new("FunctionNodeCompare")
    compare_002.name = "Compare.002"
    compare_002.data_type = 'FLOAT'
    compare_002.mode = 'ELEMENT'
    compare_002.operation = 'GREATER_THAN'

    #node Delete Geometry.001
    delete_geometry_001 = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_001.name = "Delete Geometry.001"
    delete_geometry_001.domain = 'POINT'
    delete_geometry_001.mode = 'ALL'

    #node Delete Geometry
    delete_geometry = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.name = "Delete Geometry"
    delete_geometry.domain = 'POINT'
    delete_geometry.mode = 'ALL'

    #node Compare
    compare = cutoff_group.nodes.new("FunctionNodeCompare")
    compare.name = "Compare"
    compare.data_type = 'FLOAT'
    compare.mode = 'ELEMENT'
    compare.operation = 'GREATER_THAN'

    #node Compare.001
    compare_001 = cutoff_group.nodes.new("FunctionNodeCompare")
    compare_001.name = "Compare.001"
    compare_001.data_type = 'FLOAT'
    compare_001.mode = 'ELEMENT'
    compare_001.operation = 'GREATER_THAN'

    #node Compare.003
    compare_003 = cutoff_group.nodes.new("FunctionNodeCompare")
    compare_003.name = "Compare.003"
    compare_003.data_type = 'FLOAT'
    compare_003.mode = 'ELEMENT'
    compare_003.operation = 'LESS_THAN'

    #node Delete Geometry.004
    delete_geometry_004 = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_004.name = "Delete Geometry.004"
    delete_geometry_004.domain = 'POINT'
    delete_geometry_004.mode = 'ALL'

    #node Delete Geometry.005
    delete_geometry_005 = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_005.name = "Delete Geometry.005"
    delete_geometry_005.domain = 'POINT'
    delete_geometry_005.mode = 'ALL'

    #node Compare.004
    compare_004 = cutoff_group.nodes.new("FunctionNodeCompare")
    compare_004.name = "Compare.004"
    compare_004.data_type = 'FLOAT'
    compare_004.mode = 'ELEMENT'
    compare_004.operation = 'LESS_THAN'

    #node Compare.005
    compare_005 = cutoff_group.nodes.new("FunctionNodeCompare")
    compare_005.name = "Compare.005"
    compare_005.data_type = 'FLOAT'
    compare_005.mode = 'ELEMENT'
    compare_005.operation = 'LESS_THAN'

    #node Delete Geometry.002
    delete_geometry_002 = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_002.name = "Delete Geometry.002"
    delete_geometry_002.domain = 'POINT'
    delete_geometry_002.mode = 'ALL'

    #node Delete Geometry.003
    delete_geometry_003 = cutoff_group.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_003.name = "Delete Geometry.003"
    delete_geometry_003.domain = 'POINT'
    delete_geometry_003.mode = 'ALL'

    #node Group Input
    group_input = cutoff_group.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    #node Separate XYZ
    separate_xyz = cutoff_group.nodes.new("ShaderNodeSeparateXYZ")
    separate_xyz.name = "Separate XYZ"

    #node Position
    position = cutoff_group.nodes.new("GeometryNodeInputPosition")
    position.name = "Position"





    #Set locations
    group_output.location = (872.9868774414062, 0.0)
    compare_002.location = (97.686279296875, 64.63043212890625)
    delete_geometry_001.location = (680.3541870117188, 271.48876953125)
    delete_geometry.location = (682.9868774414062, 430.9881591796875)
    compare.location = (108.29290771484375, 386.19647216796875)
    compare_001.location = (103.46710205078125, 225.82261657714844)
    compare_003.location = (88.3687744140625, -430.9881591796875)
    delete_geometry_004.location = (671.03662109375, -224.12982177734375)
    delete_geometry_005.location = (673.6693115234375, -64.63041687011719)
    compare_004.location = (98.97540283203125, -109.42210388183594)
    compare_005.location = (94.14959716796875, -269.79595947265625)
    delete_geometry_002.location = (673.7863159179688, 107.57318115234375)
    delete_geometry_003.location = (664.46875, -388.04541015625)
    group_input.location = (-882.9868774414062, 0.0)
    separate_xyz.location = (-112.14068603515625, 334.73724365234375)
    position.location = (-486.3282775878906, 285.5691833496094)

    #Set dimensions
    group_output.width, group_output.height = 140.0, 100.0
    compare_002.width, compare_002.height = 140.0, 100.0
    delete_geometry_001.width, delete_geometry_001.height = 140.0, 100.0
    delete_geometry.width, delete_geometry.height = 140.0, 100.0
    compare.width, compare.height = 140.0, 100.0
    compare_001.width, compare_001.height = 140.0, 100.0397
    compare_003.width, compare_003.height = 140.0, 100.0
    delete_geometry_004.width, delete_geometry_004.height = 140.0, 100.0
    delete_geometry_005.width, delete_geometry_005.height = 140.0, 100.0
    compare_004.width, compare_004.height = 140.0, 100.0
    compare_005.width, compare_005.height = 140.0, 100.0
    delete_geometry_002.width, delete_geometry_002.height = 140.0, 100.0
    delete_geometry_003.width, delete_geometry_003.height = 153.6519775390625, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    separate_xyz.width, separate_xyz.height = 140.0, 100.0
    position.width, position.height = 140.0, 100.0

    #initialize cutoff_group links
    #compare_002.Result -> delete_geometry_002.Selection
    cutoff_group.links.new(compare_002.outputs[0], delete_geometry_002.inputs[1])
    #delete_geometry_002.Geometry -> delete_geometry_005.Geometry
    cutoff_group.links.new(delete_geometry_002.outputs[0], delete_geometry_005.inputs[0])
    #separate_xyz.Y -> compare_001.A
    cutoff_group.links.new(separate_xyz.outputs[1], compare_001.inputs[0])
    #separate_xyz.X -> compare.A
    cutoff_group.links.new(separate_xyz.outputs[0], compare.inputs[0])
    #compare_003.Result -> delete_geometry_003.Selection
    cutoff_group.links.new(compare_003.outputs[0], delete_geometry_003.inputs[1])
    #compare.Result -> delete_geometry.Selection
    cutoff_group.links.new(compare.outputs[0], delete_geometry.inputs[1])
    #compare_005.Result -> delete_geometry_004.Selection
    cutoff_group.links.new(compare_005.outputs[0], delete_geometry_004.inputs[1])
    #delete_geometry_004.Geometry -> delete_geometry_003.Geometry
    cutoff_group.links.new(delete_geometry_004.outputs[0], delete_geometry_003.inputs[0])
    #separate_xyz.Z -> compare_002.A
    cutoff_group.links.new(separate_xyz.outputs[2], compare_002.inputs[0])
    #delete_geometry_005.Geometry -> delete_geometry_004.Geometry
    cutoff_group.links.new(delete_geometry_005.outputs[0], delete_geometry_004.inputs[0])
    #compare_004.Result -> delete_geometry_005.Selection
    cutoff_group.links.new(compare_004.outputs[0], delete_geometry_005.inputs[1])
    #compare_001.Result -> delete_geometry_001.Selection
    cutoff_group.links.new(compare_001.outputs[0], delete_geometry_001.inputs[1])
    #delete_geometry.Geometry -> delete_geometry_001.Geometry
    cutoff_group.links.new(delete_geometry.outputs[0], delete_geometry_001.inputs[0])
    #delete_geometry_001.Geometry -> delete_geometry_002.Geometry
    cutoff_group.links.new(delete_geometry_001.outputs[0], delete_geometry_002.inputs[0])
    #group_input.Geometry -> delete_geometry.Geometry
    cutoff_group.links.new(group_input.outputs[0], delete_geometry.inputs[0])
    #delete_geometry_003.Geometry -> group_output.Geometry
    cutoff_group.links.new(delete_geometry_003.outputs[0], group_output.inputs[0])
    #group_input.+x -> compare.B
    cutoff_group.links.new(group_input.outputs[1], compare.inputs[1])
    #group_input.+y -> compare_001.B
    cutoff_group.links.new(group_input.outputs[2], compare_001.inputs[1])
    #group_input.+z -> compare_002.B
    cutoff_group.links.new(group_input.outputs[3], compare_002.inputs[1])
    #group_input.-x -> compare_004.B
    cutoff_group.links.new(group_input.outputs[4], compare_004.inputs[1])
    #group_input.-y -> compare_005.B
    cutoff_group.links.new(group_input.outputs[5], compare_005.inputs[1])
    #group_input.-z -> compare_003.B
    cutoff_group.links.new(group_input.outputs[6], compare_003.inputs[1])
    #separate_xyz.X -> compare_004.A
    cutoff_group.links.new(separate_xyz.outputs[0], compare_004.inputs[0])
    #separate_xyz.Y -> compare_005.A
    cutoff_group.links.new(separate_xyz.outputs[1], compare_005.inputs[0])
    #separate_xyz.Z -> compare_003.A
    cutoff_group.links.new(separate_xyz.outputs[2], compare_003.inputs[0])
    #position.Position -> separate_xyz.Vector
    cutoff_group.links.new(position.outputs[0], separate_xyz.inputs[0])
    return cutoff_group
#initialize supercell node group
def supercell_node_group(atoms):
    cutoff_group=cutoff_group_node_group()
    supercell = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "supercell")

    supercell.color_tag = 'NONE'
    supercell.description = ""
    supercell.default_group_node_width = 140
    

    supercell.is_modifier = True

    #supercell interface
    #Socket Geometry
    geometry_socket_2 = supercell.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_2.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_3 = supercell.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_3.attribute_domain = 'POINT'

    #Socket repeat_x
    repeat_x_socket = supercell.interface.new_socket(name = "repeat_x", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_x_socket.default_value = 1
    repeat_x_socket.min_value = 1
    repeat_x_socket.max_value = 10000
    repeat_x_socket.subtype = 'NONE'
    repeat_x_socket.attribute_domain = 'POINT'

    #Socket repeat_y
    repeat_y_socket = supercell.interface.new_socket(name = "repeat_y", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_y_socket.default_value = 1
    repeat_y_socket.min_value = 1
    repeat_y_socket.max_value = 10000
    repeat_y_socket.subtype = 'NONE'
    repeat_y_socket.attribute_domain = 'POINT'

    #Socket repeat_z
    repeat_z_socket = supercell.interface.new_socket(name = "repeat_z", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_z_socket.default_value = 1
    repeat_z_socket.min_value = 1
    repeat_z_socket.max_value = 10000
    repeat_z_socket.subtype = 'NONE'
    repeat_z_socket.attribute_domain = 'POINT'

    #Socket Offset_x
    offset_x_socket = supercell.interface.new_socket(name = "Offset_x", in_out='INPUT', socket_type = 'NodeSocketInt')
    offset_x_socket.default_value = 0
    offset_x_socket.min_value = -2147483648
    offset_x_socket.max_value = 2147483647
    offset_x_socket.subtype = 'NONE'
    offset_x_socket.attribute_domain = 'POINT'

    #Socket Offset_y
    offset_y_socket = supercell.interface.new_socket(name = "Offset_y", in_out='INPUT', socket_type = 'NodeSocketInt')
    offset_y_socket.default_value = 0
    offset_y_socket.min_value = -2147483648
    offset_y_socket.max_value = 2147483647
    offset_y_socket.subtype = 'NONE'
    offset_y_socket.attribute_domain = 'POINT'

    #Socket global
    global_socket = supercell.interface.new_socket(name = "global", in_out='INPUT', socket_type = 'NodeSocketBool')
    global_socket.default_value = False
    global_socket.attribute_domain = 'POINT'

    #Socket cutoff_vectors
    cutoff_vectors_socket = supercell.interface.new_socket(name = "cutoff_vectors", in_out='INPUT', socket_type = 'NodeSocketBool')
    cutoff_vectors_socket.default_value = False
    cutoff_vectors_socket.attribute_domain = 'POINT'

    

    #Socket +x
    _x_socket_2 = supercell.interface.new_socket(name = "+x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket_2.default_value = 100.0
    _x_socket_2.min_value = -10000.0
    _x_socket_2.max_value = 10000.0
    _x_socket_2.subtype = 'NONE'
    _x_socket_2.attribute_domain = 'POINT'

    #Socket +y
    _y_socket_2 = supercell.interface.new_socket(name = "+y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket_2.default_value = 100.0
    _y_socket_2.min_value = -10000.0
    _y_socket_2.max_value = 10000.0
    _y_socket_2.subtype = 'NONE'
    _y_socket_2.attribute_domain = 'POINT'

    #Socket +z
    _z_socket_2 = supercell.interface.new_socket(name = "+z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket_2.default_value = 100.0
    _z_socket_2.min_value = -10000.0
    _z_socket_2.max_value = 10000.0
    _z_socket_2.subtype = 'NONE'
    _z_socket_2.attribute_domain = 'POINT'

    #Socket -x
    _x_socket_3 = supercell.interface.new_socket(name = "-x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket_3.default_value = -100.0
    _x_socket_3.min_value = -10000.0
    _x_socket_3.max_value = 10000.0
    _x_socket_3.subtype = 'NONE'
    _x_socket_3.attribute_domain = 'POINT'

    #Socket -y
    _y_socket_3 = supercell.interface.new_socket(name = "-y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket_3.default_value = -100.0
    _y_socket_3.min_value = -10000.0
    _y_socket_3.max_value = 10000.0
    _y_socket_3.subtype = 'NONE'
    _y_socket_3.attribute_domain = 'POINT'

    #Socket -z
    _z_socket_3 = supercell.interface.new_socket(name = "-z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket_3.default_value = -100.0
    _z_socket_3.min_value = -10000.0
    _z_socket_3.max_value = 10000.0
    _z_socket_3.subtype = 'NONE'
    _z_socket_3.attribute_domain = 'POINT'


    #initialize supercell nodes
    #node Frame.002
    frame_002 = supercell.nodes.new("NodeFrame")
    frame_002.label = "cutoff"
    frame_002.name = "Frame.002"
    frame_002.label_size = 20
    frame_002.shrink = True

    #node Frame.001
    frame_001 = supercell.nodes.new("NodeFrame")
    frame_001.label = "grid"
    frame_001.name = "Frame.001"
    frame_001.label_size = 20
    frame_001.shrink = True

    #node Frame
    frame = supercell.nodes.new("NodeFrame")
    frame.label = "offset"
    frame.name = "Frame"
    frame.label_size = 20
    frame.shrink = True

    #node Instance on Points
    instance_on_points = supercell.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points.name = "Instance on Points"
    #Selection
    instance_on_points.inputs[1].default_value = True
    #Pick Instance
    instance_on_points.inputs[3].default_value = False
    #Instance Index
    instance_on_points.inputs[4].default_value = 0
    #Rotation
    instance_on_points.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Group Output
    group_output_1 = supercell.nodes.new("NodeGroupOutput")
    group_output_1.name = "Group Output"
    group_output_1.is_active_output = True

    #node Group Input
    group_input_1 = supercell.nodes.new("NodeGroupInput")
    group_input_1.name = "Group Input"

    #node Realize Instances
    realize_instances_beforevectorcutoff = supercell.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_beforevectorcutoff.name = "Realize Instances"
    #Selection
    realize_instances_beforevectorcutoff.inputs[1].default_value = True
    #Realize All
    realize_instances_beforevectorcutoff.inputs[2].default_value = True
    #Depth
    realize_instances_beforevectorcutoff.inputs[3].default_value = 0

    #node Group.001
    cutoff_group_check = supercell.nodes.new("GeometryNodeGroup")
    cutoff_group_check.name = "Group.001"
    cutoff_group_check.node_tree = cutoff_group

    #node Mesh Line.002
    mesh_line_002 = supercell.nodes.new("GeometryNodeMeshLine")
    mesh_line_002.name = "Mesh Line.002"
    mesh_line_002.count_mode = 'TOTAL'
    mesh_line_002.mode = 'OFFSET'

    #node Mesh Line
    mesh_line = supercell.nodes.new("GeometryNodeMeshLine")
    mesh_line.name = "Mesh Line"
    mesh_line.hide = True
    mesh_line.count_mode = 'TOTAL'
    mesh_line.mode = 'OFFSET'

    #node Mesh Line.001
    mesh_line_001 = supercell.nodes.new("GeometryNodeMeshLine")
    mesh_line_001.name = "Mesh Line.001"
    mesh_line_001.hide = True
    mesh_line_001.count_mode = 'TOTAL'
    mesh_line_001.mode = 'OFFSET'

    #node Instance on Points.001
    instance_on_points_001 = supercell.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_001.name = "Instance on Points.001"
    instance_on_points_001.hide = True
    #Selection
    instance_on_points_001.inputs[1].default_value = True
    #Pick Instance
    instance_on_points_001.inputs[3].default_value = False
    #Instance Index
    instance_on_points_001.inputs[4].default_value = 0
    #Rotation
    instance_on_points_001.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points_001.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Instance on Points.002
    instance_on_points_002 = supercell.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_002.name = "Instance on Points.002"
    instance_on_points_002.hide = True
    #Selection
    instance_on_points_002.inputs[1].default_value = True
    #Pick Instance
    instance_on_points_002.inputs[3].default_value = False
    #Instance Index
    instance_on_points_002.inputs[4].default_value = 0
    #Rotation
    instance_on_points_002.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points_002.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Vector Math.001
    vector_math_001 = supercell.nodes.new("ShaderNodeVectorMath")
    vector_math_001.name = "Vector Math.001"
    vector_math_001.hide = True
    vector_math_001.operation = 'SCALE'

    #node Vector Math
    vector_math = supercell.nodes.new("ShaderNodeVectorMath")
    vector_math.name = "Vector Math"
    vector_math.hide = True
    vector_math.operation = 'SCALE'

    #node Vector Math.002
    vector_math_002 = supercell.nodes.new("ShaderNodeVectorMath")
    vector_math_002.name = "Vector Math.002"
    vector_math_002.hide = True
    vector_math_002.operation = 'SCALE'

    #node Math.002
    math_002 = supercell.nodes.new("ShaderNodeMath")
    math_002.name = "Math.002"
    math_002.hide = True
    math_002.operation = 'MULTIPLY'
    math_002.use_clamp = False
    #Value_001
    math_002.inputs[1].default_value = -1.0

    #node Math
    math = supercell.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.hide = True
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = -1.0

    #node Math.001
    math_001 = supercell.nodes.new("ShaderNodeMath")
    math_001.name = "Math.001"
    math_001.hide = True
    math_001.operation = 'MULTIPLY'
    math_001.use_clamp = False
    #Value_001
    math_001.inputs[1].default_value = -1.0

    #node Reroute
    reroute = supercell.nodes.new("NodeReroute")
    reroute.name = "Reroute"
    reroute.socket_idname = "NodeSocketInt"
    #node Reroute.001
    reroute_001 = supercell.nodes.new("NodeReroute")
    reroute_001.name = "Reroute.001"
    reroute_001.socket_idname = "NodeSocketInt"
    #node Group Input.001
    group_input_001 = supercell.nodes.new("NodeGroupInput")
    group_input_001.name = "Group Input.001"

    #node Reroute.004
    reroute_004 = supercell.nodes.new("NodeReroute")
    reroute_004.name = "Reroute.004"
    reroute_004.socket_idname = "NodeSocketInt"
    #node Reroute.002
    reroute_002 = supercell.nodes.new("NodeReroute")
    reroute_002.name = "Reroute.002"
    reroute_002.socket_idname = "NodeSocketInt"
    #node Reroute.003
    reroute_003 = supercell.nodes.new("NodeReroute")
    reroute_003.name = "Reroute.003"
    reroute_003.socket_idname = "NodeSocketInt"
    #node Reroute.007
    reroute_007 = supercell.nodes.new("NodeReroute")
    reroute_007.name = "Reroute.007"
    reroute_007.socket_idname = "NodeSocketVector"
    #node Reroute.005
    reroute_005 = supercell.nodes.new("NodeReroute")
    reroute_005.name = "Reroute.005"
    reroute_005.socket_idname = "NodeSocketVector"
    #node Reroute.006
    reroute_006 = supercell.nodes.new("NodeReroute")
    reroute_006.name = "Reroute.006"
    reroute_006.socket_idname = "NodeSocketVector"
    #node Vector.002
    vector_z = supercell.nodes.new("FunctionNodeInputVector")
    vector_z.label = "zvec"
    vector_z.name = "Vector.002"
    vector_z.hide = True
    vector_z.vector = atoms.cell[2]

    #node Vector
    vector_x = supercell.nodes.new("FunctionNodeInputVector")
    vector_x.label = "xvec"
    vector_x.name = "Vector"
    vector_x.hide = True
    vector_x.vector = atoms.cell[0]

    #node Vector.001
    vector_y = supercell.nodes.new("FunctionNodeInputVector")
    vector_y.label = "yvec"
    vector_y.hide = True
    vector_y.name = "Vector.001"
    vector_y.vector = atoms.cell[1]

    #node Switch
    switch = supercell.nodes.new("GeometryNodeSwitch")
    switch.label = "repeat_X"
    switch.name = "Switch"
    switch.input_type = 'INT'
    #True
    switch.inputs[2].default_value = 1

    #node Switch.001
    switch_001 = supercell.nodes.new("GeometryNodeSwitch")
    switch_001.label = "repeat_Y"
    switch_001.name = "Switch.001"
    switch_001.input_type = 'INT'
    #True
    switch_001.inputs[2].default_value = 1

    #node Switch.002
    switch_002 = supercell.nodes.new("GeometryNodeSwitch")
    switch_002.label = "repeat_Z"
    switch_002.name = "Switch.002"
    switch_002.input_type = 'INT'
    #True
    switch_002.inputs[2].default_value = 1

    #node Switch.003
    switch_003 = supercell.nodes.new("GeometryNodeSwitch")
    switch_003.label = "offset_X"
    switch_003.name = "Switch.003"
    switch_003.input_type = 'INT'
    #True
    switch_003.inputs[2].default_value = 0

    #node Switch.004
    switch_004 = supercell.nodes.new("GeometryNodeSwitch")
    switch_004.label = "offset_Y"
    switch_004.name = "Switch.004"
    switch_004.input_type = 'INT'
    #True
    switch_004.inputs[2].default_value = 1

    #node Switch.005
    switch_005 = supercell.nodes.new("GeometryNodeSwitch")
    switch_005.label = "offset_Z"
    switch_005.name = "Switch.005"
    switch_005.input_type = 'INT'
    #True
    switch_005.inputs[2].default_value = 0

    #node Group Input.002
    group_input_002 = supercell.nodes.new("NodeGroupInput")
    group_input_002.name = "Group Input.002"

    #node Frame.003
    frame_003 = supercell.nodes.new("NodeFrame")
    frame_003.label = "modify global value for cell repetition"
    frame_003.name = "Frame.003"
    frame_003.use_custom_color = True
    frame_003.color = (0.41161254048347473, 0.0, 0.6079999804496765)
    frame_003.label_size = 20
    frame_003.shrink = True

    #node Frame.004
    frame_004 = supercell.nodes.new("NodeFrame")
    frame_004.label = "global value for offset"
    frame_004.name = "Frame.004"
    frame_004.use_custom_color = True
    frame_004.color = (0.3529382050037384, 0.0235294159501791, 0.5098000168800354)
    frame_004.label_size = 20
    frame_004.shrink = True

    #node Frame.005
    frame_005 = supercell.nodes.new("NodeFrame")
    frame_005.name = "Frame.005"
    frame_005.use_custom_color = True
    frame_005.color = (0.08440613746643066, 0.6079999804496765, 0.0)
    frame_005.label_size = 20
    frame_005.shrink = True

    #node Realize Instances.001
    realize_instances_001 = supercell.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_001.name = "Realize Instances.001"
    #Selection
    realize_instances_001.inputs[1].default_value = True
    #Realize All
    realize_instances_001.inputs[2].default_value = True
    #Depth
    realize_instances_001.inputs[3].default_value = 0

    

    #node Join Geometry
    join_geometry = supercell.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"

    #node Merge by Distance
    merge_by_distance = supercell.nodes.new("GeometryNodeMergeByDistance")
    merge_by_distance.name = "Merge by Distance"
    merge_by_distance.mode = 'ALL'
    #Selection
    merge_by_distance.inputs[1].default_value = True
    #Distance
    merge_by_distance.inputs[2].default_value = 0.0010000000474974513

    #node Reroute.008
    reroute_element = supercell.nodes.new("NodeReroute")
    reroute_element.name = "Reroute.008"
    reroute_element.socket_idname = "NodeSocketGeometry"
    #node Switch.006
    switch_vectorcutoff = supercell.nodes.new("GeometryNodeSwitch")
    switch_vectorcutoff.name = "Switch.006"
    switch_vectorcutoff.input_type = 'GEOMETRY'
    #Switch
    switch_vectorcutoff.inputs[0].default_value = False

    #node Group Input.004
    group_input_vectorcutoff = supercell.nodes.new("NodeGroupInput")
    group_input_vectorcutoff.name = "Group Input.004"

    
    #node Named Attribute
    element_attribute = supercell.nodes.new("GeometryNodeInputNamedAttribute")
    element_attribute.name = "Named Attribute"
    element_attribute.data_type = 'FLOAT'
    #Name
    element_attribute.inputs[0].default_value = "element"

    
    #initialize supercell_001 links
    #compare.Result -> delete_geometry.Selection

    # named attribute into compare node
    # compare into delete
    # delete into switch[2]
    # input_el into switch[0]
    # reroute_008 into switch[1]

    #after n == 1
    #old_delete into delete_geometry
    #old_switch into switch[1]
    # if n len
    #switch into realize istnance_001?

    group_input_element=supercell.nodes.new("NodeGroupInput")
    group_input_element.name = "Group Input Element"
    group_input_element.label = "Element"
    old_switch = reroute_element
    for n,number in enumerate(set(atoms.get_atomic_numbers())):
        sym=chemical_symbols[number]
        #Socket cutoff_H
        cutoff_element_socket = supercell.interface.new_socket(name = f"cutoff_{sym}", in_out='INPUT', socket_type = 'NodeSocketBool')
        cutoff_element_socket.default_value = False
        cutoff_element_socket.attribute_domain = 'POINT'
        

        is_element = supercell.nodes.new("FunctionNodeCompare")
        is_element.label = f"is_element_{sym}"
        is_element.name = "Compare"
        is_element.data_type = 'INT'
        is_element.mode = 'ELEMENT'
        is_element.operation = 'EQUAL'
        is_element.inputs[3].default_value = number


        

        #node Delete Geometry
        delete_geometry_element = supercell.nodes.new("GeometryNodeDeleteGeometry")
        delete_geometry_element.name = f"Delete Geometry_{sym}"
        delete_geometry_element.domain = 'POINT'
        delete_geometry_element.mode = 'ALL'
        
        
        #node switch
        switch_element=supercell.nodes.new("GeometryNodeSwitch")
        switch_element.name = f"Switch_{sym}"
        switch_element.label = f"{sym}"
        switch_element.input_type = 'GEOMETRY'
        switch_element.inputs[0].default_value = False
        
        supercell.links.new(is_element.outputs[0], delete_geometry_element.inputs[1])
        supercell.links.new(element_attribute.outputs[0], is_element.inputs[2])
        
        supercell.links.new(group_input_element.outputs[14+n], switch_element.inputs[0])
        supercell.links.new(delete_geometry_element.outputs[0], switch_element.inputs[2])
        supercell.links.new(old_switch.outputs[0],switch_element.inputs[1])
        supercell.links.new(old_switch.outputs[0],delete_geometry_element.inputs[0])

        #supercell.links.new(old_delete.outputs[0], switch_element.inputs[2])
        if n == len(set(atoms.get_atomic_numbers()))-1:
            supercell.links.new(switch_element.outputs[0], realize_instances_beforevectorcutoff.inputs[0])
        old_switch=switch_element
        is_element.width, is_element.height = 140.0, 100.0
        switch_element.width, switch_element.height = 140.0, 100.0
        is_element.location = (400, 1000-200*n)
        delete_geometry_element.location = (700, 1000-200*n)
        switch_element.location = (1000, 1000-200*n)
    element_attribute.location = (100, 600)
    group_input_element.location = (-100, 1000)
    supercell.links.new(realize_instances_beforevectorcutoff.outputs[0], switch_vectorcutoff.inputs[1])
    supercell.links.new(realize_instances_beforevectorcutoff.outputs[0], cutoff_group_check.inputs[0])
    supercell.links.new(group_input_vectorcutoff.outputs[7], switch_vectorcutoff.inputs[0])


    #Set parents
    realize_instances_beforevectorcutoff.parent = frame_002
    cutoff_group_check.parent = frame_002
    mesh_line_002.parent = frame_001
    mesh_line.parent = frame_001
    mesh_line_001.parent = frame_001
    instance_on_points_001.parent = frame_001
    instance_on_points_002.parent = frame_001
    vector_math_001.parent = frame
    vector_math.parent = frame
    vector_math_002.parent = frame
    math_002.parent = frame
    math.parent = frame
    math_001.parent = frame
    reroute.parent = frame_003
    vector_z.parent = frame_005
    vector_x.parent = frame_005
    vector_y.parent = frame_005
    switch.parent = frame_003
    switch_001.parent = frame_003
    switch_002.parent = frame_003
    switch_003.parent = frame_004
    switch_004.parent = frame_004
    switch_005.parent = frame_004

    #Set locations
    frame_002.location = (1262.3333740234375, 551.5)
    frame_001.location = (13.666666984558105, 316.3333435058594)
    frame.location = (-920.3333129882812, -263.6666564941406)
    instance_on_points.location = (637.8240356445312, 127.73196411132812)
    group_output_1.location = (2785.444580078125, 506.9361267089844)
    group_input_1.location = (376.5082092285156, -185.8363494873047)
    realize_instances_beforevectorcutoff.location = (29.2001953125, -39.510223388671875)
    cutoff_group_check.location = (192.3704833984375, -60.467041015625)
    mesh_line_002.location = (28.935646057128906, -338.6833801269531)
    mesh_line.location = (53.452327728271484, -44.30804443359375)
    mesh_line_001.location = (62.39809036254883, -179.75738525390625)
    instance_on_points_001.location = (256.01470947265625, -122.9981689453125)
    instance_on_points_002.location = (317.39739990234375, -290.7457580566406)
    vector_math_001.location = (337.35418701171875, -96.81759643554688)
    vector_math.location = (337.35418701171875, -45.817596435546875)
    vector_math_002.location = (337.35418701171875, -147.81759643554688)
    math_002.location = (28.9537353515625, -146.12448120117188)
    math.location = (28.9537353515625, -44.124481201171875)
    math_001.location = (28.9537353515625, -95.12448120117188)
    reroute.location = (33.8333740234375, -449.462158203125)
    reroute_001.location = (-1353.8707275390625, -18.929697036743164)
    group_input_001.location = (-1792.6153564453125, 362.8712158203125)
    reroute_004.location = (-674.0416259765625, 124.515869140625)
    reroute_002.location = (-654.8232421875, 293.1017761230469)
    reroute_003.location = (-589.20556640625, -89.482421875)
    reroute_007.location = (-652.3507690429688, -134.8733367919922)
    reroute_005.location = (-565.2334594726562, 253.76507568359375)
    reroute_006.location = (-613.8910522460938, 75.26302337646484)
    vector_z.location = (50.76611328125, -115.57769775390625)
    vector_x.location = (28.9681396484375, -33.61053466796875)
    vector_y.location = (42.9547119140625, -79.21051025390625)
    switch.location = (144.450927734375, -39.410369873046875)
    switch_001.location = (138.5106201171875, -191.1127471923828)
    switch_002.location = (139.989990234375, -345.19842529296875)
    switch_003.location = (32.6162109375, -39.38666534423828)
    switch_004.location = (29.1060791015625, -203.0687255859375)
    switch_005.location = (35.7261962890625, -375.9211730957031)
    group_input_002.location = (-1581.3948974609375, -258.635009765625)
    frame_003.location = (-1376.6806640625, 436.8333435058594)
    frame_004.location = (-1324.3333740234375, -127.83333587646484)
    frame_005.location = (-1077.6666259765625, 801.1666870117188)
    realize_instances_001.location = (2425.327880859375, 437.6210021972656)
    join_geometry.location = (2189.292724609375, 361.0550842285156)
    merge_by_distance.location = (2605.523193359375, 482.2424011230469)
    reroute_element.location = (909.9058227539062, 92.05680084228516)
    switch_vectorcutoff.location = (1843.2039794921875, 658.1240234375)
    group_input_vectorcutoff.location = (1098.405029296875, 399.9941101074219)

    #Set dimensions
    frame_002.width, frame_002.height = 361.333251953125, 310.5
    frame_001.width, frame_001.height = 464.39862060546875, 504.66668701171875
    frame.width, frame.height = 506.6666564941406, 200.83334350585938
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    group_output_1.width, group_output_1.height = 140.0, 100.0
    group_input_1.width, group_input_1.height = 140.0, 100.0
    realize_instances_beforevectorcutoff.width, realize_instances_beforevectorcutoff.height = 140.0, 100.0
    cutoff_group_check.width, cutoff_group_check.height = 140.0, 100.0
    mesh_line_002.width, mesh_line_002.height = 211.4036865234375, 100.0
    mesh_line.width, mesh_line.height = 211.4036865234375, 100.0
    mesh_line_001.width, mesh_line_001.height = 186.6614227294922, 100.0
    instance_on_points_001.width, instance_on_points_001.height = 140.0, 100.0
    instance_on_points_002.width, instance_on_points_002.height = 117.73193359375, 100.0
    vector_math_001.width, vector_math_001.height = 140.0, 100.0
    vector_math.width, vector_math.height = 140.0, 100.0
    vector_math_002.width, vector_math_002.height = 140.0, 100.0
    math_002.width, math_002.height = 148.916259765625, 100.0
    math.width, math.height = 148.916259765625, 100.0
    math_001.width, math_001.height = 148.916259765625, 100.0
    reroute.width, reroute.height = 14.5, 100.0
    reroute_001.width, reroute_001.height = 14.5, 100.0
    group_input_001.width, group_input_001.height = 140.0, 100.0
    reroute_004.width, reroute_004.height = 14.5, 100.0
    reroute_002.width, reroute_002.height = 14.5, 100.0
    reroute_003.width, reroute_003.height = 14.5, 100.0
    reroute_007.width, reroute_007.height = 14.5, 100.0
    reroute_005.width, reroute_005.height = 14.5, 100.0
    reroute_006.width, reroute_006.height = 14.5, 100.0
    vector_z.width, vector_z.height = 140.0, 100.0
    vector_x.width, vector_x.height = 140.0, 100.0
    vector_y.width, vector_y.height = 140.0, 100.0
    switch.width, switch.height = 140.0, 100.0
    switch_001.width, switch_001.height = 140.0, 100.0
    switch_002.width, switch_002.height = 140.0, 100.0
    switch_003.width, switch_003.height = 140.0, 100.0
    switch_004.width, switch_004.height = 140.0, 100.0
    switch_005.width, switch_005.height = 140.0, 100.0
    group_input_002.width, group_input_002.height = 140.0, 100.0
    frame_003.width, frame_003.height = 313.6806640625, 511.16668701171875
    frame_004.width, frame_004.height = 204.666748046875, 541.8333740234375
    frame_005.width, frame_005.height = 219.99993896484375, 222.8333740234375
    realize_instances_001.width, realize_instances_001.height = 140.0, 100.0
    element_attribute.width, element_attribute.height = 140.0, 100.0
    join_geometry.width, join_geometry.height = 140.0, 100.0
    merge_by_distance.width, merge_by_distance.height = 140.0, 100.0
    reroute_element.width, reroute_element.height = 14.5, 100.0
    switch_vectorcutoff.width, switch_vectorcutoff.height = 140.0, 100.0
    group_input_vectorcutoff.width, group_input_vectorcutoff.height = 140.0, 100.0

    #initialize supercell links
    #group_input_1.Geometry -> instance_on_points.Instance
    supercell.links.new(group_input_1.outputs[0], instance_on_points.inputs[2])
    #reroute_005.Output -> mesh_line.Offset
    supercell.links.new(reroute_005.outputs[0], mesh_line.inputs[3])
    #mesh_line.Mesh -> instance_on_points_001.Points
    supercell.links.new(mesh_line.outputs[0], instance_on_points_001.inputs[0])
    #reroute_006.Output -> mesh_line_001.Offset
    supercell.links.new(reroute_006.outputs[0], mesh_line_001.inputs[3])
    #reroute_007.Output -> mesh_line_002.Offset
    supercell.links.new(reroute_007.outputs[0], mesh_line_002.inputs[3])
    #mesh_line_001.Mesh -> instance_on_points_001.Instance
    supercell.links.new(mesh_line_001.outputs[0], instance_on_points_001.inputs[2])
    #instance_on_points_001.Instances -> instance_on_points_002.Points
    supercell.links.new(instance_on_points_001.outputs[0], instance_on_points_002.inputs[0])
    #mesh_line_002.Mesh -> instance_on_points_002.Instance
    supercell.links.new(mesh_line_002.outputs[0], instance_on_points_002.inputs[2])
    #instance_on_points_002.Instances -> instance_on_points.Points
    supercell.links.new(instance_on_points_002.outputs[0], instance_on_points.inputs[0])
    #math.Value -> vector_math.Scale
    supercell.links.new(math.outputs[0], vector_math.inputs[3])
    #vector_math.Vector -> mesh_line.Start Location
    supercell.links.new(vector_math.outputs[0], mesh_line.inputs[2])
    #math_001.Value -> vector_math_001.Scale
    supercell.links.new(math_001.outputs[0], vector_math_001.inputs[3])
    #vector_math_001.Vector -> mesh_line_001.Start Location
    supercell.links.new(vector_math_001.outputs[0], mesh_line_001.inputs[2])
    #math_002.Value -> vector_math_002.Scale
    supercell.links.new(math_002.outputs[0], vector_math_002.inputs[3])
    #vector_math_002.Vector -> mesh_line_002.Start Location
    supercell.links.new(vector_math_002.outputs[0], mesh_line_002.inputs[2])
    #realize_instances.Geometry -> group_001.Geometry
    #supercell.links.new(realize_instances_beforevectorcutoff.outputs[0], cutoff_group_check.inputs[0])
    #reroute_002.Output -> mesh_line.Count
    supercell.links.new(reroute_002.outputs[0], mesh_line.inputs[0])
    #reroute_004.Output -> mesh_line_001.Count
    supercell.links.new(reroute_004.outputs[0], mesh_line_001.inputs[0])
    #reroute_003.Output -> mesh_line_002.Count
    supercell.links.new(reroute_003.outputs[0], mesh_line_002.inputs[0])
    #switch_003.Output -> math.Value
    supercell.links.new(switch_003.outputs[0], math.inputs[0])
    #switch_004.Output -> math_001.Value
    supercell.links.new(switch_004.outputs[0], math_001.inputs[0])
    #merge_by_distance.Geometry -> group_output_1.Geometry
    supercell.links.new(merge_by_distance.outputs[0], group_output_1.inputs[0])
    #group_input_001.Offset_x -> reroute.Input
    supercell.links.new(group_input_001.outputs[4], reroute.inputs[0])
    #group_input_001.Offset_y -> reroute_001.Input
    supercell.links.new(group_input_001.outputs[5], reroute_001.inputs[0])
    #switch_002.Output -> reroute_003.Input
    supercell.links.new(switch_002.outputs[0], reroute_003.inputs[0])
    #switch_001.Output -> reroute_004.Input
    supercell.links.new(switch_001.outputs[0], reroute_004.inputs[0])
    #vector_001.Vector -> reroute_006.Input
    supercell.links.new(vector_y.outputs[0], reroute_006.inputs[0])
    #vector_002.Vector -> reroute_007.Input
    supercell.links.new(vector_z.outputs[0], reroute_007.inputs[0])
    #group_input_001.global -> switch.Switch
    supercell.links.new(group_input_001.outputs[6], switch.inputs[0])
    #group_input_001.repeat_x -> switch.False
    supercell.links.new(group_input_001.outputs[1], switch.inputs[1])
    #switch.Output -> reroute_002.Input
    supercell.links.new(switch.outputs[0], reroute_002.inputs[0])
    #group_input_001.repeat_y -> switch_001.False
    supercell.links.new(group_input_001.outputs[2], switch_001.inputs[1])
    #group_input_001.repeat_z -> switch_002.False
    supercell.links.new(group_input_001.outputs[3], switch_002.inputs[1])
    #group_input_001.global -> switch_001.Switch
    supercell.links.new(group_input_001.outputs[6], switch_001.inputs[0])
    #group_input_001.global -> switch_002.Switch
    supercell.links.new(group_input_001.outputs[6], switch_002.inputs[0])
    #reroute.Output -> switch_003.False
    supercell.links.new(reroute.outputs[0], switch_003.inputs[1])
    #reroute_001.Output -> switch_004.False
    supercell.links.new(reroute_001.outputs[0], switch_004.inputs[1])
    #switch_005.Output -> math_002.Value
    supercell.links.new(switch_005.outputs[0], math_002.inputs[0])
    #reroute_001.Output -> switch_005.False
    supercell.links.new(reroute_001.outputs[0], switch_005.inputs[1])
    #group_input_002.global -> switch_003.Switch
    supercell.links.new(group_input_002.outputs[6], switch_003.inputs[0])
    #group_input_002.global -> switch_004.Switch
    supercell.links.new(group_input_002.outputs[6], switch_004.inputs[0])
    #group_input_002.global -> switch_005.Switch
    supercell.links.new(group_input_002.outputs[6], switch_005.inputs[0])
    #reroute_005.Output -> vector_math.Vector
    supercell.links.new(reroute_005.outputs[0], vector_math.inputs[0])
    #reroute_006.Output -> vector_math_001.Vector
    supercell.links.new(reroute_006.outputs[0], vector_math_001.inputs[0])
    #reroute_007.Output -> vector_math_002.Vector
    supercell.links.new(reroute_007.outputs[0], vector_math_002.inputs[0])
    #vector.Vector -> reroute_005.Input
    supercell.links.new(vector_x.outputs[0], reroute_005.inputs[0])
    #join_geometry.Geometry -> realize_instances_001.Geometry
    supercell.links.new(join_geometry.outputs[0], realize_instances_001.inputs[0])
    #realize_instances_001.Geometry -> merge_by_distance.Geometry
    supercell.links.new(realize_instances_001.outputs[0], merge_by_distance.inputs[0])
    #instance_on_points.Instances -> reroute_008.Input
    supercell.links.new(instance_on_points.outputs[0], reroute_element.inputs[0])
    #group_001.Geometry -> switch_006.True
    supercell.links.new(cutoff_group_check.outputs[0], switch_vectorcutoff.inputs[2])
    #reroute_008.Output -> realize_instances.Geometry
    #supercell.links.new(reroute_element.outputs[0], realize_instances.inputs[0])
    #reroute_008.Output -> switch_006.False
  #  supercell.links.new(reroute_element.outputs[0], switch_vectorcutoff.inputs[1])
    #group_input_004.+x -> group_001.+x
    supercell.links.new(group_input_vectorcutoff.outputs[8], cutoff_group_check.inputs[1])
    #group_input_004.+y -> group_001.+y
    supercell.links.new(group_input_vectorcutoff.outputs[9], cutoff_group_check.inputs[2])
    #group_input_004.+z -> group_001.+z
    supercell.links.new(group_input_vectorcutoff.outputs[10], cutoff_group_check.inputs[3])
    #group_input_004.-x -> group_001.-x
    supercell.links.new(group_input_vectorcutoff.outputs[11], cutoff_group_check.inputs[4])
    #group_input_004.-y -> group_001.-y
    supercell.links.new(group_input_vectorcutoff.outputs[12], cutoff_group_check.inputs[5])
    #group_input_004.-z -> group_001.-z
    supercell.links.new(group_input_vectorcutoff.outputs[13], cutoff_group_check.inputs[6])
    
    #switch_006.Output -> join_geometry.Geometry
    supercell.links.new(switch_vectorcutoff.outputs[0], join_geometry.inputs[0])


    #initialize supercell_001 nodes
    #node Compare
    is_element = supercell.nodes.new("FunctionNodeCompare")
    is_element.name = "Compare"
    is_element.hide = True
    is_element.data_type = 'INT'
    is_element.mode = 'ELEMENT'
    is_element.operation = 'EQUAL'
    #A_INT
    is_element.inputs[2].default_value = 0
    #B_INT
    is_element.inputs[3].default_value = 1

    #node Delete Geometry
    delete_geometry = supercell.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.name = "Delete Geometry"
    delete_geometry.hide = True
    delete_geometry.domain = 'POINT'
    delete_geometry.mode = 'ALL'

    #node Switch.007
    switch_007 = supercell.nodes.new("GeometryNodeSwitch")
    switch_007.name = "Switch.007"
    switch_007.hide = True
    switch_007.input_type = 'GEOMETRY'

    
    return supercell

def supercell_atoms_node_group(atoms):
    cutoff_group=cutoff_group_node_group()
    supercell_atoms = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "supercell_atoms")
    
    supercell_atoms.color_tag = 'NONE'
    supercell_atoms.description = ""
    supercell_atoms.default_group_node_width = 140
    

    supercell_atoms.is_modifier = True

    #supercell_atoms interface
    #Socket Geometry
    geometry_socket_2 = supercell_atoms.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_2.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_3 = supercell_atoms.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_3.attribute_domain = 'POINT'

    #Socket repeat_x
    repeat_x_socket = supercell_atoms.interface.new_socket(name = "repeat_x", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_x_socket.default_value = 1
    repeat_x_socket.min_value = 1
    repeat_x_socket.max_value = 10000
    repeat_x_socket.subtype = 'NONE'
    repeat_x_socket.attribute_domain = 'POINT'

    #Socket repeat_y
    repeat_y_socket = supercell_atoms.interface.new_socket(name = "repeat_y", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_y_socket.default_value = 1
    repeat_y_socket.min_value = 1
    repeat_y_socket.max_value = 10000
    repeat_y_socket.subtype = 'NONE'
    repeat_y_socket.attribute_domain = 'POINT'

    #Socket repeat_z
    repeat_z_socket = supercell_atoms.interface.new_socket(name = "repeat_z", in_out='INPUT', socket_type = 'NodeSocketInt')
    repeat_z_socket.default_value = 1
    repeat_z_socket.min_value = 1
    repeat_z_socket.max_value = 10000
    repeat_z_socket.subtype = 'NONE'
    repeat_z_socket.attribute_domain = 'POINT'

    #Socket Offset_x
    offset_x_socket = supercell_atoms.interface.new_socket(name = "Offset_x", in_out='INPUT', socket_type = 'NodeSocketInt')
    offset_x_socket.default_value = 0
    offset_x_socket.min_value = -2147483648
    offset_x_socket.max_value = 2147483647
    offset_x_socket.subtype = 'NONE'
    offset_x_socket.attribute_domain = 'POINT'

    #Socket Offset_y
    offset_y_socket = supercell_atoms.interface.new_socket(name = "Offset_y", in_out='INPUT', socket_type = 'NodeSocketInt')
    offset_y_socket.default_value = 0
    offset_y_socket.min_value = -2147483648
    offset_y_socket.max_value = 2147483647
    offset_y_socket.subtype = 'NONE'
    offset_y_socket.attribute_domain = 'POINT'

    #Socket global
    global_socket = supercell_atoms.interface.new_socket(name = "global", in_out='INPUT', socket_type = 'NodeSocketBool')
    global_socket.default_value = False
    global_socket.attribute_domain = 'POINT'

    #Socket cutoff_vectors
    cutoff_vectors_socket = supercell_atoms.interface.new_socket(name = "cutoff_vectors", in_out='INPUT', socket_type = 'NodeSocketBool')
    cutoff_vectors_socket.default_value = False
    cutoff_vectors_socket.attribute_domain = 'POINT'

    #Socket +x
    _x_socket_2 = supercell_atoms.interface.new_socket(name = "+x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket_2.default_value = 100.0
    _x_socket_2.min_value = -10000.0
    _x_socket_2.max_value = 10000.0
    _x_socket_2.subtype = 'NONE'
    _x_socket_2.attribute_domain = 'POINT'

    #Socket +y
    _y_socket_2 = supercell_atoms.interface.new_socket(name = "+y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket_2.default_value = 100.0
    _y_socket_2.min_value = -10000.0
    _y_socket_2.max_value = 10000.0
    _y_socket_2.subtype = 'NONE'
    _y_socket_2.attribute_domain = 'POINT'

    #Socket +z
    _z_socket_2 = supercell_atoms.interface.new_socket(name = "+z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket_2.default_value = 100.0
    _z_socket_2.min_value = -10000.0
    _z_socket_2.max_value = 10000.0
    _z_socket_2.subtype = 'NONE'
    _z_socket_2.attribute_domain = 'POINT'

    #Socket -x
    _x_socket_3 = supercell_atoms.interface.new_socket(name = "-x", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _x_socket_3.default_value = -100.0
    _x_socket_3.min_value = -10000.0
    _x_socket_3.max_value = 10000.0
    _x_socket_3.subtype = 'NONE'
    _x_socket_3.attribute_domain = 'POINT'

    #Socket -y
    _y_socket_3 = supercell_atoms.interface.new_socket(name = "-y", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _y_socket_3.default_value = -100.0
    _y_socket_3.min_value = -10000.0
    _y_socket_3.max_value = 10000.0
    _y_socket_3.subtype = 'NONE'
    _y_socket_3.attribute_domain = 'POINT'

    #Socket -z
    _z_socket_3 = supercell_atoms.interface.new_socket(name = "-z", in_out='INPUT', socket_type = 'NodeSocketFloat')
    _z_socket_3.default_value = -100.0
    _z_socket_3.min_value = -10000.0
    _z_socket_3.max_value = 10000.0
    _z_socket_3.subtype = 'NONE'
    _z_socket_3.attribute_domain = 'POINT'


    #initialize supercell_atoms nodes
    #node Frame.002
    frame_002 = supercell_atoms.nodes.new("NodeFrame")
    frame_002.label = "cutoff"
    frame_002.name = "Frame.002"
    frame_002.label_size = 20
    frame_002.shrink = True

    #node Frame.001
    frame_001 = supercell_atoms.nodes.new("NodeFrame")
    frame_001.label = "grid"
    frame_001.name = "Frame.001"
    frame_001.label_size = 20
    frame_001.shrink = True

    #node Frame
    frame = supercell_atoms.nodes.new("NodeFrame")
    frame.label = "offset"
    frame.name = "Frame"
    frame.label_size = 20
    frame.shrink = True

    #node Instance on Points
    instance_on_points = supercell_atoms.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points.name = "Instance on Points"
    #Selection
    instance_on_points.inputs[1].default_value = True
    #Pick Instance
    instance_on_points.inputs[3].default_value = False
    #Instance Index
    instance_on_points.inputs[4].default_value = 0
    #Rotation
    instance_on_points.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Group Output
    group_output_1 = supercell_atoms.nodes.new("NodeGroupOutput")
    group_output_1.name = "Group Output"
    group_output_1.is_active_output = True

    #node Group Input
    group_input_1 = supercell_atoms.nodes.new("NodeGroupInput")
    group_input_1.name = "Group Input"

    #node Group.001
    group_001 = supercell_atoms.nodes.new("GeometryNodeGroup")
    group_001.name = "Group.001"
    group_001.node_tree = cutoff_group

    #node Mesh Line.002
    mesh_line_002 = supercell_atoms.nodes.new("GeometryNodeMeshLine")
    mesh_line_002.name = "Mesh Line.002"
    mesh_line_002.count_mode = 'TOTAL'
    mesh_line_002.mode = 'OFFSET'

    #node Mesh Line
    mesh_line = supercell_atoms.nodes.new("GeometryNodeMeshLine")
    mesh_line.name = "Mesh Line"
    mesh_line.hide = True
    mesh_line.count_mode = 'TOTAL'
    mesh_line.mode = 'OFFSET'

    #node Mesh Line.001
    mesh_line_001 = supercell_atoms.nodes.new("GeometryNodeMeshLine")
    mesh_line_001.name = "Mesh Line.001"
    mesh_line_001.hide = True
    mesh_line_001.count_mode = 'TOTAL'
    mesh_line_001.mode = 'OFFSET'

    #node Instance on Points.001
    instance_on_points_001 = supercell_atoms.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_001.name = "Instance on Points.001"
    instance_on_points_001.hide = True
    #Selection
    instance_on_points_001.inputs[1].default_value = True
    #Pick Instance
    instance_on_points_001.inputs[3].default_value = False
    #Instance Index
    instance_on_points_001.inputs[4].default_value = 0
    #Rotation
    instance_on_points_001.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points_001.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Instance on Points.002
    instance_on_points_002 = supercell_atoms.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_002.name = "Instance on Points.002"
    instance_on_points_002.hide = True
    #Selection
    instance_on_points_002.inputs[1].default_value = True
    #Pick Instance
    instance_on_points_002.inputs[3].default_value = False
    #Instance Index
    instance_on_points_002.inputs[4].default_value = 0
    #Rotation
    instance_on_points_002.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Scale
    instance_on_points_002.inputs[6].default_value = (1.0, 1.0, 1.0)

    #node Vector Math.001
    vector_math_001 = supercell_atoms.nodes.new("ShaderNodeVectorMath")
    vector_math_001.name = "Vector Math.001"
    vector_math_001.hide = True
    vector_math_001.operation = 'SCALE'

    #node Vector Math
    vector_math = supercell_atoms.nodes.new("ShaderNodeVectorMath")
    vector_math.name = "Vector Math"
    vector_math.hide = True
    vector_math.operation = 'SCALE'

    #node Vector Math.002
    vector_math_002 = supercell_atoms.nodes.new("ShaderNodeVectorMath")
    vector_math_002.name = "Vector Math.002"
    vector_math_002.hide = True
    vector_math_002.operation = 'SCALE'

    #node Math.002
    math_002 = supercell_atoms.nodes.new("ShaderNodeMath")
    math_002.name = "Math.002"
    math_002.hide = True
    math_002.operation = 'MULTIPLY'
    math_002.use_clamp = False
    #Value_001
    math_002.inputs[1].default_value = -1.0

    #node Math
    math = supercell_atoms.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.hide = True
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = -1.0

    #node Math.001
    math_001 = supercell_atoms.nodes.new("ShaderNodeMath")
    math_001.name = "Math.001"
    math_001.hide = True
    math_001.operation = 'MULTIPLY'
    math_001.use_clamp = False
    #Value_001
    math_001.inputs[1].default_value = -1.0

    #node Reroute
    reroute = supercell_atoms.nodes.new("NodeReroute")
    reroute.name = "Reroute"
    reroute.socket_idname = "NodeSocketInt"
    #node Reroute.001
    reroute_001 = supercell_atoms.nodes.new("NodeReroute")
    reroute_001.name = "Reroute.001"
    reroute_001.socket_idname = "NodeSocketInt"
    #node Group Input.001
    group_input_001 = supercell_atoms.nodes.new("NodeGroupInput")
    group_input_001.name = "Group Input.001"

    #node Reroute.004
    reroute_004 = supercell_atoms.nodes.new("NodeReroute")
    reroute_004.name = "Reroute.004"
    reroute_004.socket_idname = "NodeSocketInt"
    #node Reroute.002
    reroute_002 = supercell_atoms.nodes.new("NodeReroute")
    reroute_002.name = "Reroute.002"
    reroute_002.socket_idname = "NodeSocketInt"
    #node Reroute.003
    reroute_003 = supercell_atoms.nodes.new("NodeReroute")
    reroute_003.name = "Reroute.003"
    reroute_003.socket_idname = "NodeSocketInt"
    #node Reroute.007
    reroute_007 = supercell_atoms.nodes.new("NodeReroute")
    reroute_007.name = "Reroute.007"
    reroute_007.socket_idname = "NodeSocketVector"
    #node Reroute.005
    reroute_005 = supercell_atoms.nodes.new("NodeReroute")
    reroute_005.name = "Reroute.005"
    reroute_005.socket_idname = "NodeSocketVector"
    #node Reroute.006
    reroute_006 = supercell_atoms.nodes.new("NodeReroute")
    reroute_006.name = "Reroute.006"
    reroute_006.socket_idname = "NodeSocketVector"
    #node Vector.002
    vector_z = supercell_atoms.nodes.new("FunctionNodeInputVector")
    vector_z.label = "zvec"
    vector_z.name = "Vector.002"
    vector_z.hide = True
    vector_z.vector = atoms.cell[2]

    #node Vector
    vector_x = supercell_atoms.nodes.new("FunctionNodeInputVector")
    vector_x.label = "xvec"
    vector_x.name = "Vector"
    vector_x.hide = True
    vector_x.vector = atoms.cell[0]

    #node Vector.001
    vector_y = supercell_atoms.nodes.new("FunctionNodeInputVector")
    vector_y.label = "yvec"
    vector_y.name = "Vector.001"
    vector_y.hide = True
    vector_y.vector = atoms.cell[1]

    #node Switch
    switch = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch.label = "repeat_X"
    switch.name = "Switch"
    switch.input_type = 'INT'
    #True
    switch.inputs[2].default_value = 3

    #node Switch.001
    switch_001 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_001.label = "repeat_Y"
    switch_001.name = "Switch.001"
    switch_001.input_type = 'INT'
    #True
    switch_001.inputs[2].default_value = 3

    #node Switch.002
    switch_002 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_002.label = "repeat_Z"
    switch_002.name = "Switch.002"
    switch_002.input_type = 'INT'
    #True
    switch_002.inputs[2].default_value = 1

    #node Switch.003
    switch_003 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_003.label = "offset_X"
    switch_003.name = "Switch.003"
    switch_003.input_type = 'INT'
    #True
    switch_003.inputs[2].default_value = 0

    #node Switch.004
    switch_004 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_004.label = "offset_Y"
    switch_004.name = "Switch.004"
    switch_004.input_type = 'INT'
    #True
    switch_004.inputs[2].default_value = 1

    #node Switch.005
    switch_005 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_005.label = "offset_Z"
    switch_005.name = "Switch.005"
    switch_005.input_type = 'INT'
    #True
    switch_005.inputs[2].default_value = 0

    #node Group Input.002
    group_input_002 = supercell_atoms.nodes.new("NodeGroupInput")
    group_input_002.name = "Group Input.002"

    #node Frame.003
    frame_003 = supercell_atoms.nodes.new("NodeFrame")
    frame_003.label = "modify global value for cell repetition"
    frame_003.name = "Frame.003"
    frame_003.use_custom_color = True
    frame_003.color = (0.41161254048347473, 0.0, 0.6079999804496765)
    frame_003.label_size = 20
    frame_003.shrink = True

    #node Frame.004
    frame_004 = supercell_atoms.nodes.new("NodeFrame")
    frame_004.label = "global value for offset"
    frame_004.name = "Frame.004"
    frame_004.use_custom_color = True
    frame_004.color = (0.3529382050037384, 0.0235294159501791, 0.5098000168800354)
    frame_004.label_size = 20
    frame_004.shrink = True

    #node Frame.005
    frame_005 = supercell_atoms.nodes.new("NodeFrame")
    frame_005.name = "Frame.005"
    frame_005.use_custom_color = True
    frame_005.color = (0.08440613746643066, 0.6079999804496765, 0.0)
    frame_005.label_size = 20
    frame_005.shrink = True

    #node Join Geometry
    join_geometry = supercell_atoms.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"

    #node Merge by Distance
    merge_by_distance = supercell_atoms.nodes.new("GeometryNodeMergeByDistance")
    merge_by_distance.name = "Merge by Distance"
    merge_by_distance.mode = 'ALL'
    #Selection
    merge_by_distance.inputs[1].default_value = True
    #Distance
    merge_by_distance.inputs[2].default_value = 0.0010000000474974513

    #node Reroute.008
    reroute_008 = supercell_atoms.nodes.new("NodeReroute")
    reroute_008.name = "Reroute.008"
    reroute_008.socket_idname = "NodeSocketGeometry"

    #node Switch.006
    switch_006 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_006.name = "Switch.006"
    switch_006.input_type = 'GEOMETRY'

    #node Group Input.004
    group_input_004 = supercell_atoms.nodes.new("NodeGroupInput")
    group_input_004.name = "Group Input.004"

    #node Compare.004
    compare_004_1 = supercell_atoms.nodes.new("FunctionNodeCompare")
    compare_004_1.name = "Compare.004"
    compare_004_1.hide = True
    compare_004_1.data_type = 'INT'
    compare_004_1.mode = 'ELEMENT'
    compare_004_1.operation = 'EQUAL'
    #A_INT
    compare_004_1.inputs[2].default_value = 0
    #B_INT
    compare_004_1.inputs[3].default_value = 1

    #node Delete Geometry
    delete_geometry_1 = supercell_atoms.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_1.name = "Delete Geometry"
    delete_geometry_1.hide = True
    delete_geometry_1.domain = 'POINT'
    delete_geometry_1.mode = 'ALL'
    #Selection
    delete_geometry_1.inputs[1].default_value = True

    #node Switch.007
    switch_007 = supercell_atoms.nodes.new("GeometryNodeSwitch")
    switch_007.name = "Switch.007"
    switch_007.hide = True
    switch_007.input_type = 'GEOMETRY'
    #Switch
    switch_007.inputs[0].default_value = False




    #Set parents
    group_001.parent = frame_002
    mesh_line_002.parent = frame_001
    mesh_line.parent = frame_001
    mesh_line_001.parent = frame_001
    instance_on_points_001.parent = frame_001
    instance_on_points_002.parent = frame_001
    vector_math_001.parent = frame
    vector_math.parent = frame
    vector_math_002.parent = frame
    math_002.parent = frame
    math.parent = frame
    math_001.parent = frame
    reroute.parent = frame_003
    vector_z.parent = frame_005
    vector_x.parent = frame_005
    vector_y.parent = frame_005
    switch.parent = frame_003
    switch_001.parent = frame_003
    switch_002.parent = frame_003
    switch_003.parent = frame_004
    switch_004.parent = frame_004
    switch_005.parent = frame_004

    #Set locations
    frame_002.location = (1456.0, 441.0)
    frame_001.location = (13.0, 317.0)
    frame.location = (-921.0, -263.0)
    instance_on_points.location = (637.8240356445312, 127.73196411132812)
    group_output_1.location = (2785.444580078125, 506.9361267089844)
    group_input_1.location = (376.5082092285156, -185.8363494873047)
    group_001.location = (29.771484375, -40.318511962890625)
    mesh_line_002.location = (29.602313995361328, -339.35003662109375)
    mesh_line.location = (54.118995666503906, -44.974700927734375)
    mesh_line_001.location = (63.06475830078125, -180.42404174804688)
    instance_on_points_001.location = (256.6813659667969, -123.66482543945312)
    instance_on_points_002.location = (318.0640563964844, -291.41241455078125)
    vector_math_001.location = (338.0208740234375, -97.4842529296875)
    vector_math.location = (338.0208740234375, -46.4842529296875)
    vector_math_002.location = (338.0208740234375, -148.4842529296875)
    math_002.location = (29.62042236328125, -146.7911376953125)
    math.location = (29.62042236328125, -44.7911376953125)
    math_001.location = (29.62042236328125, -95.7911376953125)
    reroute.location = (35.0, -449.6288146972656)
    reroute_001.location = (-1353.8707275390625, -18.929697036743164)
    group_input_001.location = (-1792.6153564453125, 362.8712158203125)
    reroute_004.location = (-674.0416259765625, 124.515869140625)
    reroute_002.location = (-654.8232421875, 293.1017761230469)
    reroute_003.location = (-589.20556640625, -89.482421875)
    reroute_007.location = (-652.3507690429688, -134.8733367919922)
    reroute_005.location = (-565.2334594726562, 253.76507568359375)
    reroute_006.location = (-613.8910522460938, 75.26302337646484)
    vector_z.location = (52.0994873046875, -117.4110107421875)
    vector_x.location = (30.301513671875, -35.44384765625)
    vector_y.location = (44.2880859375, -81.0438232421875)
    switch.location = (145.6175537109375, -39.5770263671875)
    switch_001.location = (139.67724609375, -191.27940368652344)
    switch_002.location = (141.1566162109375, -345.3650817871094)
    switch_003.location = (33.2828369140625, -40.220001220703125)
    switch_004.location = (29.772705078125, -203.90206909179688)
    switch_005.location = (36.392822265625, -376.7545166015625)
    group_input_002.location = (-1581.3948974609375, -258.635009765625)
    frame_003.location = (-1377.8472900390625, 437.0)
    frame_004.location = (-1325.0, -127.0)
    frame_005.location = (-1079.0, 803.0)
    join_geometry.location = (2189.292724609375, 361.0550842285156)
    merge_by_distance.location = (2605.523193359375, 482.2424011230469)
    reroute_008.location = (909.9058227539062, 92.05680084228516)
    switch_006.location = (1843.2039794921875, 658.1240234375)
    group_input_004.location = (1087.9091796875, 835.0838012695312)
    compare_004_1.location = (0.0, 0.0)
    delete_geometry_1.location = (0.0, 0.0)
    switch_007.location = (0.0, 0.0)

    #Set dimensions
    frame_002.width, frame_002.height = 200.0, 298.0
    frame_001.width, frame_001.height = 465.73193359375, 512.0
    frame.width, frame.height = 508.0, 203.0
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    group_output_1.width, group_output_1.height = 140.0, 100.0
    group_input_1.width, group_input_1.height = 140.0, 100.0
    group_001.width, group_001.height = 140.0, 100.0
    mesh_line_002.width, mesh_line_002.height = 211.4036865234375, 100.0
    mesh_line.width, mesh_line.height = 211.4036865234375, 100.0
    mesh_line_001.width, mesh_line_001.height = 186.6614227294922, 100.0
    instance_on_points_001.width, instance_on_points_001.height = 140.0, 100.0
    instance_on_points_002.width, instance_on_points_002.height = 117.73193359375, 100.0
    vector_math_001.width, vector_math_001.height = 140.0, 100.0
    vector_math.width, vector_math.height = 140.0, 100.0
    vector_math_002.width, vector_math_002.height = 140.0, 100.0
    math_002.width, math_002.height = 148.916259765625, 100.0
    math.width, math.height = 148.916259765625, 100.0
    math_001.width, math_001.height = 148.916259765625, 100.0
    reroute.width, reroute.height = 10.0, 100.0
    reroute_001.width, reroute_001.height = 10.0, 100.0
    group_input_001.width, group_input_001.height = 140.0, 100.0
    reroute_004.width, reroute_004.height = 10.0, 100.0
    reroute_002.width, reroute_002.height = 10.0, 100.0
    reroute_003.width, reroute_003.height = 10.0, 100.0
    reroute_007.width, reroute_007.height = 10.0, 100.0
    reroute_005.width, reroute_005.height = 10.0, 100.0
    reroute_006.width, reroute_006.height = 10.0, 100.0
    vector_z.width, vector_z.height = 140.0, 100.0
    vector_x.width, vector_x.height = 140.0, 100.0
    vector_y.width, vector_y.height = 140.0, 100.0
    switch.width, switch.height = 140.0, 100.0
    switch_001.width, switch_001.height = 140.0, 100.0
    switch_002.width, switch_002.height = 140.0, 100.0
    switch_003.width, switch_003.height = 140.0, 100.0
    switch_004.width, switch_004.height = 140.0, 100.0
    switch_005.width, switch_005.height = 140.0, 100.0
    group_input_002.width, group_input_002.height = 140.0, 100.0
    frame_003.width, frame_003.height = 315.8472900390625, 518.0
    frame_004.width, frame_004.height = 206.0, 550.0
    frame_005.width, frame_005.height = 222.0, 230.0
    join_geometry.width, join_geometry.height = 140.0, 100.0
    merge_by_distance.width, merge_by_distance.height = 140.0, 100.0
    reroute_008.width, reroute_008.height = 10.0, 100.0
    switch_006.width, switch_006.height = 140.0, 100.0
    group_input_004.width, group_input_004.height = 140.0, 100.0
    compare_004_1.width, compare_004_1.height = 140.0, 100.0
    delete_geometry_1.width, delete_geometry_1.height = 140.0, 100.0
    switch_007.width, switch_007.height = 140.0, 100.0

    #initialize supercell_atoms links
    #group_input_004.cutoff_vectors -> switch_006.Switch
    supercell_atoms.links.new(group_input_004.outputs[7], switch_006.inputs[0])
    #group_input_1.Geometry -> instance_on_points.Instance
    supercell_atoms.links.new(group_input_1.outputs[0], instance_on_points.inputs[2])
    #reroute_005.Output -> mesh_line.Offset
    supercell_atoms.links.new(reroute_005.outputs[0], mesh_line.inputs[3])
    #mesh_line.Mesh -> instance_on_points_001.Points
    supercell_atoms.links.new(mesh_line.outputs[0], instance_on_points_001.inputs[0])
    #reroute_006.Output -> mesh_line_001.Offset
    supercell_atoms.links.new(reroute_006.outputs[0], mesh_line_001.inputs[3])
    #reroute_007.Output -> mesh_line_002.Offset
    supercell_atoms.links.new(reroute_007.outputs[0], mesh_line_002.inputs[3])
    #mesh_line_001.Mesh -> instance_on_points_001.Instance
    supercell_atoms.links.new(mesh_line_001.outputs[0], instance_on_points_001.inputs[2])
    #instance_on_points_001.Instances -> instance_on_points_002.Points
    supercell_atoms.links.new(instance_on_points_001.outputs[0], instance_on_points_002.inputs[0])
    #mesh_line_002.Mesh -> instance_on_points_002.Instance
    supercell_atoms.links.new(mesh_line_002.outputs[0], instance_on_points_002.inputs[2])
    #instance_on_points_002.Instances -> instance_on_points.Points
    supercell_atoms.links.new(instance_on_points_002.outputs[0], instance_on_points.inputs[0])
    #math.Value -> vector_math.Scale
    supercell_atoms.links.new(math.outputs[0], vector_math.inputs[3])
    #vector_math.Vector -> mesh_line.Start Location
    supercell_atoms.links.new(vector_math.outputs[0], mesh_line.inputs[2])
    #math_001.Value -> vector_math_001.Scale
    supercell_atoms.links.new(math_001.outputs[0], vector_math_001.inputs[3])
    #vector_math_001.Vector -> mesh_line_001.Start Location
    supercell_atoms.links.new(vector_math_001.outputs[0], mesh_line_001.inputs[2])
    #math_002.Value -> vector_math_002.Scale
    supercell_atoms.links.new(math_002.outputs[0], vector_math_002.inputs[3])
    #vector_math_002.Vector -> mesh_line_002.Start Location
    supercell_atoms.links.new(vector_math_002.outputs[0], mesh_line_002.inputs[2])
    #reroute_002.Output -> mesh_line.Count
    supercell_atoms.links.new(reroute_002.outputs[0], mesh_line.inputs[0])
    #reroute_004.Output -> mesh_line_001.Count
    supercell_atoms.links.new(reroute_004.outputs[0], mesh_line_001.inputs[0])
    #reroute_003.Output -> mesh_line_002.Count
    supercell_atoms.links.new(reroute_003.outputs[0], mesh_line_002.inputs[0])
    #switch_003.Output -> math.Value
    supercell_atoms.links.new(switch_003.outputs[0], math.inputs[0])
    #switch_004.Output -> math_001.Value
    supercell_atoms.links.new(switch_004.outputs[0], math_001.inputs[0])
    #merge_by_distance.Geometry -> group_output_1.Geometry
    supercell_atoms.links.new(merge_by_distance.outputs[0], group_output_1.inputs[0])
    #group_input_001.Offset_x -> reroute.Input
    supercell_atoms.links.new(group_input_001.outputs[4], reroute.inputs[0])
    #group_input_001.Offset_y -> reroute_001.Input
    supercell_atoms.links.new(group_input_001.outputs[5], reroute_001.inputs[0])
    #switch_002.Output -> reroute_003.Input
    supercell_atoms.links.new(switch_002.outputs[0], reroute_003.inputs[0])
    #switch_001.Output -> reroute_004.Input
    supercell_atoms.links.new(switch_001.outputs[0], reroute_004.inputs[0])
    #vector_001.Vector -> reroute_006.Input
    supercell_atoms.links.new(vector_y.outputs[0], reroute_006.inputs[0])
    #vector_002.Vector -> reroute_007.Input
    supercell_atoms.links.new(vector_z.outputs[0], reroute_007.inputs[0])
    #group_input_001.global -> switch.Switch
    supercell_atoms.links.new(group_input_001.outputs[6], switch.inputs[0])
    #group_input_001.repeat_x -> switch.False
    supercell_atoms.links.new(group_input_001.outputs[1], switch.inputs[1])
    #switch.Output -> reroute_002.Input
    supercell_atoms.links.new(switch.outputs[0], reroute_002.inputs[0])
    #group_input_001.repeat_y -> switch_001.False
    supercell_atoms.links.new(group_input_001.outputs[2], switch_001.inputs[1])
    #group_input_001.repeat_z -> switch_002.False
    supercell_atoms.links.new(group_input_001.outputs[3], switch_002.inputs[1])
    #group_input_001.global -> switch_001.Switch
    supercell_atoms.links.new(group_input_001.outputs[6], switch_001.inputs[0])
    #group_input_001.global -> switch_002.Switch
    supercell_atoms.links.new(group_input_001.outputs[6], switch_002.inputs[0])
    #reroute.Output -> switch_003.False
    supercell_atoms.links.new(reroute.outputs[0], switch_003.inputs[1])
    #reroute_001.Output -> switch_004.False
    supercell_atoms.links.new(reroute_001.outputs[0], switch_004.inputs[1])
    #switch_005.Output -> math_002.Value
    supercell_atoms.links.new(switch_005.outputs[0], math_002.inputs[0])
    #group_input_002.global -> switch_003.Switch
    supercell_atoms.links.new(group_input_002.outputs[6], switch_003.inputs[0])
    #group_input_002.global -> switch_004.Switch
    supercell_atoms.links.new(group_input_002.outputs[6], switch_004.inputs[0])
    #group_input_002.global -> switch_005.Switch
    supercell_atoms.links.new(group_input_002.outputs[6], switch_005.inputs[0])
    #reroute_005.Output -> vector_math.Vector
    supercell_atoms.links.new(reroute_005.outputs[0], vector_math.inputs[0])
    #reroute_006.Output -> vector_math_001.Vector
    supercell_atoms.links.new(reroute_006.outputs[0], vector_math_001.inputs[0])
    #reroute_007.Output -> vector_math_002.Vector
    supercell_atoms.links.new(reroute_007.outputs[0], vector_math_002.inputs[0])
    #vector.Vector -> reroute_005.Input
    supercell_atoms.links.new(vector_x.outputs[0], reroute_005.inputs[0])
    #join_geometry.Geometry -> merge_by_distance.Geometry
    supercell_atoms.links.new(join_geometry.outputs[0], merge_by_distance.inputs[0])
    #instance_on_points.Instances -> reroute_008.Input
    supercell_atoms.links.new(instance_on_points.outputs[0], reroute_008.inputs[0])
    #group_001.Geometry -> switch_006.True
    supercell_atoms.links.new(group_001.outputs[0], switch_006.inputs[2])
    #group_input_004.+x -> group_001.+x
    supercell_atoms.links.new(group_input_004.outputs[8], group_001.inputs[1])
    #group_input_004.+y -> group_001.+y
    supercell_atoms.links.new(group_input_004.outputs[9], group_001.inputs[2])
    #group_input_004.+z -> group_001.+z
    supercell_atoms.links.new(group_input_004.outputs[10], group_001.inputs[3])
    #group_input_004.-x -> group_001.-x
    supercell_atoms.links.new(group_input_004.outputs[11], group_001.inputs[4])
    #group_input_004.-y -> group_001.-y
    supercell_atoms.links.new(group_input_004.outputs[12], group_001.inputs[5])
    #group_input_004.-z -> group_001.-z
    supercell_atoms.links.new(group_input_004.outputs[13], group_001.inputs[6])
    #reroute_008.Output -> switch_006.False
    supercell_atoms.links.new(reroute_008.outputs[0], switch_006.inputs[1])
    #reroute_008.Output -> group_001.Geometry
    supercell_atoms.links.new(reroute_008.outputs[0], group_001.inputs[0])
    #switch_006.Output -> join_geometry.Geometry
    supercell_atoms.links.new(switch_006.outputs[0], join_geometry.inputs[0])
    return supercell_atoms




def make_supercell(list_of_objects, atoms,modifier='GeometryNodes',representation='nodes'):

    if any(atoms.pbc):
        if len(atoms.cell) == 2:
            atoms.cell.append([0,0,50])
    else:
        'print()'
        return
    if representation == 'nodes':
        supercell=supercell_node_group(atoms)
    else:
        supercell=supercell_atoms_node_group(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    for obj in list_of_objects:
        if obj == list_of_objects[0]:
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.modifier_add(type='NODES')
            bpy.context.object.modifiers[modifier].node_group = supercell

        obj.select_set(True)
        
    if len(list_of_objects) > 1:
         bpy.ops.object.make_links_data(type='MODIFIERS')
    return True
