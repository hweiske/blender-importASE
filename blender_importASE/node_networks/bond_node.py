import bpy, mathutils

#initialize bonds node group
def bonds_geometry_node_group():
    if 'BONDS_GEOMETRY' in bpy.data.node_groups:
        bonds = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "BONDS_GEOMETRY")
    else:
        return

    bonds.color_tag = 'NONE'
    bonds.description = ""
    bonds.default_group_node_width = 140
    

    bonds.is_modifier = True

    #bonds interface
    #Socket Geometry
    geometry_socket = bonds.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    #Socket distance
    distance_socket = bonds.interface.new_socket(name = "distance", in_out='INPUT', socket_type = 'NodeSocketFloat')
    distance_socket.default_value = 1.700006127357483
    distance_socket.min_value = -10000.0
    distance_socket.max_value = 10000.0
    distance_socket.subtype = 'NONE'
    distance_socket.attribute_domain = 'POINT'
    distance_socket.description = "arbitrary. Distance minus atomsize (bounding box)"

    #Socket bond radius
    bond_radius_socket = bonds.interface.new_socket(name = "bond radius", in_out='INPUT', socket_type = 'NodeSocketFloat')
    bond_radius_socket.default_value = 0.020000040531158447
    bond_radius_socket.min_value = 0.0
    bond_radius_socket.max_value = 3.4028234663852886e+38
    bond_radius_socket.subtype = 'DISTANCE'
    bond_radius_socket.attribute_domain = 'POINT'

    #Socket bonded_collection
    bonded_collection_socket = bonds.interface.new_socket(name = "bonded_collection", in_out='INPUT', socket_type = 'NodeSocketCollection')
    bonded_collection_socket.attribute_domain = 'POINT'


    #initialize bonds nodes
    #node Frame.005
    frame_005 = bonds.nodes.new("NodeFrame")
    frame_005.label = "000111222"
    frame_005.name = "Frame.005"
    frame_005.label_size = 20
    frame_005.shrink = True

    #node Frame.007
    frame_007 = bonds.nodes.new("NodeFrame")
    frame_007.label = "01230123"
    frame_007.name = "Frame.007"
    frame_007.label_size = 20
    frame_007.shrink = True

    #node Frame.009
    frame_009 = bonds.nodes.new("NodeFrame")
    frame_009.label = "create all points for i,j nested loop"
    frame_009.name = "Frame.009"
    frame_009.label_size = 20
    frame_009.shrink = True

    #node Frame.019
    frame_019 = bonds.nodes.new("NodeFrame")
    frame_019.label = "store positions, colors"
    frame_019.name = "Frame.019"
    frame_019.label_size = 20
    frame_019.shrink = True

    #node Frame.017
    frame_017 = bonds.nodes.new("NodeFrame")
    frame_017.label = "create_loop"
    frame_017.name = "Frame.017"
    frame_017.label_size = 20
    frame_017.shrink = True

    #node Frame
    frame = bonds.nodes.new("NodeFrame")
    frame.label = "instances-to-points"
    frame.name = "Frame"
    frame.label_size = 20
    frame.shrink = True

    #node Frame.002
    frame_002 = bonds.nodes.new("NodeFrame")
    frame_002.label = "back to original position"
    frame_002.name = "Frame.002"
    frame_002.label_size = 20
    frame_002.shrink = True

    #node Frame.003
    frame_003 = bonds.nodes.new("NodeFrame")
    frame_003.label = "calculate bounding box"
    frame_003.name = "Frame.003"
    frame_003.label_size = 20
    frame_003.shrink = True

    #node Frame.013
    frame_013 = bonds.nodes.new("NodeFrame")
    frame_013.label = "set positions"
    frame_013.name = "Frame.013"
    frame_013.label_size = 20
    frame_013.shrink = True

    #node Frame.010
    frame_010 = bonds.nodes.new("NodeFrame")
    frame_010.label = "sample size"
    frame_010.name = "Frame.010"
    frame_010.use_custom_color = True
    frame_010.color = (0.085651695728302, 0.6079999804496765, 0.03432903438806534)
    frame_010.label_size = 20
    frame_010.shrink = True

    #node Frame.006
    frame_006 = bonds.nodes.new("NodeFrame")
    frame_006.label = "get_atomradius"
    frame_006.name = "Frame.006"
    frame_006.use_custom_color = True
    frame_006.color = (0.6079999804496765, 0.6079999804496765, 0.6079999804496765)
    frame_006.label_size = 20
    frame_006.shrink = True

    #node Frame.004
    frame_004 = bonds.nodes.new("NodeFrame")
    frame_004.label = "all indices"
    frame_004.name = "Frame.004"
    frame_004.use_custom_color = True
    frame_004.color = (0.6079999804496765, 0.6079999804496765, 0.6079999804496765)
    frame_004.label_size = 20
    frame_004.shrink = True

    #node Frame.001
    frame_001 = bonds.nodes.new("NodeFrame")
    frame_001.label = "curve profile"
    frame_001.name = "Frame.001"
    frame_001.use_custom_color = True
    frame_001.color = (0.5113166570663452, 0.41898250579833984, 0.6079999804496765)
    frame_001.label_size = 20
    frame_001.shrink = True

    #node Frame.008
    frame_008 = bonds.nodes.new("NodeFrame")
    frame_008.label = "instance curves"
    frame_008.name = "Frame.008"
    frame_008.use_custom_color = True
    frame_008.color = (0.4196130037307739, 0.349015474319458, 0.4980354309082031)
    frame_008.label_size = 20
    frame_008.shrink = True

    #node Frame.018
    frame_018 = bonds.nodes.new("NodeFrame")
    frame_018.label = "CUTOFF"
    frame_018.name = "Frame.018"
    frame_018.label_size = 20
    frame_018.shrink = True

    #node Reroute.006
    reroute_006 = bonds.nodes.new("NodeReroute")
    reroute_006.name = "Reroute.006"
    reroute_006.socket_idname = "NodeSocketGeometry"
    #node Reroute.014
    reroute_014 = bonds.nodes.new("NodeReroute")
    reroute_014.name = "Reroute.014"
    reroute_014.socket_idname = "NodeSocketColor"
    #node Reroute.004
    reroute_004 = bonds.nodes.new("NodeReroute")
    reroute_004.name = "Reroute.004"
    reroute_004.socket_idname = "NodeSocketVector"
    #node Reroute.020
    reroute_020 = bonds.nodes.new("NodeReroute")
    reroute_020.name = "Reroute.020"
    reroute_020.socket_idname = "NodeSocketFloatDistance"
    #node Reroute.007
    reroute_007 = bonds.nodes.new("NodeReroute")
    reroute_007.name = "Reroute.007"
    reroute_007.socket_idname = "NodeSocketVector"
    #node Reroute.027
    reroute_027 = bonds.nodes.new("NodeReroute")
    reroute_027.name = "Reroute.027"
    reroute_027.socket_idname = "NodeSocketFloat"
    #node Delete Geometry
    delete_geometry = bonds.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.name = "Delete Geometry"
    delete_geometry.domain = 'POINT'
    delete_geometry.mode = 'ALL'

    #node Domain Size.001
    domain_size_001 = bonds.nodes.new("GeometryNodeAttributeDomainSize")
    domain_size_001.name = "Domain Size.001"
    domain_size_001.component = 'POINTCLOUD'

    #node Math.013
    math_013 = bonds.nodes.new("ShaderNodeMath")
    math_013.name = "Math.013"
    math_013.operation = 'POWER'
    math_013.use_clamp = False
    #Value_001
    math_013.inputs[1].default_value = 2.0

    #node Points.001
    points_001 = bonds.nodes.new("GeometryNodePoints")
    points_001.name = "Points.001"
    #Position
    points_001.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Radius
    points_001.inputs[2].default_value = 0.10000000149011612

    #node Store Named Attribute
    store_named_attribute = bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute.name = "Store Named Attribute"
    store_named_attribute.hide = True
    store_named_attribute.data_type = 'FLOAT_VECTOR'
    store_named_attribute.domain = 'POINT'
    #Selection
    store_named_attribute.inputs[1].default_value = True
    #Name
    store_named_attribute.inputs[2].default_value = "start_pos"

    #node Store Named Attribute.001
    store_named_attribute_001 = bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_001.name = "Store Named Attribute.001"
    store_named_attribute_001.hide = True
    store_named_attribute_001.data_type = 'FLOAT_VECTOR'
    store_named_attribute_001.domain = 'POINT'
    #Selection
    store_named_attribute_001.inputs[1].default_value = True
    #Name
    store_named_attribute_001.inputs[2].default_value = "end_pos"

    #node Set Position
    set_position = bonds.nodes.new("GeometryNodeSetPosition")
    set_position.name = "Set Position"
    #Selection
    set_position.inputs[1].default_value = True
    #Offset
    set_position.inputs[3].default_value = (0.0, 0.0, 0.0)

    #node Join Geometry
    join_geometry = bonds.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"

    #node Group Output
    group_output = bonds.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"
    group_output.is_active_output = True

    #node Set Shade Smooth
    set_shade_smooth = bonds.nodes.new("GeometryNodeSetShadeSmooth")
    set_shade_smooth.name = "Set Shade Smooth"
    set_shade_smooth.domain = 'FACE'
    #Selection
    set_shade_smooth.inputs[1].default_value = True
    #Shade Smooth
    set_shade_smooth.inputs[2].default_value = True

    #node Math.011
    math_011 = bonds.nodes.new("ShaderNodeMath")
    math_011.name = "Math.011"
    math_011.operation = 'DIVIDE'
    math_011.use_clamp = False

    #node Math.012
    math_012 = bonds.nodes.new("ShaderNodeMath")
    math_012.name = "Math.012"
    math_012.operation = 'MODULO'
    math_012.use_clamp = False
    math_012.inputs[2].hide = True

    #node Index.001
    index_001 = bonds.nodes.new("GeometryNodeInputIndex")
    index_001.name = "Index.001"

    #node Realize Instances.001
    realize_instances_001 = bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_001.name = "Realize Instances.001"
    #Selection
    realize_instances_001.inputs[1].default_value = True
    #Realize All
    realize_instances_001.inputs[2].default_value = True
    #Depth
    realize_instances_001.inputs[3].default_value = 0

    #node Math.004
    math_004 = bonds.nodes.new("ShaderNodeMath")
    math_004.name = "Math.004"
    math_004.operation = 'FLOOR'
    math_004.use_clamp = False

    #node Store Named Attribute.002
    store_named_attribute_002 = bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_002.name = "Store Named Attribute.002"
    store_named_attribute_002.hide = True
    store_named_attribute_002.data_type = 'FLOAT'
    store_named_attribute_002.domain = 'POINT'
    #Selection
    store_named_attribute_002.inputs[1].default_value = True
    #Name
    store_named_attribute_002.inputs[2].default_value = "start_size"

    #node Reroute
    reroute = bonds.nodes.new("NodeReroute")
    reroute.name = "Reroute"
    reroute.socket_idname = "NodeSocketFloat"
    #node Reroute.001
    reroute_001 = bonds.nodes.new("NodeReroute")
    reroute_001.name = "Reroute.001"
    reroute_001.socket_idname = "NodeSocketFloat"
    #node Store Named Attribute.003
    store_named_attribute_003 = bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_003.name = "Store Named Attribute.003"
    store_named_attribute_003.hide = True
    store_named_attribute_003.data_type = 'FLOAT'
    store_named_attribute_003.domain = 'POINT'
    #Selection
    store_named_attribute_003.inputs[1].default_value = True
    #Name
    store_named_attribute_003.inputs[2].default_value = "end_size"

    #node Group Input.002
    group_input_002 = bonds.nodes.new("NodeGroupInput")
    group_input_002.name = "Group Input.002"

    #node Group Input.001
    group_input_001 = bonds.nodes.new("NodeGroupInput")
    group_input_001.label = "bond_radius"
    group_input_001.name = "Group Input.001"

    #node Collection Info
    collection_info = bonds.nodes.new("GeometryNodeCollectionInfo")
    collection_info.name = "Collection Info"
    collection_info.transform_space = 'RELATIVE'
    #Separate Children
    collection_info.inputs[1].default_value = True
    #Reset Children
    collection_info.inputs[2].default_value = False

    #node Translate Instances
    translate_instances = bonds.nodes.new("GeometryNodeTranslateInstances")
    translate_instances.name = "Translate Instances"
    #Selection
    translate_instances.inputs[1].default_value = True
    #Local Space
    translate_instances.inputs[3].default_value = True

    #node Named Attribute.001
    named_attribute_001 = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_001.name = "Named Attribute.001"
    named_attribute_001.data_type = 'FLOAT'
    #Name
    named_attribute_001.inputs[0].default_value = "length"

    #node Combine XYZ
    combine_xyz = bonds.nodes.new("ShaderNodeCombineXYZ")
    combine_xyz.name = "Combine XYZ"

    #node Vector Math.001
    vector_math_001 = bonds.nodes.new("ShaderNodeVectorMath")
    vector_math_001.name = "Vector Math.001"
    vector_math_001.operation = 'SCALE'
    #Scale
    vector_math_001.inputs[3].default_value = 0.5

    #node Bounding Box
    bounding_box = bonds.nodes.new("GeometryNodeBoundBox")
    bounding_box.name = "Bounding Box"

    #node Realize Instances.003
    realize_instances_003 = bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_003.name = "Realize Instances.003"
    #Selection
    realize_instances_003.inputs[1].default_value = True
    #Realize All
    realize_instances_003.inputs[2].default_value = True
    #Depth
    realize_instances_003.inputs[3].default_value = 0

    #node Math.003
    math_003 = bonds.nodes.new("ShaderNodeMath")
    math_003.name = "Math.003"
    math_003.operation = 'MODULO'
    math_003.use_clamp = False
    #Value_001
    math_003.inputs[1].default_value = 12.0

    #node Compare
    compare = bonds.nodes.new("FunctionNodeCompare")
    compare.name = "Compare"
    compare.data_type = 'INT'
    compare.mode = 'ELEMENT'
    compare.operation = 'EQUAL'
    #B_INT
    compare.inputs[3].default_value = 1

    #node Index
    index = bonds.nodes.new("GeometryNodeInputIndex")
    index.name = "Index"

    #node Separate Geometry
    separate_geometry = bonds.nodes.new("GeometryNodeSeparateGeometry")
    separate_geometry.name = "Separate Geometry"
    separate_geometry.domain = 'EDGE'

    #node Mesh to Curve
    mesh_to_curve = bonds.nodes.new("GeometryNodeMeshToCurve")
    mesh_to_curve.name = "Mesh to Curve"
    #Selection
    mesh_to_curve.inputs[1].default_value = True

    #node Store Named Attribute.005
    store_named_attribute_005 = bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_005.name = "Store Named Attribute.005"
    store_named_attribute_005.data_type = 'FLOAT'
    store_named_attribute_005.domain = 'POINT'
    #Selection
    store_named_attribute_005.inputs[1].default_value = True
    #Name
    store_named_attribute_005.inputs[2].default_value = "length"

    #node Curve to Points
    curve_to_points = bonds.nodes.new("GeometryNodeCurveToPoints")
    curve_to_points.name = "Curve to Points"
    curve_to_points.mode = 'COUNT'
    #Count
    curve_to_points.inputs[1].default_value = 1

    #node Instance on Points.001
    instance_on_points_001 = bonds.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_001.name = "Instance on Points.001"
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

    #node Points
    points = bonds.nodes.new("GeometryNodePoints")
    points.name = "Points"
    points.hide = True
    #Count
    points.inputs[0].default_value = 1
    #Position
    points.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Radius
    points.inputs[2].default_value = 0.3999999761581421

    #node Spline Length
    spline_length = bonds.nodes.new("GeometryNodeSplineLength")
    spline_length.name = "Spline Length"

    #node Switch
    switch = bonds.nodes.new("GeometryNodeSwitch")
    switch.name = "Switch"
    switch.input_type = 'VECTOR'

    #node Index.003
    index_003 = bonds.nodes.new("GeometryNodeInputIndex")
    index_003.name = "Index.003"

    #node Math.014
    math_014 = bonds.nodes.new("ShaderNodeMath")
    math_014.name = "Math.014"
    math_014.operation = 'MODULO'
    math_014.use_clamp = False
    #Value_001
    math_014.inputs[1].default_value = 2.0

    #node Named Attribute
    named_attribute = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute.name = "Named Attribute"
    named_attribute.data_type = 'FLOAT_VECTOR'
    #Name
    named_attribute.inputs[0].default_value = "start_pos"

    #node Named Attribute.008
    named_attribute_008 = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_008.name = "Named Attribute.008"
    named_attribute_008.data_type = 'FLOAT_VECTOR'
    #Name
    named_attribute_008.inputs[0].default_value = "end_pos"

    #node Named Attribute.009
    named_attribute_009 = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_009.name = "Named Attribute.009"
    named_attribute_009.hide = True
    named_attribute_009.data_type = 'FLOAT'
    #Name
    named_attribute_009.inputs[0].default_value = "start_size"

    #node Named Attribute.011
    named_attribute_011 = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_011.name = "Named Attribute.011"
    named_attribute_011.hide = True
    named_attribute_011.data_type = 'FLOAT'
    #Name
    named_attribute_011.inputs[0].default_value = "end_size"

    #node Reroute.005
    reroute_005 = bonds.nodes.new("NodeReroute")
    reroute_005.name = "Reroute.005"
    reroute_005.socket_idname = "NodeSocketGeometry"
    #node Sample Index.004
    sample_index_004 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_004.name = "Sample Index.004"
    sample_index_004.clamp = False
    sample_index_004.data_type = 'FLOAT_VECTOR'
    sample_index_004.domain = 'POINT'
    #Value
    sample_index_004.inputs[1].default_value = (0.0, 0.0, 0.0)

    #node Sample Index.005
    sample_index_005 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_005.name = "Sample Index.005"
    sample_index_005.clamp = False
    sample_index_005.data_type = 'FLOAT_VECTOR'
    sample_index_005.domain = 'POINT'
    #Value
    sample_index_005.inputs[1].default_value = (0.0, 0.0, 0.0)

    #node Reroute.016
    reroute_016 = bonds.nodes.new("NodeReroute")
    reroute_016.name = "Reroute.016"
    reroute_016.socket_idname = "NodeSocketFloat"
    #node Reroute.017
    reroute_017 = bonds.nodes.new("NodeReroute")
    reroute_017.name = "Reroute.017"
    reroute_017.socket_idname = "NodeSocketFloat"
    #node Sample Index.006
    sample_index_006 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_006.name = "Sample Index.006"
    sample_index_006.clamp = False
    sample_index_006.data_type = 'FLOAT'
    sample_index_006.domain = 'POINT'

    #node Reroute.010
    reroute_010 = bonds.nodes.new("NodeReroute")
    reroute_010.name = "Reroute.010"
    reroute_010.socket_idname = "NodeSocketFloat"
    #node Reroute.011
    reroute_011 = bonds.nodes.new("NodeReroute")
    reroute_011.name = "Reroute.011"
    reroute_011.socket_idname = "NodeSocketFloat"
    #node Sample Index.007
    sample_index_007 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_007.name = "Sample Index.007"
    sample_index_007.clamp = False
    sample_index_007.data_type = 'FLOAT'
    sample_index_007.domain = 'POINT'

    #node Named Attribute.002
    named_attribute_002 = bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_002.name = "Named Attribute.002"
    named_attribute_002.data_type = 'FLOAT'
    #Name
    named_attribute_002.inputs[0].default_value = "length"

    #node Reroute.018
    reroute_018 = bonds.nodes.new("NodeReroute")
    reroute_018.name = "Reroute.018"
    reroute_018.socket_idname = "NodeSocketGeometry"
    #node Position
    position = bonds.nodes.new("GeometryNodeInputPosition")
    position.name = "Position"

    #node Sample Index
    sample_index = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index.name = "Sample Index"
    sample_index.clamp = False
    sample_index.data_type = 'FLOAT_VECTOR'
    sample_index.domain = 'POINT'

    #node Sample Index.001
    sample_index_001 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_001.name = "Sample Index.001"
    sample_index_001.clamp = False
    sample_index_001.data_type = 'FLOAT_VECTOR'
    sample_index_001.domain = 'POINT'

    #node Sample Index.003
    sample_index_003 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_003.name = "Sample Index.003"
    sample_index_003.clamp = False
    sample_index_003.data_type = 'FLOAT_VECTOR'
    sample_index_003.domain = 'POINT'
    #Value
    sample_index_003.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Index
    sample_index_003.inputs[2].default_value = 0

    #node Reroute.009
    reroute_009 = bonds.nodes.new("NodeReroute")
    reroute_009.name = "Reroute.009"
    reroute_009.socket_idname = "NodeSocketFloat"
    #node Reroute.008
    reroute_008 = bonds.nodes.new("NodeReroute")
    reroute_008.name = "Reroute.008"
    reroute_008.socket_idname = "NodeSocketFloat"
    #node Reroute.013
    reroute_013 = bonds.nodes.new("NodeReroute")
    reroute_013.name = "Reroute.013"
    reroute_013.socket_idname = "NodeSocketFloat"
    #node Reroute.003
    reroute_003 = bonds.nodes.new("NodeReroute")
    reroute_003.name = "Reroute.003"
    reroute_003.socket_idname = "NodeSocketGeometry"
    #node Sample Index.002
    sample_index_002 = bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_002.name = "Sample Index.002"
    sample_index_002.clamp = False
    sample_index_002.data_type = 'FLOAT_VECTOR'
    sample_index_002.domain = 'POINT'
    #Value
    sample_index_002.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Index
    sample_index_002.inputs[2].default_value = 0

    #node Reroute.002
    reroute_002 = bonds.nodes.new("NodeReroute")
    reroute_002.name = "Reroute.002"
    reroute_002.socket_idname = "NodeSocketGeometry"
    #node Reroute.012
    reroute_012 = bonds.nodes.new("NodeReroute")
    reroute_012.name = "Reroute.012"
    reroute_012.socket_idname = "NodeSocketFloat"
    #node Set Material.006
    set_material_006 = bonds.nodes.new("GeometryNodeSetMaterial")
    set_material_006.name = "Set Material.006"
    #Selection
    set_material_006.inputs[1].default_value = True
    if "BONDS_MAT" in bpy.data.materials:
        set_material_006.inputs[2].default_value = bpy.data.materials["BONDS_MAT"]

    #node Switch.005
    switch_005 = bonds.nodes.new("GeometryNodeSwitch")
    switch_005.name = "Switch.005"
    switch_005.input_type = 'GEOMETRY'
    #Switch
    switch_005.inputs[0].default_value = False

    #node Curve to Mesh
    curve_to_mesh = bonds.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.name = "Curve to Mesh"
    #Fill Caps
    curve_to_mesh.inputs[2].default_value = False

    #node Reroute.021
    reroute_021 = bonds.nodes.new("NodeReroute")
    reroute_021.name = "Reroute.021"
    reroute_021.socket_idname = "NodeSocketGeometry"
    #node Curve Circle
    curve_circle = bonds.nodes.new("GeometryNodeCurvePrimitiveCircle")
    curve_circle.name = "Curve Circle"
    curve_circle.hide = True
    curve_circle.mode = 'RADIUS'
    #Resolution
    curve_circle.inputs[0].default_value = 32
    #Radius
    curve_circle.inputs[4].default_value = 1.0

    #node Curve to Mesh.001
    curve_to_mesh_001 = bonds.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh_001.name = "Curve to Mesh.001"
    #Fill Caps
    curve_to_mesh_001.inputs[2].default_value = False

    #node Group Input.003
    group_input_003 = bonds.nodes.new("NodeGroupInput")
    group_input_003.name = "Group Input.003"
    group_input_003.hide = True

    #node Curve Circle.001
    curve_circle_001 = bonds.nodes.new("GeometryNodeCurvePrimitiveCircle")
    curve_circle_001.name = "Curve Circle.001"
    curve_circle_001.mode = 'RADIUS'
    #Resolution
    curve_circle_001.inputs[0].default_value = 16

    #node Set Curve Radius
    set_curve_radius = bonds.nodes.new("GeometryNodeSetCurveRadius")
    set_curve_radius.name = "Set Curve Radius"
    #Selection
    set_curve_radius.inputs[1].default_value = True
    #Radius
    set_curve_radius.inputs[2].default_value = 0.004999999888241291

    #node Reroute.015
    reroute_015 = bonds.nodes.new("NodeReroute")
    reroute_015.name = "Reroute.015"
    reroute_015.socket_idname = "NodeSocketVector"
    #node Realize Instances
    realize_instances = bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances.name = "Realize Instances"
    #Selection
    realize_instances.inputs[1].default_value = True
    #Realize All
    realize_instances.inputs[2].default_value = True
    #Depth
    realize_instances.inputs[3].default_value = 0

    #node Curve Line
    curve_line = bonds.nodes.new("GeometryNodeCurvePrimitiveLine")
    curve_line.name = "Curve Line"
    curve_line.hide = True
    curve_line.mode = 'POINTS'
    #Start
    curve_line.inputs[0].default_value = (0.0, 0.0, 0.0)
    #End
    curve_line.inputs[1].default_value = (0.0, 0.0, 1.0)

    #node Instance on Points
    instance_on_points = bonds.nodes.new("GeometryNodeInstanceOnPoints")
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

    #node Compare.005
    compare_005 = bonds.nodes.new("FunctionNodeCompare")
    compare_005.name = "Compare.005"
    compare_005.data_type = 'FLOAT'
    compare_005.mode = 'ELEMENT'
    compare_005.operation = 'GREATER_THAN'

    #node Boolean Math
    boolean_math = bonds.nodes.new("FunctionNodeBooleanMath")
    boolean_math.name = "Boolean Math"
    boolean_math.operation = 'OR'

    #node Vector Math
    vector_math = bonds.nodes.new("ShaderNodeVectorMath")
    vector_math.name = "Vector Math"
    vector_math.operation = 'DISTANCE'

    #node Math.002
    math_002 = bonds.nodes.new("ShaderNodeMath")
    math_002.name = "Math.002"
    math_002.operation = 'SUBTRACT'
    math_002.use_clamp = False

    #node Compare.006
    compare_006 = bonds.nodes.new("FunctionNodeCompare")
    compare_006.name = "Compare.006"
    compare_006.data_type = 'INT'
    compare_006.mode = 'ELEMENT'
    compare_006.operation = 'EQUAL'

    #node Math
    math = bonds.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = 1.0

    #node Math.001
    math_001 = bonds.nodes.new("ShaderNodeMath")
    math_001.name = "Math.001"
    math_001.operation = 'ADD'
    math_001.use_clamp = False

    #node Math.006
    math_006 = bonds.nodes.new("ShaderNodeMath")
    math_006.name = "Math.006"
    math_006.operation = 'MULTIPLY'
    math_006.use_clamp = False
    #Value_001
    math_006.inputs[1].default_value = 1.0000001192092896




    #Set parents
    domain_size_001.parent = frame_009
    math_013.parent = frame_009
    points_001.parent = frame_009
    store_named_attribute.parent = frame_019
    store_named_attribute_001.parent = frame_019
    math_011.parent = frame_017
    math_012.parent = frame_017
    index_001.parent = frame_017
    realize_instances_001.parent = frame
    math_004.parent = frame_017
    store_named_attribute_002.parent = frame_019
    store_named_attribute_003.parent = frame_019
    named_attribute_001.parent = frame_002
    combine_xyz.parent = frame_002
    vector_math_001.parent = frame_002
    bounding_box.parent = frame_003
    realize_instances_003.parent = frame_003
    math_003.parent = frame_003
    compare.parent = frame_003
    index.parent = frame_003
    separate_geometry.parent = frame_003
    mesh_to_curve.parent = frame_003
    store_named_attribute_005.parent = frame_003
    curve_to_points.parent = frame_003
    instance_on_points_001.parent = frame_003
    points.parent = frame_003
    spline_length.parent = frame_003
    switch.parent = frame_013
    index_003.parent = frame_013
    math_014.parent = frame_013
    named_attribute.parent = frame_013
    named_attribute_008.parent = frame_013
    named_attribute_009.parent = frame_010
    named_attribute_011.parent = frame_010
    reroute_005.parent = frame_006
    sample_index_004.parent = frame_006
    sample_index_005.parent = frame_006
    reroute_016.parent = frame_006
    reroute_017.parent = frame_006
    sample_index_006.parent = frame_006
    reroute_010.parent = frame_006
    reroute_011.parent = frame_006
    sample_index_007.parent = frame_006
    named_attribute_002.parent = frame_006
    reroute_018.parent = frame_006
    position.parent = frame_004
    sample_index.parent = frame_004
    sample_index_001.parent = frame_004
    sample_index_003.parent = frame_004
    reroute_009.parent = frame_004
    reroute_008.parent = frame_004
    reroute_013.parent = frame_004
    reroute_003.parent = frame_004
    sample_index_002.parent = frame_004
    reroute_002.parent = frame_004
    reroute_012.parent = frame_004
    set_material_006.parent = frame_001
    switch_005.parent = frame_001
    curve_to_mesh.parent = frame_001
    reroute_021.parent = frame_001
    curve_circle.parent = frame_001
    curve_to_mesh_001.parent = frame_001
    group_input_003.parent = frame_001
    curve_circle_001.parent = frame_001
    set_curve_radius.parent = frame_001
    reroute_015.parent = frame_008
    realize_instances.parent = frame_008
    curve_line.parent = frame_008
    instance_on_points.parent = frame_008
    compare_005.parent = frame_018
    boolean_math.parent = frame_018
    vector_math.parent = frame_018
    math_002.parent = frame_018
    compare_006.parent = frame_018
    math.parent = frame_018
    math_001.parent = frame_018
    math_006.parent = frame_018

    #Set locations
    frame_005.location = (-3140.510498046875, -1904.2864990234375)
    frame_007.location = (-3135.767578125, -2116.841552734375)
    frame_009.location = (-3181.666748046875, -1455.8333740234375)
    frame_019.location = (-1459.0, -1335.6666259765625)
    frame_017.location = (-3037.0, -1758.5)
    frame.location = (-4235.0, -1564.5)
    frame_002.location = (-5349.66650390625, -2105.833251953125)
    frame_003.location = (-6175.66650390625, -1492.5)
    frame_013.location = (-821.6666870117188, -1891.1666259765625)
    frame_010.location = (-2375.666748046875, -1963.6666259765625)
    frame_006.location = (-3847.938232421875, -2869.833251953125)
    frame_004.location = (-3824.549560546875, -2313.833251953125)
    frame_001.location = (1074.3333740234375, -927.1666870117188)
    frame_008.location = (-271.6666564941406, -1293.8333740234375)
    frame_018.location = (-2172.333251953125, -1689.1666259765625)
    reroute_006.location = (2803.35205078125, -1283.796875)
    reroute_014.location = (-215.2777099609375, -1474.3909912109375)
    reroute_004.location = (-2375.80859375, -1578.4334716796875)
    reroute_020.location = (-5490.40625, -1094.441162109375)
    reroute_007.location = (-2400.44189453125, -1548.634765625)
    reroute_027.location = (-140.0753173828125, -1507.0865478515625)
    delete_geometry.location = (-851.4942016601562, -1383.24267578125)
    domain_size_001.location = (29.085693359375, -128.757080078125)
    math_013.location = (240.842041015625, -39.6812744140625)
    points_001.location = (443.87060546875, -83.1571044921875)
    store_named_attribute.location = (30.86669921875, -44.3043212890625)
    store_named_attribute_001.location = (34.7431640625, -76.582275390625)
    set_position.location = (450.714111328125, -1319.991943359375)
    join_geometry.location = (2937.53759765625, -1224.6448974609375)
    group_output.location = (3741.489990234375, -1239.6837158203125)
    set_shade_smooth.location = (3344.00537109375, -1165.5565185546875)
    math_011.location = (28.739501953125, -39.683349609375)
    math_012.location = (178.139404296875, -349.00146484375)
    index_001.location = (30.088134765625, -228.886474609375)
    realize_instances_001.location = (28.9775390625, -39.3519287109375)
    math_004.location = (195.13427734375, -102.0716552734375)
    store_named_attribute_002.location = (29.869384765625, -116.4759521484375)
    reroute.location = (-2375.56982421875, -1597.0758056640625)
    reroute_001.location = (-2340.423095703125, -1622.2178955078125)
    store_named_attribute_003.location = (29.3232421875, -158.4775390625)
    group_input_002.location = (-2435.63623046875, -1784.1600341796875)
    group_input_001.location = (-6890.75732421875, -1679.4515380859375)
    collection_info.location = (-6461.70654296875, -1858.63037109375)
    translate_instances.location = (-4389.6572265625, -1704.3221435546875)
    named_attribute_001.location = (28.955078125, -89.60107421875)
    combine_xyz.location = (634.04541015625, -71.77880859375)
    vector_math_001.location = (794.76708984375, -39.33642578125)
    bounding_box.location = (69.5498046875, -192.764892578125)
    realize_instances_003.location = (264.73095703125, -209.09619140625)
    math_003.location = (172.32470703125, -317.693603515625)
    compare.location = (360.6064453125, -322.438720703125)
    index.location = (28.9951171875, -343.0909423828125)
    separate_geometry.location = (545.8583984375, -232.15087890625)
    mesh_to_curve.location = (587.80908203125, -92.4969482421875)
    store_named_attribute_005.location = (1019.888671875, -77.61767578125)
    curve_to_points.location = (1220.6123046875, -39.664794921875)
    instance_on_points_001.location = (1561.91162109375, -60.552001953125)
    points.location = (1415.45361328125, -160.5032958984375)
    spline_length.location = (552.41748046875, -374.9117431640625)
    switch.location = (420.64202880859375, -144.48779296875)
    index_003.location = (29.05926513671875, -122.7100830078125)
    math_014.location = (225.81976318359375, -39.693603515625)
    named_attribute.location = (225.99456787109375, -211.3751220703125)
    named_attribute_008.location = (222.96746826171875, -365.8988037109375)
    named_attribute_009.location = (29.7978515625, -44.5545654296875)
    named_attribute_011.location = (28.990966796875, -82.6453857421875)
    reroute_005.location = (33.833251953125, -346.431884765625)
    sample_index_004.location = (986.665771484375, -39.6650390625)
    sample_index_005.location = (984.484130859375, -287.500244140625)
    reroute_016.location = (389.083251953125, -377.364013671875)
    reroute_017.location = (469.48779296875, -164.764892578125)
    sample_index_006.location = (984.484130859375, -287.500244140625)
    reroute_010.location = (946.107666015625, -167.074462890625)
    reroute_011.location = (934.743896484375, -368.2568359375)
    sample_index_007.location = (986.665771484375, -39.6650390625)
    named_attribute_002.location = (592.584228515625, -395.543701171875)
    reroute_018.location = (845.536376953125, -346.107421875)
    position.location = (642.7607421875, -401.120361328125)
    sample_index.location = (986.665771484375, -39.6015625)
    sample_index_001.location = (984.484130859375, -287.436767578125)
    sample_index_003.location = (984.484130859375, -287.436767578125)
    reroute_009.location = (946.107666015625, -167.010986328125)
    reroute_008.location = (934.743896484375, -368.193359375)
    reroute_013.location = (469.48779296875, -164.701416015625)
    reroute_003.location = (845.536376953125, -346.0439453125)
    sample_index_002.location = (986.665771484375, -39.6015625)
    reroute_002.location = (33.833251953125, -346.368408203125)
    reroute_012.location = (389.083251953125, -377.300537109375)
    set_material_006.location = (1025.8441162109375, -326.07916259765625)
    switch_005.location = (805.7330322265625, -242.70440673828125)
    curve_to_mesh.location = (585.5738525390625, -349.49139404296875)
    reroute_021.location = (376.4156494140625, -339.73040771484375)
    curve_circle.location = (427.6524658203125, -418.33843994140625)
    curve_to_mesh_001.location = (590.5220947265625, -220.72076416015625)
    group_input_003.location = (29.2069091796875, -130.22064208984375)
    curve_circle_001.location = (240.4832763671875, -39.6986083984375)
    set_curve_radius.location = (419.3665771484375, -290.44024658203125)
    reroute_015.location = (609.2327880859375, -123.7274169921875)
    realize_instances.location = (368.2862854003906, -39.2337646484375)
    curve_line.location = (28.674957275390625, -177.2742919921875)
    instance_on_points.location = (198.92251586914062, -54.5875244140625)
    compare_005.location = (720.173828125, -82.3499755859375)
    boolean_math.location = (982.4013671875, -170.4898681640625)
    vector_math.location = (229.5009765625, -39.308837890625)
    math_002.location = (425.8466796875, -42.927978515625)
    compare_006.location = (393.8809814453125, -480.2308349609375)
    math.location = (474.9549560546875, -204.21923828125)
    math_001.location = (29.0693359375, -282.761474609375)
    math_006.location = (193.4967041015625, -284.9173583984375)

    #Set dimensions
    frame_005.width, frame_005.height = 169.232177734375, 50.3642578125
    frame_007.width, frame_007.height = 198.0, 47.47119140625
    frame_009.width, frame_009.height = 612.666748046875, 285.1666259765625
    frame_019.width, frame_019.height = 204.0, 211.5
    frame_017.width, frame_017.height = 364.0, 519.166748046875
    frame.width, frame.height = 198.0, 202.5
    frame_002.width, frame_002.height = 964.0, 234.5
    frame_003.width, frame_003.height = 1724.7099609375, 492.5
    frame_013.width, frame_013.height = 589.3333740234375, 511.1666259765625
    frame_010.width, frame_010.height = 198.666748046875, 135.5001220703125
    frame_006.width, frame_006.height = 1155.60498046875, 540.5
    frame_004.width, frame_004.height = 1155.549560546875, 501.166748046875
    frame_001.width, frame_001.height = 1194.6666259765625, 491.83331298828125
    frame_008.width, frame_008.height = 643.066162109375, 399.833251953125
    frame_018.width, frame_018.height = 1151.333251953125, 650.5001220703125
    reroute_006.width, reroute_006.height = 14.5, 100.0
    reroute_014.width, reroute_014.height = 14.5, 100.0
    reroute_004.width, reroute_004.height = 14.5, 100.0
    reroute_020.width, reroute_020.height = 14.5, 100.0
    reroute_007.width, reroute_007.height = 14.5, 100.0
    reroute_027.width, reroute_027.height = 14.5, 100.0
    delete_geometry.width, delete_geometry.height = 140.0, 100.0
    domain_size_001.width, domain_size_001.height = 140.0, 100.0
    math_013.width, math_013.height = 140.0, 100.0
    points_001.width, points_001.height = 140.0, 100.0
    store_named_attribute.width, store_named_attribute.height = 140.0, 100.0
    store_named_attribute_001.width, store_named_attribute_001.height = 140.0, 100.0
    set_position.width, set_position.height = 126.54090881347656, 100.0
    join_geometry.width, join_geometry.height = 140.0, 100.0
    group_output.width, group_output.height = 140.0, 100.0
    set_shade_smooth.width, set_shade_smooth.height = 140.0, 100.0
    math_011.width, math_011.height = 140.0, 100.0
    math_012.width, math_012.height = 140.0, 100.0
    index_001.width, index_001.height = 140.0, 100.0
    realize_instances_001.width, realize_instances_001.height = 140.0, 100.0
    math_004.width, math_004.height = 140.0, 100.0
    store_named_attribute_002.width, store_named_attribute_002.height = 140.0, 100.0
    reroute.width, reroute.height = 14.5, 100.0
    reroute_001.width, reroute_001.height = 14.5, 100.0
    store_named_attribute_003.width, store_named_attribute_003.height = 140.0, 100.0
    group_input_002.width, group_input_002.height = 140.0, 100.0
    group_input_001.width, group_input_001.height = 140.0, 100.0
    collection_info.width, collection_info.height = 140.0, 100.0
    translate_instances.width, translate_instances.height = 140.0, 100.0
    named_attribute_001.width, named_attribute_001.height = 140.0, 100.0
    combine_xyz.width, combine_xyz.height = 140.0, 100.0
    vector_math_001.width, vector_math_001.height = 140.0, 100.0
    bounding_box.width, bounding_box.height = 115.14681243896484, 100.0
    realize_instances_003.width, realize_instances_003.height = 140.0, 100.0
    math_003.width, math_003.height = 140.0, 100.0
    compare.width, compare.height = 140.0, 100.0
    index.width, index.height = 140.0, 100.0
    separate_geometry.width, separate_geometry.height = 140.0, 100.0
    mesh_to_curve.width, mesh_to_curve.height = 140.0, 100.0
    store_named_attribute_005.width, store_named_attribute_005.height = 140.0, 100.0
    curve_to_points.width, curve_to_points.height = 140.0, 100.0
    instance_on_points_001.width, instance_on_points_001.height = 134.04345703125, 100.0
    points.width, points.height = 140.0, 100.0
    spline_length.width, spline_length.height = 140.0, 100.0
    switch.width, switch.height = 140.0, 100.0
    index_003.width, index_003.height = 140.0, 100.0
    math_014.width, math_014.height = 140.0, 100.0
    named_attribute.width, named_attribute.height = 140.0, 100.0
    named_attribute_008.width, named_attribute_008.height = 140.0, 100.0
    named_attribute_009.width, named_attribute_009.height = 140.0, 100.0
    named_attribute_011.width, named_attribute_011.height = 140.0, 100.0
    reroute_005.width, reroute_005.height = 14.5, 100.0
    sample_index_004.width, sample_index_004.height = 140.0, 100.0
    sample_index_005.width, sample_index_005.height = 140.0, 100.0
    reroute_016.width, reroute_016.height = 14.5, 100.0
    reroute_017.width, reroute_017.height = 14.5, 100.0
    sample_index_006.width, sample_index_006.height = 140.0, 100.0
    reroute_010.width, reroute_010.height = 14.5, 100.0
    reroute_011.width, reroute_011.height = 14.5, 100.0
    sample_index_007.width, sample_index_007.height = 140.0, 100.0
    named_attribute_002.width, named_attribute_002.height = 140.0, 100.0
    reroute_018.width, reroute_018.height = 14.5, 100.0
    position.width, position.height = 140.0, 100.0
    sample_index.width, sample_index.height = 140.0, 100.0
    sample_index_001.width, sample_index_001.height = 140.0, 100.0
    sample_index_003.width, sample_index_003.height = 140.0, 100.0
    reroute_009.width, reroute_009.height = 14.5, 100.0
    reroute_008.width, reroute_008.height = 14.5, 100.0
    reroute_013.width, reroute_013.height = 14.5, 100.0
    reroute_003.width, reroute_003.height = 14.5, 100.0
    sample_index_002.width, sample_index_002.height = 140.0, 100.0
    reroute_002.width, reroute_002.height = 14.5, 100.0
    reroute_012.width, reroute_012.height = 14.5, 100.0
    set_material_006.width, set_material_006.height = 140.0, 100.0
    switch_005.width, switch_005.height = 140.0, 100.0
    curve_to_mesh.width, curve_to_mesh.height = 140.0, 100.0
    reroute_021.width, reroute_021.height = 14.5, 100.0
    curve_circle.width, curve_circle.height = 140.0, 100.0
    curve_to_mesh_001.width, curve_to_mesh_001.height = 140.0, 100.0
    group_input_003.width, group_input_003.height = 140.0, 100.0
    curve_circle_001.width, curve_circle_001.height = 140.0, 100.0
    set_curve_radius.width, set_curve_radius.height = 140.0, 100.0
    reroute_015.width, reroute_015.height = 14.5, 100.0
    realize_instances.width, realize_instances.height = 140.0, 100.0
    curve_line.width, curve_line.height = 140.0, 100.0
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    compare_005.width, compare_005.height = 151.46044921875, 100.0
    boolean_math.width, boolean_math.height = 140.0, 100.0
    vector_math.width, vector_math.height = 140.0, 100.0
    math_002.width, math_002.height = 100.0, 100.0
    compare_006.width, compare_006.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0
    math_001.width, math_001.height = 140.0, 100.0
    math_006.width, math_006.height = 140.0, 100.0

    #initialize bonds links
    #set_shade_smooth.Geometry -> group_output.Geometry
    bonds.links.new(set_shade_smooth.outputs[0], group_output.inputs[0])
    #named_attribute.Attribute -> switch.False
    bonds.links.new(named_attribute.outputs[0], switch.inputs[1])
    #compare_006.Result -> boolean_math.Boolean
    bonds.links.new(compare_006.outputs[0], boolean_math.inputs[1])
    #store_named_attribute_002.Geometry -> store_named_attribute_003.Geometry
    bonds.links.new(store_named_attribute_002.outputs[0], store_named_attribute_003.inputs[0])
    #points_001.Points -> store_named_attribute.Geometry
    bonds.links.new(points_001.outputs[0], store_named_attribute.inputs[0])
    #store_named_attribute.Geometry -> store_named_attribute_001.Geometry
    bonds.links.new(store_named_attribute.outputs[0], store_named_attribute_001.inputs[0])
    #curve_circle.Curve -> curve_to_mesh.Profile Curve
    bonds.links.new(curve_circle.outputs[0], curve_to_mesh.inputs[1])
    #sample_index.Value -> reroute_007.Input
    bonds.links.new(sample_index.outputs[0], reroute_007.inputs[0])
    #sample_index_001.Value -> reroute_004.Input
    bonds.links.new(sample_index_001.outputs[0], reroute_004.inputs[0])
    #reroute_007.Output -> store_named_attribute.Value
    bonds.links.new(reroute_007.outputs[0], store_named_attribute.inputs[3])
    #reroute_004.Output -> store_named_attribute_001.Value
    bonds.links.new(reroute_004.outputs[0], store_named_attribute_001.inputs[3])
    #named_attribute_008.Attribute -> switch.True
    bonds.links.new(named_attribute_008.outputs[0], switch.inputs[2])
    #curve_line.Curve -> instance_on_points.Instance
    bonds.links.new(curve_line.outputs[0], instance_on_points.inputs[2])
    #instance_on_points.Instances -> realize_instances.Geometry
    bonds.links.new(instance_on_points.outputs[0], realize_instances.inputs[0])
    #math_014.Value -> switch.Switch
    bonds.links.new(math_014.outputs[0], switch.inputs[0])
    #math_012.Value -> reroute_012.Input
    bonds.links.new(math_012.outputs[0], reroute_012.inputs[0])
    #index_003.Index -> math_014.Value
    bonds.links.new(index_003.outputs[0], math_014.inputs[0])
    #math_004.Value -> reroute_013.Input
    bonds.links.new(math_004.outputs[0], reroute_013.inputs[0])
    #switch.Output -> reroute_015.Input
    bonds.links.new(switch.outputs[0], reroute_015.inputs[0])
    #position.Position -> sample_index.Value
    bonds.links.new(position.outputs[0], sample_index.inputs[1])
    #reroute_015.Output -> set_position.Position
    bonds.links.new(reroute_015.outputs[0], set_position.inputs[2])
    #reroute_009.Output -> sample_index.Index
    bonds.links.new(reroute_009.outputs[0], sample_index.inputs[2])
    #reroute_008.Output -> sample_index_001.Index
    bonds.links.new(reroute_008.outputs[0], sample_index_001.inputs[2])
    #reroute_003.Output -> sample_index.Geometry
    bonds.links.new(reroute_003.outputs[0], sample_index.inputs[0])
    #switch_005.Output -> set_material_006.Geometry
    bonds.links.new(switch_005.outputs[0], set_material_006.inputs[0])
    #reroute_012.Output -> reroute_008.Input
    bonds.links.new(reroute_012.outputs[0], reroute_008.inputs[0])
    #reroute_013.Output -> reroute_009.Input
    bonds.links.new(reroute_013.outputs[0], reroute_009.inputs[0])
    #domain_size_001.Point Count -> math_013.Value
    bonds.links.new(domain_size_001.outputs[0], math_013.inputs[0])
    #math_011.Value -> compare_006.A
    bonds.links.new(math_011.outputs[0], compare_006.inputs[2])
    #domain_size_001.Point Count -> math_011.Value
    bonds.links.new(domain_size_001.outputs[0], math_011.inputs[1])
    #index_001.Index -> math_011.Value
    bonds.links.new(index_001.outputs[0], math_011.inputs[0])
    #math_011.Value -> math_004.Value
    bonds.links.new(math_011.outputs[0], math_004.inputs[0])
    #index_001.Index -> math_012.Value
    bonds.links.new(index_001.outputs[0], math_012.inputs[0])
    #domain_size_001.Point Count -> math_012.Value
    bonds.links.new(domain_size_001.outputs[0], math_012.inputs[1])
    #math_012.Value -> compare_006.B
    bonds.links.new(math_012.outputs[0], compare_006.inputs[3])
    #reroute_003.Output -> sample_index_001.Geometry
    bonds.links.new(reroute_003.outputs[0], sample_index_001.inputs[0])
    #math_013.Value -> points_001.Count
    bonds.links.new(math_013.outputs[0], points_001.inputs[0])
    #reroute_002.Output -> reroute_003.Input
    bonds.links.new(reroute_002.outputs[0], reroute_003.inputs[0])
    #set_material_006.Geometry -> reroute_006.Input
    bonds.links.new(set_material_006.outputs[0], reroute_006.inputs[0])
    #store_named_attribute_001.Geometry -> store_named_attribute_002.Geometry
    bonds.links.new(store_named_attribute_001.outputs[0], store_named_attribute_002.inputs[0])
    #position.Position -> sample_index_001.Value
    bonds.links.new(position.outputs[0], sample_index_001.inputs[1])
    #realize_instances.Geometry -> set_position.Geometry
    bonds.links.new(realize_instances.outputs[0], set_position.inputs[0])
    #reroute_021.Output -> set_curve_radius.Curve
    bonds.links.new(reroute_021.outputs[0], set_curve_radius.inputs[0])
    #set_curve_radius.Curve -> curve_to_mesh.Curve
    bonds.links.new(set_curve_radius.outputs[0], curve_to_mesh.inputs[0])
    #curve_to_mesh.Mesh -> switch_005.True
    bonds.links.new(curve_to_mesh.outputs[0], switch_005.inputs[2])
    #reroute_021.Output -> curve_to_mesh_001.Curve
    bonds.links.new(reroute_021.outputs[0], curve_to_mesh_001.inputs[0])
    #curve_circle_001.Curve -> curve_to_mesh_001.Profile Curve
    bonds.links.new(curve_circle_001.outputs[0], curve_to_mesh_001.inputs[1])
    #curve_to_mesh_001.Mesh -> switch_005.False
    bonds.links.new(curve_to_mesh_001.outputs[0], switch_005.inputs[1])
    #math_002.Value -> compare_005.A
    bonds.links.new(math_002.outputs[0], compare_005.inputs[0])
    #reroute_004.Output -> vector_math.Vector
    bonds.links.new(reroute_004.outputs[0], vector_math.inputs[1])
    #reroute_007.Output -> vector_math.Vector
    bonds.links.new(reroute_007.outputs[0], vector_math.inputs[0])
    #math.Value -> compare_005.B
    bonds.links.new(math.outputs[0], compare_005.inputs[1])
    #group_input_003.bond radius -> curve_circle_001.Radius
    bonds.links.new(group_input_003.outputs[1], curve_circle_001.inputs[4])
    #boolean_math.Boolean -> delete_geometry.Selection
    bonds.links.new(boolean_math.outputs[0], delete_geometry.inputs[1])
    #delete_geometry.Geometry -> instance_on_points.Points
    bonds.links.new(delete_geometry.outputs[0], instance_on_points.inputs[0])
    #compare_005.Result -> boolean_math.Boolean
    bonds.links.new(compare_005.outputs[0], boolean_math.inputs[0])
    #store_named_attribute_003.Geometry -> delete_geometry.Geometry
    bonds.links.new(store_named_attribute_003.outputs[0], delete_geometry.inputs[0])
    #vector_math.Value -> math_002.Value
    bonds.links.new(vector_math.outputs[1], math_002.inputs[0])
    #realize_instances_001.Geometry -> domain_size_001.Geometry
    bonds.links.new(realize_instances_001.outputs[0], domain_size_001.inputs[0])
    #realize_instances_001.Geometry -> reroute_002.Input
    bonds.links.new(realize_instances_001.outputs[0], reroute_002.inputs[0])
    #join_geometry.Geometry -> set_shade_smooth.Geometry
    bonds.links.new(join_geometry.outputs[0], set_shade_smooth.inputs[0])
    #set_position.Geometry -> reroute_021.Input
    bonds.links.new(set_position.outputs[0], reroute_021.inputs[0])
    #reroute_006.Output -> join_geometry.Geometry
    bonds.links.new(reroute_006.outputs[0], join_geometry.inputs[0])
    #group_input_001.bonded_collection -> collection_info.Collection
    bonds.links.new(group_input_001.outputs[2], collection_info.inputs[0])
    #collection_info.Instances -> bounding_box.Geometry
    bonds.links.new(collection_info.outputs[0], bounding_box.inputs[0])
    #bounding_box.Bounding Box -> realize_instances_003.Geometry
    bonds.links.new(bounding_box.outputs[0], realize_instances_003.inputs[0])
    #index.Index -> math_003.Value
    bonds.links.new(index.outputs[0], math_003.inputs[0])
    #math_003.Value -> compare.A
    bonds.links.new(math_003.outputs[0], compare.inputs[0])
    #math_003.Value -> compare.A
    bonds.links.new(math_003.outputs[0], compare.inputs[2])
    #realize_instances_003.Geometry -> separate_geometry.Geometry
    bonds.links.new(realize_instances_003.outputs[0], separate_geometry.inputs[0])
    #compare.Result -> separate_geometry.Selection
    bonds.links.new(compare.outputs[0], separate_geometry.inputs[1])
    #separate_geometry.Selection -> mesh_to_curve.Mesh
    bonds.links.new(separate_geometry.outputs[0], mesh_to_curve.inputs[0])
    #mesh_to_curve.Curve -> store_named_attribute_005.Geometry
    bonds.links.new(mesh_to_curve.outputs[0], store_named_attribute_005.inputs[0])
    #spline_length.Length -> store_named_attribute_005.Value
    bonds.links.new(spline_length.outputs[0], store_named_attribute_005.inputs[3])
    #reroute_010.Output -> sample_index_004.Index
    bonds.links.new(reroute_010.outputs[0], sample_index_004.inputs[2])
    #reroute_011.Output -> sample_index_005.Index
    bonds.links.new(reroute_011.outputs[0], sample_index_005.inputs[2])
    #reroute_018.Output -> sample_index_004.Geometry
    bonds.links.new(reroute_018.outputs[0], sample_index_004.inputs[0])
    #reroute_016.Output -> reroute_011.Input
    bonds.links.new(reroute_016.outputs[0], reroute_011.inputs[0])
    #reroute_017.Output -> reroute_010.Input
    bonds.links.new(reroute_017.outputs[0], reroute_010.inputs[0])
    #reroute_018.Output -> sample_index_005.Geometry
    bonds.links.new(reroute_018.outputs[0], sample_index_005.inputs[0])
    #reroute_005.Output -> reroute_018.Input
    bonds.links.new(reroute_005.outputs[0], reroute_018.inputs[0])
    #realize_instances_001.Geometry -> reroute_005.Input
    bonds.links.new(realize_instances_001.outputs[0], reroute_005.inputs[0])
    #reroute_013.Output -> reroute_017.Input
    bonds.links.new(reroute_013.outputs[0], reroute_017.inputs[0])
    #reroute_012.Output -> reroute_016.Input
    bonds.links.new(reroute_012.outputs[0], reroute_016.inputs[0])
    #named_attribute_002.Attribute -> sample_index_006.Value
    bonds.links.new(named_attribute_002.outputs[0], sample_index_006.inputs[1])
    #reroute.Output -> store_named_attribute_002.Value
    bonds.links.new(reroute.outputs[0], store_named_attribute_002.inputs[3])
    #reroute_001.Output -> store_named_attribute_003.Value
    bonds.links.new(reroute_001.outputs[0], store_named_attribute_003.inputs[3])
    #sample_index_007.Value -> reroute.Input
    bonds.links.new(sample_index_007.outputs[0], reroute.inputs[0])
    #sample_index_006.Value -> reroute_001.Input
    bonds.links.new(sample_index_006.outputs[0], reroute_001.inputs[0])
    #reroute_011.Output -> sample_index_006.Index
    bonds.links.new(reroute_011.outputs[0], sample_index_006.inputs[2])
    #reroute_010.Output -> sample_index_007.Index
    bonds.links.new(reroute_010.outputs[0], sample_index_007.inputs[2])
    #named_attribute_002.Attribute -> sample_index_007.Value
    bonds.links.new(named_attribute_002.outputs[0], sample_index_007.inputs[1])
    #reroute_018.Output -> sample_index_007.Geometry
    bonds.links.new(reroute_018.outputs[0], sample_index_007.inputs[0])
    #reroute_018.Output -> sample_index_006.Geometry
    bonds.links.new(reroute_018.outputs[0], sample_index_006.inputs[0])
    #math_001.Value -> math_006.Value
    bonds.links.new(math_001.outputs[0], math_006.inputs[0])
    #group_input_002.distance -> math.Value
    bonds.links.new(group_input_002.outputs[0], math.inputs[0])
    #curve_to_points.Points -> instance_on_points_001.Points
    bonds.links.new(curve_to_points.outputs[0], instance_on_points_001.inputs[0])
    #points.Points -> instance_on_points_001.Instance
    bonds.links.new(points.outputs[0], instance_on_points_001.inputs[2])
    #store_named_attribute_005.Geometry -> curve_to_points.Curve
    bonds.links.new(store_named_attribute_005.outputs[0], curve_to_points.inputs[0])
    #instance_on_points_001.Instances -> translate_instances.Instances
    bonds.links.new(instance_on_points_001.outputs[0], translate_instances.inputs[0])
    #named_attribute_001.Attribute -> combine_xyz.X
    bonds.links.new(named_attribute_001.outputs[0], combine_xyz.inputs[0])
    #named_attribute_001.Attribute -> combine_xyz.Y
    bonds.links.new(named_attribute_001.outputs[0], combine_xyz.inputs[1])
    #named_attribute_001.Attribute -> combine_xyz.Z
    bonds.links.new(named_attribute_001.outputs[0], combine_xyz.inputs[2])
    #vector_math_001.Vector -> translate_instances.Translation
    bonds.links.new(vector_math_001.outputs[0], translate_instances.inputs[2])
    #combine_xyz.Vector -> vector_math_001.Vector
    bonds.links.new(combine_xyz.outputs[0], vector_math_001.inputs[0])
    #translate_instances.Instances -> realize_instances_001.Geometry
    bonds.links.new(translate_instances.outputs[0], realize_instances_001.inputs[0])
    #math_006.Value -> math_002.Value
    bonds.links.new(math_006.outputs[0], math_002.inputs[1])
    #named_attribute_011.Attribute -> math_001.Value
    bonds.links.new(named_attribute_011.outputs[0], math_001.inputs[1])
    #named_attribute_009.Attribute -> math_001.Value
    bonds.links.new(named_attribute_009.outputs[0], math_001.inputs[0])
    return bonds


def bonds_node_group():

    bonds = mat.node_tree
    #start with a clean node tree
    for node in bonds.nodes:
        bonds.nodes.remove(node)
    bonds.color_tag = 'NONE'
    bonds.description = ""
    bonds.default_group_node_width = 140
    

    #bonds interface

    #initialize bonds nodes
    #node Material Output
    material_output = bonds.nodes.new("ShaderNodeOutputMaterial")
    material_output.name = "Material Output"
    material_output.is_active_output = True
    material_output.target = 'ALL'
    #Displacement
    material_output.inputs[2].default_value = (0.0, 0.0, 0.0)
    #Thickness
    material_output.inputs[3].default_value = 0.0

    #node Attribute.001
    attribute_001 = bonds.nodes.new("ShaderNodeAttribute")
    attribute_001.name = "Attribute.001"
    attribute_001.attribute_name = "EMISSION_CURVE"
    attribute_001.attribute_type = 'GEOMETRY'

    #node Math
    math = bonds.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = 0.10000000149011612

    #node Attribute
    attribute = bonds.nodes.new("ShaderNodeAttribute")
    attribute.name = "Attribute"
    attribute.attribute_name = "COLOR_CURVE"
    attribute.attribute_type = 'GEOMETRY'

    #node Principled BSDF
    principled_bsdf = bonds.nodes.new("ShaderNodeBsdfPrincipled")
    principled_bsdf.name = "Principled BSDF"
    principled_bsdf.distribution = 'GGX'
    principled_bsdf.subsurface_method = 'RANDOM_WALK_SKIN'
    #Base Color
    principled_bsdf.inputs[0].default_value = (0.800000011920929, 0.800000011920929, 0.800000011920929, 1.0)
    #Metallic
    principled_bsdf.inputs[1].default_value = 0.0
    #Roughness
    principled_bsdf.inputs[2].default_value = 0.5
    #IOR
    principled_bsdf.inputs[3].default_value = 1.4500000476837158
    #Alpha
    principled_bsdf.inputs[4].default_value = 1.0
    #Normal
    principled_bsdf.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Diffuse Roughness
    principled_bsdf.inputs[7].default_value = 0.0
    #Subsurface Weight
    principled_bsdf.inputs[8].default_value = 0.0
    #Subsurface Radius
    principled_bsdf.inputs[9].default_value = (1.0, 0.20000000298023224, 0.10000000149011612)
    #Subsurface Scale
    principled_bsdf.inputs[10].default_value = 0.05000000074505806
    #Subsurface IOR
    principled_bsdf.inputs[11].default_value = 1.399999976158142
    #Subsurface Anisotropy
    principled_bsdf.inputs[12].default_value = 0.0
    #Specular IOR Level
    principled_bsdf.inputs[13].default_value = 0.5
    #Specular Tint
    principled_bsdf.inputs[14].default_value = (1.0, 1.0, 1.0, 1.0)
    #Anisotropic
    principled_bsdf.inputs[15].default_value = 0.0
    #Anisotropic Rotation
    principled_bsdf.inputs[16].default_value = 0.0
    #Tangent
    principled_bsdf.inputs[17].default_value = (0.0, 0.0, 0.0)
    #Transmission Weight
    principled_bsdf.inputs[18].default_value = 0.05740181356668472
    #Coat Weight
    principled_bsdf.inputs[19].default_value = 0.0
    #Coat Roughness
    principled_bsdf.inputs[20].default_value = 0.029999999329447746
    #Coat IOR
    principled_bsdf.inputs[21].default_value = 1.5
    #Coat Tint
    principled_bsdf.inputs[22].default_value = (1.0, 1.0, 1.0, 1.0)
    #Coat Normal
    principled_bsdf.inputs[23].default_value = (0.0, 0.0, 0.0)
    #Sheen Weight
    principled_bsdf.inputs[24].default_value = 0.0
    #Sheen Roughness
    principled_bsdf.inputs[25].default_value = 0.5
    #Sheen Tint
    principled_bsdf.inputs[26].default_value = (1.0, 1.0, 1.0, 1.0)
    #Thin Film Thickness
    principled_bsdf.inputs[29].default_value = 0.0
    #Thin Film IOR
    principled_bsdf.inputs[30].default_value = 1.3300000429153442


    #Set locations
    material_output.location = (344.0, 300.0)
    attribute_001.location = (-559.7167358398438, -72.7958755493164)
    math.location = (-292.2492370605469, -216.01373291015625)
    attribute.location = (-572.2391967773438, 203.6114959716797)
    principled_bsdf.location = (-37.92226028442383, 256.66119384765625)

    #Set dimensions
    material_output.width, material_output.height = 140.0, 100.0
    attribute_001.width, attribute_001.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0
    attribute.width, attribute.height = 140.0, 100.0
    principled_bsdf.width, principled_bsdf.height = 240.0, 100.0

    #initialize bonds links
    #principled_bsdf.BSDF -> material_output.Surface
    bonds.links.new(principled_bsdf.outputs[0], material_output.inputs[0])
    #attribute.Color -> principled_bsdf.Emission Color
    bonds.links.new(attribute.outputs[0], principled_bsdf.inputs[27])
    #attribute_001.Fac -> math.Value
    bonds.links.new(attribute_001.outputs[2], math.inputs[0])
    #math.Value -> principled_bsdf.Emission Strength
    bonds.links.new(math.outputs[0], principled_bsdf.inputs[28])
    return bonds


#initialize BONDS_MAT node group
def bonds_node_group(mat):

    bonds = mat.node_tree
    #start with a clean node tree
    for node in bonds.nodes:
        bonds.nodes.remove(node)
    bonds.color_tag = 'NONE'
    bonds.description = ""
    bonds.default_group_node_width = 140
    

    #bonds interface

    #initialize bonds nodes
    #node Material Output
    material_output = bonds.nodes.new("ShaderNodeOutputMaterial")
    material_output.name = "Material Output"
    material_output.is_active_output = True
    material_output.target = 'ALL'
    #Displacement
    material_output.inputs[2].default_value = (0.0, 0.0, 0.0)
    #Thickness
    material_output.inputs[3].default_value = 0.0

    #node Attribute
    attribute = bonds.nodes.new("ShaderNodeAttribute")
    attribute.name = "Attribute"
    attribute.attribute_name = "COLOR_CURVE"
    attribute.attribute_type = 'GEOMETRY'

    #node Principled BSDF
    principled_bsdf = bonds.nodes.new("ShaderNodeBsdfPrincipled")
    principled_bsdf.name = "Principled BSDF"
    principled_bsdf.distribution = 'GGX'
    principled_bsdf.subsurface_method = 'RANDOM_WALK_SKIN'
    #Base Color
    principled_bsdf.inputs[0].default_value = (0.800000011920929, 0.800000011920929, 0.800000011920929, 1.0)
    #Metallic
    principled_bsdf.inputs[1].default_value = 0.0
    #Roughness
    principled_bsdf.inputs[2].default_value = 0.5
    #IOR
    principled_bsdf.inputs[3].default_value = 1.4500000476837158
    #Alpha
    principled_bsdf.inputs[4].default_value = 1.0
    #Normal
    principled_bsdf.inputs[5].default_value = (0.0, 0.0, 0.0)
    #Diffuse Roughness
    principled_bsdf.inputs[7].default_value = 0.0
    #Subsurface Weight
    principled_bsdf.inputs[8].default_value = 0.0
    #Subsurface Radius
    principled_bsdf.inputs[9].default_value = (1.0, 0.20000000298023224, 0.10000000149011612)
    #Subsurface Scale
    principled_bsdf.inputs[10].default_value = 0.05000000074505806
    #Subsurface IOR
    principled_bsdf.inputs[11].default_value = 1.399999976158142
    #Subsurface Anisotropy
    principled_bsdf.inputs[12].default_value = 0.0
    #Specular IOR Level
    principled_bsdf.inputs[13].default_value = 0.5
    #Specular Tint
    principled_bsdf.inputs[14].default_value = (1.0, 1.0, 1.0, 1.0)
    #Anisotropic
    principled_bsdf.inputs[15].default_value = 0.0
    #Anisotropic Rotation
    principled_bsdf.inputs[16].default_value = 0.0
    #Tangent
    principled_bsdf.inputs[17].default_value = (0.0, 0.0, 0.0)
    #Transmission Weight
    principled_bsdf.inputs[18].default_value = 0.05740181356668472
    #Coat Weight
    principled_bsdf.inputs[19].default_value = 0.0
    #Coat Roughness
    principled_bsdf.inputs[20].default_value = 0.029999999329447746
    #Coat IOR
    principled_bsdf.inputs[21].default_value = 1.5
    #Coat Tint
    principled_bsdf.inputs[22].default_value = (1.0, 1.0, 1.0, 1.0)
    #Coat Normal
    principled_bsdf.inputs[23].default_value = (0.0, 0.0, 0.0)
    #Sheen Weight
    principled_bsdf.inputs[24].default_value = 0.0
    #Sheen Roughness
    principled_bsdf.inputs[25].default_value = 0.5
    #Sheen Tint
    principled_bsdf.inputs[26].default_value = (1.0, 1.0, 1.0, 1.0)
    #Emission Strength
    principled_bsdf.inputs[28].default_value = 0.009999999776482582
    #Thin Film Thickness
    principled_bsdf.inputs[29].default_value = 0.0
    #Thin Film IOR
    principled_bsdf.inputs[30].default_value = 1.3300000429153442


    #Set locations
    material_output.location = (344.0, 300.0)
    attribute.location = (-455.671142578125, 84.5844955444336)
    principled_bsdf.location = (-113.27938079833984, 240.1623992919922)

    #Set dimensions
    material_output.width, material_output.height = 140.0, 100.0
    attribute.width, attribute.height = 140.0, 100.0
    principled_bsdf.width, principled_bsdf.height = 240.0, 100.0

    #initialize bonds links
    #principled_bsdf.BSDF -> material_output.Surface
    bonds.links.new(principled_bsdf.outputs[0], material_output.inputs[0])
    #attribute.Color -> principled_bsdf.Emission Color
    bonds.links.new(attribute.outputs[0], principled_bsdf.inputs[27])
    return bonds



def make_bonds():
    
    bpy.primitive_cube_add(size=1, location=(0, 0, 0))
    bonds_obj = bpy.context.object
    bonds_obj.name = "bonds_object"
    if "BONDS_MAT" not in bpy.data.materials:
        mat = bpy.data.materials.new(name = "BONDS_MAT")
        mat.use_nodes = True
    else:
        mat = bpy.data.materials["BONDS_MAT"]
    bonds = bonds_node_group(mat)
    bonds_geometry_node_group()
    return bonds_obj