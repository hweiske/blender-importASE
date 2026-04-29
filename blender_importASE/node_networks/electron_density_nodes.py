import bpy
    
#initialize visualize_edensity node group
def visualize_edensity_node_group(): #from node2python
    visualize_edensity= bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "visualize_edensity")

    #initialize visualize_edensity nodes
    #node Frame
    frame = visualize_edensity.nodes.new("NodeFrame")
    frame.label = "compare -direction"

    #node Frame.001
    frame_001 = visualize_edensity.nodes.new("NodeFrame")
    frame_001.label = "compare +direction"

    #node Volume to Mesh.001
    volume_to_mesh_001 = visualize_edensity.nodes.new("GeometryNodeVolumeToMesh")
    volume_to_mesh_001.resolution_mode = 'GRID'
    #Voxel Size
    volume_to_mesh_001.inputs[1].default_value = 0.30000001192092896
    #Voxel Amount
    volume_to_mesh_001.inputs[2].default_value = 64.0
    #Adaptivity
    volume_to_mesh_001.inputs[4].default_value = 0.0

    #node Join Geometry
    join_geometry = visualize_edensity.nodes.new("GeometryNodeJoinGeometry")

    #visualize_edensity outputs
    #output Geometry
    #visualize_edensity.outputs.new('NodeSocketGeometry', "Geometry")
    #visualize_edensity.outputs[0].attribute_domain = 'POINT'


    #node Group Output
    group_output = visualize_edensity.nodes.new("NodeGroupOutput")

    #node Set Shade Smooth
    set_shade_smooth = visualize_edensity.nodes.new("GeometryNodeSetShadeSmooth")
    #Selection
    set_shade_smooth.inputs[1].default_value = True
    #Shade Smooth
    set_shade_smooth.inputs[2].default_value = True

    #node Delete Geometry.005
    delete_geometry_005 = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_005.domain = 'POINT'
    delete_geometry_005.mode = 'ALL'

    #node Delete Geometry.004
    delete_geometry_004 = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_004.domain = 'POINT'
    delete_geometry_004.mode = 'ALL'

    #node Delete Geometry.003
    delete_geometry_003 = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_003.domain = 'POINT'
    delete_geometry_003.mode = 'ALL'

    #node Compare.005
    compare_005 = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare_005.data_type = 'FLOAT'
    compare_005.operation = 'LESS_THAN'
    compare_005.mode = 'ELEMENT'
    #A_INT
    compare_005.inputs[2].default_value = 0
    #B_INT
    compare_005.inputs[3].default_value = 0
    #A_VEC3
    compare_005.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare_005.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare_005.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare_005.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare_005.inputs[8].default_value = ""
    #B_STR
    compare_005.inputs[9].default_value = ""
    #C
    compare_005.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare_005.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare_005.inputs[12].default_value = 0.0010000000474974513

    #node Compare.003
    compare_003 = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare_003.data_type = 'FLOAT'
    compare_003.operation = 'LESS_THAN'
    compare_003.mode = 'ELEMENT'
    #A_INT
    compare_003.inputs[2].default_value = 0
    #B_INT
    compare_003.inputs[3].default_value = 0
    #A_VEC3
    compare_003.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare_003.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare_003.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare_003.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare_003.inputs[8].default_value = ""
    #B_STR
    compare_003.inputs[9].default_value = ""
    #C
    compare_003.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare_003.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare_003.inputs[12].default_value = 0.0010000000474974513

    #node Compare.004
    compare_004 = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare_004.data_type = 'FLOAT'
    compare_004.operation = 'LESS_THAN'
    compare_004.mode = 'ELEMENT'
    #A_INT
    compare_004.inputs[2].default_value = 0
    #B_INT
    compare_004.inputs[3].default_value = 0
    #A_VEC3
    compare_004.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare_004.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare_004.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare_004.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare_004.inputs[8].default_value = ""
    #B_STR
    compare_004.inputs[9].default_value = ""
    #C
    compare_004.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare_004.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare_004.inputs[12].default_value = 0.0010000000474974513

    #node Delete Geometry
    delete_geometry = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.domain = 'POINT'
    delete_geometry.mode = 'ALL'

    #node Compare
    compare = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare.data_type = 'FLOAT'
    compare.operation = 'LESS_THAN'
    compare.mode = 'ELEMENT'
    #A_INT
    compare.inputs[2].default_value = 0
    #B_INT
    compare.inputs[3].default_value = 0
    #A_VEC3
    compare.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare.inputs[8].default_value = ""
    #B_STR
    compare.inputs[9].default_value = ""
    #C
    compare.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare.inputs[12].default_value = 0.0010000000474974513

    #node Compare.001
    compare_001 = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare_001.data_type = 'FLOAT'
    compare_001.operation = 'LESS_THAN'
    compare_001.mode = 'ELEMENT'
    #A_INT
    compare_001.inputs[2].default_value = 0
    #B_INT
    compare_001.inputs[3].default_value = 0
    #A_VEC3
    compare_001.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare_001.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare_001.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare_001.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare_001.inputs[8].default_value = ""
    #B_STR
    compare_001.inputs[9].default_value = ""
    #C
    compare_001.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare_001.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare_001.inputs[12].default_value = 0.0010000000474974513

    #node Delete Geometry.001
    delete_geometry_001 = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_001.domain = 'POINT'
    delete_geometry_001.mode = 'ALL'

    #node Compare.002
    compare_002 = visualize_edensity.nodes.new("FunctionNodeCompare")
    compare_002.data_type = 'FLOAT'
    compare_002.operation = 'LESS_THAN'
    compare_002.mode = 'ELEMENT'
    #A_INT
    compare_002.inputs[2].default_value = 0
    #B_INT
    compare_002.inputs[3].default_value = 0
    #A_VEC3
    compare_002.inputs[4].default_value = (0.0, 0.0, 0.0)
    #B_VEC3
    compare_002.inputs[5].default_value = (0.0, 0.0, 0.0)
    #A_COL
    compare_002.inputs[6].default_value = (0.0, 0.0, 0.0, 0.0)
    #B_COL
    compare_002.inputs[7].default_value = (0.0, 0.0, 0.0, 0.0)
    #A_STR
    compare_002.inputs[8].default_value = ""
    #B_STR
    compare_002.inputs[9].default_value = ""
    #C
    compare_002.inputs[10].default_value = 0.8999999761581421
    #Angle
    compare_002.inputs[11].default_value = 0.08726649731397629
    #Epsilon
    compare_002.inputs[12].default_value = 0.0010000000474974513

    #node Delete Geometry.002
    delete_geometry_002 = visualize_edensity.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry_002.domain = 'POINT'
    delete_geometry_002.mode = 'ALL'

    #node Separate XYZ.001
    separate_xyz_001 = visualize_edensity.nodes.new("ShaderNodeSeparateXYZ")

    #node Math.002
    math_002 = visualize_edensity.nodes.new("ShaderNodeMath")
    math_002.operation = 'SUBTRACT'
    #Value_002
    math_002.inputs[2].default_value = 0.5

    #node Math.001
    math_001 = visualize_edensity.nodes.new("ShaderNodeMath")
    math_001.operation = 'SUBTRACT'
    #Value_002
    math_001.inputs[2].default_value = 0.5

    #node Math.003
    math_003 = visualize_edensity.nodes.new("ShaderNodeMath")
    math_003.operation = 'SUBTRACT'
    #Value_002
    math_003.inputs[2].default_value = 0.5

    #node Bounding Box
    bounding_box = visualize_edensity.nodes.new("GeometryNodeBoundBox")

    #node Position
    position = visualize_edensity.nodes.new("GeometryNodeInputPosition")

    #node Separate XYZ
    separate_xyz = visualize_edensity.nodes.new("ShaderNodeSeparateXYZ")

    #visualize_edensity.interface.items_tree
    #input Geometry
    #visualize_edensity.interface.items_tree.new('NodeSocketGeometry', "Geometry")
    #visualize_edensity.interface.items_tree[0].attribute_domain = 'POINT'
    visualize_edensity.interface.new_socket('Geometry',in_out='OUTPUT',socket_type='NodeSocketGeometry')
    visualize_edensity.interface.new_socket('Geometry',in_out='INPUT',socket_type='NodeSocketGeometry')
    visualize_edensity.interface.new_socket('isovalue',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff X',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff Y',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff Z',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff -X',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff -Y',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('cutoff -Z',in_out='INPUT',socket_type='NodeSocketFloat')
    visualize_edensity.interface.new_socket('+ material',in_out='INPUT',socket_type='NodeSocketMaterial')
    visualize_edensity.interface.new_socket('- material',in_out='INPUT',socket_type='NodeSocketMaterial')

    #input isovalue

    visualize_edensity.interface.items_tree[2].default_value = 0.03
    visualize_edensity.interface.items_tree[2].min_value = 0
    visualize_edensity.interface.items_tree[2].max_value = 100
    visualize_edensity.interface.items_tree[2].attribute_domain = 'POINT'

    #input cutoff X
    visualize_edensity.interface.items_tree[3].default_value = 0.0
    visualize_edensity.interface.items_tree[3].min_value = 0.0
    visualize_edensity.interface.items_tree[3].max_value = 10000.0
    visualize_edensity.interface.items_tree[3].attribute_domain = 'POINT'

    #input cutoff Y
    visualize_edensity.interface.items_tree[4].default_value = 0.0
    visualize_edensity.interface.items_tree[4].min_value = 0.0
    visualize_edensity.interface.items_tree[4].max_value = 100.0
    visualize_edensity.interface.items_tree[4].attribute_domain = 'POINT'

    #input cutoff Z
    visualize_edensity.interface.items_tree[5].default_value = 0.0
    visualize_edensity.interface.items_tree[5].min_value = 0.0
    visualize_edensity.interface.items_tree[5].max_value = 100.0
    visualize_edensity.interface.items_tree[5].attribute_domain = 'POINT'

   #input cutoff -X
    visualize_edensity.interface.items_tree[6].default_value = 0.0
    visualize_edensity.interface.items_tree[6].min_value = 0.0
    visualize_edensity.interface.items_tree[6].max_value = 1000.0
    visualize_edensity.interface.items_tree[6].attribute_domain = 'POINT'

    #input cutoff -Y
    visualize_edensity.interface.items_tree[7].default_value = 0.0
    visualize_edensity.interface.items_tree[7].min_value = 0.0
    visualize_edensity.interface.items_tree[7].max_value = 10000.0
    visualize_edensity.interface.items_tree[7].attribute_domain = 'POINT'

    #input cutoff -Z
    visualize_edensity.interface.items_tree[8].default_value = 0.0
    visualize_edensity.interface.items_tree[8].min_value = 0.0
    visualize_edensity.interface.items_tree[8].max_value = 10000.0
    visualize_edensity.interface.items_tree[8].attribute_domain = 'POINT'

    
    matp = newShader("+ material",  0, 0, 1)
    bpy.context.active_object.data.materials.append(matp)
    matm = newShader("- material",  1, 0, 0)
    bpy.context.active_object.data.materials.append(matm)



    #input material +
    visualize_edensity.interface.items_tree[9].attribute_domain = 'POINT'
    visualize_edensity.interface.items_tree[9].default_value = matp

    #input material -
    visualize_edensity.interface.items_tree[10].attribute_domain = 'POINT'
    visualize_edensity.interface.items_tree[10].default_value = matm


    #node Group Input
    group_input = visualize_edensity.nodes.new("NodeGroupInput")

    #node Set Material
    set_material = visualize_edensity.nodes.new("GeometryNodeSetMaterial")
    #Selection
    set_material.inputs[1].default_value = True

    #node Set Material.001
    set_material_001 = visualize_edensity.nodes.new("GeometryNodeSetMaterial")
    #Selection
    set_material_001.inputs[1].default_value = True

    #node Volume to Mesh
    volume_to_mesh = visualize_edensity.nodes.new("GeometryNodeVolumeToMesh")
    volume_to_mesh.resolution_mode = 'GRID'
    #Voxel Size
    volume_to_mesh.inputs[1].default_value = 0.30000001192092896
    #Voxel Amount
    volume_to_mesh.inputs[2].default_value = 64.0
    #Adaptivity
    volume_to_mesh.inputs[4].default_value = 0.0

    #node Math
    math = visualize_edensity.nodes.new("ShaderNodeMath")
    math.operation = 'MULTIPLY'
    #Value_001
    math.inputs[1].default_value = -1.0
    #Value_002
    math.inputs[2].default_value = 0.5

    #Set parents
    delete_geometry_005.parent = frame
    delete_geometry_004.parent = frame
    delete_geometry_003.parent = frame
    compare_005.parent = frame
    compare_003.parent = frame
    compare_004.parent = frame
    delete_geometry.parent = frame_001
    compare.parent = frame_001
    compare_001.parent = frame_001
    delete_geometry_001.parent = frame_001
    compare_002.parent = frame_001
    delete_geometry_002.parent = frame_001
    separate_xyz_001.parent = frame
    math_002.parent = frame
    math_001.parent = frame
    math_003.parent = frame
    bounding_box.parent = frame

    #Set locations
    frame.location = (0.0, 0.0)
    frame_001.location = (0.0, 0.0)
    volume_to_mesh_001.location = (-25.607208251953125, -156.8311004638672)
    join_geometry.location = (621.0, 93.33197021484375)
    group_output.location = (1577.626953125, 109.21952056884766)
    set_shade_smooth.location = (818.9706420898438, 147.1034393310547)
    delete_geometry_005.location = (1216.2142333984375, -746.0068359375)
    delete_geometry_004.location = (1208.96435546875, -860.6949462890625)
    delete_geometry_003.location = (1212.4752197265625, -663.9649658203125)
    compare_005.location = (1218.83837890625, -791.5283203125)
    compare_003.location = (1218.8431396484375, -625.2783203125)
    compare_004.location = (1203.624267578125, -700.189453125)
    delete_geometry.location = (1104.0548095703125, -143.5847625732422)
    compare.location = (1099.885009765625, -95.4963150024414)
    compare_001.location = (1095.203857421875, -179.80917358398438)
    delete_geometry_001.location = (1100.77197265625, -220.94252014160156)
    compare_002.location = (1110.41796875, -271.1480712890625)
    delete_geometry_002.location = (1103.7412109375, -332.3167419433594)
    separate_xyz_001.location = (872.7716674804688, -815.6917114257812)
    math_002.location = (1031.317138671875, -759.0149536132812)
    math_001.location = (1030.524658203125, -801.5042114257812)
    math_003.location = (1038.015380859375, -726.0418701171875)
    bounding_box.location = (878.5214233398438, -682.5440063476562)
    position.location = (413.29278564453125, -416.83447265625)
    separate_xyz.location = (561.9790649414062, -364.5952453613281)
    group_input.location = (-674.58447265625, -288.4441833496094)
    set_material.location = (400.64227294921875, 113.72370910644531)
    set_material_001.location = (392.53057861328125, -39.666114807128906)
    volume_to_mesh.location = (-284.9970397949219, 83.57960510253906)
    math.location = (-183.65554809570312, -147.9720916748047)

    #Set dimensions
    frame.width, frame.height = 546.0, 336.0
    frame_001.width, frame_001.height = 215.0, 337.0
    volume_to_mesh_001.width, volume_to_mesh_001.height = 170.0, 100.0
    join_geometry.width, join_geometry.height = 140.0, 100.0
    group_output.width, group_output.height = 140.0, 100.0
    set_shade_smooth.width, set_shade_smooth.height = 140.0, 100.0
    delete_geometry_005.width, delete_geometry_005.height = 140.0, 100.0
    delete_geometry_004.width, delete_geometry_004.height = 140.0, 100.0
    delete_geometry_003.width, delete_geometry_003.height = 140.0, 100.0
    compare_005.width, compare_005.height = 140.0, 100.0
    compare_003.width, compare_003.height = 140.0, 100.0
    compare_004.width, compare_004.height = 140.0, 100.0
    delete_geometry.width, delete_geometry.height = 140.0, 100.0
    compare.width, compare.height = 140.0, 100.0
    compare_001.width, compare_001.height = 140.0, 100.0
    delete_geometry_001.width, delete_geometry_001.height = 140.0, 100.0
    compare_002.width, compare_002.height = 140.0, 100.0
    delete_geometry_002.width, delete_geometry_002.height = 140.0, 100.0
    separate_xyz_001.width, separate_xyz_001.height = 140.0, 100.0
    math_002.width, math_002.height = 140.0, 100.0
    math_001.width, math_001.height = 140.0, 100.0
    math_003.width, math_003.height = 140.0, 100.0
    bounding_box.width, bounding_box.height = 140.0, 100.0
    position.width, position.height = 140.0, 100.0
    separate_xyz.width, separate_xyz.height = 140.0, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    set_material.width, set_material.height = 140.0, 100.0
    set_material_001.width, set_material_001.height = 140.0, 100.0
    volume_to_mesh.width, volume_to_mesh.height = 170.0, 100.0
    math.width, math.height = 140.0, 100.0

    #initialize visualize_edensity links
    #group_input.Geometry -> volume_to_mesh.Volume
    visualize_edensity.links.new(group_input.outputs[0], volume_to_mesh.inputs[0])
    #group_input.isovalue -> math.Value
    visualize_edensity.links.new(group_input.outputs[1], math.inputs[0])
    #group_input.isovalue -> volume_to_mesh.Threshold
    visualize_edensity.links.new(group_input.outputs[1], volume_to_mesh.inputs[3])
    #math.Value -> volume_to_mesh_001.Threshold
    visualize_edensity.links.new(math.outputs[0], volume_to_mesh_001.inputs[3])
    #group_input.Geometry -> volume_to_mesh_001.Volume
    visualize_edensity.links.new(group_input.outputs[0], volume_to_mesh_001.inputs[0])
    #volume_to_mesh.Mesh -> set_material.Geometry
    visualize_edensity.links.new(volume_to_mesh.outputs[0], set_material.inputs[0])
    #set_material.Geometry -> join_geometry.Geometry
    visualize_edensity.links.new(set_material.outputs[0], join_geometry.inputs[0])
    #set_material_001.Geometry -> join_geometry.Geometry
    visualize_edensity.links.new(set_material_001.outputs[0], join_geometry.inputs[0])
    #volume_to_mesh_001.Mesh -> set_material_001.Geometry
    visualize_edensity.links.new(volume_to_mesh_001.outputs[0], set_material_001.inputs[0])
    #join_geometry.Geometry -> set_shade_smooth.Geometry
    visualize_edensity.links.new(join_geometry.outputs[0], set_shade_smooth.inputs[0])
    #separate_xyz.X -> compare.A
    visualize_edensity.links.new(separate_xyz.outputs[0], compare.inputs[0])
    #position.Position -> separate_xyz.Vector
    visualize_edensity.links.new(position.outputs[0], separate_xyz.inputs[0])
    #compare.Result -> delete_geometry.Selection
    visualize_edensity.links.new(compare.outputs[0], delete_geometry.inputs[1])
    #compare_001.Result -> delete_geometry_001.Selection
    visualize_edensity.links.new(compare_001.outputs[0], delete_geometry_001.inputs[1])
    #delete_geometry.Geometry -> delete_geometry_001.Geometry
    visualize_edensity.links.new(delete_geometry.outputs[0], delete_geometry_001.inputs[0])
    #separate_xyz.Y -> compare_001.A
    visualize_edensity.links.new(separate_xyz.outputs[1], compare_001.inputs[0])
    #compare_002.Result -> delete_geometry_002.Selection
    visualize_edensity.links.new(compare_002.outputs[0], delete_geometry_002.inputs[1])
    #delete_geometry_001.Geometry -> delete_geometry_002.Geometry
    visualize_edensity.links.new(delete_geometry_001.outputs[0], delete_geometry_002.inputs[0])
    #separate_xyz.Z -> compare_002.A
    visualize_edensity.links.new(separate_xyz.outputs[2], compare_002.inputs[0])
    #set_shade_smooth.Geometry -> delete_geometry.Geometry
    visualize_edensity.links.new(set_shade_smooth.outputs[0], delete_geometry.inputs[0])
    #group_input.cutoff X -> compare.B
    visualize_edensity.links.new(group_input.outputs[2], compare.inputs[1])
    #group_input.cutoff Y -> compare_001.B
    visualize_edensity.links.new(group_input.outputs[3], compare_001.inputs[1])
    #group_input.cutoff Z -> compare_002.B
    visualize_edensity.links.new(group_input.outputs[4], compare_002.inputs[1])
    #compare_003.Result -> delete_geometry_003.Selection
    visualize_edensity.links.new(compare_003.outputs[0], delete_geometry_003.inputs[1])
    #compare_004.Result -> delete_geometry_005.Selection
    visualize_edensity.links.new(compare_004.outputs[0], delete_geometry_005.inputs[1])
    #delete_geometry_003.Geometry -> delete_geometry_005.Geometry
    visualize_edensity.links.new(delete_geometry_003.outputs[0], delete_geometry_005.inputs[0])
    #compare_005.Result -> delete_geometry_004.Selection
    visualize_edensity.links.new(compare_005.outputs[0], delete_geometry_004.inputs[1])
    #delete_geometry_005.Geometry -> delete_geometry_004.Geometry
    visualize_edensity.links.new(delete_geometry_005.outputs[0], delete_geometry_004.inputs[0])
    #delete_geometry_002.Geometry -> delete_geometry_003.Geometry
    visualize_edensity.links.new(delete_geometry_002.outputs[0], delete_geometry_003.inputs[0])
    #math_001.Value -> compare_005.A
    visualize_edensity.links.new(math_001.outputs[0], compare_005.inputs[0])
    #delete_geometry_004.Geometry -> group_output.Geometry
    visualize_edensity.links.new(delete_geometry_004.outputs[0], group_output.inputs[0])
    #set_shade_smooth.Geometry -> bounding_box.Geometry
    visualize_edensity.links.new(set_shade_smooth.outputs[0], bounding_box.inputs[0])
    #bounding_box.Max -> separate_xyz_001.Vector
    visualize_edensity.links.new(bounding_box.outputs[2], separate_xyz_001.inputs[0])
    #separate_xyz_001.Z -> math_001.Value
    visualize_edensity.links.new(separate_xyz_001.outputs[2], math_001.inputs[0])
    #group_input.cutoff -Z -> compare_005.B
    visualize_edensity.links.new(group_input.outputs[7], compare_005.inputs[1])
    #separate_xyz.Z -> math_001.Value
    visualize_edensity.links.new(separate_xyz.outputs[2], math_001.inputs[1])
    #separate_xyz.Y -> math_002.Value
    visualize_edensity.links.new(separate_xyz.outputs[1], math_002.inputs[1])
    #separate_xyz_001.Y -> math_002.Value
    visualize_edensity.links.new(separate_xyz_001.outputs[1], math_002.inputs[0])
    #separate_xyz.X -> math_003.Value
    visualize_edensity.links.new(separate_xyz.outputs[0], math_003.inputs[1])
    #separate_xyz_001.X -> math_003.Value
    visualize_edensity.links.new(separate_xyz_001.outputs[0], math_003.inputs[0])
    #math_003.Value -> compare_003.A
    visualize_edensity.links.new(math_003.outputs[0], compare_003.inputs[0])
    #math_002.Value -> compare_004.A
    visualize_edensity.links.new(math_002.outputs[0], compare_004.inputs[0])
    #group_input.cutoff -X -> compare_003.B
    visualize_edensity.links.new(group_input.outputs[5], compare_003.inputs[1])
    #group_input.cutoff -Y -> compare_004.B
    visualize_edensity.links.new(group_input.outputs[6], compare_004.inputs[1])
    #group_input.material + -> set_material.Material
    visualize_edensity.links.new(group_input.outputs[8], set_material.inputs[2])
    #group_input.material - -> set_material_001.Material
    visualize_edensity.links.new(group_input.outputs[9], set_material_001.inputs[2])
    return visualize_edensity

def newMaterial(id):
    
    mat = bpy.data.materials.get(id)
    
    if mat is None:
        mat = bpy.data.materials.new(name=id)
    
    mat.use_nodes = True
    
  #  if mat.node_tree:
  #      mat.node_tree.links.clear()
  #      mat.node_tree.nodes.clear()
    
    return mat
    
def newShader(id, r, g, b):
    
    mat = newMaterial(id)
    
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.new(type='ShaderNodeBsdfPrincipled')
    s=nodes["Principled BSDF"]
    s.inputs[0].default_value=(r,g,b,1) #color
    s.inputs[1].default_value=0 #metal
    s.inputs[2].default_value=0.6 #roughness
    s.inputs[3].default_value=1.45 #IOR
    s.inputs[4].default_value=0.6 #alpha
    mat.node_tree.links.clear()
    links.new(nodes['Principled BSDF'].outputs['BSDF'],nodes['Material Output'].inputs['Surface']) 
    return mat

