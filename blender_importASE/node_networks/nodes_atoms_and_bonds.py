import bpy
from ..utils import atomcolors
from ase.data import covalent_radii, chemical_symbols, colors

import bpy, mathutils

def read_structure(atoms,name, animate=True):
    if animate:
        trajectory=atoms
        atoms=trajectory[0]
    vertices=atoms.get_positions()
    object_name=name
    mesh = bpy.data.meshes.new(name=object_name)
    obj = bpy.data.objects.new(name=object_name, object_data=mesh)
    bpy.context.collection.objects.link(obj)
    # Create the mesh from the vertex list
    mesh.from_pydata(vertices, [], [])  # No edges or faces
    if "element" not in mesh.attributes:
        mesh.attributes.new(name="element", type='FLOAT', domain='POINT')
    if "atom_radius" not in mesh.attributes:
        mesh.attributes.new(name="atom_radius", type='FLOAT', domain='POINT')

    element = mesh.attributes["element"].data
    rad = mesh.attributes["atom_radius"].data

    for i, value in enumerate(element):
        atom=atoms[i]
        value.value = atom.number  # Example: setting index value
        # Update the mesh   
    for i, value in enumerate(rad):
        atom=atoms[i]
        rad=covalent_radii[atom.number]
        value.value = rad 
    mesh.update()
    vertx=obj.data.vertices
    #doesnt work yet
    if animate:
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
       # bpy.data.scenes['Scene'].animall_properties.key_point_location = True
        vertx=obj.data.vertices
        for n,frame in enumerate(trajectory):
            #bpy.data.scenes['Scene'].frame_current=n
            for nv,v in enumerate(vertx):
                v.co=frame.positions[nv]
        
                v.keyframe_insert(data_path="co", frame=n)    
        #    bpy.context.view_layer.update()
        #    bpy.ops.object.mode_set(mode='EDIT')
        #    bpy.ops.view3d.insert_keyframe_animall()
        #    bpy.ops.object.mode_set(mode='OBJECT')

    return(obj, mesh)

#initialize set_atoms node group
def set_atoms_node_group():
    set_atoms = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "set_atoms")

    set_atoms.color_tag = 'NONE'
    set_atoms.description = ""
    set_atoms.default_group_node_width = 140
    


    #set_atoms interface
    #Socket Instances
    instances_socket = set_atoms.interface.new_socket(name = "Instances", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    instances_socket.attribute_domain = 'POINT'

    #Socket Points
    points_socket = set_atoms.interface.new_socket(name = "Points", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    points_socket.attribute_domain = 'POINT'

    #Socket Selection
    selection_socket = set_atoms.interface.new_socket(name = "Selection", in_out='INPUT', socket_type = 'NodeSocketBool')
    selection_socket.default_value = True
    selection_socket.attribute_domain = 'POINT'
    selection_socket.hide_value = True

    #Socket Scale
    scale_socket = set_atoms.interface.new_socket(name = "Scale", in_out='INPUT', socket_type = 'NodeSocketVector')
    scale_socket.default_value = (1.0, 1.0, 1.0)
    scale_socket.min_value = -3.4028234663852886e+38
    scale_socket.max_value = 3.4028234663852886e+38
    scale_socket.subtype = 'XYZ'
    scale_socket.attribute_domain = 'POINT'

    #Socket resolution
    resolution_socket = set_atoms.interface.new_socket(name = "resolution", in_out='INPUT', socket_type = 'NodeSocketInt')
    resolution_socket.default_value = 32
    resolution_socket.min_value = 3
    resolution_socket.max_value = 1024
    resolution_socket.subtype = 'NONE'
    resolution_socket.attribute_domain = 'POINT'


    #initialize set_atoms nodes
    #node Group Output
    group_output = set_atoms.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"
    group_output.is_active_output = True

    #node Group Input
    group_input = set_atoms.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    #node Instance on Points
    instance_on_points = set_atoms.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points.name = "Instance on Points"
    #Pick Instance
    instance_on_points.inputs[3].default_value = False
    #Instance Index
    instance_on_points.inputs[4].default_value = 0
    #Rotation
    instance_on_points.inputs[5].default_value = (0.0, 0.0, 0.0)

    #node UV Sphere
    uv_sphere = set_atoms.nodes.new("GeometryNodeMeshUVSphere")
    uv_sphere.name = "UV Sphere"
    #Radius
    uv_sphere.inputs[2].default_value = 0.5

    #node Shade Smooth
    shade_smooth = set_atoms.nodes.new("GeometryNodeSetShadeSmooth")
    shade_smooth.name = "Shade Smooth"
    shade_smooth.domain = 'FACE'
    #Selection
    shade_smooth.inputs[1].default_value = True
    #Shade Smooth
    shade_smooth.inputs[2].default_value = True

    #node Math
    math = set_atoms.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'MULTIPLY'
    math.use_clamp = False
    #Value_001
    math.inputs[1].default_value = 2.0





    #Set locations
    group_output.location = (799.6827392578125, -45.560829162597656)
    group_input.location = (-281.32958984375, -62.15126037597656)
    instance_on_points.location = (388.9354248046875, -31.57483673095703)
    uv_sphere.location = (118.87770080566406, -174.6345672607422)
    shade_smooth.location = (649.55029296875, -45.560829162597656)
    math.location = (-61.0, -204.416015625)

    #Set dimensions
    group_output.width, group_output.height = 140.0, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    uv_sphere.width, uv_sphere.height = 140.0, 100.0
    shade_smooth.width, shade_smooth.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0

    #initialize set_atoms links
    #group_input.Points -> instance_on_points.Points
    set_atoms.links.new(group_input.outputs[0], instance_on_points.inputs[0])
    #group_input.Scale -> instance_on_points.Scale
    set_atoms.links.new(group_input.outputs[2], instance_on_points.inputs[6])
    #instance_on_points.Instances -> shade_smooth.Geometry
    set_atoms.links.new(instance_on_points.outputs[0], shade_smooth.inputs[0])
    #shade_smooth.Geometry -> group_output.Instances
    set_atoms.links.new(shade_smooth.outputs[0], group_output.inputs[0])
    #uv_sphere.Mesh -> instance_on_points.Instance
    set_atoms.links.new(uv_sphere.outputs[0], instance_on_points.inputs[2])
    #group_input.Selection -> instance_on_points.Selection
    set_atoms.links.new(group_input.outputs[1], instance_on_points.inputs[1])
    #group_input.resolution -> math.Value
    set_atoms.links.new(group_input.outputs[3], math.inputs[0])
    #math.Value -> uv_sphere.Segments
    set_atoms.links.new(math.outputs[0], uv_sphere.inputs[0])
    #group_input.resolution -> uv_sphere.Rings
    set_atoms.links.new(group_input.outputs[3], uv_sphere.inputs[1])
    return set_atoms

#initialize atoms_from_verts node group
def atoms_and_bonds(obj, atoms, modifier='GeometryNodes',bondmat=None):
    atomcolor=atomcolors()
    
    atoms_and_bonds = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = f"atoms_and_bonds_{atoms.get_chemical_formula()}")

    atoms_and_bonds.color_tag = 'NONE'
    atoms_and_bonds.description = ""
    atoms_and_bonds.default_group_node_width = 140
    

    atoms_and_bonds.is_modifier = True

    #atoms_from_verts interface
    #Socket Geometry
    geometry_socket = atoms_and_bonds.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_1 = atoms_and_bonds.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_1.attribute_domain = 'POINT'

    
    #Socket Bond distance
    bond_distance = atoms_and_bonds.interface.new_socket(name = "bond_distance", in_out='INPUT', socket_type = 'NodeSocketFloat')
    bond_distance.default_value = 1
    bond_distance.min_value = -10000.0
    bond_distance.max_value = 10000.0
    bond_distance.subtype = 'NONE'
    bond_distance.attribute_domain = 'POINT'

    #Socket bond radius
    radius_socket = atoms_and_bonds.interface.new_socket(name = "bond_radius", in_out='INPUT', socket_type = 'NodeSocketFloat')
    radius_socket.default_value = 0.075
    radius_socket.min_value = 0
    radius_socket.max_value = 3.4028234663852886e+38
    radius_socket.subtype = 'DISTANCE'
    radius_socket.attribute_domain = 'POINT'

    #Socket resolution
    resolution = atoms_and_bonds.interface.new_socket(name = "RESOLUTION", in_out='INPUT', socket_type = 'NodeSocketFloat')
    resolution.default_value = 16
    resolution.min_value = 3
    resolution.max_value = 64
    resolution.socket_type = 'NodeSocketInt'
    resolution.attribute_domain = 'POINT'

    #initialize atoms_from_verts nodes
    #node Group Input
    group_input_at_atoms = atoms_and_bonds.nodes.new("NodeGroupInput")
    group_input_at_atoms.name = "Group Input"

    

    #node Named Attribute
    radius_attribute = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    radius_attribute.label = "radius_attribute"
    radius_attribute.name = "Named Attribute"
    radius_attribute.data_type = 'FLOAT'
    #Name
    radius_attribute.inputs[0].default_value = "atom_radius"

    #node Named Attribute.001
    element_attribute = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    element_attribute.label = "element_attribute"
    element_attribute.name = "Named Attribute.001"
    element_attribute.data_type = 'INT'
    #Name
    element_attribute.inputs[0].default_value = "element"

    #node Join Geometry
    join_geometry_atoms = atoms_and_bonds.nodes.new("GeometryNodeJoinGeometry")
    join_geometry_atoms.name = "Join Geometry"
    
    #node Compare Elements
    numbers=list(set(atoms.get_atomic_numbers()))

    #setup up color attribute
    #color attribute
    color_attribute = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    color_attribute.label = "color_attribute"
    color_attribute.name = "Named Attribute"
    color_attribute.data_type = 'FLOAT_COLOR'
    color_attribute.inputs[2].default_value = "color"
    

    atoms_and_bonds.links.new(group_input_at_atoms.outputs[0], color_attribute.inputs[0])
        

    for n,number in enumerate(numbers):

        sym=chemical_symbols[number]
        compare = atoms_and_bonds.nodes.new("FunctionNodeCompare")
        compare.label = f"is_element_{number}"
        compare.name = "Compare"
        compare.data_type = 'INT'
        compare.mode = 'ELEMENT'
        compare.operation = 'EQUAL'
        compare.inputs[3].default_value = number
        compare.location = (-1200.4930419921875, -100+200*n)
        #node Group
        group = atoms_and_bonds.nodes.new("GeometryNodeGroup")
        group.label = f"set_atoms_{number}"
        group.name = "Group"
        group.node_tree = bpy.data.node_groups['set_atoms']
        group.location = (-458.2143249511719, -100+200*n)
        compare.width, compare.height = 140.0, 100.0
        group.width, group.height = 140.0, 100.0

        #node color
        color = atoms_and_bonds.nodes.new("FunctionNodeInputColor")
        color.name = sym
        if sym in atomcolor.color_dict:
            color.value = list(atomcolor.color_dict[sym]) + [1]
        else:
            color.value = list(colors.jmol_colors[number]) + [1]
        color.location = (-750, -800+200*n)
        color.width, color.height = 140.0, 100.0

       
        


        #node Set Material
        set_material = atoms_and_bonds.nodes.new("GeometryNodeSetMaterial")
        set_material.name = f"Set Material - {sym}"
        #Selection
        set_material.inputs[1].default_value = True
        if sym in bpy.data.materials:
            set_material.inputs[2].default_value = bpy.data.materials[sym]
        set_material.location = (-268.2143249511719, -100+140*n)
        set_material.width, set_material.height = 140.0, 100.0
        
        atoms_and_bonds.links.new(compare.outputs[0], group.inputs[1])
        atoms_and_bonds.links.new(color_attribute.outputs[0], group.inputs[0])
        atoms_and_bonds.links.new(radius_attribute.outputs[0], group.inputs[2])
        atoms_and_bonds.links.new(element_attribute.outputs[0], compare.inputs[2])
        atoms_and_bonds.links.new(group.outputs[0], set_material.inputs[0])
        atoms_and_bonds.links.new(group_input_at_atoms.outputs[3], group.inputs[3])
        atoms_and_bonds.links.new(set_material.outputs[0], join_geometry_atoms.inputs[0])

        
        
        #switches
        if n > 0:
            if n < len(numbers):
                switch_color = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
                switch_color.name = f"Switch Color {sym}-{chemical_symbols[numbers[n-1]]}"
                switch_color.label = f"Switch Color {sym}-{chemical_symbols[numbers[n-1]]}"
                switch_color.inputs[0].default_value = False
                switch_color.input_type = 'RGBA'
                atoms_and_bonds.links.new(compare.outputs[0], switch_color.inputs[0]) #for every n 
            if n == 1:
                old=atoms_and_bonds.nodes[chemical_symbols[numbers[n-1]]]
                atoms_and_bonds.links.new(old.outputs[0], switch_color.inputs[1]) #old is color n-1
                atoms_and_bonds.links.new(color.outputs[0], switch_color.inputs[2]) #color is n
            else:
                old=atoms_and_bonds.nodes[f"Switch Color {chemical_symbols[numbers[n-1]]}-{chemical_symbols[numbers[n-2]]}"]
                atoms_and_bonds.links.new(old.outputs[0], switch_color.inputs[1]) #old is switch_color n-1-2
                atoms_and_bonds.links.new(color.outputs[0], switch_color.inputs[2]) #color is n
            switch_color.location = (-500, -1000+200*n)
        
    if len(numbers) > 1:    
        atoms_and_bonds.links.new(switch_color.outputs[0], color_attribute.inputs[3])
    else:
        atoms_and_bonds.links.new(color.outputs[0], color_attribute.inputs[3])


    


    #Set locations
    group_input_at_atoms.location = (-1500 , 63.91105270385742)
    radius_attribute.location = (-1200, -300)
    element_attribute.location = (-1500, -66.59099578857422)
    join_geometry_atoms.location = (0, 1.9488945007324219)
    color_attribute.location = (-600, 0)
    #Set dimensions
    group_input_at_atoms.width, group_input_at_atoms.height = 140.0, 100.0
    radius_attribute.width, radius_attribute.height = 140.0, 100.0
    element_attribute.width, element_attribute.height = 140.0, 100.0
    
    join_geometry_atoms.width, join_geometry_atoms.height = 140.0, 100.0
    


###################  BONDS ################################

    

    #initialize atoms_from_verts nodes
    #node Frame.008
    frame_008 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_008.label = "instance curves"
    frame_008.name = "Frame.008"
    frame_008.label_size = 20
    frame_008.shrink = True

    #node Frame.001
    frame_001 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_001.label = "curve profile"
    frame_001.name = "Frame.001"
    frame_001.label_size = 20
    frame_001.shrink = True

    #node Frame.016
    frame_016 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_016.label = "color_curve"
    frame_016.name = "Frame.016"
    frame_016.label_size = 20
    frame_016.shrink = True

    #node Frame.005
    frame_005 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_005.label = "000111222"
    frame_005.name = "Frame.005"
    frame_005.label_size = 20
    frame_005.shrink = True

    #node Frame.007
    frame_007 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_007.label = "01230123"
    frame_007.name = "Frame.007"
    frame_007.label_size = 20
    frame_007.shrink = True

    #node Frame.017
    frame_017 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_017.label = "create_loop"
    frame_017.name = "Frame.017"
    frame_017.label_size = 20
    frame_017.shrink = True

    #node Frame.013
    frame_013 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_013.label = "set positions"
    frame_013.name = "Frame.013"
    frame_013.label_size = 20
    frame_013.shrink = True

    #node Frame.018
    frame_018 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_018.label = "CUTOFF"
    frame_018.name = "Frame.018"
    frame_018.label_size = 20
    frame_018.shrink = True

    #node Frame.014
    frame_014 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_014.label = "sample_rad"
    frame_014.name = "Frame.014"
    frame_014.label_size = 20
    frame_014.shrink = True

    #node Frame.012
    frame_012 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_012.label = "sample colors"
    frame_012.name = "Frame.012"
    frame_012.label_size = 20
    frame_012.shrink = True

    #node Frame.019
    frame_019 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_019.label = "store positions, colors"
    frame_019.name = "Frame.019"
    frame_019.label_size = 20
    frame_019.shrink = True

    #node Frame.009
    frame_009 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_009.label = "create all points for i,j nested loop"
    frame_009.name = "Frame.009"
    frame_009.label_size = 20
    frame_009.shrink = True

    #node Frame.004
    frame_004 = atoms_and_bonds.nodes.new("NodeFrame")
    frame_004.label = "all indices"
    frame_004.name = "Frame.004"
    frame_004.label_size = 20
    frame_004.shrink = True

    #node Reroute.006
    reroute_006 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_006.name = "Reroute.006"
    reroute_006.socket_idname = "NodeSocketGeometry"
    #node Reroute.014
    reroute_014 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_014.name = "Reroute.014"
    reroute_014.socket_idname = "NodeSocketColor"
    #node Reroute.015
    reroute_015 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_015.name = "Reroute.015"
    reroute_015.socket_idname = "NodeSocketVector"
    #node Reroute.004
    reroute_004 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_004.name = "Reroute.004"
    reroute_004.socket_idname = "NodeSocketVector"
    #node Reroute.019
    reroute_019 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_019.name = "Reroute.019"
    reroute_019.socket_idname = "NodeSocketFloat"
    #node Set Curve Radius
    set_curve_radius = atoms_and_bonds.nodes.new("GeometryNodeSetCurveRadius")
    set_curve_radius.name = "Set Curve Radius"
    #Selection
    set_curve_radius.inputs[1].default_value = True
    #Radius
    set_curve_radius.inputs[2].default_value = 0.005

    #node Set Material.006
    set_material_006 = atoms_and_bonds.nodes.new("GeometryNodeSetMaterial")
    set_material_006.name = "Set Material.006"
    #Selection
    set_material_006.inputs[1].default_value = True
    set_material_006.inputs[2].default_value = bondmat

    #node Switch.005
    switch_005 = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
    switch_005.name = "Switch.005"
    switch_005.input_type = 'GEOMETRY'
    #Switch
    switch_005.inputs[0].default_value = False

    #node Curve to Mesh
    curve_to_mesh = atoms_and_bonds.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.name = "Curve to Mesh"
    #Fill Caps
    curve_to_mesh.inputs[2].default_value = False

    #node Reroute.021
    reroute_021 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_021.name = "Reroute.021"
    reroute_021.socket_idname = "NodeSocketGeometry"
    #node Curve Circle
    curve_circle = atoms_and_bonds.nodes.new("GeometryNodeCurvePrimitiveCircle")
    curve_circle.name = "Curve Circle"
    curve_circle.hide = True
    curve_circle.mode = 'RADIUS'
    #Resolution
    curve_circle.inputs[0].default_value = 32
    #Radius
    curve_circle.inputs[4].default_value = 1.0

    #node Curve to Mesh.001
    curve_to_mesh_001 = atoms_and_bonds.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh_001.name = "Curve to Mesh.001"
    #Fill Caps
    curve_to_mesh_001.inputs[2].default_value = False

    #node Reroute.007
    reroute_007 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_007.name = "Reroute.007"
    reroute_007.socket_idname = "NodeSocketVector"
    #node Switch.001
    switch_001 = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
    switch_001.name = "Switch.001"
    switch_001.input_type = 'RGBA'

    #node Math.005
    math_005 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_005.name = "Math.005"
    math_005.operation = 'MODULO'
    math_005.use_clamp = False
    #Value_001
    math_005.inputs[1].default_value = 2.0

    #node Index.002
    index_002 = atoms_and_bonds.nodes.new("GeometryNodeInputIndex")
    index_002.name = "Index.002"

    #node Set Position
    set_position = atoms_and_bonds.nodes.new("GeometryNodeSetPosition")
    set_position.name = "Set Position"
    #Selection
    set_position.inputs[1].default_value = True
    #Offset
    set_position.inputs[3].default_value = (0.0, 0.0, 0.0)

    #node Store Named Attribute.005
    store_named_attribute_005 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_005.name = "Store Named Attribute.005"
    store_named_attribute_005.data_type = 'FLOAT_COLOR'
    store_named_attribute_005.domain = 'POINT'
    #Selection
    store_named_attribute_005.inputs[1].default_value = True
    #Name
    store_named_attribute_005.inputs[2].default_value = "COLOR_CURVE"

    #node Reroute.017
    reroute_017 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_017.name = "Reroute.017"
    reroute_017.socket_idname = "NodeSocketColor"
    #node Reroute.028
    reroute_028 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_028.name = "Reroute.028"
    reroute_028.socket_idname = "NodeSocketColor"
    #node Math.004
    math_004 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_004.name = "Math.004"
    math_004.operation = 'FLOOR'
    math_004.use_clamp = False

    #node Math.011
    math_011 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_011.name = "Math.011"
    math_011.operation = 'DIVIDE'
    math_011.use_clamp = False

    #node Math.012
    math_012 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_012.name = "Math.012"
    math_012.operation = 'MODULO'
    math_012.use_clamp = False
    math_012.inputs[2].hide = True

    #node Index.001
    index_001 = atoms_and_bonds.nodes.new("GeometryNodeInputIndex")
    index_001.name = "Index.001"

    #node Switch
    switch = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
    switch.name = "Switch"
    switch.input_type = 'VECTOR'

    #node Index.003
    index_003 = atoms_and_bonds.nodes.new("GeometryNodeInputIndex")
    index_003.name = "Index.003"

    #node Math.014
    math_014 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_014.name = "Math.014"
    math_014.operation = 'MODULO'
    math_014.use_clamp = False
    #Value_001
    math_014.inputs[1].default_value = 2.0

    #node Named Attribute
    radius_attribute = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    radius_attribute.name = "Named Attribute"
    radius_attribute.data_type = 'FLOAT_VECTOR'
    #Name
    radius_attribute.inputs[0].default_value = "start_pos"

    #node Named Attribute.008
    named_attribute_008 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_008.name = "Named Attribute.008"
    named_attribute_008.data_type = 'FLOAT_VECTOR'
    #Name
    named_attribute_008.inputs[0].default_value = "end_pos"

    #node Realize Instances
    realize_instances = atoms_and_bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances.name = "Realize Instances"
    #Selection
    realize_instances.inputs[1].default_value = True
    #Realize All
    realize_instances.inputs[2].default_value = True
    #Depth
    realize_instances.inputs[3].default_value = 0

    #node Reroute.023
    reroute_023 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_023.name = "Reroute.023"
    reroute_023.socket_idname = "NodeSocketFloatDistance"
    #node Group Input.003
    group_input_003 = atoms_and_bonds.nodes.new("NodeGroupInput")
    group_input_003.name = "Group Input.003"
    group_input_003.hide = True

    #node Join Geometry.001
    join_geometry_001 = atoms_and_bonds.nodes.new("GeometryNodeJoinGeometry")
    join_geometry_001.name = "Join Geometry.001"

    #node Curve Circle.001
    curve_circle_001 = atoms_and_bonds.nodes.new("GeometryNodeCurvePrimitiveCircle")
    curve_circle_001.name = "Curve Circle.001"
    curve_circle_001.mode = 'RADIUS'
    #Resolution
    curve_circle_001.inputs[0].default_value = 32

    #node Delete Geometry
    delete_geometry = atoms_and_bonds.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.name = "Delete Geometry"
    delete_geometry.domain = 'POINT'
    delete_geometry.mode = 'ALL'

    #node Curve Line
    curve_line = atoms_and_bonds.nodes.new("GeometryNodeCurvePrimitiveLine")
    curve_line.name = "Curve Line"
    curve_line.hide = True
    curve_line.mode = 'POINTS'
    #Start
    curve_line.inputs[0].default_value = (0.0, 0.0, 0.0)
    #End
    curve_line.inputs[1].default_value = (0.0, 0.0, 1.0)

    #node Instance on Points
    instance_on_points = atoms_and_bonds.nodes.new("GeometryNodeInstanceOnPoints")
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
    compare_005 = atoms_and_bonds.nodes.new("FunctionNodeCompare")
    compare_005.name = "Compare.005"
    compare_005.data_type = 'FLOAT'
    compare_005.mode = 'ELEMENT'
    compare_005.operation = 'GREATER_THAN'

    #node Compare.006
    compare_006 = atoms_and_bonds.nodes.new("FunctionNodeCompare")
    compare_006.name = "Compare.006"
    compare_006.data_type = 'INT'
    compare_006.mode = 'ELEMENT'
    compare_006.operation = 'EQUAL'

    #node Boolean Math
    boolean_math = atoms_and_bonds.nodes.new("FunctionNodeBooleanMath")
    boolean_math.name = "Boolean Math"
    boolean_math.operation = 'OR'

    #node Sample Index.008
    sample_index_008 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_008.name = "Sample Index.008"
    sample_index_008.clamp = False
    sample_index_008.data_type = 'FLOAT'
    sample_index_008.domain = 'POINT'

    #node Sample Index.009
    sample_index_009 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_009.name = "Sample Index.009"
    sample_index_009.clamp = False
    sample_index_009.data_type = 'FLOAT'
    sample_index_009.domain = 'POINT'

    #node Reroute.030
    reroute_030 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_030.name = "Reroute.030"
    reroute_030.socket_idname = "NodeSocketFloat"
    #node Reroute.031
    reroute_031 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_031.name = "Reroute.031"
    reroute_031.socket_idname = "NodeSocketFloat"
    #node Reroute.018
    reroute_018 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_018.name = "Reroute.018"
    reroute_018.socket_idname = "NodeSocketGeometry"
    #node Named Attribute.010
    named_attribute_010 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_010.name = "Named Attribute.010"
    named_attribute_010.data_type = 'FLOAT'
    #Name
    named_attribute_010.inputs[0].default_value = "atom_radius"

    #node Sample Index.004
    sample_index_004 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_004.name = "Sample Index.004"
    sample_index_004.clamp = False
    sample_index_004.data_type = 'FLOAT_COLOR'
    sample_index_004.domain = 'POINT'

    #node Sample Index.005
    sample_index_005 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_005.name = "Sample Index.005"
    sample_index_005.clamp = False
    sample_index_005.data_type = 'FLOAT_COLOR'
    sample_index_005.domain = 'POINT'

    #node Reroute.016
    reroute_016 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_016.name = "Reroute.016"
    reroute_016.socket_idname = "NodeSocketGeometry"
    #node Named Attribute.002
    named_attribute_002 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_002.name = "Named Attribute.002"
    named_attribute_002.data_type = 'FLOAT_COLOR'
    #Name
    named_attribute_002.inputs[0].default_value = "color"

    #node Reroute.025
    reroute_025 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_025.name = "Reroute.025"
    reroute_025.socket_idname = "NodeSocketFloat"
    #node Reroute.026
    reroute_026 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_026.name = "Reroute.026"
    reroute_026.socket_idname = "NodeSocketFloat"
    #node Group Input.002
    group_input_002 = atoms_and_bonds.nodes.new("NodeGroupInput")
    group_input_002.name = "Group Input.002"

    #node Store Named Attribute.001
    store_named_attribute_001 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_001.name = "Store Named Attribute.001"
    store_named_attribute_001.hide = True
    store_named_attribute_001.data_type = 'FLOAT_VECTOR'
    store_named_attribute_001.domain = 'POINT'
    #Selection
    store_named_attribute_001.inputs[1].default_value = True
    #Name
    store_named_attribute_001.inputs[2].default_value = "end_pos"

    #node Store Named Attribute
    store_named_attribute = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute.name = "Store Named Attribute"
    store_named_attribute.hide = True
    store_named_attribute.data_type = 'FLOAT_VECTOR'
    store_named_attribute.domain = 'POINT'
    #Selection
    store_named_attribute.inputs[1].default_value = True
    #Name
    store_named_attribute.inputs[2].default_value = "start_pos"

    #node Store Named Attribute.003
    store_named_attribute_003 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_003.name = "Store Named Attribute.003"
    store_named_attribute_003.hide = True
    store_named_attribute_003.data_type = 'FLOAT_COLOR'
    store_named_attribute_003.domain = 'POINT'
    #Selection
    store_named_attribute_003.inputs[1].default_value = True
    #Name
    store_named_attribute_003.inputs[2].default_value = "end_col"

    #node Store Named Attribute.002
    store_named_attribute_002 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_002.name = "Store Named Attribute.002"
    store_named_attribute_002.hide = True
    store_named_attribute_002.data_type = 'FLOAT_COLOR'
    store_named_attribute_002.domain = 'POINT'
    #Selection
    store_named_attribute_002.inputs[1].default_value = True
    #Name
    store_named_attribute_002.inputs[2].default_value = "start_col"

    #node Store Named Attribute.006
    store_named_attribute_006 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_006.name = "Store Named Attribute.006"
    store_named_attribute_006.hide = True
    store_named_attribute_006.data_type = 'FLOAT'
    store_named_attribute_006.domain = 'POINT'
    #Selection
    store_named_attribute_006.inputs[1].default_value = True
    #Name
    store_named_attribute_006.inputs[2].default_value = "start_rad"

    #node Store Named Attribute.004
    store_named_attribute_004 = atoms_and_bonds.nodes.new("GeometryNodeStoreNamedAttribute")
    store_named_attribute_004.name = "Store Named Attribute.004"
    store_named_attribute_004.hide = True
    store_named_attribute_004.data_type = 'FLOAT'
    store_named_attribute_004.domain = 'POINT'
    #Selection
    store_named_attribute_004.inputs[1].default_value = True
    #Name
    store_named_attribute_004.inputs[2].default_value = "end_rad"

    #node Group Output
    group_output = atoms_and_bonds.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"
    group_output.is_active_output = True

    #node Vector Math
    vector_math = atoms_and_bonds.nodes.new("ShaderNodeVectorMath")
    vector_math.name = "Vector Math"
    vector_math.operation = 'DISTANCE'

    #node Math.002
    math_002 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_002.name = "Math.002"
    math_002.operation = 'ADD'
    math_002.use_clamp = False
    #Value_001
    math_002.inputs[1].default_value = 0.0

    #node Math.001
    math_001 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_001.name = "Math.001"
    math_001.operation = 'MULTIPLY'
    math_001.use_clamp = False

    #node Reroute.001
    reroute_001 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_001.name = "Reroute.001"
    reroute_001.socket_idname = "NodeSocketGeometry"
    #node Points
    points = atoms_and_bonds.nodes.new("GeometryNodePoints")
    points.name = "Points"
    #Count
    points.inputs[0].default_value = 1
    #Position
    points.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Radius
    points.inputs[2].default_value = 0.10000000149011612

    #node Instance on Points.002
    instance_on_points_002 = atoms_and_bonds.nodes.new("GeometryNodeInstanceOnPoints")
    instance_on_points_002.name = "Instance on Points.002"
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

    #node Named Attribute.003
    named_attribute_003 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_003.name = "Named Attribute.003"
    named_attribute_003.hide = True
    named_attribute_003.data_type = 'FLOAT_COLOR'
    #Name
    named_attribute_003.inputs[0].default_value = "start_col"

    #node Named Attribute.004
    named_attribute_004 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_004.name = "Named Attribute.004"
    named_attribute_004.hide = True
    named_attribute_004.data_type = 'FLOAT_COLOR'
    #Name
    named_attribute_004.inputs[0].default_value = "end_col"

    #node Domain Size.001
    domain_size_001 = atoms_and_bonds.nodes.new("GeometryNodeAttributeDomainSize")
    domain_size_001.name = "Domain Size.001"
    domain_size_001.component = 'POINTCLOUD'

    #node Math.013
    math_013 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_013.name = "Math.013"
    math_013.operation = 'POWER'
    math_013.use_clamp = False
    #Value_001
    math_013.inputs[1].default_value = 2.0

    #node Points.001
    points_001 = atoms_and_bonds.nodes.new("GeometryNodePoints")
    points_001.name = "Points.001"
    #Position
    points_001.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Radius
    points_001.inputs[2].default_value = 0.10000000149011612

    #node Realize Instances.001
    realize_instances_001 = atoms_and_bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_001.name = "Realize Instances.001"
    #Selection
    realize_instances_001.inputs[1].default_value = True
    #Realize All
    realize_instances_001.inputs[2].default_value = True
    #Depth
    realize_instances_001.inputs[3].default_value = 0

    #node Position
    position = atoms_and_bonds.nodes.new("GeometryNodeInputPosition")
    position.name = "Position"

    #node Reroute.002
    reroute_002 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_002.name = "Reroute.002"
    reroute_002.socket_idname = "NodeSocketGeometry"
    #node Sample Index
    sample_index = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index.name = "Sample Index"
    sample_index.clamp = False
    sample_index.data_type = 'FLOAT_VECTOR'
    sample_index.domain = 'POINT'

    #node Sample Index.001
    sample_index_001 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_001.name = "Sample Index.001"
    sample_index_001.clamp = False
    sample_index_001.data_type = 'FLOAT_VECTOR'
    sample_index_001.domain = 'POINT'

    #node Sample Index.002
    sample_index_002 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_002.name = "Sample Index.002"
    sample_index_002.clamp = False
    sample_index_002.data_type = 'FLOAT_VECTOR'
    sample_index_002.domain = 'POINT'
    #Value
    sample_index_002.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Index
    sample_index_002.inputs[2].default_value = 0

    #node Sample Index.003
    sample_index_003 = atoms_and_bonds.nodes.new("GeometryNodeSampleIndex")
    sample_index_003.name = "Sample Index.003"
    sample_index_003.clamp = False
    sample_index_003.data_type = 'FLOAT_VECTOR'
    sample_index_003.domain = 'POINT'
    #Value
    sample_index_003.inputs[1].default_value = (0.0, 0.0, 0.0)
    #Index
    sample_index_003.inputs[2].default_value = 0

    #node Reroute.003
    reroute_003 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_003.name = "Reroute.003"
    reroute_003.socket_idname = "NodeSocketGeometry"
    #node Reroute.009
    reroute_009 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_009.name = "Reroute.009"
    reroute_009.socket_idname = "NodeSocketFloat"
    #node Reroute.008
    reroute_008 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_008.name = "Reroute.008"
    reroute_008.socket_idname = "NodeSocketFloat"
    #node Reroute.012
    reroute_012 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_012.name = "Reroute.012"
    reroute_012.socket_idname = "NodeSocketFloat"
    #node Reroute.013
    reroute_013 = atoms_and_bonds.nodes.new("NodeReroute")
    reroute_013.name = "Reroute.013"
    reroute_013.socket_idname = "NodeSocketFloat"
    

    #node Merge by Distance
    merge_by_distance = atoms_and_bonds.nodes.new("GeometryNodeMergeByDistance")
    merge_by_distance.name = "Merge by Distance"
    merge_by_distance.mode = 'ALL'
    #Selection
    merge_by_distance.inputs[1].default_value = True
    #Distance
    merge_by_distance.inputs[2].default_value = 0.0010000020265579224

    #node Named Attribute.009
    named_attribute_009 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_009.name = "Named Attribute.009"
    named_attribute_009.data_type = 'FLOAT'
    #Name
    named_attribute_009.inputs[0].default_value = "start_rad"

    #node Named Attribute.011
    named_attribute_011 = atoms_and_bonds.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_011.name = "Named Attribute.011"
    named_attribute_011.hide = True
    named_attribute_011.data_type = 'FLOAT'
    #Name
    named_attribute_011.inputs[0].default_value = "end_rad"

    #node Math.003
    math_003 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_003.name = "Math.003"
    math_003.operation = 'MULTIPLY'
    math_003.use_clamp = False

    #node Compare.003
    compare_003 = atoms_and_bonds.nodes.new("FunctionNodeCompare")
    compare_003.name = "Compare.003"
    compare_003.data_type = 'FLOAT'
    compare_003.mode = 'ELEMENT'
    compare_003.operation = 'LESS_THAN'
    #B
    compare_003.inputs[1].default_value = 0.5

    #node Switch.003
    switch_003 = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
    switch_003.name = "Switch.003"
    switch_003.input_type = 'FLOAT'
    #False
    switch_003.inputs[1].default_value = 2.0
    #True
    switch_003.inputs[2].default_value = 1.0

    #node Compare.004
    compare_004 = atoms_and_bonds.nodes.new("FunctionNodeCompare")
    compare_004.name = "Compare.004"
    compare_004.data_type = 'FLOAT'
    compare_004.mode = 'ELEMENT'
    compare_004.operation = 'LESS_THAN'
    #B
    compare_004.inputs[1].default_value = 0.5

    #node Switch.004
    switch_004 = atoms_and_bonds.nodes.new("GeometryNodeSwitch")
    switch_004.name = "Switch.004"
    switch_004.input_type = 'FLOAT'
    #False
    switch_004.inputs[1].default_value = 2.0
    #True
    switch_004.inputs[2].default_value = 1.0

    #node Math
    math = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = 'ADD'
    math.use_clamp = False

    #node Math.006
    math_006 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_006.name = "Math.006"
    math_006.operation = 'MULTIPLY'
    math_006.use_clamp = False
    math_006.inputs[0].default_value = 1
    math_006.inputs[1].default_value = 1

    #node Math.007
    math_007 = atoms_and_bonds.nodes.new("ShaderNodeMath")
    math_007.name = "Math.007"
    math_007.operation = 'MULTIPLY'
    math_007.use_clamp = False
    math_007.inputs[0].default_value = 1
    math_007.inputs[1].default_value = 1


    #node Realize Instances.002
    realize_instances_002 = atoms_and_bonds.nodes.new("GeometryNodeRealizeInstances")
    realize_instances_002.name = "Realize Instances.002"
    #Selection
    realize_instances_002.inputs[1].default_value = True
    #Realize All
    realize_instances_002.inputs[2].default_value = True
    #Depth
    realize_instances_002.inputs[3].default_value = 0




    #Set parents
    reroute_015.parent = frame_008
    set_curve_radius.parent = frame_001
    set_material_006.parent = frame_001
    switch_005.parent = frame_001
    curve_to_mesh.parent = frame_001
    reroute_021.parent = frame_001
    curve_circle.parent = frame_001
    curve_to_mesh_001.parent = frame_001
    switch_001.parent = frame_016
    math_005.parent = frame_016
    index_002.parent = frame_016
    math_004.parent = frame_017
    math_011.parent = frame_017
    math_012.parent = frame_017
    index_001.parent = frame_017
    switch.parent = frame_013
    index_003.parent = frame_013
    math_014.parent = frame_013
    radius_attribute.parent = frame_013
    named_attribute_008.parent = frame_013
    realize_instances.parent = frame_008
    group_input_003.parent = frame_001
    curve_circle_001.parent = frame_001
    curve_line.parent = frame_008
    instance_on_points.parent = frame_008
    compare_005.parent = frame_018
    compare_006.parent = frame_018
    boolean_math.parent = frame_018
    sample_index_008.parent = frame_014
    sample_index_009.parent = frame_014
    reroute_030.parent = frame_014
    reroute_031.parent = frame_014
    reroute_018.parent = frame_014
    named_attribute_010.parent = frame_014
    sample_index_004.parent = frame_012
    sample_index_005.parent = frame_012
    reroute_016.parent = frame_012
    named_attribute_002.parent = frame_012
    reroute_025.parent = frame_012
    reroute_026.parent = frame_012
    store_named_attribute_001.parent = frame_019
    store_named_attribute.parent = frame_019
    store_named_attribute_003.parent = frame_019
    store_named_attribute_002.parent = frame_019
    store_named_attribute_006.parent = frame_019
    store_named_attribute_004.parent = frame_019
    vector_math.parent = frame_018
    math_002.parent = frame_018
    math_001.parent = frame_018
    named_attribute_003.parent = frame_016
    named_attribute_004.parent = frame_016
    domain_size_001.parent = frame_009
    math_013.parent = frame_009
    points_001.parent = frame_009
    position.parent = frame_004
    reroute_002.parent = frame_004
    sample_index.parent = frame_004
    sample_index_001.parent = frame_004
    sample_index_002.parent = frame_004
    sample_index_003.parent = frame_004
    reroute_003.parent = frame_004
    reroute_009.parent = frame_004
    reroute_008.parent = frame_004
    reroute_012.parent = frame_004
    reroute_013.parent = frame_004
    named_attribute_009.parent = frame_018
    named_attribute_011.parent = frame_018
    math_003.parent = frame_018
    switch_004.parent = frame_018
    math.parent = frame_018
    math_006.parent = frame_018
    math_007.parent = frame_018

    #Set locations
    frame_008.location = (-507.6666564941406, -1299.8333740234375)
    frame_001.location = (1305.0, -927.1666870117188)
    frame_016.location = (-943.6666870117188, -2357.166748046875)
    frame_005.location = (-3140.510498046875, -1904.2864990234375)
    frame_007.location = (-3135.767578125, -2116.841552734375)
    frame_017.location = (-3269.666748046875, -1758.5)
    frame_013.location = (-1041.0, -1779.1666259765625)
    frame_018.location = (-2862.333251953125, -1689.1666259765625)
    frame_014.location = (-2835.0, -3081.833251953125)
    frame_012.location = (-2519.0, -2315.166748046875)
    frame_019.location = (-1691.0, -1335.6666259765625)
    frame_009.location = (-3414.333251953125, -1455.8333740234375)
    frame_004.location = (-4057.217041015625, -2313.833251953125)
    reroute_006.location = (2685.3515625, -1283.796875)
    reroute_014.location = (-215.2777099609375, -1474.3909912109375)
    reroute_015.location = (609.5299072265625, -123.8023681640625)
    reroute_004.location = (-2608.47509765625, -1578.4334716796875)
    reroute_019.location = (1033.3616943359375, -1128.0413818359375)
    set_curve_radius.location = (420.4346923828125, -280.97491455078125)
    set_material_006.location = (1025.79833984375, -326.07916259765625)
    switch_005.location = (805.6875, -242.70440673828125)
    curve_to_mesh.location = (585.5281982421875, -349.49139404296875)
    reroute_021.location = (376.3699951171875, -339.73040771484375)
    curve_circle.location = (427.6068115234375, -418.33843994140625)
    curve_to_mesh_001.location = (590.4764404296875, -220.72076416015625)
    reroute_007.location = (-2633.1083984375, -1548.634765625)
    switch_001.location = (567.9396362304688, -39.8232421875)
    math_005.location = (346.04608154296875, -61.42822265625)
    index_002.location = (29.14959716796875, -200.949951171875)
    set_position.location = (332.7141418457031, -1319.991943359375)
    store_named_attribute_005.location = (719.0224609375, -1300.6259765625)
    reroute_017.location = (-1709.025146484375, -1579.417724609375)
    reroute_028.location = (-1700.907470703125, -1621.119384765625)
    math_004.location = (195.134521484375, -102.0716552734375)
    math_011.location = (28.73974609375, -39.683349609375)
    math_012.location = (178.1396484375, -349.00146484375)
    index_001.location = (30.08837890625, -228.886474609375)
    switch.location = (420.275634765625, -144.3544921875)
    index_003.location = (28.69287109375, -122.5767822265625)
    math_014.location = (225.453369140625, -39.560302734375)
    radius_attribute.location = (225.628173828125, -211.2418212890625)
    named_attribute_008.location = (222.60107421875, -365.7655029296875)
    realize_instances.location = (368.5834045410156, -39.3087158203125)
    reroute_023.location = (1182.2025146484375, -1182.3656005859375)
    group_input_003.location = (29.1612548828125, -130.22064208984375)
    join_geometry_001.location = (3597.325439453125, -158.23992919921875)
    curve_circle_001.location = (240.4376220703125, -39.6986083984375)
    delete_geometry.location = (-969.494140625, -1383.24267578125)
    curve_line.location = (28.972076416015625, -177.3492431640625)
    instance_on_points.location = (199.21963500976562, -54.6624755859375)
    compare_005.location = (1292.1741943359375, -82.3499755859375)
    compare_006.location = (997.9556884765625, -455.2247314453125)
    boolean_math.location = (1554.4017333984375, -170.4898681640625)
    sample_index_008.location = (515.2841796875, -39.532470703125)
    sample_index_009.location = (513.1025390625, -287.367431640625)
    reroute_030.location = (444.43798828125, -216.706298828125)
    reroute_031.location = (460.463134765625, -462.0830078125)
    reroute_018.location = (285.13818359375, -282.128173828125)
    named_attribute_010.location = (29.1884765625, -253.385498046875)
    sample_index_004.location = (274.17578125, -39.73193359375)
    sample_index_005.location = (271.994140625, -287.56689453125)
    reroute_016.location = (44.02978515625, -282.32763671875)
    named_attribute_002.location = (29.20751953125, -373.29541015625)
    reroute_025.location = (203.32958984375, -216.90576171875)
    reroute_026.location = (219.354736328125, -462.282470703125)
    group_input_002.location = (-2488.740478515625, -1668.6180419921875)
    store_named_attribute_001.location = (33.205078125, -76.582275390625)
    store_named_attribute.location = (30.2001953125, -44.3043212890625)
    store_named_attribute_003.location = (31.30126953125, -159.7994384765625)
    store_named_attribute_002.location = (29.202880859375, -116.4759521484375)
    store_named_attribute_006.location = (30.53125, -202.8370361328125)
    store_named_attribute_004.location = (32.62939453125, -246.1605224609375)
    group_output.location = (4869.4736328125, 111.9588623046875)
    vector_math.location = (686.83447265625, -39.308837890625)
    math_002.location = (857.7833251953125, -41.7724609375)
    math_001.location = (634.37109375, -175.628173828125)
    reroute_001.location = (-4828.6337890625, -1449.430908203125)
    points.location = (-5156.953125, -1573.485595703125)
    instance_on_points_002.location = (-4744.6337890625, -1556.60009765625)
    named_attribute_003.location = (348.67474365234375, -237.367919921875)
    named_attribute_004.location = (340.82489013671875, -277.072265625)
    domain_size_001.location = (29.085693359375, -128.757080078125)
    math_013.location = (240.842041015625, -39.6812744140625)
    points_001.location = (443.87060546875, -83.1571044921875)
    realize_instances_001.location = (-4544.4560546875, -1595.501953125)
    position.location = (642.7607421875, -401.120361328125)
    reroute_002.location = (33.833251953125, -346.368408203125)
    sample_index.location = (986.665771484375, -39.6015625)
    sample_index_001.location = (984.484130859375, -287.436767578125)
    sample_index_002.location = (986.665771484375, -39.6015625)
    sample_index_003.location = (984.484130859375, -287.436767578125)
    reroute_003.location = (845.536376953125, -346.0439453125)
    reroute_009.location = (946.107666015625, -167.010986328125)
    reroute_008.location = (934.743896484375, -368.193359375)
    reroute_012.location = (389.083251953125, -377.300537109375)
    reroute_013.location = (469.48779296875, -164.701416015625)
    merge_by_distance.location = (4291.6298828125, 160.23275756835938)
    named_attribute_009.location = (38.109375, -267.3577880859375)
    named_attribute_011.location = (29.117431640625, -441.4185791015625)
    math_003.location = (659.7685546875, -378.9822998046875)
    compare_003.location = (-2606.014892578125, -1813.4000244140625)
    switch_003.location = (-2401.684326171875, -1835.3531494140625)
    compare_004.location = (-2592.161865234375, -1989.02490234375)
    switch_004.location = (457.185302734375, -320.656005859375)
    math.location = (1019.6214599609375, -204.21923828125)
    math_006.location = (814.3336181640625, -204.5140380859375)
    math_007.location = (839.6668701171875, -366.2752685546875)
    realize_instances_002.location = (4685.3271484375, 180.95066833496094)

    #Set dimensions
    frame_008.width, frame_008.height = 643.36328125, 399.833251953125
    frame_001.width, frame_001.height = 1194.666748046875, 491.83331298828125
    frame_016.width, frame_016.height = 736.6666870117188, 330.0
    frame_005.width, frame_005.height = 169.232177734375, 50.3642578125
    frame_007.width, frame_007.height = 198.0, 47.47119140625
    frame_017.width, frame_017.height = 364.0, 519.166748046875
    frame_013.width, frame_013.height = 589.3333740234375, 510.5001220703125
    frame_018.width, frame_018.height = 1723.333251953125, 625.8333740234375
    frame_014.width, frame_014.height = 684.0, 501.166748046875
    frame_012.width, frame_012.height = 443.333251953125, 518.5
    frame_019.width, frame_019.height = 202.0, 299.5
    frame_009.width, frame_009.height = 612.66650390625, 285.1666259765625
    frame_004.width, frame_004.height = 1155.55029296875, 501.166748046875
    reroute_006.width, reroute_006.height = 14.5, 100.0
    reroute_014.width, reroute_014.height = 14.5, 100.0
    reroute_015.width, reroute_015.height = 14.5, 100.0
    reroute_004.width, reroute_004.height = 14.5, 100.0
    reroute_019.width, reroute_019.height = 14.5, 100.0
    set_curve_radius.width, set_curve_radius.height = 140.0, 100.0
    set_material_006.width, set_material_006.height = 140.0, 100.0
    switch_005.width, switch_005.height = 140.0, 100.0
    curve_to_mesh.width, curve_to_mesh.height = 140.0, 100.0
    reroute_021.width, reroute_021.height = 14.5, 100.0
    curve_circle.width, curve_circle.height = 140.0, 100.0
    curve_to_mesh_001.width, curve_to_mesh_001.height = 140.0, 100.0
    reroute_007.width, reroute_007.height = 14.5, 100.0
    switch_001.width, switch_001.height = 140.0, 100.0
    math_005.width, math_005.height = 140.0, 100.0
    index_002.width, index_002.height = 140.0, 100.0
    set_position.width, set_position.height = 140.0, 100.0
    store_named_attribute_005.width, store_named_attribute_005.height = 140.0, 100.0
    reroute_017.width, reroute_017.height = 14.5, 100.0
    reroute_028.width, reroute_028.height = 14.5, 100.0
    math_004.width, math_004.height = 140.0, 100.0
    math_011.width, math_011.height = 140.0, 100.0
    math_012.width, math_012.height = 140.0, 100.0
    index_001.width, index_001.height = 140.0, 100.0
    switch.width, switch.height = 140.0, 100.0
    index_003.width, index_003.height = 140.0, 100.0
    math_014.width, math_014.height = 140.0, 100.0
    radius_attribute.width, radius_attribute.height = 140.0, 100.0
    named_attribute_008.width, named_attribute_008.height = 140.0, 100.0
    realize_instances.width, realize_instances.height = 140.0, 100.0
    reroute_023.width, reroute_023.height = 14.5, 100.0
    group_input_003.width, group_input_003.height = 140.0, 100.0
    join_geometry_001.width, join_geometry_001.height = 140.0, 100.0
    curve_circle_001.width, curve_circle_001.height = 140.0, 100.0
    delete_geometry.width, delete_geometry.height = 140.0, 100.0
    curve_line.width, curve_line.height = 140.0, 100.0
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    compare_005.width, compare_005.height = 151.46044921875, 100.0
    compare_006.width, compare_006.height = 140.0, 100.0
    boolean_math.width, boolean_math.height = 140.0, 100.0
    sample_index_008.width, sample_index_008.height = 140.0, 100.0
    sample_index_009.width, sample_index_009.height = 140.0, 100.0
    reroute_030.width, reroute_030.height = 14.5, 100.0
    reroute_031.width, reroute_031.height = 14.5, 100.0
    reroute_018.width, reroute_018.height = 14.5, 100.0
    named_attribute_010.width, named_attribute_010.height = 140.0, 100.0
    sample_index_004.width, sample_index_004.height = 140.0, 100.0
    sample_index_005.width, sample_index_005.height = 140.0, 100.0
    reroute_016.width, reroute_016.height = 14.5, 100.0
    named_attribute_002.width, named_attribute_002.height = 140.0, 100.0
    reroute_025.width, reroute_025.height = 14.5, 100.0
    reroute_026.width, reroute_026.height = 14.5, 100.0
    group_input_002.width, group_input_002.height = 140.0, 100.0
    store_named_attribute_001.width, store_named_attribute_001.height = 140.0, 100.0
    store_named_attribute.width, store_named_attribute.height = 140.0, 100.0
    store_named_attribute_003.width, store_named_attribute_003.height = 140.0, 100.0
    store_named_attribute_002.width, store_named_attribute_002.height = 140.0, 100.0
    store_named_attribute_006.width, store_named_attribute_006.height = 140.0, 100.0
    store_named_attribute_004.width, store_named_attribute_004.height = 140.0, 100.0
    group_output.width, group_output.height = 140.0, 100.0
    vector_math.width, vector_math.height = 140.0, 100.0
    math_002.width, math_002.height = 100.0, 100.0
    math_001.width, math_001.height = 140.0, 100.0
    reroute_001.width, reroute_001.height = 14.5, 100.0
    points.width, points.height = 140.0, 100.0
    instance_on_points_002.width, instance_on_points_002.height = 140.0, 100.0
    named_attribute_003.width, named_attribute_003.height = 140.0, 100.0
    named_attribute_004.width, named_attribute_004.height = 140.0, 100.0
    domain_size_001.width, domain_size_001.height = 140.0, 100.0
    math_013.width, math_013.height = 140.0, 100.0
    points_001.width, points_001.height = 140.0, 100.0
    realize_instances_001.width, realize_instances_001.height = 140.0, 100.0
    position.width, position.height = 140.0, 100.0
    reroute_002.width, reroute_002.height = 14.5, 100.0
    sample_index.width, sample_index.height = 140.0, 100.0
    sample_index_001.width, sample_index_001.height = 140.0, 100.0
    sample_index_002.width, sample_index_002.height = 140.0, 100.0
    sample_index_003.width, sample_index_003.height = 140.0, 100.0
    reroute_003.width, reroute_003.height = 14.5, 100.0
    reroute_009.width, reroute_009.height = 14.5, 100.0
    reroute_008.width, reroute_008.height = 14.5, 100.0
    reroute_012.width, reroute_012.height = 14.5, 100.0
    reroute_013.width, reroute_013.height = 14.5, 100.0
    merge_by_distance.width, merge_by_distance.height = 140.0, 100.0
    named_attribute_009.width, named_attribute_009.height = 158.470458984375, 100.0
    named_attribute_011.width, named_attribute_011.height = 140.0, 100.0
    math_003.width, math_003.height = 140.0, 100.0
    compare_003.width, compare_003.height = 140.0, 100.0
    switch_003.width, switch_003.height = 140.0, 100.0
    compare_004.width, compare_004.height = 140.0, 100.0
    switch_004.width, switch_004.height = 140.0, 100.0
    math.width, math.height = 140.0, 100.0
    math_006.width, math_006.height = 140.0, 100.0
    math_007.width, math_007.height = 140.0, 100.0
    realize_instances_002.width, realize_instances_002.height = 140.0, 100.0

    #initialize atoms_from_verts links
    #realize_instances_002.Geometry -> group_output.Geometry
    atoms_and_bonds.links.new(realize_instances_002.outputs[0], group_output.inputs[0])
    #named_attribute.Attribute -> switch.False
    atoms_and_bonds.links.new(radius_attribute.outputs[0], switch.inputs[1])
    #compare_006.Result -> boolean_math.Boolean
    atoms_and_bonds.links.new(compare_006.outputs[0], boolean_math.inputs[1])
    #store_named_attribute_002.Geometry -> store_named_attribute_003.Geometry
    atoms_and_bonds.links.new(store_named_attribute_002.outputs[0], store_named_attribute_003.inputs[0])
    #named_attribute_002.Attribute -> sample_index_005.Value
    atoms_and_bonds.links.new(named_attribute_002.outputs[0], sample_index_005.inputs[1])
    #points_001.Points -> store_named_attribute.Geometry
    atoms_and_bonds.links.new(points_001.outputs[0], store_named_attribute.inputs[0])
    #store_named_attribute.Geometry -> store_named_attribute_001.Geometry
    atoms_and_bonds.links.new(store_named_attribute.outputs[0], store_named_attribute_001.inputs[0])
    #curve_circle.Curve -> curve_to_mesh.Profile Curve
    atoms_and_bonds.links.new(curve_circle.outputs[0], curve_to_mesh.inputs[1])
    #sample_index.Value -> reroute_007.Input
    atoms_and_bonds.links.new(sample_index.outputs[0], reroute_007.inputs[0])
    #sample_index_001.Value -> reroute_004.Input
    atoms_and_bonds.links.new(sample_index_001.outputs[0], reroute_004.inputs[0])
    #reroute_007.Output -> store_named_attribute.Value
    atoms_and_bonds.links.new(reroute_007.outputs[0], store_named_attribute.inputs[3])
    #reroute_004.Output -> store_named_attribute_001.Value
    atoms_and_bonds.links.new(reroute_004.outputs[0], store_named_attribute_001.inputs[3])
    #named_attribute_008.Attribute -> switch.True
    atoms_and_bonds.links.new(named_attribute_008.outputs[0], switch.inputs[2])
    #reroute_003.Output -> reroute_016.Input
    atoms_and_bonds.links.new(reroute_003.outputs[0], reroute_016.inputs[0])
    #curve_line.Curve -> instance_on_points.Instance
    atoms_and_bonds.links.new(curve_line.outputs[0], instance_on_points.inputs[2])
    #instance_on_points.Instances -> realize_instances.Geometry
    atoms_and_bonds.links.new(instance_on_points.outputs[0], realize_instances.inputs[0])
    #math_014.Value -> switch.Switch
    atoms_and_bonds.links.new(math_014.outputs[0], switch.inputs[0])
    #math_012.Value -> reroute_012.Input
    atoms_and_bonds.links.new(math_012.outputs[0], reroute_012.inputs[0])
    #named_attribute_003.Attribute -> switch_001.False
    atoms_and_bonds.links.new(named_attribute_003.outputs[0], switch_001.inputs[1])
    #index_003.Index -> math_014.Value
    atoms_and_bonds.links.new(index_003.outputs[0], math_014.inputs[0])
    #math_004.Value -> reroute_013.Input
    atoms_and_bonds.links.new(math_004.outputs[0], reroute_013.inputs[0])
    #switch_001.Output -> reroute_014.Input
    atoms_and_bonds.links.new(switch_001.outputs[0], reroute_014.inputs[0])
    #switch.Output -> reroute_015.Input
    atoms_and_bonds.links.new(switch.outputs[0], reroute_015.inputs[0])
    #named_attribute_004.Attribute -> switch_001.True
    atoms_and_bonds.links.new(named_attribute_004.outputs[0], switch_001.inputs[2])
    #position.Position -> sample_index.Value
    atoms_and_bonds.links.new(position.outputs[0], sample_index.inputs[1])
    #sample_index_004.Value -> reroute_017.Input
    atoms_and_bonds.links.new(sample_index_004.outputs[0], reroute_017.inputs[0])
    #reroute_015.Output -> set_position.Position
    atoms_and_bonds.links.new(reroute_015.outputs[0], set_position.inputs[2])
    #sample_index_005.Value -> reroute_028.Input
    atoms_and_bonds.links.new(sample_index_005.outputs[0], reroute_028.inputs[0])
    #named_attribute_002.Attribute -> sample_index_004.Value
    atoms_and_bonds.links.new(named_attribute_002.outputs[0], sample_index_004.inputs[1])
    #reroute_009.Output -> sample_index.Index
    atoms_and_bonds.links.new(reroute_009.outputs[0], sample_index.inputs[2])
    #reroute_008.Output -> sample_index_001.Index
    atoms_and_bonds.links.new(reroute_008.outputs[0], sample_index_001.inputs[2])
    #reroute_003.Output -> sample_index.Geometry
    atoms_and_bonds.links.new(reroute_003.outputs[0], sample_index.inputs[0])
    #reroute_014.Output -> store_named_attribute_005.Value
    atoms_and_bonds.links.new(reroute_014.outputs[0], store_named_attribute_005.inputs[3])
    #switch_005.Output -> set_material_006.Geometry
    atoms_and_bonds.links.new(switch_005.outputs[0], set_material_006.inputs[0])
    #reroute_016.Output -> sample_index_005.Geometry
    atoms_and_bonds.links.new(reroute_016.outputs[0], sample_index_005.inputs[0])
    #reroute_012.Output -> reroute_008.Input
    atoms_and_bonds.links.new(reroute_012.outputs[0], reroute_008.inputs[0])
    #reroute_026.Output -> sample_index_005.Index
    atoms_and_bonds.links.new(reroute_026.outputs[0], sample_index_005.inputs[2])
    #reroute_013.Output -> reroute_009.Input
    atoms_and_bonds.links.new(reroute_013.outputs[0], reroute_009.inputs[0])
    #reroute_025.Output -> sample_index_004.Index
    atoms_and_bonds.links.new(reroute_025.outputs[0], sample_index_004.inputs[2])
    #reroute_017.Output -> store_named_attribute_002.Value
    atoms_and_bonds.links.new(reroute_017.outputs[0], store_named_attribute_002.inputs[3])
    #reroute_016.Output -> sample_index_004.Geometry
    atoms_and_bonds.links.new(reroute_016.outputs[0], sample_index_004.inputs[0])
    #domain_size_001.Point Count -> math_013.Value
    atoms_and_bonds.links.new(domain_size_001.outputs[0], math_013.inputs[0])
    #index_002.Index -> math_005.Value
    atoms_and_bonds.links.new(index_002.outputs[0], math_005.inputs[0])
    #math_005.Value -> switch_001.Switch
    atoms_and_bonds.links.new(math_005.outputs[0], switch_001.inputs[0])
    #math_011.Value -> compare_006.A
    atoms_and_bonds.links.new(math_011.outputs[0], compare_006.inputs[2])
    #domain_size_001.Point Count -> math_011.Value
    atoms_and_bonds.links.new(domain_size_001.outputs[0], math_011.inputs[1])
    #index_001.Index -> math_011.Value
    atoms_and_bonds.links.new(index_001.outputs[0], math_011.inputs[0])
    #math_011.Value -> math_004.Value
    atoms_and_bonds.links.new(math_011.outputs[0], math_004.inputs[0])
    #index_001.Index -> math_012.Value
    atoms_and_bonds.links.new(index_001.outputs[0], math_012.inputs[0])
    #domain_size_001.Point Count -> math_012.Value
    atoms_and_bonds.links.new(domain_size_001.outputs[0], math_012.inputs[1])
    #reroute_006.Output -> join_geometry_001.Geometry
    atoms_and_bonds.links.new(reroute_006.outputs[0], join_geometry_001.inputs[0])
    #math_012.Value -> compare_006.B
    atoms_and_bonds.links.new(math_012.outputs[0], compare_006.inputs[3])
    #reroute_003.Output -> sample_index_001.Geometry
    atoms_and_bonds.links.new(reroute_003.outputs[0], sample_index_001.inputs[0])
    #math_013.Value -> points_001.Count
    atoms_and_bonds.links.new(math_013.outputs[0], points_001.inputs[0])
    #reroute_002.Output -> reroute_003.Input
    atoms_and_bonds.links.new(reroute_002.outputs[0], reroute_003.inputs[0])
    #set_material_006.Geometry -> reroute_006.Input
    atoms_and_bonds.links.new(set_material_006.outputs[0], reroute_006.inputs[0])
    #store_named_attribute_001.Geometry -> store_named_attribute_002.Geometry
    atoms_and_bonds.links.new(store_named_attribute_001.outputs[0], store_named_attribute_002.inputs[0])
    #reroute_028.Output -> store_named_attribute_003.Value
    atoms_and_bonds.links.new(reroute_028.outputs[0], store_named_attribute_003.inputs[3])
    #position.Position -> sample_index_001.Value
    atoms_and_bonds.links.new(position.outputs[0], sample_index_001.inputs[1])
    #realize_instances.Geometry -> set_position.Geometry
    atoms_and_bonds.links.new(realize_instances.outputs[0], set_position.inputs[0])
    #set_position.Geometry -> store_named_attribute_005.Geometry
    atoms_and_bonds.links.new(set_position.outputs[0], store_named_attribute_005.inputs[0])
    #reroute_021.Output -> set_curve_radius.Curve
    atoms_and_bonds.links.new(reroute_021.outputs[0], set_curve_radius.inputs[0])
    #set_curve_radius.Curve -> curve_to_mesh.Curve
    atoms_and_bonds.links.new(set_curve_radius.outputs[0], curve_to_mesh.inputs[0])
    #curve_to_mesh.Mesh -> switch_005.True
    atoms_and_bonds.links.new(curve_to_mesh.outputs[0], switch_005.inputs[2])
    #store_named_attribute_005.Geometry -> reroute_021.Input
    atoms_and_bonds.links.new(store_named_attribute_005.outputs[0], reroute_021.inputs[0])
    #reroute_021.Output -> curve_to_mesh_001.Curve
    atoms_and_bonds.links.new(reroute_021.outputs[0], curve_to_mesh_001.inputs[0])
    #curve_circle_001.Curve -> curve_to_mesh_001.Profile Curve
    atoms_and_bonds.links.new(curve_circle_001.outputs[0], curve_to_mesh_001.inputs[1])
    #curve_to_mesh_001.Mesh -> switch_005.False
    atoms_and_bonds.links.new(curve_to_mesh_001.outputs[0], switch_005.inputs[1])
    #math_002.Value -> compare_005.A
    atoms_and_bonds.links.new(math_002.outputs[0], compare_005.inputs[0])
    #reroute_004.Output -> vector_math.Vector
    atoms_and_bonds.links.new(reroute_004.outputs[0], vector_math.inputs[1])
    #reroute_007.Output -> vector_math.Vector
    atoms_and_bonds.links.new(reroute_007.outputs[0], vector_math.inputs[0])
    #reroute_009.Output -> reroute_025.Input
    atoms_and_bonds.links.new(reroute_009.outputs[0], reroute_025.inputs[0])
    #reroute_008.Output -> reroute_026.Input
    atoms_and_bonds.links.new(reroute_008.outputs[0], reroute_026.inputs[0])
    #reroute_001.Output -> instance_on_points_002.Points
    atoms_and_bonds.links.new(reroute_001.outputs[0], instance_on_points_002.inputs[0])
    #points.Points -> instance_on_points_002.Instance
    atoms_and_bonds.links.new(points.outputs[0], instance_on_points_002.inputs[2])
    #instance_on_points_002.Instances -> realize_instances_001.Geometry
    atoms_and_bonds.links.new(instance_on_points_002.outputs[0], realize_instances_001.inputs[0])
    #######################put points to output ####
    atoms_and_bonds.links.new(instance_on_points_002.outputs[0], join_geometry_001.inputs[0])
    #group_input_003.Radius -> curve_circle_001.Radius
    atoms_and_bonds.links.new(group_input_003.outputs[2], curve_circle_001.inputs[4])
    #boolean_math.Boolean -> delete_geometry.Selection
    atoms_and_bonds.links.new(boolean_math.outputs[0], delete_geometry.inputs[1])
    #delete_geometry.Geometry -> instance_on_points.Points
    atoms_and_bonds.links.new(delete_geometry.outputs[0], instance_on_points.inputs[0])
    #compare_005.Result -> boolean_math.Boolean
    atoms_and_bonds.links.new(compare_005.outputs[0], boolean_math.inputs[0])
    #store_named_attribute_006.Geometry -> store_named_attribute_004.Geometry
    atoms_and_bonds.links.new(store_named_attribute_006.outputs[0], store_named_attribute_004.inputs[0])
    #store_named_attribute_003.Geometry -> store_named_attribute_006.Geometry
    atoms_and_bonds.links.new(store_named_attribute_003.outputs[0], store_named_attribute_006.inputs[0])
    #store_named_attribute_004.Geometry -> delete_geometry.Geometry
    atoms_and_bonds.links.new(store_named_attribute_004.outputs[0], delete_geometry.inputs[0])
    #named_attribute_010.Attribute -> sample_index_009.Value
    atoms_and_bonds.links.new(named_attribute_010.outputs[0], sample_index_009.inputs[1])
    #named_attribute_010.Attribute -> sample_index_008.Value
    atoms_and_bonds.links.new(named_attribute_010.outputs[0], sample_index_008.inputs[1])
    #reroute_018.Output -> sample_index_009.Geometry
    atoms_and_bonds.links.new(reroute_018.outputs[0], sample_index_009.inputs[0])
    #reroute_031.Output -> sample_index_009.Index
    atoms_and_bonds.links.new(reroute_031.outputs[0], sample_index_009.inputs[2])
    #reroute_030.Output -> sample_index_008.Index
    atoms_and_bonds.links.new(reroute_030.outputs[0], sample_index_008.inputs[2])
    #reroute_018.Output -> sample_index_008.Geometry
    atoms_and_bonds.links.new(reroute_018.outputs[0], sample_index_008.inputs[0])
    #reroute_016.Output -> reroute_018.Input
    atoms_and_bonds.links.new(reroute_016.outputs[0], reroute_018.inputs[0])
    #reroute_025.Output -> reroute_030.Input
    atoms_and_bonds.links.new(reroute_025.outputs[0], reroute_030.inputs[0])
    #reroute_026.Output -> reroute_031.Input
    atoms_and_bonds.links.new(reroute_026.outputs[0], reroute_031.inputs[0])
    #vector_math.Value -> math_002.Value
    atoms_and_bonds.links.new(vector_math.outputs[1], math_002.inputs[0])
    #realize_instances_001.Geometry -> domain_size_001.Geometry
    atoms_and_bonds.links.new(realize_instances_001.outputs[0], domain_size_001.inputs[0])
    #realize_instances_001.Geometry -> reroute_002.Input
    atoms_and_bonds.links.new(realize_instances_001.outputs[0], reroute_002.inputs[0])

    #removed shade smooth
    atoms_and_bonds.links.new(merge_by_distance.outputs[0], realize_instances_002.inputs[0])
                              
                              
    #join_geometry_001.Geometry -> merge_by_distance.Geometry
    atoms_and_bonds.links.new(join_geometry_001.outputs[0], merge_by_distance.inputs[0])
    #sample_index_008.Value -> store_named_attribute_006.Value
    atoms_and_bonds.links.new(sample_index_008.outputs[0], store_named_attribute_006.inputs[3])
    #sample_index_009.Value -> store_named_attribute_004.Value
    atoms_and_bonds.links.new(sample_index_009.outputs[0], store_named_attribute_004.inputs[3])
    #named_attribute_009.Attribute -> compare_003.A
    atoms_and_bonds.links.new(named_attribute_009.outputs[0], compare_003.inputs[0])
    #compare_003.Result -> switch_003.Switch
    atoms_and_bonds.links.new(compare_003.outputs[0], switch_003.inputs[0])
    #switch_003.Output -> math_001.Value
    atoms_and_bonds.links.new(switch_003.outputs[0], math_001.inputs[0])
    #compare_004.Result -> switch_004.Switch
    atoms_and_bonds.links.new(compare_004.outputs[0], switch_004.inputs[0])
    #named_attribute_011.Attribute -> compare_004.A
    atoms_and_bonds.links.new(named_attribute_011.outputs[0], compare_004.inputs[0])
    #switch_004.Output -> math_003.Value
    atoms_and_bonds.links.new(switch_004.outputs[0], math_003.inputs[0])
    #named_attribute_009.Attribute -> math_001.Value
    atoms_and_bonds.links.new(named_attribute_009.outputs[0], math_001.inputs[1])
    #named_attribute_011.Attribute -> math_003.Value
    atoms_and_bonds.links.new(named_attribute_011.outputs[0], math_003.inputs[1])
    #math_006.Value -> math.Value
    atoms_and_bonds.links.new(math_006.outputs[0], math.inputs[0])
    #math_007.Value -> math.Value
    atoms_and_bonds.links.new(math_007.outputs[0], math.inputs[1])
    #math_001.Value -> math_006.Value
    atoms_and_bonds.links.new(math_001.outputs[0], math_006.inputs[0])
    #group_input_002.B -> math_006.Value
    atoms_and_bonds.links.new(group_input_002.outputs[1], math_006.inputs[1])
    #math_003.Value -> math_007.Value
    atoms_and_bonds.links.new(math_003.outputs[0], math_007.inputs[0])
    #group_input_002.B -> math_007.Value
    atoms_and_bonds.links.new(group_input_002.outputs[1], math_007.inputs[1])
    #math.Value -> compare_005.B
    atoms_and_bonds.links.new(math.outputs[0], compare_005.inputs[1])

    #reroute_001.Output -> join_geometry_001.Geometry
    atoms_and_bonds.links.new(reroute_001.outputs[0], join_geometry_001.inputs[0])

    atoms_and_bonds.links.new(color_attribute.outputs[0], reroute_001.inputs[0])

    atoms_and_bonds.links.new(join_geometry_atoms.outputs[0], join_geometry_001.inputs[0])

    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='NODES')
    obj.modifiers[modifier].node_group = atoms_and_bonds
    #attach materials to atoms  
    for number in set(atoms.get_atomic_numbers()):
        sym = chemical_symbols[number]
        obj.data.materials.append(bpy.data.materials[sym])
    obj.data.materials.append(bondmat)
    return  atoms_and_bonds
