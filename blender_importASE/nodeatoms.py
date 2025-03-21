import bpy, mathutils
from ase.data import covalent_radii, chemical_symbols

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
    if "radius" not in mesh.attributes:
        mesh.attributes.new(name="radius", type='FLOAT', domain='POINT')

    element = mesh.attributes["element"].data
    rad = mesh.attributes["radius"].data

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
        bpy.data.scenes['Scene'].animall_properties.key_point_location = True
        vertx=obj.data.vertices
        for n,frame in enumerate(trajectory):
            bpy.data.scenes['Scene'].frame_current=n
            for nv,v in enumerate(vertx):
                v.co=frame.positions[nv]
            
            bpy.context.view_layer.update()
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.view3d.insert_keyframe_animall()
            bpy.ops.object.mode_set(mode='OBJECT')

    return(obj, mesh)
#initialize set_atoms node group
def set_atoms_node_group():
    if 'set_atoms' in bpy.data.node_groups:
        return
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
    #Segments
    uv_sphere.inputs[0].default_value = 32
    #Rings
    uv_sphere.inputs[1].default_value = 16
    #Radius
    uv_sphere.inputs[2].default_value = 1.0


    shade= set_atoms.nodes.new("GeometryNodeSetShadeSmooth")
    shade.name = "Shade Smooth"
    shade.inputs[2].default_value = True
    shade.location = (350, 0)
    shade.width, shade.height = 140.0, 100.0

    



    #Set locations
    group_output.location = (500.13250732421875, 0.0)
    group_input.location = (-435.13250732421875, 0.0)
    instance_on_points.location = (235.13250732421875, 30.57642364501953)
    uv_sphere.location = (-241.16090393066406, -142.24588012695312)

    #Set dimensions
    group_output.width, group_output.height = 140.0, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    instance_on_points.width, instance_on_points.height = 140.0, 100.0
    uv_sphere.width, uv_sphere.height = 140.0, 100.0

    #initialize set_atoms links
    #group_input.Points -> instance_on_points.Points
    set_atoms.links.new(group_input.outputs[0], instance_on_points.inputs[0])
    #group_input.Scale -> instance_on_points.Scale
    set_atoms.links.new(group_input.outputs[2], instance_on_points.inputs[6])
    #instance_on_points.Instances -> group_output.Instances
    set_atoms.links.new(instance_on_points.outputs[0],shade.inputs[0])
    set_atoms.links.new(shade.outputs[0], group_output.inputs[0])
    #uv_sphere.Mesh -> instance_on_points.Instance
    set_atoms.links.new(uv_sphere.outputs[0], instance_on_points.inputs[2])
    #group_input.Selection -> instance_on_points.Selection
    set_atoms.links.new(group_input.outputs[1], instance_on_points.inputs[1])
    return set_atoms


#initialize atoms_from_verts node group
def atoms_from_verts_node_group(atoms,name,animate=True):
    if animate:
        trajectory=atoms
        atoms=trajectory[0]
    
    atoms_from_verts = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = f"atoms_from_verts_{atoms.get_chemical_formula()}")

    atoms_from_verts.color_tag = 'NONE'
    atoms_from_verts.description = ""
    atoms_from_verts.default_group_node_width = 140
    

    atoms_from_verts.is_modifier = True

    #atoms_from_verts interface
    #Socket Geometry
    geometry_socket = atoms_from_verts.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_1 = atoms_from_verts.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_1.attribute_domain = 'POINT'


    #initialize atoms_from_verts nodes
    #node Group Input
    group_input_1 = atoms_from_verts.nodes.new("NodeGroupInput")
    group_input_1.name = "Group Input"

    #node Group Output
    group_output_1 = atoms_from_verts.nodes.new("NodeGroupOutput")
    group_output_1.name = "Group Output"
    group_output_1.is_active_output = True

    

    #node Named Attribute
    named_attribute = atoms_from_verts.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute.label = "radius_attribute"
    named_attribute.name = "Named Attribute"
    named_attribute.data_type = 'FLOAT'
    #Name
    named_attribute.inputs[0].default_value = "radius"

    #node Named Attribute.001
    named_attribute_001 = atoms_from_verts.nodes.new("GeometryNodeInputNamedAttribute")
    named_attribute_001.label = "element_attribute"
    named_attribute_001.name = "Named Attribute.001"
    named_attribute_001.data_type = 'INT'
    #Name
    named_attribute_001.inputs[0].default_value = "element"

    #node Join Geometry
    join_geometry = atoms_from_verts.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"
    
    #node Compare Elements
    for n,number in enumerate(set(atoms.get_atomic_numbers())):
        sym=chemical_symbols[number]
        compare = atoms_from_verts.nodes.new("FunctionNodeCompare")
        compare.label = f"is_element_{number}"
        compare.name = "Compare"
        compare.data_type = 'INT'
        compare.mode = 'ELEMENT'
        compare.operation = 'EQUAL'
        compare.inputs[3].default_value = number
        compare.location = (-834.4930419921875, -100+200*n)
        #node Group
        group = atoms_from_verts.nodes.new("GeometryNodeGroup")
        group.label = f"set_atoms_{number}"
        group.name = "Group"
        group.node_tree = bpy.data.node_groups['set_atoms']
        group.location = (-458.2143249511719, -100+200*n)
        compare.width, compare.height = 140.0, 100.0
        group.width, group.height = 140.0, 100.0

        #node Set Material
        set_material = atoms_from_verts.nodes.new("GeometryNodeSetMaterial")
        set_material.name = f"Set Material - {sym}"
        #Selection
        set_material.inputs[1].default_value = True
        if sym in bpy.data.materials:
            set_material.inputs[2].default_value = bpy.data.materials[sym]
        set_material.location = (-268.2143249511719, -100+140*n)
        set_material.width, set_material.height = 140.0, 100.0
        
        atoms_from_verts.links.new(compare.outputs[0], group.inputs[1])
        atoms_from_verts.links.new(group_input_1.outputs[0], group.inputs[0])
        atoms_from_verts.links.new(named_attribute.outputs[0], group.inputs[2])
        atoms_from_verts.links.new(named_attribute_001.outputs[0], compare.inputs[2])
        atoms_from_verts.links.new(group.outputs[0], set_material.inputs[0])
        atoms_from_verts.links.new(set_material.outputs[0], join_geometry.inputs[0])
        
    
    
    #node Frame
    frame = atoms_from_verts.nodes.new("NodeFrame")
    frame.label = " "
    frame.name = "Frame"
    frame.label_size = 20
    frame.shrink = True

    atoms_from_verts.links.new(join_geometry.outputs[0], group_output_1.inputs[0])



    #Set locations
    group_input_1.location = (-650, 63.91105270385742)
    group_output_1.location = (111.78567504882812, 0.4488945007324219)
    named_attribute.location = (-650, -112.09099578857422)
    named_attribute_001.location = (-1024.4930419921875, -66.59099578857422)
    join_geometry.location = (-70.63062286376953, 1.9488945007324219)
    frame.location = (712.62744140625, -409.5729064941406)

    #Set dimensions
    group_input_1.width, group_input_1.height = 140.0, 100.0
    group_output_1.width, group_output_1.height = 140.0, 100.0
    named_attribute.width, named_attribute.height = 140.0, 100.0
    named_attribute_001.width, named_attribute_001.height = 140.0, 100.0
    
    join_geometry.width, join_geometry.height = 140.0, 100.0
    frame.width, frame.height = 150.0, 100.0
    if animate:
        obj,mesh=read_structure(trajectory,name,animate=animate)
    else:
        obj,mesh=read_structure(atoms,name,animate=animate)
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='NODES')
    obj.modifiers['GeometryNodes'].node_group = atoms_from_verts
    return atoms_from_verts
