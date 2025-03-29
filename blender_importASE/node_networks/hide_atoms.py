from ase.data import chemical_symbols
import bpy
def hide_node_group(atoms):
    hide_atoms = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "hide atoms")

    hide_atoms.color_tag = 'NONE'
    hide_atoms.description = ""
    hide_atoms.default_group_node_width = 140
    

    hide_atoms.is_modifier = True

    
    group_input=hide_atoms.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"
    
    #node Group Output
    group_output = hide_atoms.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"
    group_output.is_active_output = True

    element_attribute = hide_atoms.nodes.new("GeometryNodeInputNamedAttribute")
    element_attribute.name = "Attribute"
    element_attribute.data_type = 'INT'
    element_attribute.inputs[0].default_value = "element"
#Socket Geometry
    geometry_socket = hide_atoms.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'
    #Socket Geometry
    geometry_socket = hide_atoms.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    old_switch = group_input
    for n,number in enumerate(set(atoms.get_atomic_numbers())):

        sym=chemical_symbols[number]
        #Socket cutoff_H
        cutoff_element_socket = hide_atoms.interface.new_socket(name = f"cutoff_{sym}", in_out='INPUT', socket_type = 'NodeSocketBool')
        cutoff_element_socket.default_value = False
        cutoff_element_socket.attribute_domain = 'POINT'
        
    
        is_element = hide_atoms.nodes.new("FunctionNodeCompare")
        is_element.label = f"is_element_{sym}"
        is_element.name = "Compare"
        is_element.data_type = 'INT'
        is_element.mode = 'ELEMENT'
        is_element.operation = 'EQUAL'
        is_element.inputs[3].default_value = number


        

        #node Delete Geometry
        delete_geometry_element = hide_atoms.nodes.new("GeometryNodeDeleteGeometry")
        delete_geometry_element.name = f"Delete Geometry_{sym}"
        delete_geometry_element.domain = 'POINT'
        delete_geometry_element.mode = 'ALL'
        
        
        #node switch
        switch_element=hide_atoms.nodes.new("GeometryNodeSwitch")
        switch_element.name = f"Switch_{sym}"
        switch_element.label = f"{sym}"
        switch_element.input_type = 'GEOMETRY'
        switch_element.inputs[0].default_value = False
        
        hide_atoms.links.new(is_element.outputs[0], delete_geometry_element.inputs[1])
        hide_atoms.links.new(element_attribute.outputs[0], is_element.inputs[2])
        
        hide_atoms.links.new(group_input.outputs[1+n], switch_element.inputs[0])
        hide_atoms.links.new(delete_geometry_element.outputs[0], switch_element.inputs[2])

        hide_atoms.links.new(old_switch.outputs[0],switch_element.inputs[1])
        hide_atoms.links.new(old_switch.outputs[0],delete_geometry_element.inputs[0])

        #supercell.links.new(old_delete.outputs[0], switch_element.inputs[2])
        if n == len(set(atoms.get_atomic_numbers()))-1:
            hide_atoms.links.new(switch_element.outputs[0], group_output.inputs[0])
        old_switch=switch_element
        ## old_delete=delete_geometry_element
        is_element.width, is_element.height = 140.0, 100.0
        switch_element.width, switch_element.height = 140.0, 100.0
        is_element.location = (500, 1000-200*n)
        delete_geometry_element.location = (800, 1000-200*n)
        switch_element.location = (1100, 1000-200*n)
    element_attribute.location = (200, 800)
    group_input.location = (0, 1000)
    group_output.location = (1400, 1000)
    #hide_atoms.links.new(realize_instances_beforevectorcutoff.outputs[0], switch_vectorcutoff.inputs[1])
    #hide_atoms.links.new(group_input_vectorcutoff.outputs[7], switch_vectorcutoff.inputs[0])


    return hide_atoms


def hide_atoms(obj, atoms,modifier='GeometryNodes'):
     
    hide=hide_node_group(atoms)
    
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='NODES')
    obj.modifiers[modifier].node_group = hide
    return
