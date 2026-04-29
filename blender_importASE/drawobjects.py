import bpy
import ase
from ase.data import covalent_radii, vdw_alvarez
import numpy as np
import ase.neighborlist


def draw_atoms(atoms, scale=1,resolution=16, representation="Balls'n'Sticks"):
    cnt = 0
    list_of_atoms=[]
    # bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0, 0, 0), segments=resolution, ring_count=resolution)
    if bpy.app.version[1] == 0: #use_auto_smoot dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    for n, atom in enumerate(atoms):
        ob = sphere.copy()
        ob.data = sphere.data.copy()
        ob.location = atom.position
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].name = atom.symbol
        if representation == "Balls'n'Sticks":
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [covalent_radii[
                                                                                               atom.number]  * scale, ] * 3
        elif representation == 'Licorice':
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [0.1] * 3
        elif representation == 'VDW':
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [vdw_alvarez.vdw_radii[
                                                                                               atom.number]] * 3
        # sprint(bpy.data.node_groups)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
            bpy.data.materials[atom.symbol])
        cnt += 1
        #full_object_name = bpy.utils.object.full_name(ob)
        #print(full_object_name)
        #list_of_atoms.append(full_object_name)
        #print(ob)
        list_of_atoms.append(ob)
        bpy.ops.object.transform_apply(location=False,rotation=True,scale=True)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_sphere'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_atoms


def draw_bonds(atoms,resolution=16):
    list_of_bonds=[]
    nl = ase.neighborlist.NeighborList([covalent_radii[atomic_number] * 0.9 for atomic_number in atoms.numbers],
                                       self_interaction=False, bothways=True)
    nl.update(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='bonds')
    except Exception:
        None
    # bpy.ops.surface.primitive_nurbs_surface_cylinder_add(radius=1.0, enter_editmode=False, align='WORLD',
    # location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bond = create_half_bond(resolution=resolution)
    cnt = 0
    for atom in atoms:
        if nl.get_neighbors(atom.index)[0].size > 0:
            neighbors, offsets = nl.get_neighbors(atom.index)
            for neighbor, offset in zip(neighbors, offsets):
                displacements = [0.5 * (atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)),
                                 0.5 * (atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor])]
                for n, displacement in enumerate(displacements):
                    if n == 0:
                        location = atom.position + (displacement / 2)
                    # else:
                    # location=atoms[neighbor].position + (displacement/2)
                    distance = atoms.get_distance(atom.index, neighbor, mic=True) / 2
                    ob = bond.copy()
                    ob.data = bond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
                    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}'
                    bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
                        bpy.data.materials[f'{atom.symbol}-bond'])
                    ob.location = location
                    ob.scale = (1, 1, 1)
                    # if n == 0:
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacement[1], displacement[0])
                    theta = np.arccos(displacement[2] / distance)
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    list_of_bonds.append(ob)
                    bpy.ops.object.transform_apply(location=False,rotation=True,scale=True)
                    break
                cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_bond'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_bonds,nl


def draw_bonds_new(atoms,resolution=16):
    #print("Using Longbond mech")
    list_of_bonds=[]
    bondlengths=[] # for animation
    nl = ase.neighborlist.NeighborList([covalent_radii[atomic_number] * 0.9 for atomic_number in atoms.numbers],
                                       self_interaction=False, bothways=False)
    nl.update(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='bonds')
    except Exception:
        pass
    # create half bond
    hbond = create_half_bond(resolution=resolution)
    bond = create_full_bond(resolution=resolution)
    cnt = 0
    cell = atoms.get_cell()
    # make a list of all bonds
    for atom in atoms:
        if nl.get_neighbors(atom.index)[0].size > 0:
            neighbors, offsets = nl.get_neighbors(atom.index)
            for neighbor, offset in zip(neighbors, offsets):
                neighbor_position = atoms.positions[neighbor] + offset.dot(atoms.cell)
                # Hier liegt das problem: die Zelle ist nicht 1x1x1
                #print(f'PBC: {atoms.pbc}')
                if atoms.pbc.all():
                    #print(atoms.pbc, 'need to check for unit cell')
                    is_same_unit_cell = is_inside_cell(neighbor_position, cell)
                    if is_same_unit_cell:
                        dis=round(atoms.get_distance(atom.index,neighbor,mic=False),3)
                        dismin=round(atoms.get_distance(atom.index,neighbor,mic=True),3)
                        if dis != dismin:
                            #print(dis,dismin)
                            is_same_unit_cell = False
                else:
                    is_same_unit_cell = True
                if is_same_unit_cell:
                    # print("Longbond")
                    # Draw one complex bond between the atoms, assign two materials if necessary
                    displacements = [0.5 * (atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)),
                                     0.5 * (atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor])]
                    location = atom.position + (displacements[0] / 2)
                    distance = atoms.get_distance(atom.index, neighbor, mic=True)
                    ob = bond.copy()
                    ob.data = bond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
                    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}+{atoms[neighbor].symbol}{neighbor}-{atom.symbol}{atom.index}'
                    # This still needs changing for colorbond support
                    assign_to_longbond(ob, f'{atom.symbol}-bond', f'{atoms[neighbor].symbol}-bond')
                    ob.location = location
                    ob.scale = (1, 1, 1)
                    # if n == 0:
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / (distance / 2))
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    list_of_bonds.append(ob)#anim
                    bondlengths.append('long')#anim
                    bpy.ops.object.transform_apply(location=False,rotation=True,scale=True) 
                else:
                    # print("Shortbond")
                    # create two bon fragments on either end
                    # atom to neighbor
                    displacements = [0.5 * (atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)),
                                     0.5 * (atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor])]
                    location = atom.position + (displacements[0] / 2)
                    distance = atoms.get_distance(atom.index, neighbor, mic=True) / 2
                    ob = hbond.copy()
                    ob.data = hbond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
                    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}'
                    bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
                        bpy.data.materials[f'{atom.symbol}-bond'])
                    ob.location = location
                    ob.scale = (1, 1, 1)
                    # if n == 0:
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / distance)
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    list_of_bonds.append(ob) #anim
                    bondlengths.append('short1') #anim
                    # neighbor to atom
                    displacements = [0.5 * (atom.position - atoms.positions[neighbor] - np.dot(offset, atoms.cell)),
                                     0.5 * (atoms.positions[neighbor] + np.dot(offset, atoms.cell) + atom.position)]
                    location = atoms.positions[neighbor] + (displacements[0] / 2)
                    distance = atoms.get_distance(atom.index, neighbor, mic=True) / 2
                    ob = hbond.copy()
                    ob.data = hbond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
                    ob.name = f'{atoms[neighbor].symbol}{neighbor}-{atom.symbol}{atom.index}'
                    # This still needs changing for colorbond support
                    bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
                        bpy.data.materials[f'{atom.symbol}-bond'])
                    ob.location = location
                    ob.scale = (1, 1, 1)
                    # if n == 0:
                    ob.dimensions = (0.2, 0.2, distance)
                    phi = np.arctan2(displacements[0][1], displacements[0][0])
                    theta = np.arccos(displacements[0][2] / distance)
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    list_of_bonds.append(ob)#anim
                    bondlengths.append('short2')
                    bpy.ops.object.transform_apply(location=False,rotation=True,scale=True)
            cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_bond'].select_set(True)
    bpy.data.objects['ref_bondx2'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_bonds,nl,bondlengths


def draw_unit_cell(atoms):
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='cell')
    except Exception:
        None

    # SETUP MATERIAL
    matu = bpy.data.materials.new(name='unit_cell')
    matu.use_nodes = True
    tu = matu.node_tree
    su = tu.nodes['Principled BSDF']
    COL = (0.1, 0.1, 0.1, 1)
    su.inputs[0].default_value = COL
    bpy.ops.mesh.primitive_cylinder_add(vertices=16)
    if bpy.app.version[1] == 0: #use_auto_smoot dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    cell = bpy.context.object
    cell.name = 'ref_cell'
    X = [0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    Y = [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    Z = [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]
    for n in range(1, len(X)):
        pos1 = np.array([X[n - 1], Y[n - 1], Z[n - 1]])
        pos2 = np.array([X[n], Y[n], Z[n]])
        location1 = np.dot(pos1, atoms.cell)
        location2 = np.dot(pos2, atoms.cell)
        #print(n, location1, location2, pos1, pos2)
        displacement = location2 - location1
        distance = np.linalg.norm(location1 - location2)
        ob = cell.copy()
        ob.data = cell.data.copy()
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
        ob.name = 'unitcell-cylinder'
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
            bpy.data.materials['unit_cell'])
        ob.location = location1 + (displacement / 2)
        ob.scale = (1, 1, 1)
        ob.dimensions = (0.1, 0.1, distance)
        phi = np.arctan2(displacement[1], displacement[0])
        #print(displacement, distance)
        theta = np.arccos(displacement[2] / distance)
        ob.rotation_euler[1] = theta
        ob.rotation_euler[2] = phi
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_cell'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.join()
    bpy.ops.object.select_all(action='DESELECT')
    return None


def create_half_bond(resolution=16):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_cylinder_add(vertices=resolution)
    if bpy.app.version[1] == 0: #use_auto_smoot dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    bond = bpy.context.object
    bond.name = 'ref_bond'
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="FACE")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    # Inset upper face:
    bond.data.polygons[14].select = True
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.inset(thickness=0.1)
    bpy.ops.mesh.select_all(action='DESELECT')
    # Inset lower face:
    bpy.ops.object.mode_set(mode='OBJECT')
    bond.data.polygons[17].select = True
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.inset(thickness=0.1)
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    return bond


def create_full_bond(resolution=16):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_cylinder_add(vertices=resolution)
    if bpy.app.version[1] == 0: #use_auto_smoot dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    bond = bpy.context.object
    bond.name = 'ref_bondx2'
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="FACE")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    # Inset upper face:
    bond.data.polygons[14].select = True
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value": (0, 0, 2)})
    bpy.ops.mesh.inset(thickness=0.1)
    bpy.ops.mesh.select_all(action='DESELECT')
    # Inset lower face:
    bpy.ops.object.mode_set(mode='OBJECT')
    bond.data.polygons[17].select = True
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.inset(thickness=0.1)
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    return bond


def is_inside_cell(pos, cell):
    # Get the atom's position
    atom_pos = pos

    # Invert the cell matrix to get the transformation matrix
    inv_cell = np.linalg.inv(cell)

    # Calculate fractional coordinates of the atom
    fractional_coords = np.dot(atom_pos, inv_cell)

    # Check if the fractional coordinates are within [0, 1] in all dimensions
    return all(0 <= coord <= 1 for coord in fractional_coords)


def assign_to_longbond(bond, mat1, mat2):
    mat1_r = bpy.data.materials.get(mat1)
    mat2_r = bpy.data.materials.get(mat2)
    bond.data.materials.append(mat1_r)
    bond.data.materials.append(mat2_r)
    list_bottom = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                   61, 62, 63, 64, 65]
    list_top = [14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                43, 44, 45, 46, 47, 48, 49]
    for face in list_bottom:
        bond.data.polygons[face].material_index = 0
    for face in list_top:
        bond.data.polygons[face].material_index = 1

