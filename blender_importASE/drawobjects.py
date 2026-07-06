import bpy
import bmesh
from mathutils import Vector
import ase
from ase.data import covalent_radii
from .utils import get_vdw_radius
import numpy as np
import ase.neighborlist
from .utils import atomcolors


def draw_atoms(atoms, scale=1,resolution=16, representation="Balls'n'Sticks"):
    cnt = 0
    list_of_atoms=[]
    # bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0, 0, 0), segments=resolution, ring_count=resolution)
    if bpy.app.version < (4, 1, 0): #use_auto_smooth dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    sphere.name = bpy.context.object.name
    for n, atom in enumerate(atoms):
        ob = sphere.copy()
        ob.data = sphere.data.copy()
        ob.location = atom.position
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].name = atom.symbol
        size = bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale
        if representation == "Balls'n'Sticks":
            size = [covalent_radii[atom.number]  * scale, ] * 3
        if representation == "bonds_fromnodes":
            size = [covalent_radii[atom.number]  * scale, ] * 3
        elif representation == 'Licorice':
            size = [0.1] * 3
        elif representation == 'VDW':
            # get_vdw_radius also covers Z > 103, where the Alvarez table
            # would raise an IndexError
            size = [get_vdw_radius(atom.number)] * 3
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = size
        # sprint(bpy.data.node_groups)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(
            bpy.data.materials[atom.symbol])
        cnt += 1
        #full_object_name = bpy.utils.object.full_name(ob)
        #print(full_object_name)
        #list_of_atoms.append(full_object_name)
        #print(ob)
        list_of_atoms.append(ob)
        bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects[sphere.name].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_atoms


def draw_bonds(atoms,resolution=16):
    list_of_bonds=[]
    bondlengths=[] # for animation
    nl = ase.neighborlist.NeighborList([covalent_radii[atomic_number] * 0.9 for atomic_number in atoms.numbers],
                                       self_interaction=False, bothways=True)
    nl.update(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    try:
        scene = bpy.context.scene
        bond_collection = bpy.data.collections.new("bonds")
        scene.collection.children.link(bond_collection)
    except Exception:
        print("Group bonds was not created, proceeding without group")
    # bpy.ops.surface.primitive_nurbs_surface_cylinder_add(radius=1.0, enter_editmode=False, align='WORLD',
    # location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bond = create_bond(resolution=resolution, half_bond=True)
    cnt = 0
    for atom in atoms:
        if nl.get_neighbors(atom.index)[0].size > 0:
            neighbors, offsets = nl.get_neighbors(atom.index)
            for neighbor, offset in zip(neighbors, offsets):
                make_shortbond(atom, atoms, offset, neighbor, bond, list_of_bonds, bondlengths,
                                   type='short1')
                make_shortbond(atom, atoms, offset, neighbor, bond, list_of_bonds, bondlengths,
                                   type='short2')
                cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bond.select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_bonds,nl


def draw_longbonds(atoms,resolution=16, colorbonds=False):
    cb = atomcolors()
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
    hbond = create_bond(resolution=resolution, half_bond=True)
    bond = create_bond(resolution=resolution, half_bond=False)
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
                    make_longbond(atom, atoms, offset, neighbor, bond, list_of_bonds, bondlengths, cb, colorbonds=colorbonds)
                else:
                    make_shortbond(atom, atoms, offset, neighbor, hbond, list_of_bonds, bondlengths,
                                   type='short1', colorbonds=colorbonds)
                    make_shortbond(atom, atoms, offset, neighbor, hbond, list_of_bonds, bondlengths,
                                   type='short2', colorbonds=colorbonds)
            cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_bond'].select_set(True)
    bpy.data.objects['ref_bondx2'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return list_of_bonds,nl,bondlengths

def make_longbond(atom, atoms, offset, neighbor, bond, list_of_bonds, bondlengths, cb, colorbonds=False):
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
    create_bond_mat(cb, ob, atom.symbol, atoms[neighbor].symbol, smooth=True, colorbond=colorbonds)
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
    bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)

def make_shortbond(atom, atoms, offset, neighbor, bond, list_of_bonds, bondlengths, type='short1', colorbonds=False):
    if type == "short1":
        displacements = [0.5 * (atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)),
                                        0.5 * (atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor])]
        location = atom.position + (displacements[0] / 2)
    else:
        displacements = [0.5 * (atom.position - atoms.positions[neighbor] - np.dot(offset, atoms.cell)),
                                     0.5 * (atoms.positions[neighbor] + np.dot(offset, atoms.cell) + atom.position)]
        location = atoms.positions[neighbor] + (displacements[0] / 2)
    distance = atoms.get_distance(atom.index, neighbor, mic=True) / 2
    ob = bond.copy()
    ob.data = bond.data.copy()
    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}'
    if colorbonds:
        ob.data.materials.append(bpy.data.materials[f'{atom.symbol}-bond'])
    else:
        ob.data.materials.append(bpy.data.materials['Gray-bond'])
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
    if bpy.app.version < (4, 1, 0): #use_auto_smooth dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    cell = bpy.context.object
    cell.name = 'ref_cell'
    X = [0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    Y = [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    Z = [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]
    cylinders = []
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
        cylinders.append(ob)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_cell'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    # join only the cell cylinders; select_all('SELECT') would also join the
    # atoms and bonds, and join() needs an active object in background mode
    for ob in cylinders:
        ob.select_set(True)
    bpy.context.view_layer.objects.active = cylinders[-1]
    bpy.ops.object.join()
    bpy.ops.object.select_all(action='DESELECT')
    return None

def create_bond(resolution=16, half_bond=False):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_cylinder_add(vertices=resolution)
    if bpy.app.version < (4, 1, 0): #use_auto_smooth dropped after 4.0
        bpy.ops.object.shade_smooth(use_auto_smooth=True)
    else:
        bpy.ops.object.shade_smooth()
    bond = bpy.context.object
    if half_bond:
        bond.name = 'ref_bond'
    else:
        bond.name = 'ref_bondx2'
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="FACE")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type='FACE')
    bm = bmesh.from_edit_mesh(bond.data)
    z = Vector((0.0, 0.0, 1.0))
    tol = 1e-3
    for f in bm.faces:
        n = f.normal.normalized() if f.normal.length > tol else f.normal
        dot = n.dot(z)
        # Select faces whose normals are nearly +Z or -Z
        if abs(dot - 1.0) <= tol:
            f.select = True
        else:
            f.select = False
    bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value": (0, 0, 2)})
    if not half_bond:
        bpy.ops.transform.translate(value=(0, 0, -0.3))
        bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value": (0, 0, 0.3)})
        bpy.ops.transform.resize(value=(1.3, 1.3, 1))
    bpy.ops.mesh.inset(thickness=0.1)
    for f in bm.faces:
        n = f.normal.normalized() if f.normal.length > tol else f.normal
        dot = n.dot(z)
        # Select faces whose normals are nearly +Z or -Z
        if abs(dot + 1.0) <= tol:
            f.select = True
        else:
            f.select = False
    bpy.ops.transform.translate(value=(0, 0, 0.3))
    bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value": (0, 0, -0.3)})
    bpy.ops.transform.resize(value=(1.3, 1.3, 1))
    bpy.ops.mesh.inset(thickness=0.1)
    bmesh.update_edit_mesh(bond.data, loop_triangles=False, destructive=False)
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


def create_bond_mat(cb, bond, atom_1, atom_2, smooth, colorbond=False):
    if colorbond:
        mat = cb.create_bondmat(atom_1, atom_2, smooth, name=f'{atom_1}-{atom_2}-bond')
        bond.data.materials.append(mat)
    else:
        mat = bpy.data.materials['Gray-bond']
        bond.data.materials.append(mat)
