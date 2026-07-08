"""Export features.

- export_xyz: write the active nodes-representation object back to a real
  .xyz file (vertex positions + the stored 'element' attribute).
- export_3dprint: export a 3D_print-representation collection as STL
  files - the atoms joined per element, the bonds (with their icosphere
  joints) as one file, and simple resin supports - zipped into a single
  archive for the slicer.
"""
import os
import tempfile
import zipfile
import math

import bpy
import numpy as np
from mathutils import Vector
from ase import Atoms
from ase.data import chemical_symbols
import ase.io


def export_xyz(obj, filepath):
    """Write vertices + 'element' attribute of a nodes-representation
    object to xyz, in world coordinates."""
    mesh = obj.data
    if 'element' not in mesh.attributes:
        raise ValueError(
            f"'{obj.name}' has no 'element' attribute - select the structure "
            "object of a nodes-representation import")
    n = len(mesh.vertices)
    numbers = np.empty(n)
    mesh.attributes['element'].data.foreach_get('value', numbers)
    co = np.empty(n * 3)
    mesh.vertices.foreach_get('co', co)
    co = co.reshape(-1, 3)
    world = np.array(obj.matrix_world.to_3x3())
    translation = np.array(obj.matrix_world.translation)
    positions = co @ world.T + translation
    atoms = Atoms(numbers=numbers.round().astype(int), positions=positions)
    ase.io.write(filepath, atoms, format='xyz')
    return len(atoms)


def build_supports(atom_objects, collection, base_radius=0.25, tip_radius=0.1,
                   support_layer=0.8, plate_thickness=0.6, plate_holes=True,
                   plate_margin=1.0, pillar_length=2.0, segments=12):
    """Resin supports: a base plate below the model and tapered pillars up
    to every atom that would print as an island.

    Working bottom-up, an atom counts as supported when it is bonded
    (covalent radii * 1.2) to an already-supported atom at most
    'support_layer' higher - everything else gets a pillar, so floating
    fragments and upward-only branches are supported too.

    A pillar whose vertical path would run through another atom sphere
    starts on top of the highest blocking sphere instead of the plate
    (support stacking). With plate_holes, a regular grid of square holes
    is punched into the plate to save material; pillars that would land
    on a hole are moved onto material.

    base_radius/tip_radius: pillar radius at the plate and at the atom
    contact point. Returns the created object.
    """
    from ase.data import covalent_radii, atomic_numbers

    centers, sphere_radii, bond_radii = [], [], []
    for ob in atom_objects:
        corners = [ob.matrix_world @ Vector(c) for c in ob.bound_box]
        origin = np.array(ob.matrix_world.translation)
        centers.append(origin)
        zs = [v.z for v in corners]
        sphere_radii.append((max(zs) - min(zs)) / 2)
        symbol = ob.name.split('.')[0]
        bond_radii.append(covalent_radii[atomic_numbers.get(symbol, 6)])
    centers = np.array(centers)
    sphere_radii = np.array(sphere_radii)
    bond_radii = np.array(bond_radii)
    n = len(centers)

    # bond connectivity from distances (matches the drawn bonds closely)
    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=2)
    cutoff = 1.2 * (bond_radii[:, None] + bond_radii[None, :])
    bonded = (distances < cutoff) & ~np.eye(n, dtype=bool)

    # bottom-up island analysis: pillar every atom without a supported,
    # not-too-much-higher bonded neighbor
    order = np.argsort(centers[:, 2])
    is_supported = np.zeros(n, dtype=bool)
    pillar_atoms = []
    for i in order:
        holds = [j for j in np.nonzero(bonded[i])[0]
                 if is_supported[j] and centers[j, 2] <= centers[i, 2] + support_layer]
        if not holds:
            pillar_atoms.append(i)
        is_supported[i] = True

    zfloor = float((centers[:, 2] - sphere_radii).min())
    plate_top = zfloor - pillar_length

    # pillar start points: the plate, or the top of the highest atom sphere
    # blocking the vertical path (support stacking)
    plate_pillars, stacked_pillars = [], []
    for i in pillar_atoms:
        x, y = centers[i, 0], centers[i, 1]
        target_z = centers[i, 2] - sphere_radii[i]
        start_z = None
        for j in range(n):
            if j == i:
                continue
            d_xy = math.hypot(centers[j, 0] - x, centers[j, 1] - y)
            if d_xy >= sphere_radii[j]:
                continue
            top_at_column = centers[j, 2] + math.sqrt(sphere_radii[j]**2 - d_xy**2)
            if top_at_column < target_z and (start_z is None or top_at_column > start_z):
                start_z = top_at_column
        if start_z is None:
            plate_pillars.append((x, y, target_z))
        else:
            stacked_pillars.append((x, y, target_z, start_z - 0.1))  # embed base

    verts, faces = [], []

    def add_ring(x, y, z, radius):
        start = len(verts)
        for k in range(segments):
            a = 2 * math.pi * k / segments
            verts.append((x + radius * math.cos(a), y + radius * math.sin(a), z))
        return start

    def add_pillar(base_x, base_y, base_z, tip_x, tip_y, tip_z):
        bottom = add_ring(base_x, base_y, base_z, base_radius)
        top = add_ring(tip_x, tip_y, tip_z + 0.1, tip_radius)  # embed tip
        for k in range(segments):
            m = (k + 1) % segments
            faces.append([bottom + k, bottom + m, top + m, top + k])
        faces.append([top + k for k in range(segments)][::-1])
        faces.append([bottom + k for k in range(segments)])

    # plate layout: a solid slab, optionally with a regular grid of square
    # holes punched in to save material
    xs = [p[0] for p in plate_pillars] or [c[0] for c in centers]
    ys = [p[1] for p in plate_pillars] or [c[1] for c in centers]
    x0, x1 = min(xs) - plate_margin, max(xs) + plate_margin
    y0, y1 = min(ys) - plate_margin, max(ys) + plate_margin
    hole = max(0.6, 2 * base_radius)   # hole edge length
    pitch = 2.2 * hole                 # hole spacing (majority stays material)

    def hole_grid():
        """Regularly spaced square holes, kept clear of the plate border."""
        holes_x, holes_y = [], []
        hx = x0 + pitch - hole / 2
        while hx + hole < x1 - (pitch - hole):
            holes_x.append(hx)
            hx += pitch
        hy = y0 + pitch - hole / 2
        while hy + hole < y1 - (pitch - hole):
            holes_y.append(hy)
            hy += pitch
        return holes_x, holes_y

    def push_off_holes(x, y, holes_x, holes_y):
        """Move a pillar base off any hole so it lands on material."""
        margin = base_radius + 0.05
        for hx in holes_x:
            for hy in holes_y:
                if hx - margin < x < hx + hole + margin and hy - margin < y < hy + hole + margin:
                    # push out through the nearest hole edge
                    moves = [(hx - margin - x, 0), (hx + hole + margin - x, 0),
                             (0, hy - margin - y), (0, hy + hole + margin - y)]
                    dx, dy = min(moves, key=lambda m: abs(m[0]) + abs(m[1]))
                    return x + dx, y + dy
        return x, y

    def add_box(bx0, bx1, by0, by1, bz0, bz1):
        if bx1 <= bx0 or by1 <= by0:
            return
        base = len(verts)
        verts.extend([(bx0, by0, bz0), (bx1, by0, bz0), (bx1, by1, bz0), (bx0, by1, bz0),
                      (bx0, by0, bz1), (bx1, by0, bz1), (bx1, by1, bz1), (bx0, by1, bz1)])
        faces.extend([[base, base + 1, base + 2, base + 3][::-1],
                      [base + 4, base + 5, base + 6, base + 7],
                      [base, base + 1, base + 5, base + 4],
                      [base + 1, base + 2, base + 6, base + 5],
                      [base + 2, base + 3, base + 7, base + 6],
                      [base + 3, base, base + 4, base + 7]])

    z0, z1 = plate_top - plate_thickness, plate_top
    if plate_pillars:
        if plate_holes:
            holes_x, holes_y = hole_grid()
            for x, y, target_z in plate_pillars:
                bx, by = push_off_holes(x, y, holes_x, holes_y)
                add_pillar(bx, by, plate_top, x, y, target_z)
            # full-width strips between the hole rows...
            y_edges = [y0] + [e for hy in holes_y for e in (hy, hy + hole)] + [y1]
            for yA, yB in zip(y_edges[0::2], y_edges[1::2]):
                add_box(x0, x1, yA, yB, z0, z1)
            # ...and the pieces between holes within each hole row
            x_edges = [x0] + [e for hx in holes_x for e in (hx, hx + hole)] + [x1]
            for hy in holes_y:
                for xA, xB in zip(x_edges[0::2], x_edges[1::2]):
                    add_box(xA, xB, hy, hy + hole, z0, z1)
        else:
            for x, y, target_z in plate_pillars:
                add_pillar(x, y, plate_top, x, y, target_z)
            add_box(x0, x1, y0, y1, z0, z1)

    for x, y, target_z, start_z in stacked_pillars:
        add_pillar(x, y, start_z, x, y, target_z)

    mesh = bpy.data.meshes.new('supports')
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    obj = bpy.data.objects.new('supports', mesh)
    collection.objects.link(obj)
    obj['ase_auto_supports'] = True
    return obj


def collect_print_objects(context):
    """The atom groups, bonds, and support objects of the active object's
    import collection. Returns (collection, groups, bonds, supports)."""
    active = context.active_object
    if active is None or not active.users_collection:
        raise ValueError('select an object of the imported structure first')
    collection = active.users_collection[0]
    # operate on the top-level import collection if the active object sits
    # in the 'atoms' sub-collection
    for coll in bpy.data.collections:
        if collection.name in [c.name for c in coll.children]:
            collection = coll
            break

    groups = {}
    bonds, supports = [], []
    for ob in collection.all_objects:
        if ob.type != 'MESH':
            continue
        base = ob.name.split('.')[0]
        if base == 'bonds_object':
            bonds.append(ob)
        elif 'support' in ob.name.lower():
            supports.append(ob)
        elif base in chemical_symbols:
            groups.setdefault(base, []).append(ob)
    if not groups:
        raise ValueError(
            f"no atom objects found in collection '{collection.name}' - "
            "use the 3D print representation")
    return collection, groups, bonds, supports


def rebuild_supports(context, **params):
    """(Re)generate the automatic supports of the active structure with the
    given build_supports parameters. Auto-generated supports (marked with
    the ase_auto_supports property) are replaced; user-made support
    objects are left alone."""
    collection, groups, bonds, supports = collect_print_objects(context)
    kept = []
    for ob in supports:
        if ob.get('ase_auto_supports'):
            bpy.data.objects.remove(ob, do_unlink=True)
        else:
            kept.append(ob)
    atom_objects = [ob for objs in groups.values() for ob in objs]
    kept.append(build_supports(atom_objects, collection, **params))
    return kept


def export_3dprint(context, filepath, generate_supports=True,
                   base_radius=0.25, tip_radius=0.1, support_layer=0.8,
                   plate_thickness=0.6, plate_holes=True):
    """Export the collection of the active object as per-element STLs,
    the bonds, and supports, zipped into `filepath`."""
    collection, groups, bonds, supports = collect_print_objects(context)

    user_supports = [ob for ob in supports if not ob.get('ase_auto_supports')]
    if generate_supports and not user_supports:
        # regenerate the auto supports with the current parameters
        supports = rebuild_supports(context, base_radius=base_radius,
                                    tip_radius=tip_radius,
                                    support_layer=support_layer,
                                    plate_thickness=plate_thickness,
                                    plate_holes=plate_holes)
    elif user_supports:
        supports = user_supports

    tmpdir = tempfile.mkdtemp(prefix='ase_3dprint_')

    def write_stl(objects, name):
        bpy.ops.object.select_all(action='DESELECT')
        for ob in objects:
            ob.select_set(True)
        context.view_layer.objects.active = objects[0]
        stl_path = os.path.join(tmpdir, name)
        bpy.ops.wm.stl_export(filepath=stl_path, export_selected_objects=True,
                              apply_modifiers=True)
        return stl_path

    stl_files = []
    for element in sorted(groups):
        stl_files.append(write_stl(groups[element], f'atoms_{element}.stl'))
    if bonds:
        stl_files.append(write_stl(bonds, 'bonds.stl'))
    if supports:
        stl_files.append(write_stl(supports, 'supports.stl'))

    with zipfile.ZipFile(filepath, 'w', zipfile.ZIP_DEFLATED) as archive:
        for stl in stl_files:
            archive.write(stl, os.path.basename(stl))
    for stl in stl_files:
        os.remove(stl)
    os.rmdir(tmpdir)
    return [os.path.basename(s) for s in stl_files]
