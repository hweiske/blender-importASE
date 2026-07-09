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

    A pillar never runs through another atom sphere (they would fuse on
    the print): its base is offset sideways so the tapered pillar routes
    around any blocking atom up to the target atom's underside. With
    plate_holes, a regular grid of square holes is punched into the plate
    to save material; pillars that would land on a hole are moved onto
    material.

    base_radius/tip_radius: pillar radius at the plate and at the atom
    contact point. support_layer: minimum height a bonded/touching
    neighbor must sit below an atom to hold it up - larger values treat
    shallower bonds as overhangs and add more pillars. Returns the
    created object.
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

    # bond / contact connectivity (bonds match the drawn ones; 'touching'
    # also catches non-bonded spheres one rests on)
    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=2)
    bonded = (distances < 1.2 * (bond_radii[:, None] + bond_radii[None, :]))
    touching = distances < (sphere_radii[:, None] + sphere_radii[None, :] + 0.3)
    holds = (bonded | touching) & ~np.eye(n, dtype=bool)

    # bottom-up island analysis: an atom is grounded only if it rests on a
    # grounded neighbor that is genuinely BELOW it (support_layer). A
    # neighbor at the same height or above never counts - those atoms are
    # printed before their anchor and need their own pillar.
    order = np.argsort(centers[:, 2])
    grounded = np.zeros(n, dtype=bool)
    pillar_atoms = []
    for i in order:
        below = centers[:, 2] < centers[i, 2] - support_layer
        if np.any(below & grounded & holds[i]):
            grounded[i] = True
        else:
            pillar_atoms.append(i)
            grounded[i] = True

    zfloor = float((centers[:, 2] - sphere_radii).min())
    plate_top = zfloor - pillar_length

    def seg_point_dists(a, b, pts):
        ab = b - a
        denom = float(ab @ ab)
        if denom < 1e-9:
            return np.linalg.norm(pts - a, axis=1)
        t = np.clip((pts - a) @ ab / denom, 0.0, 1.0)
        return np.linalg.norm(pts - (a + t[:, None] * ab), axis=1)

    def clear_base(i, tip):
        """A base xy whose pillar segment (plate -> tip) clears every other
        atom sphere; offset sideways (bypass) when the vertical line is
        blocked. Falls back to the least-bad offset if nothing fully clears."""
        need = sphere_radii + base_radius + 0.15
        need[i] = -1.0  # ignore the target atom itself

        def min_clearance(bx, by):
            a = np.array([bx, by, plate_top])
            return float((seg_point_dists(a, tip, centers) - need).min())

        if min_clearance(tip[0], tip[1]) > 0:
            return tip[0], tip[1]
        best, best_c = (tip[0], tip[1]), -1e18
        max_off = 3 * float(sphere_radii.max()) + 2.0
        r = 0.3
        while r <= max_off:
            for k in range(12):
                ang = 2 * math.pi * k / 12
                bx, by = tip[0] + r * math.cos(ang), tip[1] + r * math.sin(ang)
                c = min_clearance(bx, by)
                if c > 0:
                    return bx, by
                if c > best_c:
                    best, best_c = (bx, by), c
            r += 0.3
        return best

    # pillar = (base_x, base_y, tip_x, tip_y, tip_z), base routed to clear atoms
    pillars = []
    for i in pillar_atoms:
        tip = np.array([centers[i, 0], centers[i, 1], centers[i, 2] - sphere_radii[i]])
        bx, by = clear_base(i, tip)
        pillars.append((bx, by, tip[0], tip[1], tip[2]))

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
    xs = [p[0] for p in pillars] or [c[0] for c in centers]
    ys = [p[1] for p in pillars] or [c[1] for c in centers]
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
    if pillars:
        if plate_holes:
            holes_x, holes_y = hole_grid()
            for bx, by, tx, ty, tz in pillars:
                px, py = push_off_holes(bx, by, holes_x, holes_y)
                add_pillar(px, py, plate_top, tx, ty, tz)
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
            for bx, by, tx, ty, tz in pillars:
                add_pillar(bx, by, plate_top, tx, ty, tz)
            add_box(x0, x1, y0, y1, z0, z1)

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
