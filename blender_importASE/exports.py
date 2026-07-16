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
from mathutils import Matrix, Vector
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
                   plate_margin=1.0, pillar_length=2.0, segments=12,
                   bond_radius=0.15, bond_objects=None):
    """Resin supports: a base plate below the model and tapered pillars up
    to every atom that would print as an island.

    Working bottom-up, an atom counts as supported when it is bonded
    (covalent radii * 1.2) to an already-supported atom at most
    'support_layer' higher - everything else gets a pillar, so floating
    fragments and upward-only branches are supported too.

    A pillar never runs through another atom sphere or bond cylinder
    (they would fuse on the print): its base is offset sideways so the
    tapered pillar routes around any obstacle up to the target atom's
    underside. With plate_holes, a regular grid of square holes is
    punched into the plate to save material; pillars that would land on
    a hole are moved onto material.

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

    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=2)

    # A BVH of the actually drawn bond meshes is the ground truth for both
    # connectivity and obstacle avoidance - no covalent-radius guess can
    # match the node tree's bonds across ionic (no Na-Na) and cluster
    # (B-B) cases at once.
    from mathutils.bvhtree import BVHTree
    bond_bvh = None
    if bond_objects:
        deps = bpy.context.evaluated_depsgraph_get()
        bverts, bpolys = [], []
        for bo in bond_objects:
            ev = bo.evaluated_get(deps)
            me = ev.to_mesh()
            mw = bo.matrix_world
            off = len(bverts)
            bverts.extend([mw @ v.co for v in me.vertices])
            bpolys.extend([[off + vi for vi in p.vertices] for p in me.polygons])
            ev.to_mesh_clear()
        if bpolys:
            bond_bvh = BVHTree.FromPolygons(bverts, bpolys)

    def bond_between(a, b):
        """True if a drawn bond tube runs between atoms a and b (sampled
        points along the center line lie on/inside the bond geometry)."""
        pa, pb = Vector(centers[a]), Vector(centers[b])
        hits = 0
        for t in (0.35, 0.5, 0.65):
            loc, _, _, d = bond_bvh.find_nearest(pa.lerp(pb, t), bond_radius + 0.4)
            if loc is not None and d < bond_radius + 0.25:
                hits += 1
        return hits >= 2

    if bond_bvh is not None:
        bonded = np.zeros((n, n), dtype=bool)
        for a in range(n):
            for b in range(a + 1, n):
                if distances[a, b] < 3.6 and bond_between(a, b):
                    bonded[a, b] = bonded[b, a] = True
    else:  # fallback when no bond geometry is available
        bonded = distances < 1.2 * (bond_radii[:, None] + bond_radii[None, :])
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
        """Find a base xy so the pillar (plate -> tip) does not skewer any
        other atom. Prefers a base that also clears every drawn bond; if
        none does, accepts one that at least clears the atoms (a support
        grazing a bond just fuses with it on the print). Returns None when
        no base avoids the atoms - better to drop the pillar than punch it
        through an atom."""
        need = sphere_radii + base_radius + 0.15
        need[i] = -1.0  # ignore the target atom itself
        target = centers[i]
        skip = sphere_radii[i] + 0.15  # ignore samples inside the target atom

        def clearances(bx, by):
            """(atom_clear, bond_clear): signed clearances of the pillar
            segment to the nearest other atom and to the drawn bonds."""
            p0 = np.array([bx, by, plate_top])
            atom_c = float((seg_point_dists(p0, tip, centers) - need).min())
            bond_c = 1e9
            if bond_bvh is not None:
                steps = max(2, int(float(np.linalg.norm(tip - p0)) / 0.15))
                for s in range(steps + 1):
                    pt = p0 + (tip - p0) * (s / steps)
                    if np.linalg.norm(pt - target) < skip:
                        continue
                    loc, _, _, d = bond_bvh.find_nearest(Vector(pt), base_radius + 0.3)
                    if loc is not None:
                        bond_c = min(bond_c, d - base_radius - 0.1)
            return atom_c, bond_c

        atom_c, bond_c = clearances(tip[0], tip[1])
        if atom_c > 0 and bond_c > 0:
            return tip[0], tip[1]
        # spiral outward; keep the best atom-clearing base as a fallback
        best_atom = (tip[0], tip[1], bond_c) if atom_c > 0 else None
        max_off = 5 * float(sphere_radii.max()) + 4.0
        r = 0.25
        while r <= max_off:
            for k in range(24):
                ang = 2 * math.pi * k / 24
                bx, by = tip[0] + r * math.cos(ang), tip[1] + r * math.sin(ang)
                atom_c, bond_c = clearances(bx, by)
                if atom_c > 0 and bond_c > 0:
                    return bx, by
                if atom_c > 0 and (best_atom is None or bond_c > best_atom[2]):
                    best_atom = (bx, by, bond_c)
            r += 0.25
        return best_atom[:2] if best_atom is not None else None

    # pillar = (base_x, base_y, tip_x, tip_y, tip_z); skip an atom when no
    # pillar path avoids piercing another atom (rather than skewer it)
    pillars = []
    skipped = 0
    for i in pillar_atoms:
        tip = np.array([centers[i, 0], centers[i, 1], centers[i, 2] - sphere_radii[i]])
        base = clear_base(i, tip)
        if base is None:
            skipped += 1
            continue
        pillars.append((base[0], base[1], tip[0], tip[1], tip[2]))
    if skipped:
        print(f'supports: {skipped} atom(s) left unsupported (no pillar path '
              'avoids other atoms)')

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
    kept.append(build_supports(atom_objects, collection, bond_objects=bonds, **params))
    return kept


def export_3dprint(context, filepath, generate_supports=True,
                   base_radius=0.25, tip_radius=0.1, support_layer=0.8,
                   plate_thickness=0.6, plate_holes=True, plate_gap=2.0):
    """Export the collection of the active object as per-element STLs,
    the bonds, and supports, zipped into `filepath`.

    Existing supports in the scene (whatever "Rebuild 3D-print supports"
    last produced, or any you made yourself) are exported as-is, so
    re-exporting after a rebuild picks up exactly what you see. The
    support parameters below are only used to generate supports the
    first time, when none exist yet.
    """
    collection, groups, bonds, supports = collect_print_objects(context)

    # bake any unapplied object scale into the mesh data before exporting
    # (scenes imported before the draw_atoms fix carry scaled atom objects)
    for group in groups.values():
        for ob in group:
            if tuple(ob.scale) != (1.0, 1.0, 1.0):
                ob.data.transform(Matrix.Diagonal((*ob.scale, 1.0)))
                ob.scale = (1.0, 1.0, 1.0)

    if not supports and generate_supports:
        supports = rebuild_supports(context, base_radius=base_radius,
                                    tip_radius=tip_radius,
                                    support_layer=support_layer,
                                    plate_thickness=plate_thickness,
                                    plate_holes=plate_holes,
                                    pillar_length=plate_gap)

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
