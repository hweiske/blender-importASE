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
from contextlib import contextmanager

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
                   plate_thickness=0.6, plate_holes=True,
                   plate_margin=1.0, pillar_length=2.0, segments=12,
                   bond_objects=None, pillar_material=False,
                   max_angle=math.radians(30)):
    """Resin supports: a base plate below the model and tapered pillars up
    to the underside of every atom that does not rest on another one.

    Atoms are spheres, so a bond does not hold one up: the bottom of the
    sphere still hangs free and prints as an overhang. An atom only rests
    on another when its lowest point lies inside that atom's sphere;
    every other atom gets a pillar to its lowest point, whatever its
    element and whatever it is bonded to.

    A pillar never runs through another atom sphere or bond cylinder
    (they would fuse on the print): its base is offset sideways so the
    tapered pillar routes around any obstacle up to the target atom's
    underside. An atom no pillar can reach that way is left out (and
    counted in a printed note). With plate_holes, a regular grid of
    square holes is punched into the plate to save material; pillars
    that would land on a hole are moved onto material.

    base_radius/tip_radius: pillar radius at the plate and at the atom
    contact point. max_angle: steepest tilt of a pillar from vertical
    (radians); a pillar that would have to lean further to get around
    the other atoms would itself print as an overhang, so that atom is
    left out instead.

    pillar_material: build the pillars as one object per element of the
    atom they hold up (marked 'ase_support_element', with that element's
    material), so they print in the atom's filament and blend in; the
    plate stays its own object. Returns the created objects, plate first.
    """
    centers, sphere_radii, symbols = [], [], []
    for ob in atom_objects:
        # measured on the sphere's own vertices in world space: the object
        # origin is not the center once the transforms are applied, and the
        # corners of a rotated object's bounding box span more than the
        # sphere (pillar tips would stop short of it); the plain mesh also
        # leaves out the outline shell
        co = np.empty(3 * len(ob.data.vertices))
        ob.data.vertices.foreach_get('co', co)
        mw = np.array(ob.matrix_world)
        co = co.reshape(-1, 3) @ mw[:3, :3].T + mw[:3, 3]
        lo, hi = co.min(axis=0), co.max(axis=0)
        centers.append((lo + hi) / 2)
        sphere_radii.append((hi[2] - lo[2]) / 2)
        symbols.append(ob.name.split('.')[0])
    centers = np.array(centers)
    sphere_radii = np.array(sphere_radii)
    n = len(centers)

    # the drawn bond meshes are obstacles the pillars route around
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

    # an atom rests on another when its lowest point lies inside that
    # sphere; every other atom has a free underside and needs a pillar
    lowest = centers - np.outer(sphere_radii, [0.0, 0.0, 1.0])
    inside = (np.linalg.norm(lowest[:, None, :] - centers[None, :, :], axis=2)
              < sphere_radii[None, :]) & ~np.eye(n, dtype=bool)
    pillar_atoms = [i for i in np.argsort(centers[:, 2]) if not inside[i].any()]

    zfloor = float((centers[:, 2] - sphere_radii).min())
    plate_top = zfloor - pillar_length

    def seg_point_dists(a, b, pts):
        ab = b - a
        denom = float(ab @ ab)
        if denom < 1e-9:
            return np.linalg.norm(pts - a, axis=1)
        t = np.clip((pts - a) @ ab / denom, 0.0, 1.0)
        return np.linalg.norm(pts - (a + t[:, None] * ab), axis=1)

    def max_offset(tip):
        """Furthest a pillar base may sit sideways from below its tip."""
        return (tip[2] - plate_top) * math.tan(min(max_angle, math.radians(89.9)))

    def clear_base(i, tip):
        """Find a base xy so the pillar (plate -> tip) does not skewer any
        other atom. Prefers a base that also clears every drawn bond; if
        none does, accepts one that at least clears the atoms (a support
        grazing a bond just fuses with it on the print). Returns None when
        no base avoids the atoms - better to drop the pillar than punch it
        through an atom."""
        others = np.arange(n) != i  # the target atom itself is no obstacle
        target = centers[i]
        skip = sphere_radii[i] + 0.15  # ignore samples inside the target atom

        def clearances(bx, by):
            """(atom_clear, bond_clear): signed clearances of the pillar to
            the nearest other atom and to the drawn bonds, sampled along
            it with the pillar's own (tapered) radius at each point - the
            thin tip fits past a neighbor beside the atom's underside."""
            p0 = np.array([bx, by, plate_top])
            steps = max(2, int(float(np.linalg.norm(tip - p0)) / 0.1))
            t = np.linspace(0.0, 1.0, steps + 1)
            pts = p0 + np.outer(t, tip - p0)
            radius = base_radius + (tip_radius - base_radius) * t
            # only atoms near the pillar's center line need the sampled check
            near = others & (seg_point_dists(p0, tip, centers)
                             < sphere_radii + base_radius + 0.1)
            atom_c = 1e9
            if near.any():
                gaps = (np.linalg.norm(pts[:, None, :] - centers[None, near, :], axis=2)
                        - sphere_radii[None, near] - radius[:, None])
                atom_c = float(gaps.min()) - 0.05
            bond_c = 1e9
            if bond_bvh is not None:
                for pt, rad in zip(pts[::2], radius[::2]):
                    if np.linalg.norm(pt - target) < skip:
                        continue
                    loc, _, _, d = bond_bvh.find_nearest(Vector(pt), rad + 0.3)
                    if loc is not None:
                        bond_c = min(bond_c, d - rad - 0.1)
            return atom_c, bond_c

        atom_c, bond_c = clearances(tip[0], tip[1])
        if atom_c > 0 and bond_c > 0:
            return tip[0], tip[1]
        # spiral outward, no further than max_angle allows; keep the best
        # atom-clearing base as a fallback
        best_atom = (tip[0], tip[1], bond_c) if atom_c > 0 else None
        max_off = min(5 * float(sphere_radii.max()) + 4.0,
                      max_offset(tip))
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
        pillars.append((base[0], base[1], tip[0], tip[1], tip[2], symbols[i]))
    verts, faces = [], []  # the plate (and the pillars, unless split)
    pillar_meshes = {}     # element -> (verts, faces) with pillar_material

    def add_ring(verts, x, y, z, radius):
        start = len(verts)
        for k in range(segments):
            a = 2 * math.pi * k / segments
            verts.append((x + radius * math.cos(a), y + radius * math.sin(a), z))
        return start

    def add_pillar(base_x, base_y, base_z, tip_x, tip_y, tip_z, symbol):
        pv, pf = (pillar_meshes.setdefault(symbol, ([], []))
                  if pillar_material else (verts, faces))
        bottom = add_ring(pv, base_x, base_y, base_z, base_radius)
        top = add_ring(pv, tip_x, tip_y, tip_z + 0.1, tip_radius)  # embed tip
        for k in range(segments):
            m = (k + 1) % segments
            pf.append([bottom + k, bottom + m, top + m, top + k])
        pf.append([top + k for k in range(segments)][::-1])
        pf.append([bottom + k for k in range(segments)])

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

    def push_off_holes(x, y, tip, holes_x, holes_y):
        """Move a pillar base off any hole so it lands on material, through
        the nearest hole edge that keeps the pillar within max_angle;
        None when every way off the hole tilts it too far."""
        margin = base_radius + 0.05
        for hx in holes_x:
            for hy in holes_y:
                if hx - margin < x < hx + hole + margin and hy - margin < y < hy + hole + margin:
                    moves = [(hx - margin - x, 0), (hx + hole + margin - x, 0),
                             (0, hy - margin - y), (0, hy + hole + margin - y)]
                    for dx, dy in sorted(moves, key=lambda m: abs(m[0]) + abs(m[1])):
                        if (math.hypot(x + dx - tip[0], y + dy - tip[1])
                                <= max_offset(tip) + 1e-6):
                            return x + dx, y + dy
                    return None
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
            for bx, by, tx, ty, tz, symbol in pillars:
                moved = push_off_holes(bx, by, (tx, ty, tz), holes_x, holes_y)
                if moved is None:
                    skipped += 1
                    continue
                add_pillar(*moved, plate_top, tx, ty, tz, symbol)
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
            for bx, by, tx, ty, tz, symbol in pillars:
                add_pillar(bx, by, plate_top, tx, ty, tz, symbol)
            add_box(x0, x1, y0, y1, z0, z1)

    if skipped:
        print(f'supports: {skipped} atom(s) left unsupported (no pillar within '
              f'{math.degrees(max_angle):.0f} deg of vertical avoids the other atoms)')

    # remembered so the export and the rebuild operator can start from the
    # supports as they are, and tell whether new settings need a rebuild
    params = {
        'base_radius': base_radius, 'tip_radius': tip_radius,
        'plate_thickness': plate_thickness,
        'plate_holes': plate_holes, 'pillar_length': pillar_length,
        'pillar_material': pillar_material, 'max_angle': max_angle}

    def make_object(name, verts, faces, symbol=None):
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(verts, [], faces)
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        collection.objects.link(obj)
        obj['ase_auto_supports'] = True
        obj['ase_support_params'] = params
        if symbol is not None:
            obj['ase_support_element'] = symbol
            atom = atom_objects[symbols.index(symbol)]
            if atom.active_material is not None:
                mesh.materials.append(atom.active_material)
        return obj

    objects = [make_object('supports', verts, faces)]
    for symbol in sorted(pillar_meshes):
        objects.append(make_object(f'supports_{symbol}', *pillar_meshes[symbol],
                                   symbol=symbol))
    return objects


@contextmanager
def outlines_disabled(objects):
    """Switch off the 'outline' geometry-nodes modifiers of `objects` for
    the duration of the block.

    The outline is a render effect: it joins the mesh with a copy pushed
    out along the normals, so evaluating it for print geometry would give
    every atom, bond and joint a second, slightly larger shell nested
    around the real one (and inflate the sizes the supports are fitted
    to)."""
    switched = []
    for ob in objects:
        for mod in ob.modifiers:
            if (mod.type == 'NODES' and mod.node_group is not None
                    and mod.node_group.name.startswith('outline')):
                # both flags: Blender 5.2's STL exporter evaluates the
                # render visibility, older versions the viewport one
                switched.append((mod, mod.show_viewport, mod.show_render))
                mod.show_viewport = mod.show_render = False
    bpy.context.view_layer.update()
    try:
        yield
    finally:
        for mod, viewport, render in switched:
            mod.show_viewport, mod.show_render = viewport, render
        bpy.context.view_layer.update()


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


def stored_support_params(context):
    """The build_supports parameters the active structure's automatic
    supports were made with, or None (no structure, no automatic
    supports, or supports from a version that did not record them)."""
    try:
        supports = collect_print_objects(context)[3]
    except ValueError:
        return None
    for ob in supports:
        if ob.get('ase_auto_supports') and 'ase_support_params' in ob:
            return {k: v for k, v in ob['ase_support_params'].items()}
    return None


def _same_params(a, b):
    return a is not None and set(a) == set(b) and all(
        abs(float(a[k]) - float(b[k])) < 1e-6 for k in b)


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
    with outlines_disabled(atom_objects + bonds):
        kept.extend(build_supports(atom_objects, collection,
                                   bond_objects=bonds, **params))
    return kept


_STL_DTYPE = np.dtype([('normal', '<f4', 3), ('verts', '<f4', (3, 3)),
                       ('attr', '<u2')])


def _read_stl(path):
    """Header and triangle records of a binary STL."""
    with open(path, 'rb') as f:
        data = f.read()
    n = int.from_bytes(data[80:84], 'little')
    return data[:80], np.frombuffer(data[84:84 + 50 * n], dtype=_STL_DTYPE).copy()


def _write_stl(path, header, tris):
    with open(path, 'wb') as f:
        f.write(header)
        f.write(len(tris).to_bytes(4, 'little'))
        f.write(tris.tobytes())


def place_for_print(stl_paths, scale):
    """Move and scale the STL files as one rigid body: centered on the
    origin in x/y, the lowest point (the support plate) at z=0, and
    `scale` mm per Angstrom.

    STL carries no units - slicers read the numbers as mm. Left in
    Angstrom, a molecule is a few mm across and single parts are small
    enough for a slicer to offer a unit conversion on some files but not
    others; and parts far from the origin are re-centered one by one.
    The same transform on every file keeps them aligned."""
    loaded = {p: _read_stl(p) for p in stl_paths}
    verts = [t['verts'].reshape(-1, 3) for _, t in loaded.values() if len(t)]
    if not verts:
        return
    verts = np.concatenate(verts)
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    offset = np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    for path, (header, tris) in loaded.items():
        tris['verts'] = (tris['verts'] - offset) * scale
        _write_stl(path, header, tris)


def export_3dprint(context, filepath, generate_supports=True,
                   base_radius=0.25, tip_radius=0.1,
                   plate_thickness=0.6, plate_holes=True, plate_gap=2.0,
                   print_scale=10.0, pillar_material=False,
                   max_angle=math.radians(30)):
    """Export the collection of the active object as per-element STLs,
    the bonds, and supports, zipped into `filepath`.

    All files share one frame: print_scale mm per Angstrom, centered on
    the origin with the plate on z=0 (see place_for_print), so loading
    them together as one multi-part object keeps them aligned.

    With generate_supports, the automatic supports are built with the
    support parameters below: existing automatic supports made with the
    same parameters are exported as they are (exactly what you see),
    otherwise they are rebuilt. Supports you made yourself (no
    ase_auto_supports mark) are always exported as-is, and no automatic
    ones are added next to them. Without generate_supports, whatever
    supports exist are exported unchanged.

    pillar_material: the pillars go into the STL of the element they hold
    up (atoms_<element>.stl), so the slicer prints them in that atom's
    filament; only the plate stays in supports.stl.
    """
    collection, groups, bonds, supports = collect_print_objects(context)

    # bake any unapplied object scale into the mesh data before exporting
    # (scenes imported before the draw_atoms fix carry scaled atom objects)
    for group in groups.values():
        for ob in group:
            if tuple(ob.scale) != (1.0, 1.0, 1.0):
                ob.data.transform(Matrix.Diagonal((*ob.scale, 1.0)))
                ob.scale = (1.0, 1.0, 1.0)

    params = dict(base_radius=base_radius, tip_radius=tip_radius,
                  plate_thickness=plate_thickness, plate_holes=plate_holes,
                  pillar_length=plate_gap, pillar_material=pillar_material,
                  max_angle=max_angle)
    user_made = [ob for ob in supports if not ob.get('ase_auto_supports')]
    if generate_supports and not user_made and not _same_params(
            stored_support_params(context), params):
        supports = rebuild_supports(context, **params)

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
    atom_objects = [ob for objs in groups.values() for ob in objs]
    # pillars built per element (pillar_material) print with their atoms
    pillars = {}
    for ob in supports:
        if ob.get('ase_support_element') in groups:
            pillars.setdefault(ob['ase_support_element'], []).append(ob)
    plate = [ob for ob in supports
             if not any(ob in objs for objs in pillars.values())]
    with outlines_disabled(atom_objects + bonds + supports):
        for element in sorted(groups):
            stl_files.append(write_stl(groups[element] + pillars.get(element, []),
                                       f'atoms_{element}.stl'))
        if bonds:
            stl_files.append(write_stl(bonds, 'bonds.stl'))
        if plate:
            stl_files.append(write_stl(plate, 'supports.stl'))

    place_for_print(stl_files, print_scale)

    with zipfile.ZipFile(filepath, 'w', zipfile.ZIP_DEFLATED) as archive:
        for stl in stl_files:
            archive.write(stl, os.path.basename(stl))
    for stl in stl_files:
        os.remove(stl)
    os.rmdir(tmpdir)
    return [os.path.basename(s) for s in stl_files]
