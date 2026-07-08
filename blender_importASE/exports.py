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


def build_supports(atom_objects, collection, pillar_radius=0.25,
                   support_layer=0.8, plate_thickness=0.6, plate_margin=1.0,
                   tip_fraction=0.4, pillar_length=2.0, segments=12):
    """Resin supports: a base plate below the model and tapered pillars up
    to every atom that would print as an island. Working bottom-up, an
    atom counts as supported when it is bonded (covalent radii * 1.2) to
    an already-supported atom that is at most 'support_layer' higher -
    everything else gets a pillar, so floating fragments and upward-only
    branches are supported too. Returns the created object."""
    from ase.data import covalent_radii, atomic_numbers

    centers, bottoms, radii = [], [], []
    for ob in atom_objects:
        corners = [ob.matrix_world @ Vector(c) for c in ob.bound_box]
        origin = np.array(ob.matrix_world.translation)
        centers.append(origin)
        bottoms.append((origin[0], origin[1], min(v.z for v in corners)))
        symbol = ob.name.split('.')[0]
        radii.append(covalent_radii[atomic_numbers.get(symbol, 6)])
    centers = np.array(centers)
    radii = np.array(radii)
    n = len(centers)

    # bond connectivity from distances (matches the drawn bonds closely)
    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=2)
    cutoff = 1.2 * (radii[:, None] + radii[None, :])
    bonded = (distances < cutoff) & ~np.eye(n, dtype=bool)

    # bottom-up island analysis: pillar every atom without a supported,
    # not-too-much-higher bonded neighbor
    order = np.argsort(centers[:, 2])
    is_supported = np.zeros(n, dtype=bool)
    supported = []
    for i in order:
        holds = [j for j in np.nonzero(bonded[i])[0]
                 if is_supported[j] and centers[j, 2] <= centers[i, 2] + support_layer]
        if not holds:
            supported.append(bottoms[i])
        is_supported[i] = True
    zfloor = min(b[2] for b in bottoms)

    plate_top = zfloor - pillar_length
    verts, faces = [], []

    def add_ring(x, y, z, radius):
        start = len(verts)
        for i in range(segments):
            a = 2 * math.pi * i / segments
            verts.append((x + radius * math.cos(a), y + radius * math.sin(a), z))
        return start

    for x, y, z in supported:
        bottom = add_ring(x, y, plate_top, pillar_radius)
        top = add_ring(x, y, z + 0.1, pillar_radius * tip_fraction)  # embed tip
        for i in range(segments):
            j = (i + 1) % segments
            faces.append([bottom + i, bottom + j, top + j, top + i])
        # cap the tip
        faces.append([top + i for i in range(segments)][::-1])

    # base plate under all pillars
    xs = [b[0] for b in supported]
    ys = [b[1] for b in supported]
    x0, x1 = min(xs) - plate_margin, max(xs) + plate_margin
    y0, y1 = min(ys) - plate_margin, max(ys) + plate_margin
    z0, z1 = plate_top - plate_thickness, plate_top
    base = len(verts)
    verts.extend([(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                  (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)])
    faces.extend([[base, base + 1, base + 2, base + 3][::-1],
                  [base + 4, base + 5, base + 6, base + 7],
                  [base, base + 1, base + 5, base + 4],
                  [base + 1, base + 2, base + 6, base + 5],
                  [base + 2, base + 3, base + 7, base + 6],
                  [base + 3, base, base + 4, base + 7]])

    mesh = bpy.data.meshes.new('supports')
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    obj = bpy.data.objects.new('supports', mesh)
    collection.objects.link(obj)
    return obj


def export_3dprint(context, filepath, generate_supports=True,
                   pillar_radius=0.25, support_layer=0.8):
    """Export the collection of the active object as per-element STLs,
    the bonds, and supports, zipped into `filepath`."""
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

    if generate_supports and not supports:
        atom_objects = [ob for objs in groups.values() for ob in objs]
        supports = [build_supports(atom_objects, collection,
                                   pillar_radius=pillar_radius,
                                   support_layer=support_layer)]

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
