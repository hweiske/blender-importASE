import bpy
import numpy as np
from ase import Atoms
from ase.io.cube import read_cube
from ase.calculators.vasp import VaspChargeDensity
from ase.units import Bohr as BOHR_TO_ANG
#if bpy.app.version[1] < 4:         didn't work
#    import pyopenvdb as vdb
#else:
#    import openvdb as vdb
try:
    import openvdb as vdb
except ImportError:
    try:
        import pyopenvdb as vdb
    except ImportError:
        vdb = None
import os
from .node_networks.compat import set_mod_input, get_mod_input
from .node_networks.electron_density_nodes import (visualize_edensity_node_group,
                                                   density_materials, newShader)
from .utils import toggle
import os.path


def is_vasp_density(filename):
    """True for VASP volumetric files (CHGCAR, CHG, PARCHG, AECCAR*)."""
    base = os.path.basename(filename).upper()
    return any(key in base for key in ('CHGCAR', 'PARCHG', 'AECCAR')) or base.startswith('CHG')


def read_vasp_density(filepath):
    """Read a VASP charge density file once (atoms and grid).

    Robust against real-world files: transparently decompresses gzipped
    files, retries with a sanitized species line when the header carries
    POTCAR-style slashes (e.g. 'Fe/', which some ase versions cannot
    parse), and turns ase's silent empty result into a clear error.
    """
    import tempfile

    tmp_path = None
    path = filepath
    with open(filepath, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    if gzipped:
        import gzip
        with tempfile.NamedTemporaryFile(delete=False, suffix='_CHGCAR',
                                         mode='wb') as tmp:
            with gzip.open(filepath, 'rb') as src:
                tmp.write(src.read())
            tmp_path = tmp.name
        path = tmp_path

    lines = []
    try:
        density = VaspChargeDensity(path)

        if not density.atoms:
            # retry with '/' stripped from the comment and species lines
            with open(path, errors='replace') as fh:
                lines = fh.readlines()
            if len(lines) > 5 and ('/' in lines[0] or '/' in lines[5]):
                lines[0] = lines[0].replace('/', ' ')
                lines[5] = lines[5].replace('/', ' ')
                with tempfile.NamedTemporaryFile(delete=False, suffix='_CHGCAR',
                                                 mode='w') as tmp:
                    tmp.writelines(lines)
                    retry_path = tmp.name
                try:
                    density = VaspChargeDensity(retry_path)
                finally:
                    os.remove(retry_path)

        if not density.atoms:
            header = [line.rstrip() for line in lines[:8]]
            raise ValueError(
                f"could not read a structure from '{os.path.basename(filepath)}': "
                "ase could not parse the POSCAR-style header (for VASP-4 style "
                "files the first line must list the element symbols). "
                "Header read:\n  " + "\n  ".join(header))
    finally:
        if tmp_path is not None:
            os.remove(tmp_path)
    return density


def data2vol(volume, spacing, origin, filepath, modifier='GeometryNodes',
             plus_material=None, minus_material=None):
    """Turn a volumetric numpy grid into a Blender volume object with the
    visualize_edensity node group attached.

    volume:  3d array indexed along the three cell vectors
    spacing: three vectors, one grid step along each cell vector
    origin:  cartesian origin of the grid
    plus_material/minus_material: materials for the +/- isosurfaces;
    defaults to the blue/red '+ material'/'- material' pair
    """
    if vdb is None:
        raise ImportError(
            "Neither 'openvdb' nor 'pyopenvdb' is installed. "
            "Please install one of them to import volumetric density files."
        )
    GRID = vdb.FloatGrid()
    GRID.copyFromArray(np.ascontiguousarray(volume, dtype=float))
    SX = list(spacing[0]) + [0.]
    SY = list(spacing[1]) + [0.]
    SZ = list(spacing[2]) + [0.]
    GRID.transform = vdb.createLinearTransform(
        [[SX[0], SX[1], SX[2], SX[3]], [SY[0], SY[1], SY[2], SY[3]], [SZ[0], SZ[1], SZ[2], SZ[3]], [0, 0, 0, 1]])

    GRID.gridClass = vdb.GridClass.FOG_VOLUME
    GRID.name = 'density'
    TMPFILE = os.path.splitext(filepath)[0] + '_density.vdb'
    vdb.write(TMPFILE, GRID)
    _ = bpy.ops.object.volume_import(filepath=TMPFILE, location=origin)
    density_obj = bpy.context.active_object
    node = visualize_edensity_node_group()
    bpy.ops.object.modifier_add(type='NODES')
    bpy.context.object.modifiers[modifier].node_group = node
    # the isosurfaces take their material from the object's own slots
    # (slot 0 = positive, slot 1 = negative), not from a modifier input
    density_materials(density_obj, plus_material, minus_material)
    # the lattice vectors the 'offset a/b/c' sockets step along: one grid
    # spacing times the number of samples along that axis
    mod = bpy.context.object.modifiers[modifier]
    identifiers = {item.name: item.identifier for item in node.interface.items_tree
                   if getattr(item, 'in_out', None) == 'INPUT'}
    for axis, step, count in zip('abc', spacing, np.shape(volume)):
        if f'cell {axis}' in identifiers:
            set_mod_input(mod, identifiers[f'cell {axis}'],
                          [float(v) * int(count) for v in step])
    # what a later supercell rebuild tiles from: the single-cell grid and
    # its true shape (the VDB's own active-voxel box can be smaller when
    # the density is zero at the border)
    density_obj['ase_base_vdb'] = TMPFILE
    density_obj['ase_grid_shape'] = list(np.shape(volume))
    toggle(bpy.context.object, SET=False)
    return density_obj


def density_cell_vectors(density_obj):
    """The lattice vectors of a density's own grid, from the VDB transform.

    One index step along each axis is that lattice vector divided by the
    number of samples, and the transform carries it for a triclinic cell
    too - which is why this asks the grid rather than the bounding box.
    Returns None when the grid cannot be read.
    """
    if vdb is None:
        return None
    path = density_obj.get('ase_base_vdb') or density_obj.data.filepath
    path = bpy.path.abspath(path) if path else ''
    if not path or not os.path.exists(path):
        return None
    grid = vdb.read(path, 'density')
    shape = density_obj.get('ase_grid_shape')
    if shape is None:
        box = grid.evalActiveVoxelBoundingBox()
        shape = tuple(np.asarray(box[1]) - np.asarray(box[0]) + 1)
    origin = np.asarray(grid.transform.indexToWorld((0, 0, 0)))
    vectors = []
    for axis, count in enumerate(shape):
        step = np.zeros(3)
        step[axis] = 1
        vectors.append((np.asarray(grid.transform.indexToWorld(tuple(step))) - origin)
                       * int(count))
    return vectors


def upgrade_density_nodes(density_obj):
    """Point an existing density at the current visualize_edensity group.

    A modifier keeps whatever node group it was created with, so a file
    saved by an older add-on never gains what a new revision adds - the
    offsets, or the material indices - until its densities are rebuilt.
    This carries the settings over by socket name, fills in the lattice
    vectors from the grid itself and leaves the material slots as they are.
    Returns True when something was changed.
    """
    node = visualize_edensity_node_group()
    changed = False
    for mod in density_obj.modifiers:
        if (mod.type != 'NODES' or mod.node_group is None
                or not mod.node_group.name.startswith('visualize_edensity')):
            continue
        if mod.node_group is not node:
            old = {item.name: get_mod_input(mod, item.identifier)
                   for item in mod.node_group.interface.items_tree
                   if getattr(item, 'in_out', None) == 'INPUT'
                   and item.socket_type not in ('NodeSocketGeometry',
                                                'NodeSocketMaterial')}
            mod.node_group = node
            for item in node.interface.items_tree:
                if (getattr(item, 'in_out', None) == 'INPUT'
                        and item.name in old and old[item.name] is not None):
                    set_mod_input(mod, item.identifier, old[item.name])
            changed = True
        identifiers = {item.name: item.identifier for item in node.interface.items_tree
                       if getattr(item, 'in_out', None) == 'INPUT'}
        cell = None
        for axis, name in enumerate('abc'):
            key = identifiers.get(f'cell {name}')
            if key is None:
                continue
            current = get_mod_input(mod, key)
            if current is not None and any(abs(float(v)) > 1e-9 for v in current):
                continue
            if cell is None:
                cell = density_cell_vectors(density_obj)
                if cell is None:
                    break
            set_mod_input(mod, key, [float(v) for v in cell[axis]])
            changed = True
    if changed:
        density_materials(density_obj, keep_existing=True)
    return changed


def density_supercell(density_obj, repeat=(1, 1, 1)):
    """Rebuild a density volume as a supercell by tiling its grid.

    The grid of a periodic calculation is itself periodic - it holds one
    cell's worth of samples and does not repeat the far plane - so tiling
    it *is* the density of the supercell. That is what repeating the
    isosurface mesh of a single cell cannot do: each copy would still be
    capped at the cell boundary it was generated in, leaving a flat cut
    face wherever a lobe crosses it. Here the marching cubes runs over the
    whole tiled grid instead, so the surface continues through the interior
    boundaries and only the outside of the supercell is closed off.

    Always tiles from the single-cell grid the import wrote, so changing
    the repeats never compounds. Returns the written .vdb path.
    """
    if vdb is None:
        raise ImportError("openvdb is needed to rebuild a density supercell.")
    upgrade_density_nodes(density_obj)
    base = density_obj.get('ase_base_vdb')
    if base is None:
        # imported before this was recorded: the volume still points at its
        # own single-cell grid, so adopt that as the base
        base = density_obj.data.filepath
        if not base:
            raise ValueError(f"{density_obj.name} has no grid file to tile from")
        density_obj['ase_base_vdb'] = base
    base = bpy.path.abspath(base)
    if not os.path.exists(base):
        raise FileNotFoundError(f"the density grid {base} is gone - re-import the file")

    grid = vdb.read(base, 'density')
    shape = density_obj.get('ase_grid_shape')
    if shape is None:
        # an older import: fall back to the active voxel box, which is the
        # true shape unless the density is exactly zero at a border plane
        bbox = grid.evalActiveVoxelBoundingBox()
        shape = tuple(np.asarray(bbox[1]) - np.asarray(bbox[0]) + 1)
    volume = np.zeros(tuple(int(n) for n in shape), dtype=float)
    grid.copyToArray(volume, ijk=(0, 0, 0))

    repeat = tuple(max(1, int(n)) for n in repeat)
    if repeat == (1, 1, 1):
        path = base
    else:
        tiled = vdb.FloatGrid()
        tiled.copyFromArray(np.ascontiguousarray(np.tile(volume, repeat)))
        tiled.transform = grid.transform      # voxel size is unchanged
        tiled.gridClass = grid.gridClass
        tiled.name = 'density'
        stem = os.path.splitext(base)[0]
        path = f'{stem}_{repeat[0]}x{repeat[1]}x{repeat[2]}.vdb'
        vdb.write(path, tiled)

    density_obj.data.filepath = path
    density_obj['ase_density_repeat'] = list(repeat)
    density_obj.data.update_tag()
    density_obj.update_tag()
    return path


def cube2vol(filename, filepath=os.environ.get('HOME'), modifier='GeometryNodes'):
    with open(filename, 'r') as f:
        atoms = read_cube(f, read_data=True, verbose=True)
        ORIGIN = atoms['origin']
    VOLUME = atoms['data']
    SPACING = atoms['spacing']
    return data2vol(VOLUME, SPACING, ORIGIN, filename, modifier=modifier)


def is_ams_tape41(filename):
    """True for an AMS/BAND TAPE41 (binary KF) volume file.

    TAPE41 carries no informative extension - by SCM convention it is
    either literally named ``TAPE41``, renamed to ``*.t41`` (see the AMS
    docs: "renaming it to foobar.t41 will allow AMSview to read it"), or
    just kept as a ``*.TAPE41``/``*_TAPE41`` suffix (the most common case
    in practice: copying it out of a results directory next to other
    same-named files, e.g. ``restart.TAPE41``). Anything else is left to
    the normal ase.io.read dispatch.
    """
    base = os.path.basename(filename).upper()
    return base == 'TAPE41' or base.endswith('.T41') or base.endswith('TAPE41')


def _read_atoms_from_kf(kf):
    """Read the atomic geometry out of a TAPE41's Geometry section.

    Same layout ase.io.cube uses internally, just reached through
    plams' KFFile instead of the cube text format: ``nnuc`` atoms,
    ``labels`` a fixed-width concatenated string of element symbols,
    ``xyznuc`` a flat xyz array in Bohr.
    """
    nnuc = kf.read('Geometry', 'nnuc')
    raw_labels = kf.read('Geometry', 'labels')
    xyz_bohr = np.array(kf.read('Geometry', 'xyznuc')).reshape(-1, 3)
    label_len = len(raw_labels) // nnuc
    symbols = [raw_labels[i * label_len:(i + 1) * label_len].strip() for i in range(nnuc)]
    return Atoms(symbols=symbols, positions=xyz_bohr * BOHR_TO_ANG)


def read_tape41(filepath, volumes=None, kfkey='FOO'):
    """Read atoms and one or more named volumes out of an AMS TAPE41.

    volumes: specific FOO-section names to read (e.g. the NOCV pairs
    ``['dRhoNOCV=1,k=1', 'dRhoNOCV=2,k=1']``); None reads every entry
    under kfkey. Returns ``(atoms, {name: (data, spacing, origin)})``
    with data indexed like ase.io.cube (shape (nx, ny, nz)) and
    spacing/origin in Angstrom, ready for data2vol().
    """
    from scm.plams.tools.kftools import KFFile

    kf = KFFile(str(filepath))
    atoms = _read_atoms_from_kf(kf)

    nx = kf.read('Grid', 'nr of points x')
    ny = kf.read('Grid', 'nr of points y')
    nz = kf.read('Grid', 'nr of points z')
    vx = np.array(kf.read('Grid', 'x-vector')) * BOHR_TO_ANG
    vy = np.array(kf.read('Grid', 'y-vector')) * BOHR_TO_ANG
    vz = np.array(kf.read('Grid', 'z-vector')) * BOHR_TO_ANG
    origin = np.array(kf.read('Grid', 'Start_point')) * BOHR_TO_ANG

    skeleton = kf.get_skeleton()
    names = volumes if volumes else list(skeleton.get(kfkey, []))
    if not names:
        raise ValueError(f"no volumes found under KF section '{kfkey}' in {filepath}")

    result = {}
    for name in names:
        raw = kf.read(kfkey, name)
        data = np.reshape(np.array(raw), (nx, ny, nz), order='F')
        result[name] = (data, (vx, vy, vz), origin)
    return atoms, result


def tape41_import(filepath, filename, modifier='GeometryNodes', volumes=None, found=None):
    """Import one or more volumes from an AMS TAPE41 (e.g. NOCV deformation
    densities from a PEDANOCV restart) straight into Blender, without a
    separate cube-file conversion step.

    Every requested volume becomes its own volume object with the usual
    blue/red '+ material'/'- material' pair (same convention as cube2vol),
    named after its KF variable so e.g. 'dRhoNOCV=1,k=1' and
    'dRhoNOCV=2,k=1' land as two separate, independently toggleable
    objects - matching how NOCV pairs are conventionally rendered one at
    a time, not overlaid. Returns the list of density objects.

    Pass `found` (the volumes dict from an already-done read_tape41()
    call) to skip re-reading the KF file - the importer already needs the
    atoms out of it before this point, so it reads once and hands the
    volumes dict through rather than opening the file twice.
    """
    if found is None:
        _, found = read_tape41(filepath, volumes=volumes)
    objects = []
    for name, (data, spacing, origin) in found.items():
        obj = data2vol(data, spacing, origin, f'{filepath}_{name}', modifier=modifier)
        obj.name = name
        objects.append(obj)
    return objects


def chgcar2vol(filename, modifier='GeometryNodes', density=None):
    """Import a VASP charge density (CHGCAR/PARCHG/AECCAR) as volumes.

    The grid spans the unit cell, so the spacing vectors are the cell
    vectors divided by the grid shape; VaspChargeDensity already divides
    the values by the cell volume (density in e/A^3).

    Returns a list of volume objects: the total charge density, and for
    spin-polarized files a second volume with the spin difference (green
    isosurface for spin-up excess, pink for spin-down).
    """
    if density is None:
        density = VaspChargeDensity(filename)
    atoms = density.atoms[-1]
    volume = density.chg[-1]
    cell = atoms.get_cell()
    spacing = [cell[i] / volume.shape[i] for i in range(3)]
    objects = [data2vol(volume, spacing, (0.0, 0.0, 0.0), filename, modifier=modifier)]

    chgdiff = getattr(density, 'chgdiff', [])
    if len(chgdiff):
        spin_up_mat = bpy.data.materials.get('+ spin material') or newShader('+ spin material', 0.1, 0.75, 0.25)   # green
        spin_down_mat = bpy.data.materials.get('- spin material') or newShader('- spin material', 0.95, 0.35, 0.6)  # pink
        spin_obj = data2vol(chgdiff[-1], spacing, (0.0, 0.0, 0.0),
                            os.path.splitext(filename)[0] + '_spin',
                            modifier=modifier,
                            plus_material=spin_up_mat,
                            minus_material=spin_down_mat)
        objects.append(spin_obj)
    return objects
