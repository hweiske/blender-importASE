import bpy
import numpy as np
from ase.io.cube import read_cube
from ase.calculators.vasp import VaspChargeDensity
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
from .node_networks.electron_density_nodes import visualize_edensity_node_group, newShader
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
    if not os.path.isfile(TMPFILE):
        vdb.write(TMPFILE, GRID)
    _ = bpy.ops.object.volume_import(filepath=TMPFILE, location=origin)
    density_obj = bpy.context.active_object
    visualize_edensity_node_group()
    bpy.ops.object.modifier_add(type='NODES')
    node = bpy.data.node_groups["visualize_edensity"]
    bpy.context.object.modifiers[modifier].node_group = node
    if plus_material is None:
        plus_material = bpy.data.materials["+ material"]
    if minus_material is None:
        minus_material = bpy.data.materials["- material"]
    bpy.context.object.modifiers[modifier]["Socket_9"] = plus_material
    bpy.context.object.modifiers[modifier]["Socket_10"] = minus_material
    toggle(bpy.context.object, SET=False)
    return density_obj


def cube2vol(filename, filepath=os.environ.get('HOME'), modifier='GeometryNodes'):
    with open(filename, 'r') as f:
        atoms = read_cube(f, read_data=True, verbose=True)
        ORIGIN = atoms['origin']
    VOLUME = atoms['data']
    SPACING = atoms['spacing']
    return data2vol(VOLUME, SPACING, ORIGIN, filename, modifier=modifier)


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
        spin_up_mat = newShader('+ spin material', 0.1, 0.75, 0.25)   # green
        spin_down_mat = newShader('- spin material', 0.95, 0.35, 0.6)  # pink
        spin_obj = data2vol(chgdiff[-1], spacing, (0.0, 0.0, 0.0),
                            os.path.splitext(filename)[0] + '_spin',
                            modifier=modifier,
                            plus_material=spin_up_mat,
                            minus_material=spin_down_mat)
        objects.append(spin_obj)
    return objects
