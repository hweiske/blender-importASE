"""Headless smoke test for blender_importASE on Blender 5.x.

Run: blender -b --factory-startup --python smoke_test.py
"""
import sys
import traceback
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import bpy
import numpy as np
from ase import Atoms
import ase.io

print(f"### Blender {bpy.app.version_string}, Python {sys.version.split()[0]}")

# --- build test inputs -------------------------------------------------
SCRATCH = '/tmp/blender_importASE_smoketest'
os.makedirs(SCRATCH, exist_ok=True)

# small periodic crystal
si = Atoms('Si2O4', positions=[(0, 0, 0), (2.3, 2.3, 2.3), (1.2, 1.2, 1.2),
                               (3.4, 3.4, 1.2), (1.2, 3.4, 3.4), (3.4, 1.2, 3.4)],
           cell=[4.6, 4.6, 4.6], pbc=True)
ase.io.write(f'{SCRATCH}/crystal.xyz', si)

# tiny trajectory (5 frames of jiggling water)
rng = np.random.default_rng(42)
water = Atoms('OH2', positions=[(0, 0, 0), (0.76, 0.59, 0), (-0.76, 0.59, 0)])
traj = []
for i in range(5):
    im = water.copy()
    im.positions += rng.normal(0, 0.05, im.positions.shape)
    traj.append(im)
ase.io.write(f'{SCRATCH}/traj.xyz', traj)

# MO-like cube with +/- lobes and a color-density cube for the
# marching-cubes density-mesh importer
mo_cell = Atoms('OH2', positions=[(4, 4, 4), (4.76, 4.59, 4), (3.24, 4.59, 4)], cell=[8, 8, 8])
mx, my, mz = np.mgrid[0:8:24j, 0:8:24j, 0:8:24j]
mo_data = (mx - 4) * np.exp(-((mx - 4)**2 + (my - 4)**2 + (mz - 4)**2) / 2.5)
from ase.io.cube import write_cube as _write_cube
with open(f'{SCRATCH}/mo.cube', 'w') as f:
    _write_cube(f, mo_cell, data=mo_data)
with open(f'{SCRATCH}/colordens.cube', 'w') as f:
    _write_cube(f, mo_cell, data=np.exp(-((mz - 4)**2) / 8.0))

# bonded chain + charges csv for the partial-charges importer
chain = Atoms('OHC', positions=[(0, 0, 0), (0.5, 0, 0), (1.0, 0, 0)])
ase.io.write(f'{SCRATCH}/chargemol.xyz', chain)
with open(f'{SCRATCH}/charges.csv', 'w') as f:
    f.write('element,charge\nO,-0.6\nH,0.25\nC,0.35\n')

# rocksalt supercell for the coordination-polyhedra importer
from ase.build import bulk
ase.io.write(f'{SCRATCH}/nacl.extxyz', bulk('NaCl', 'rocksalt', a=5.64) * (2, 2, 2))

# density grids: gaussian blob around a water molecule (.cube and CHGCAR)
from ase.io.cube import write_cube
from ase.calculators.vasp import VaspChargeDensity
water_cell = Atoms('OH2', positions=[(2, 2, 2), (2.76, 2.59, 2), (1.24, 2.59, 2)],
                   cell=[4, 4, 4], pbc=True)
x, y, z = np.mgrid[0:4:20j, 0:4:20j, 0:4:20j]
rho = np.exp(-((x - 2)**2 + (y - 2)**2 + (z - 2)**2))
with open(f'{SCRATCH}/water.cube', 'w') as f:
    write_cube(f, water_cell, data=rho)
vcd = VaspChargeDensity(filename=None)
vcd.atoms = [water_cell]
vcd.chg = [rho]
# spin-polarized: also exercises the green/pink spin-difference volume
vcd.chgdiff = [rho * np.sign(x - 2)]
vcd.write(f'{SCRATCH}/CHGCAR')
for stale in ('water_density.vdb', 'CHGCAR_density.vdb'):
    try:
        os.remove(f'{SCRATCH}/{stale}')
    except FileNotFoundError:
        pass

# --- register addon ----------------------------------------------------
results = {}

def step(name, fn):
    try:
        fn()
        results[name] = 'OK'
        print(f"### PASS {name}")
    except Exception:
        results[name] = 'FAIL'
        print(f"### FAIL {name}")
        traceback.print_exc()

def fresh_scene():
    bpy.ops.wm.read_factory_settings(use_empty=False)
    for ob in list(bpy.data.objects):
        bpy.data.objects.remove(ob, do_unlink=True)
    for coll in ('meshes', 'materials', 'node_groups', 'collections'):
        data = getattr(bpy.data, coll)
        for item in list(data):
            data.remove(item)

import blender_importASE
step('register', blender_importASE.register)

from blender_importASE.ui import import_ase_molecule

def run_import(path, **kw):
    fresh_scene()
    import_ase_molecule(path, os.path.basename(path), **kw)
step('nodes_crystal', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation='nodes', animate=False))
step('ballsnsticks_crystal', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation="Balls'n'Sticks", long_bonds=True, unit_cell=True, animate=False))
step('ballsnsticks_nolongbond', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation="Balls'n'Sticks", long_bonds=False, animate=False))
step('licorice', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation='Licorice', long_bonds=True, animate=False))
step('vdw', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation='VDW', animate=False))
step('3D_print', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation='3D_print', animate=False))
step('trajectory_nodes', lambda: run_import(f'{SCRATCH}/traj.xyz',
     representation='nodes', animate=True))
step('trajectory_keyframes', lambda: run_import(f'{SCRATCH}/traj.xyz',
     representation="Balls'n'Sticks", overwrite=False, long_bonds=False, animate=True))
step('trajectory_longbonds', lambda: run_import(f'{SCRATCH}/traj.xyz',
     representation="Balls'n'Sticks", overwrite=False, long_bonds=True, animate=True))
step('cube_density', lambda: run_import(f'{SCRATCH}/water.cube',
     representation='nodes', animate=False, read_density=True))
step('chgcar_density', lambda: run_import(f'{SCRATCH}/CHGCAR',
     representation='nodes', animate=False, read_density=True))
def run_polyhedra():
    fresh_scene()
    from blender_importASE.polyhedra import import_polyhedra
    import_polyhedra(f'{SCRATCH}/nacl.extxyz', 'nacl.extxyz')
    obj = next(o for o in bpy.data.objects if 'polyhedra' in o.name)
    assert len(obj.data.polygons) > 0, 'no polyhedra faces generated'

step('polyhedra', run_polyhedra)

def run_density_mesh():
    from importlib import util
    if util.find_spec('skimage') is None:
        print('scikit-image not installed - skipping the actual import')
        return
    fresh_scene()
    from blender_importASE.density_mesh import import_density_mesh
    import_density_mesh(f'{SCRATCH}/mo.cube', 'mo.cube', iso_value=0.05,
                        color_filepath=f'{SCRATCH}/colordens.cube')
    obj = bpy.data.objects['mo_isomesh']
    assert len(obj.data.polygons) > 0, 'no isosurface faces generated'
    assert 'density_color' in obj.data.color_attributes, 'missing color attribute'

step('density_mesh', run_density_mesh)

def run_charges():
    fresh_scene()
    from blender_importASE.charges import import_charges
    import_charges(f'{SCRATCH}/chargemol.xyz', 'chargemol.xyz',
                   charge_filepath=f'{SCRATCH}/charges.csv')
    obj = next(o for o in bpy.data.objects if 'charges' in o.name)
    charge_vals = [d.value for d in obj.data.attributes['charge'].data]
    assert charge_vals and abs(charge_vals[0] + 0.6) < 1e-5, 'charge attribute wrong'
    slots = [s.material.name for s in obj.material_slots if s.material]
    assert 'charge_atoms' in slots and 'color_curve_charge' in slots, 'charge materials missing'

step('charges', run_charges)

def run_export_3dprint():
    import zipfile
    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/crystal.xyz', 'crystal.xyz',
                        representation='3D_print', animate=False,
                        read_density=False, outline=False, add_supercell=False)
    atom = next(o for o in bpy.data.objects if o.name.split('.')[0] in ('Si', 'O'))
    bpy.context.view_layer.objects.active = atom
    zip_path = f'{SCRATCH}/print_export.zip'
    bpy.ops.export_mesh.ase_3dprint(filepath=zip_path)
    with zipfile.ZipFile(zip_path) as z:
        names = set(z.namelist())
    assert names == {'atoms_O.stl', 'atoms_Si.stl', 'bonds.stl', 'supports.stl'}, names

step('export_3dprint', run_export_3dprint)

def run_export_xyz():
    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/crystal.xyz', 'crystal.xyz',
                        representation='nodes', animate=False,
                        read_density=False, outline=False, add_supercell=False)
    obj = bpy.data.objects['O4Si2_crystal']
    bpy.context.view_layer.objects.active = obj
    xyz_path = f'{SCRATCH}/roundtrip.xyz'
    bpy.ops.export_mesh.ase_xyz(filepath=xyz_path)
    back = ase.io.read(xyz_path)
    assert back.get_chemical_formula() == 'O4Si2', back.get_chemical_formula()

step('export_xyz', run_export_xyz)
step('operator_via_ops', lambda: (
    fresh_scene(),
    bpy.ops.import_mesh.ase(directory=SCRATCH, files=[{"name": "crystal.xyz"}]),
))
step('unregister', blender_importASE.unregister)

print('### SUMMARY')
for k, v in results.items():
    print(f"### {v:4s} {k}")
if any(v == 'FAIL' for v in results.values()):
    sys.exit(1)
