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
step('bonds_fromnodes', lambda: run_import(f'{SCRATCH}/crystal.xyz',
     representation='bonds_fromnodes', animate=False))
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
