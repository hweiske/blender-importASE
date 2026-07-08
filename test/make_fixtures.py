"""(Re)generate the test fixture set in test/fixtures/.

Each file maps to importer features that must keep working (see
smoke_test.py). Deterministic - safe to re-run; add new structures here
or drop files into test/fixtures/ directly.

Run with any python that has ase + numpy:
    python test/make_fixtures.py
"""
import os
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io.cube import write_cube
from ase.calculators.vasp import VaspChargeDensity
import ase.io  # noqa: F401 - used via ase.io.write

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
os.makedirs(FIXTURES, exist_ok=True)


def main():
    # note: some fixtures are real, user-provided data and are NOT
    # regenerated here (crystal.cif, LED_dens.cube, LED_color.cube, ...);
    # this script only maintains the synthetic minimal set

    # tiny trajectory (5 frames of jiggling water): animation paths
    rng = np.random.default_rng(42)
    water = Atoms('OH2', positions=[(0, 0, 0), (0.76, 0.59, 0), (-0.76, 0.59, 0)])
    traj = []
    for _ in range(5):
        im = water.copy()
        im.positions += rng.normal(0, 0.05, im.positions.shape)
        traj.append(im)
    ase.io.write(f'{FIXTURES}/traj.xyz', traj)

    # MO-like cube with +/- lobes and a color-density cube:
    # marching-cubes density-mesh importer
    mo_cell = Atoms('OH2', positions=[(4, 4, 4), (4.76, 4.59, 4), (3.24, 4.59, 4)],
                    cell=[8, 8, 8])
    mx, my, mz = np.mgrid[0:8:24j, 0:8:24j, 0:8:24j]
    mo_data = (mx - 4) * np.exp(-((mx - 4)**2 + (my - 4)**2 + (mz - 4)**2) / 2.5)
    with open(f'{FIXTURES}/mo.cube', 'w') as f:
        write_cube(f, mo_cell, data=mo_data)

    # bonded chain + charges csv: partial-charges importer
    chain = Atoms('OHC', positions=[(0, 0, 0), (0.5, 0, 0), (1.0, 0, 0)])
    ase.io.write(f'{FIXTURES}/chargemol.xyz', chain)
    with open(f'{FIXTURES}/charges.csv', 'w') as f:
        f.write('element,charge\nO,-0.6\nH,0.25\nC,0.35\n')

    # rocksalt supercell: coordination-polyhedra importer
    ase.io.write(f'{FIXTURES}/nacl.extxyz', bulk('NaCl', 'rocksalt', a=5.64) * (2, 2, 2))

    # density grids around a water molecule: volume density importers
    water_cell = Atoms('OH2', positions=[(2, 2, 2), (2.76, 2.59, 2), (1.24, 2.59, 2)],
                       cell=[4, 4, 4], pbc=True)
    x, y, z = np.mgrid[0:4:20j, 0:4:20j, 0:4:20j]
    rho = np.exp(-((x - 2)**2 + (y - 2)**2 + (z - 2)**2))
    with open(f'{FIXTURES}/water.cube', 'w') as f:
        write_cube(f, water_cell, data=rho)
    vcd = VaspChargeDensity(filename=None)
    vcd.atoms = [water_cell]
    vcd.chg = [rho]
    # spin-polarized: also exercises the green/pink spin-difference volume
    vcd.chgdiff = [rho * np.sign(x - 2)]
    vcd.write(f'{FIXTURES}/CHGCAR')

    print(f'fixtures written to {FIXTURES}')


if __name__ == '__main__':
    main()
