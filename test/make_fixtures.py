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
    # only written when absent: the checked-in CHGCAR is real, user-provided
    # data (spin-polarized bcc Fe with a POTCAR-style 'Fe/' species line)
    if not os.path.exists(f'{FIXTURES}/CHGCAR'):
        vcd = VaspChargeDensity(filename=None)
        vcd.atoms = [water_cell]
        vcd.chg = [rho]
        # spin-polarized: also exercises the green/pink spin-difference volume
        vcd.chgdiff = [rho * np.sign(x - 2)]
        vcd.write(f'{FIXTURES}/CHGCAR')

    # minimal synthetic AMS TAPE41: two small named volumes under [FOO], the
    # same layout a real PEDANOCV restart's NOCVdRhoPlot writes (see
    # import_cubefiles.read_tape41). Only written when plams is importable -
    # unlike the other fixtures this one needs an extra dependency, so it's
    # optional the way the addon's own plams support is optional.
    try:
        from scm.plams.tools.kftools import KFFile
        BOHR = 1.8897259886  # Angstrom -> Bohr, matching import_cubefiles' inverse
        nx = ny = nz = 6
        gx, gy, gz = np.mgrid[0:1:nx * 1j, 0:1:ny * 1j, 0:1:nz * 1j]
        lobe1 = np.exp(-((gx - 0.3) ** 2 + (gy - 0.5) ** 2 + (gz - 0.5) ** 2) / 0.02)
        lobe2 = -np.exp(-((gx - 0.7) ** 2 + (gy - 0.5) ** 2 + (gz - 0.5) ** 2) / 0.02)
        tape41_path = f'{FIXTURES}/synthetic.TAPE41'
        if os.path.exists(tape41_path):
            os.remove(tape41_path)
        kf = KFFile(tape41_path)
        positions_bohr = np.array([(2, 2, 2), (2.76, 2.59, 2), (1.24, 2.59, 2)]) * BOHR
        kf.write('Geometry', 'nnuc', 3)
        # fixed-width, equal-length slots per label (matching real AMS TAPE41
        # output) - read_tape41 divides the raw string length by nnuc, so
        # uneven padding (plain space-separated symbols) slices wrong
        kf.write('Geometry', 'labels', ''.join(s.ljust(8) for s in ['O', 'H', 'H']))
        kf.write('Geometry', 'xyznuc', positions_bohr.flatten().tolist())
        kf.write('Grid', 'nr of points x', nx)
        kf.write('Grid', 'nr of points y', ny)
        kf.write('Grid', 'nr of points z', nz)
        kf.write('Grid', 'x-vector', [4.0 * BOHR, 0.0, 0.0])
        kf.write('Grid', 'y-vector', [0.0, 4.0 * BOHR, 0.0])
        kf.write('Grid', 'z-vector', [0.0, 0.0, 4.0 * BOHR])
        kf.write('Grid', 'Start_point', [0.0, 0.0, 0.0])
        kf.write('FOO', 'dRhoNOCV=1,k=1', lobe1.flatten(order='F').tolist())
        kf.write('FOO', 'dRhoNOCV=2,k=1', lobe2.flatten(order='F').tolist())
        print(f'wrote {tape41_path}')
    except ImportError:
        print('plams not installed - skipping synthetic.TAPE41 (optional fixture)')

    print(f'fixtures written to {FIXTURES}')


if __name__ == '__main__':
    main()
