"""Headless smoke test for blender_importASE on Blender 5.x.

Run: blender -b --factory-startup --python smoke_test.py
"""
import sys
import traceback
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import bpy
import ase.io

print(f"### Blender {bpy.app.version_string}, Python {sys.version.split()[0]}")

# --- test inputs --------------------------------------------------------
# The structures live in test/fixtures/ (regenerate with make_fixtures.py;
# every file there is part of the contract - the importers must handle it).
# They are copied to a scratch directory so outputs (.vdb files etc.)
# never end up in the repository.
import shutil
FIXTURES = os.path.join(REPO, 'test', 'fixtures')
SCRATCH = '/tmp/blender_importASE_smoketest'
shutil.rmtree(SCRATCH, ignore_errors=True)
shutil.copytree(FIXTURES, SCRATCH)

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
step('nodes_crystal', lambda: run_import(f'{SCRATCH}/crystal.cif',
     representation='nodes', animate=False))
step('ballsnsticks_crystal', lambda: run_import(f'{SCRATCH}/crystal.cif',
     representation="Balls'n'Sticks", long_bonds=True, unit_cell=True, animate=False))
step('ballsnsticks_nolongbond', lambda: run_import(f'{SCRATCH}/crystal.cif',
     representation="Balls'n'Sticks", long_bonds=False, animate=False))
step('licorice', lambda: run_import(f'{SCRATCH}/crystal.cif',
     representation='Licorice', long_bonds=True, animate=False))
step('vdw', lambda: run_import(f'{SCRATCH}/crystal.cif',
     representation='VDW', animate=False))
step('3D_print', lambda: run_import(f'{SCRATCH}/crystal.cif',
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
    import_polyhedra(f'{SCRATCH}/nacl.extxyz', 'nacl.extxyz', outline=True)
    faces_obj = next(o for o in bpy.data.objects if o.name.endswith('_faces'))
    assert len(faces_obj.data.polygons) > 0, 'no polyhedra faces generated'
    assert not faces_obj.modifiers, 'polyhedra faces must stay modifier-free'
    structure = next(o for o in bpy.data.objects
                     if 'polyhedra' in o.name and o.type == 'MESH'
                     and not o.name.endswith(('_faces', '_table')))
    names = [m.node_group.name for m in structure.modifiers if m.node_group]
    assert any(n.startswith('outline') for n in names), names

step('polyhedra', run_polyhedra)

def run_density_mesh():
    from importlib import util
    if util.find_spec('skimage') is None:
        print('scikit-image not installed - skipping the actual import')
        return
    fresh_scene()
    from blender_importASE.density_mesh import import_density_mesh
    # plain +/- lobes (no color file); the colored path is covered by led_pair
    import_density_mesh(f'{SCRATCH}/mo.cube', 'mo.cube', iso_value=0.05)
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
    from ase.data import chemical_symbols
    fresh_scene()
    reference = ase.io.read(f'{SCRATCH}/crystal.cif')
    import_ase_molecule(f'{SCRATCH}/crystal.cif', 'crystal.cif',
                        representation='3D_print', animate=False,
                        read_density=False, outline=False, add_supercell=False)
    atom = next(o for o in bpy.data.objects
                if o.name.split('.')[0] in chemical_symbols and o.type == 'MESH')
    bpy.context.view_layer.objects.active = atom
    zip_path = f'{SCRATCH}/print_export.zip'
    bpy.ops.export_mesh.ase_3dprint(filepath=zip_path)
    with zipfile.ZipFile(zip_path) as z:
        names = set(z.namelist())
    expected = {f'atoms_{el}.stl' for el in set(reference.get_chemical_symbols())}
    expected |= {'bonds.stl', 'supports.stl'}
    assert names == expected, (names, expected)

step('export_3dprint', run_export_3dprint)

def run_export_xyz():
    fresh_scene()
    reference = ase.io.read(f'{SCRATCH}/crystal.cif')
    import_ase_molecule(f'{SCRATCH}/crystal.cif', 'crystal.cif',
                        representation='nodes', animate=False,
                        read_density=False, outline=False, add_supercell=False)
    obj = bpy.data.objects[reference.get_chemical_formula() + '_crystal']
    bpy.context.view_layer.objects.active = obj
    xyz_path = f'{SCRATCH}/roundtrip.xyz'
    bpy.ops.export_mesh.ase_xyz(filepath=xyz_path)
    back = ase.io.read(xyz_path)
    assert back.get_chemical_formula() == reference.get_chemical_formula(), \
        back.get_chemical_formula()

step('export_xyz', run_export_xyz)

def run_structure_sweep():
    """Every structure file in test/fixtures must import with the default
    nodes representation - drop new structures there to extend the set."""
    failures = []
    for fname in sorted(os.listdir(SCRATCH)):
        if os.path.splitext(fname)[1].lower() not in ('.xyz', '.extxyz', '.cif'):
            continue
        if fname == 'roundtrip.xyz':
            continue
        try:
            fresh_scene()
            import_ase_molecule(f'{SCRATCH}/{fname}', fname,
                                representation='nodes', animate=False,
                                read_density=False)
        except Exception as exc:
            failures.append(f'{fname}: {exc!r}')
    assert not failures, failures

step('structure_sweep', run_structure_sweep)

def run_led_pair():
    """Real density + color-density pair through the marching-cubes mesh
    importer."""
    from importlib import util
    if util.find_spec('skimage') is None:
        print('scikit-image not installed - skipping the actual import')
        return
    fresh_scene()
    from blender_importASE.density_mesh import import_density_mesh
    import_density_mesh(f'{SCRATCH}/LED_dens.cube', 'LED_dens.cube',
                        iso_value=0.05, preset='LED', import_atoms=True,
                        color_filepath=f'{SCRATCH}/LED_color.cube')
    obj = bpy.data.objects['LED_dens_isomesh']
    assert len(obj.data.polygons) > 0, 'no isosurface faces generated'
    assert obj.data.materials[0].name == 'LED material', obj.data.materials[0].name
    # the atoms came along as the nodes representation
    assert any('LED_dens' in o.name and o is not obj for o in bpy.data.objects), \
        'structure object missing'

step('led_pair', run_led_pair)
step('operator_via_ops', lambda: (
    fresh_scene(),
    bpy.ops.import_mesh.ase(directory=SCRATCH, files=[{"name": "crystal.cif"}]),
))
step('unregister', blender_importASE.unregister)

print('### SUMMARY')
for k, v in results.items():
    print(f"### {v:4s} {k}")
if any(v == 'FAIL' for v in results.values()):
    sys.exit(1)
