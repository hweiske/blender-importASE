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
def run_tape41_density():
    # an actual import, not importlib.util.find_spec: a partial/broken
    # install (e.g. a failed pip leaving a bare 'scm' namespace dir behind)
    # can make find_spec succeed while the real import still fails - the
    # same distinction check_dependency() draws in __init__.py.
    try:
        import scm.plams  # noqa: F401
    except ImportError:
        print('plams not installed - skipping the actual import')
        return
    run_import(f'{SCRATCH}/synthetic.TAPE41', representation='nodes',
              animate=False, read_density=True)
    names = {o.name for o in bpy.data.objects if o.type == 'VOLUME'}
    assert names == {'dRhoNOCV=1,k=1', 'dRhoNOCV=2,k=1'}, names

step('tape41_density', run_tape41_density)
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

def run_polyhedra_molecules():
    """The 3x3x3 molecule selection: every molecule reaching into the cell
    comes out whole, a framework still closes its boundary polyhedra, the
    margin pulls in the surrounding molecules, and 'import unit cell'
    draws the flat black cell."""
    import numpy as np
    from blender_importASE.polyhedra import build_polyhedra_atoms, import_polyhedra

    def per_carbon(new_atoms):
        """(C, H) neighbor counts of every carbon, measured on the
        imported non-periodic positions - a whole benzene is all (2, 1)."""
        pos = new_atoms.get_positions()
        sym = new_atoms.get_chemical_symbols()
        dist = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(dist, 99.0)
        return [(sum(1 for j in range(len(pos)) if sym[j] == 'C' and dist[i, j] < 1.5),
                 sum(1 for j in range(len(pos)) if sym[j] == 'H' and dist[i, j] < 1.2))
                for i in range(len(pos)) if sym[i] == 'C']

    # a benzene ring centered on the cell corner: four of its images reach
    # into the cell, and each of those arrives as a complete ring
    atoms = ase.io.read(f'{SCRATCH}/molcrystal.extxyz')
    whole, _ = build_polyhedra_atoms(atoms, complete_molecules=True)
    assert whole.get_chemical_formula() == 'C24H24', whole.get_chemical_formula()
    assert all(c == (2, 1) for c in per_carbon(whole)), sorted(set(per_carbon(whole)))
    # ... which the plain image expansion does not manage: open fragments
    cut, _ = build_polyhedra_atoms(atoms, complete_molecules=False)
    assert any(c != (2, 1) for c in per_carbon(cut)), sorted(set(per_carbon(cut)))

    # a real bromoantimonate: every alkylammonium carbon and nitrogen must
    # come out 4-coordinate and every Sb with its six Br, or a molecule was
    # cut by the cell
    real = ase.io.read(f'{SCRATCH}/crystal.cif')

    def miscoordinated(new_atoms):
        pos = new_atoms.get_positions()
        sym = np.array(new_atoms.get_chemical_symbols())
        dist = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(dist, 99.0)
        light = ((dist < 1.8) & np.isin(sym, ['C', 'N', 'H'])[None, :]).sum(1)
        bromine = ((dist < 3.2) & (sym == 'Br')[None, :]).sum(1)
        return (int((light[np.isin(sym, ['C', 'N'])] != 4).sum()),
                int((bromine[sym == 'Sb'] != 6).sum()))

    grown, faces = build_polyhedra_atoms(real, complete_molecules=True)
    assert miscoordinated(grown) == (0, 0), miscoordinated(grown)
    assert len(faces) > 0, 'no polyhedra on the bromoantimonate'
    legacy, _ = build_polyhedra_atoms(real, complete_molecules=False)
    assert miscoordinated(legacy)[0] > 0, 'expansion-only path unexpectedly complete'
    # the margin brings the neighboring molecules along
    wider, _ = build_polyhedra_atoms(real, complete_molecules=True, cell_margin=3.0)
    assert len(wider) > len(grown), (len(wider), len(grown))
    assert miscoordinated(wider) == (0, 0), miscoordinated(wider)

    # one whole copy per molecule instead of every image reaching the cell
    single, _ = build_polyhedra_atoms(real, complete_molecules=True,
                                      all_images=False)
    assert single.get_chemical_formula() == real.get_chemical_formula(), \
        single.get_chemical_formula()
    assert miscoordinated(single) == (0, 0), miscoordinated(single)

    # the regression the fixtures above cannot catch: a Br bridging the
    # hydrogens of a benzene and of its own periodic image. At a bond
    # criterion loose enough to count that 2.30 A H...Br contact (which is
    # what ASE's neighbor-list skin silently did) the ring is part of an
    # endless chain and cannot be completed; at 1.3 x (r1+r2) = 1.96 A it
    # is a molecule and the Br is a free ion.
    from blender_importASE.polyhedra import (bond_neighbors, grow_shells,
                                             _components)
    chain = ase.io.read(f'{SCRATCH}/hbond_chain.extxyz')

    def closed_components(cutoff):
        nl = bond_neighbors(chain, cutoff)
        nl.update(chain)
        neighbors = [nl.get_neighbors(i) for i in range(len(chain))]
        return [grow_shells(neighbors, [(component[0], (0, 0, 0))])[2]
                for component in _components([idx for idx, _ in neighbors])]

    assert closed_components(1.3) == [True, True], closed_components(1.3)
    assert closed_components(1.6) == [False], closed_components(1.6)
    tight, _ = build_polyhedra_atoms(chain, complete_molecules=True)
    assert tight.get_chemical_formula() == 'C6H6Br', tight.get_chemical_formula()
    assert all(c == (2, 1) for c in per_carbon(tight)), sorted(set(per_carbon(tight)))

    # a framework has no molecule to complete - its atoms in the cell are
    # grown by framework_shells bonded shells, and every boundary
    # polyhedron still closes at the default of 1
    nacl = ase.io.read(f'{SCRATCH}/nacl.extxyz')
    counts = []
    for shells in (0, 1, 2):
        frame, faces = build_polyhedra_atoms(nacl, complete_molecules=True,
                                             framework_shells=shells)
        counts.append((len(frame), len(faces)))
    assert counts[0] == (len(nacl), 16), counts[0]
    assert counts[1] == (71, 244), counts[1]   # unchanged by this rewrite
    assert counts[2][0] > counts[1][0], counts

    fresh_scene()
    import_polyhedra(f'{SCRATCH}/molcrystal.extxyz', 'molcrystal.extxyz',
                     unit_cell=True)
    cells = [o for o in bpy.data.objects if 'unitcell' in o.name]
    assert len(cells) == 1, [o.name for o in cells]
    mats = [m.name for m in bpy.data.materials if m.name.startswith('unit_cell')]
    assert mats == ['unit_cell'], mats
    mat = cells[0].data.materials[0]
    output = next(n for n in mat.node_tree.nodes if n.type == 'OUTPUT_MATERIAL')
    source = output.inputs['Surface'].links[0].from_node
    assert source.type == 'RGB', source.type
    assert tuple(source.outputs[0].default_value) == (0.0, 0.0, 0.0, 1.0), \
        tuple(source.outputs[0].default_value)
    # a molecule without a cell gets no cell object (and no 3x3x3 pass)
    fresh_scene()
    import_polyhedra(f'{SCRATCH}/water.xyz', 'water.xyz', unit_cell=True)
    assert not [o for o in bpy.data.objects if 'unitcell' in o.name], 'cell drawn without a cell'

step('polyhedra_molecules', run_polyhedra_molecules)

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
def run_element_colors():
    """Per-element colors: the default scheme, the sidebar swatch's path
    into the materials, and the sync of the 'atom_color' attribute the
    colored bonds read."""
    from blender_importASE import element_colors as ec
    from blender_importASE.utils import default_element_color
    # antimony's default is #C75D8F, stored linear like every base color
    assert max(abs(a - b) for a, b in
               zip(default_element_color('Sb'),
                   (0.571125, 0.109462, 0.274677, 1.0))) < 1e-5, \
        default_element_color('Sb')

    fresh_scene()
    ase.io.write(f'{SCRATCH}/antimony.xyz',
                 ase.Atoms('SbH', positions=[(0, 0, 0), (1.7, 0, 0)]))
    import_ase_molecule(f'{SCRATCH}/antimony.xyz', 'antimony.xyz',
                        representation='nodes', animate=False,
                        read_density=False)
    obj = next(o for o in bpy.data.objects
               if o.type == 'MESH' and 'atom_color' in o.data.attributes)
    assert ec.structure_symbols(obj) == ['H', 'Sb'], ec.structure_symbols(obj)
    assert ec.element_color_socket('Sb') is not None, 'no Sb color to draw'

    def sb_color():
        numbers = [round(d.value) for d in obj.data.attributes['element'].data]
        i = numbers.index(51)
        return tuple(obj.data.attributes['atom_color'].data[i].color)

    # imported at the default color, materials and attribute agreeing
    assert max(abs(a - b) for a, b in
               zip(sb_color(), default_element_color('Sb'))) < 1e-5, sb_color()

    # a picked color reaches the atom material, the bond materials and the
    # attribute
    ec.set_element_color('Sb', (0.1, 0.2, 0.3))
    assert sb_color()[:3] == (0.10000000149011612, 0.20000000298023224,
                              0.30000001192092896), sb_color()
    bsdf = bpy.data.materials['Sb-bond'].node_tree.nodes['Principled BSDF']
    assert tuple(bsdf.inputs[0].default_value)[:3] == sb_color()[:3], 'bond material'

    # editing the material alone (Material Properties, a script) is carried
    # over to the attribute by the depsgraph handler
    atom_bsdf = bpy.data.materials['Sb'].node_tree.nodes['Principled BSDF']
    atom_bsdf.inputs[0].default_value = (0.4, 0.5, 0.6, 1.0)
    bpy.context.view_layer.update()
    assert abs(sb_color()[0] - 0.4) < 1e-6, sb_color()

    # a second import keeps the picked color instead of resetting it
    import_ase_molecule(f'{SCRATCH}/antimony.xyz', 'antimony.xyz',
                        representation='nodes', animate=False,
                        read_density=False)
    assert abs(ec.get_element_color('Sb')[0] - 0.4) < 1e-6, ec.get_element_color('Sb')
    again = [o for o in bpy.data.objects
             if o.type == 'MESH' and 'atom_color' in o.data.attributes][-1]
    numbers = [round(d.value) for d in again.data.attributes['element'].data]
    color = tuple(again.data.attributes['atom_color'].data[numbers.index(51)].color)
    assert abs(color[0] - 0.4) < 1e-6, color

    # ... and the operators put it back
    bpy.context.view_layer.objects.active = obj
    bpy.ops.ase.reset_element_colors()
    bpy.ops.ase.sync_element_colors()
    assert max(abs(a - b) for a, b in
               zip(sb_color(), default_element_color('Sb'))) < 1e-5, sb_color()

def run_lead_defaults():
    """Lead's own entry in the color scheme: a violet, fully metallic."""
    from blender_importASE.utils import default_element_color
    assert default_element_color('Pb') == (0.2, 0.0, 0.5, 1.0), default_element_color('Pb')

    fresh_scene()
    ase.io.write(f'{SCRATCH}/lead.xyz', ase.Atoms('Pb2', positions=[(0, 0, 0), (3.5, 0, 0)]))
    import_ase_molecule(f'{SCRATCH}/lead.xyz', 'lead.xyz', representation='nodes',
                        animate=False, read_density=False)
    bsdf = next(n for n in bpy.data.materials['Pb'].node_tree.nodes
                if n.type == 'BSDF_PRINCIPLED')
    assert tuple(round(v, 4) for v in bsdf.inputs['Base Color'].default_value) \
        == (0.2, 0.0, 0.5, 1.0), tuple(bsdf.inputs['Base Color'].default_value)
    assert bsdf.inputs['Metallic'].default_value == 1.0, bsdf.inputs['Metallic'].default_value
    assert bsdf.inputs['Roughness'].default_value == 0.5, bsdf.inputs['Roughness'].default_value

step('lead_defaults', run_lead_defaults)
step('element_colors', run_element_colors)
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
