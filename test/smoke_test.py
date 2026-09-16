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
    # the faces are shaded by a Principled/Glass mix, both halves tinted by
    # the atom_color attribute
    material = faces_obj.data.materials[0]
    kinds = {n.bl_idname for n in material.node_tree.nodes}
    assert {'ShaderNodeBsdfPrincipled', 'ShaderNodeBsdfGlass', 'ShaderNodeMixShader',
            'ShaderNodeAttribute'} <= kinds, sorted(kinds)
    mix = next(n for n in material.node_tree.nodes if n.bl_idname == 'ShaderNodeMixShader')
    assert mix.inputs[0].default_value == 0.5, mix.inputs[0].default_value
    output = next(n for n in material.node_tree.nodes
                  if n.bl_idname == 'ShaderNodeOutputMaterial')
    # == not is: every lookup hands out a fresh python proxy for the node,
    # so identity comparison fails even on the same node
    assert output.inputs['Surface'].links[0].from_node == mix, \
        'the mix does not drive the surface'
    for node_type, socket in (('ShaderNodeBsdfPrincipled', 'Base Color'),
                              ('ShaderNodeBsdfGlass', 'Color')):
        node = next(n for n in material.node_tree.nodes if n.bl_idname == node_type)
        assert node.inputs[socket].is_linked, f'{node_type} is not tinted by atom_color'
    structure = next(o for o in bpy.data.objects
                     if 'polyhedra' in o.name and o.type == 'MESH'
                     and not o.name.endswith(('_faces', '_table')))
    names = [m.node_group.name for m in structure.modifiers if m.node_group]
    # the stack of the nodes representation: hide, supercell, atoms, outline
    kinds = [n.split('_')[0].split('.')[0].split(' ')[0] for n in names]
    assert kinds == ['hide', 'supercell', 'atoms', 'outline'], names
    # the faces only get the structure's supercell, never the outline
    assert [m.node_group.name for m in faces_obj.modifiers] == [names[1]], \
        [m.node_group.name for m in faces_obj.modifiers]

    def evaluated(obj, attr='vertices'):
        depsgraph = bpy.context.evaluated_depsgraph_get()
        mesh = obj.evaluated_get(depsgraph).to_mesh()
        count = len(getattr(mesh, attr))
        obj.evaluated_get(depsgraph).to_mesh_clear()
        return count

    from blender_importASE.node_networks.compat import set_mod_input
    faces_before = evaluated(faces_obj, 'polygons')
    assert faces_before == len(faces_obj.data.polygons), faces_before

    # hiding an element takes its atoms out of the structure
    hide = structure.modifiers[0]
    hide_na = next(i.identifier for i in hide.node_group.interface.items_tree
                   if i.name == 'cutoff_Na')
    atoms_before = evaluated(structure)
    set_mod_input(hide, hide_na, True)
    assert evaluated(structure) < atoms_before, 'hiding Na removed nothing'
    set_mod_input(hide, hide_na, False)

    # repeating the structure's supercell drives the faces along
    supercell = structure.modifiers[1]
    repeat_x = next(i.identifier for i in supercell.node_group.interface.items_tree
                    if i.name == 'repeat_x')
    set_mod_input(supercell, repeat_x, 2)
    # the supercell merges by distance, so polyhedra that coincide across
    # the cell boundary collapse into one: more faces, but no doubles
    import numpy as np
    depsgraph = bpy.context.evaluated_depsgraph_get()
    mesh = faces_obj.evaluated_get(depsgraph).to_mesh()
    co = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get('co', co)
    co = np.round(co.reshape(-1, 3), 3)
    polygons = [tuple(sorted(map(tuple, co[list(p.vertices)]))) for p in mesh.polygons]
    faces_obj.evaluated_get(depsgraph).to_mesh_clear()
    assert faces_before < len(polygons) <= 2 * faces_before, (len(polygons), faces_before)
    assert len(set(polygons)) == len(polygons), 'duplicate polyhedra faces'

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

def run_polyhedra_adps():
    """Thermal ellipsoids: the tensor math against the U_eq a refinement
    writes, the drawn ellipsoid size against a known tensor, and the adps
    switch and ring material in the node tree."""
    import numpy as np
    from blender_importASE.adp import atom_adps, ellipsoids, probability_scale
    from blender_importASE.polyhedra import import_polyhedra

    assert abs(probability_scale(0.5) - 1.5382) < 1e-4

    # every symmetry copy of a site keeps the site's U_eq (the trace is
    # invariant under rotation, the CIF's U_eq is not computed by us); the
    # benzoic acid is monoclinic in a nonstandard P21/n setting with riding
    # Uiso hydrogens, the neutron urea has anisotropic H on special positions
    for fixture, tolerance in (('benzoic_acid_cod7252065.cif', 1.5e-4),
                               ('urea_neutron_cod1008775.cif', None)):
        atoms = ase.io.read(f'{SCRATCH}/{fixture}', store_tags=True)
        U = atom_adps(atoms)
        kinds = atoms.arrays['spacegroup_kinds']
        _, axes, valid = ellipsoids(U)
        assert valid.all(), f'{fixture}: {(~valid).sum()} invalid tensors'
        for site, ueq in enumerate(atoms.info['_atom_site_u_iso_or_equiv']
                                   if tolerance else []):
            traces = np.trace(U[kinds == site], axis1=1, axis2=2) / 3
            assert np.allclose(traces, float(ueq), atol=tolerance), (fixture, site, traces, ueq)
        # symmetry copies share the eigenvalues
        for site in set(kinds):
            spread = np.ptp(axes[kinds == site], axis=0).max()
            assert spread < 1e-9, (fixture, site, spread)

    # one atom with an axis-aligned tensor: the evaluated mesh spans the
    # ellipsoid's semi-axes (+ the ring tube) at 50 % probability
    single = f'{SCRATCH}/adp_single.cif'
    with open(single, 'w') as fh:
        fh.write("data_single\n_cell_length_a 10\n_cell_length_b 10\n"
                 "_cell_length_c 10\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                 "_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'P 1'\n"
                 "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
                 "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
                 "C1 C 0.5 0.5 0.5\nloop_\n_atom_site_aniso_label\n"
                 "_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n"
                 "_atom_site_aniso_U_33\n_atom_site_aniso_U_12\n"
                 "_atom_site_aniso_U_13\n_atom_site_aniso_U_23\n"
                 "C1 0.01 0.04 0.09 0 0 0\n")
    fresh_scene()
    obj = import_polyhedra(single, 'adp_single.cif', adps=True, outline=False)
    from blender_importASE.controls import find_ase_modifier
    atoms_mod, _ = find_ase_modifier(obj)
    kinds = [m.node_group.name.split('_')[0].split('.')[0] for m in obj.modifiers]
    assert kinds == ['hide atoms', 'supercell', 'atoms', 'adp'], kinds
    ring = obj.modifiers['adp_rings'].node_group.interface.items_tree['ring_radius']
    depsgraph = bpy.context.evaluated_depsgraph_get()
    mesh = obj.evaluated_get(depsgraph).to_mesh()
    coords = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get('co', coords)
    coords = coords.reshape(-1, 3)
    extent = (coords.max(axis=0) - coords.min(axis=0)) / 2
    expected = 1.5382 * np.sqrt([0.01, 0.04, 0.09]) + ring.default_value
    assert np.allclose(extent, expected, rtol=0.02), (extent, expected)
    slots = [m.name for m in obj.data.materials]
    assert slots[-1] == 'adp_rings', slots
    # the materials the faces actually resolve to on the evaluated mesh (an
    # empty material slipping in at slot 0 shifts them all by one)
    face_slots = np.empty(len(mesh.polygons), dtype=np.int32)
    mesh.polygons.foreach_get('material_index', face_slots)
    used = {mesh.materials[i].name if mesh.materials[i] else None
            for i in set(face_slots.tolist())}
    assert used == {'C', 'adp_rings'}, used
    obj.evaluated_get(depsgraph).to_mesh_clear()

    # the adps switch brings back the plain sphere and drops the rings
    from blender_importASE.node_networks.compat import set_mod_input
    switch = atoms_mod.node_group.interface.items_tree['adps']
    set_mod_input(atoms_mod, switch.identifier, False)
    depsgraph = bpy.context.evaluated_depsgraph_get()
    mesh = obj.evaluated_get(depsgraph).to_mesh()
    coords = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get('co', coords)
    extent = np.ptp(coords.reshape(-1, 3), axis=0) / 2
    assert np.allclose(extent, extent[0], rtol=0.02), f'not a sphere: {extent}'
    obj.evaluated_get(depsgraph).to_mesh_clear()

    # the rings are a flat black emission
    ring_material = bpy.data.materials['adp_rings']
    ring_nodes = {n.bl_idname for n in ring_material.node_tree.nodes}
    assert ring_nodes == {'ShaderNodeEmission', 'ShaderNodeOutputMaterial'}, ring_nodes

    def ring_faces(obj):
        """(faces on the ring material, materials the faces resolve to)"""
        depsgraph = bpy.context.evaluated_depsgraph_get()
        mesh = obj.evaluated_get(depsgraph).to_mesh()
        face_slots = np.empty(len(mesh.polygons), dtype=np.int32)
        mesh.polygons.foreach_get('material_index', face_slots)
        names = [m.name if m else None for m in mesh.materials]
        used = {names[i] for i in set(face_slots.tolist())}
        count = int(sum(1 for i in face_slots if names[i] == 'adp_rings'))
        obj.evaluated_get(depsgraph).to_mesh_clear()
        return count, used

    # real structures import with the house-style outline, every face still
    # on its own element / bond / ring material; the rings pass the outline
    # as curves, so it neither shells nor doubles them
    for fixture in ('benzoic_acid_cod7252065.cif', 'urea_neutron_cod1008775.cif'):
        counts = {}
        for outline in (False, True):
            fresh_scene()
            obj = import_polyhedra(f'{SCRATCH}/{fixture}', fixture, adps=True, outline=outline)
            counts[outline], used = ring_faces(obj)
        symbols = set(ase.io.read(f'{SCRATCH}/{fixture}').get_chemical_symbols())
        assert symbols | {'adp_rings', 'outline_color'} <= used, used
        assert any(name and name.startswith('BOND') for name in used), used
        assert None not in used, used
        assert counts[True] == counts[False] > 0, counts

        # hydrogens are spheres by default; hydrogen_adps gives them rings
        atoms_mod, _ = find_ase_modifier(obj)
        switch = atoms_mod.node_group.interface.items_tree['hydrogen_adps']
        set_mod_input(atoms_mod, switch.identifier, True)
        with_h, _ = ring_faces(obj)
        assert with_h > counts[True] > 0, (with_h, counts[True])

    # SHELX: the .res SHELXL wrote into the benzoic acid CIF gives the same
    # structure and tensors (to the CIF's rounding), riding H included
    from blender_importASE.shelx import read_shelx
    from ase.spacegroup import crystal
    res = read_shelx(f'{SCRATCH}/benzoic_acid_cod7252065.res')
    cif = ase.io.read(f'{SCRATCH}/benzoic_acid_cod7252065.cif', store_tags=True)
    distance = np.linalg.norm(res.positions[:, None] - cif.positions[None], axis=2)
    match = distance.argmin(axis=1)
    assert len(res) == len(cif) == len(set(match.tolist())), (len(res), len(cif))
    assert distance.min(axis=1).max() < 2e-3
    assert res.get_chemical_symbols() == [cif.get_chemical_symbols()[i] for i in match]
    assert np.abs(atom_adps(res) - atom_adps(cif)[match]).max() < 1e-3

    # the less common parts of the format on a synthetic Cc (LATT -7):
    # long-form SFAC, a fixed coordinate (10 + p), continuation lines, a
    # label repeated in a second residue, a riding hydrogen
    with open(f'{SCRATCH}/adp_synthetic.ins', 'w') as fh:
        fh.write("TITL synthetic in Cc\nCELL 0.71073 7.1 8.3 9.2 90 101 90\n"
                 "ZERR 4 0.001 0.001 0.001 0 0.01 0\nLATT -7\nSYMM X, -Y, 1/2+Z\n"
                 "SFAC C 2.3100 20.8439 1.0200 10.2075 1.5886 0.5687 0.8650 51.6512 =\n"
                 "   0.2156 0.0033 0.0016 1.1500 0.7700 12.0110\nSFAC H\n"
                 "UNIT 8 4\nFVAR 1.0 0.6\nRESI 1 MOL\n"
                 "C1 1 0.1234 0.2345 10.34560 11.0 0.02 0.03 0.04 =\n"
                 "   0.002 -0.003 0.004\nH1 2 0.2 0.3 0.4 11.0 -1.2\n"
                 "RESI 2 MOL\nC1 1 0.6234 0.1345 0.1456 21.0 0.025 ! a comment\n"
                 "HKLF 4\nEND\nQ1 1 0.5 0.5 0.5 11.0 0.05 0.3\n")
    synthetic = read_shelx(f'{SCRATCH}/adp_synthetic.ins')
    reference = crystal(['C', 'H', 'C'], [[0.1234, 0.2345, 0.3456], [0.2, 0.3, 0.4],
                                          [0.6234, 0.1345, 0.1456]],
                        spacegroup=9, cellpar=(7.1, 8.3, 9.2, 90, 101, 90))
    assert len(synthetic) == len(reference) == 12, len(synthetic)
    assert synthetic.info['_atom_site_label'] == ['C1', 'H1', 'C1#2']
    U = atom_adps(synthetic)
    kinds = synthetic.arrays['spacegroup_kinds']
    ueq_c1 = np.trace(U[kinds == 0][0]) / 3
    assert np.isclose(np.trace(U[kinds == 1][0]) / 3, 1.2 * ueq_c1), 'riding H'
    assert np.isclose(np.trace(U[kinds == 2][0]) / 3, 0.025)
    assert ellipsoids(U)[2].all()

    # and the importer takes the .res directly
    fresh_scene()
    obj = import_polyhedra(f'{SCRATCH}/benzoic_acid_cod7252065.res',
                           'benzoic_acid_cod7252065.res', adps=True, outline=True)
    count, used = ring_faces(obj)
    assert count > 0 and {'C', 'H', 'O', 'adp_rings'} <= used, (count, used)

    # a file without displacement parameters imports as usual with adps on
    fresh_scene()
    obj = import_polyhedra(f'{SCRATCH}/nacl.extxyz', 'nacl.extxyz', adps=True)
    assert 'adp_rings' not in obj.modifiers, [m.name for m in obj.modifiers]

step('polyhedra_adps', run_polyhedra_adps)

def run_import_adps():
    """The regular importer: ellipsoids by default when the file carries an
    anisotropic table (nodes representation), nothing otherwise, and SHELX
    files read directly."""
    from blender_importASE.controls import find_ase_modifier

    def structure():
        return next(o for o in bpy.data.objects if o.type == 'MESH'
                    and find_ase_modifier(o)[0] is not None)

    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/crystal.cif', 'crystal.cif',
                        representation='nodes', animate=False, read_density=False)
    obj = structure()
    assert obj.modifiers[-1].name == 'adp_rings', [m.name for m in obj.modifiers]
    atoms_mod, names = find_ase_modifier(obj)
    assert {'adps', 'hydrogen_adps', 'adp_scale'} <= set(names), names
    assert obj.data.attributes['adp_valid'].data[0].value

    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/crystal.cif', 'crystal.cif', adps=False,
                        representation='nodes', animate=False, read_density=False)
    obj = structure()
    assert 'adp_rings' not in obj.modifiers and 'adp_valid' not in obj.data.attributes

    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/water.xyz', 'water.xyz',
                        representation='nodes', animate=False, read_density=False)
    assert 'adp_rings' not in structure().modifiers

    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/benzoic_acid_cod7252065.res', 'benzoic_acid_cod7252065.res',
                        representation='nodes', animate=False, read_density=False)
    obj = structure()
    assert len(obj.data.vertices) == 64, len(obj.data.vertices)
    assert obj.modifiers[-1].name == 'adp_rings', [m.name for m in obj.modifiers]

    # other representations draw objects of their own - ADPs do not apply
    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/crystal.cif', 'crystal.cif',
                        representation="Balls'n'Sticks", animate=False, read_density=False)
    assert not any(m.name == 'adp_rings' for o in bpy.data.objects for m in o.modifiers)

step('import_adps', run_import_adps)

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

def run_density_mesh_shells():
    """A nest of isosurfaces, each colored by the isovalue it stands for
    and more transparent the further out it sits."""
    from importlib import util
    if util.find_spec('skimage') is None:
        print('scikit-image not installed - skipping the actual import')
        return
    import numpy as np
    from blender_importASE.density_mesh import (import_density_mesh, shell_levels,
                                                COLOR_ATTRIBUTE)

    # geometric by default, so the levels stay apart on a density that
    # falls off exponentially
    levels = shell_levels(0.02, 4, limit=0.5)
    assert len(levels) == 4 and abs(levels[0] - 0.02) < 1e-9, levels
    assert abs(levels[-1] - 0.25) < 1e-9, levels          # half the data maximum
    ratios = [levels[i + 1] / levels[i] for i in range(3)]
    assert max(ratios) - min(ratios) < 1e-6, ratios
    assert shell_levels(0.02, 1, limit=0.5) == [0.02], shell_levels(0.02, 1, limit=0.5)
    linear = shell_levels(0.02, 4, spacing='LINEAR', limit=0.5)
    steps = [linear[i + 1] - linear[i] for i in range(3)]
    assert max(steps) - min(steps) < 1e-9, steps

    def shades(obj):
        attr = obj.data.color_attributes[COLOR_ATTRIBUTE]
        rgba = np.empty(len(attr.data) * 4)
        attr.data.foreach_get('color', rgba)
        rgba = rgba.reshape(-1, 4)
        return (sorted({round(float(v), 4) for v in rgba[:, 0]}),
                sorted({round(float(v), 4) for v in rgba[:, 3]}))

    fresh_scene()
    single = import_density_mesh(f'{SCRATCH}/mo.cube', 'mo.cube', iso_value=0.02,
                                 shells=1, import_atoms=False)
    colors, alphas = shades(single)
    assert (colors, alphas) == ([0.0, 1.0], [1.0]), (colors, alphas)
    assert single.data.materials[0].name == 'density_mesh material', \
        single.data.materials[0].name
    single_verts = len(single.data.vertices)

    fresh_scene()
    nest = import_density_mesh(f'{SCRATCH}/mo.cube', 'mo.cube', iso_value=0.02,
                               shells=3, import_atoms=False, layered=False)
    colors, alphas = shades(nest)
    # the color channel is the shell's own strength, so the ramp reads as
    # a color map: 0 outermost, 1 innermost, both signs on the same scale
    assert colors == [0.0, 0.5, 1.0], colors
    assert alphas == colors, (colors, alphas)
    assert len(nest.data.vertices) > single_verts, (len(nest.data.vertices), single_verts)

    # jet, and an Emission rather than a Principled BSDF: a contour map is
    # a color map, not a lit surface
    material = nest.data.materials[0]
    assert material.name == 'density_jet material', material.name
    kinds = {n.bl_idname for n in material.node_tree.nodes}
    assert 'ShaderNodeEmission' in kinds, sorted(kinds)
    assert 'ShaderNodeBsdfPrincipled' not in kinds, sorted(kinds)
    # an HSV sweep from blue to red taking the long way round the hue
    # circle, so it runs blue - cyan - green - yellow - red
    ramp = material.node_tree.nodes['Color Ramp'].color_ramp
    assert ramp.color_mode == 'HSV', ramp.color_mode
    assert ramp.hue_interpolation == 'FAR', ramp.hue_interpolation
    stops = [(round(e.position, 3), tuple(round(c, 3) for c in e.color[:3]))
             for e in ramp.elements]
    assert stops == [(0.1, (0.0, 0.0, 1.0)), (1.0, (1.0, 0.0, 0.0))], stops
    middle = tuple(round(c, 2) for c in ramp.evaluate(0.5)[:3])
    assert middle[1] > 0.5 and middle[0] < 0.5, middle      # green, not magenta

    # four nodes and nothing else: attribute, ramp, emission, output
    assert len(material.node_tree.nodes) == 4, [n.name for n in material.node_tree.nodes]
    output = next(n for n in material.node_tree.nodes
                  if n.bl_idname == 'ShaderNodeOutputMaterial')
    assert output.inputs['Surface'].links[0].from_node.bl_idname \
        == 'ShaderNodeEmission', 'the surface is not the emission'

    # every shell winds outwards - marching cubes winds by the gradient,
    # which flips between the two signs
    mesh = nest.data
    positions = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get('co', positions)
    positions = positions.reshape(-1, 3)
    attr = mesh.color_attributes[COLOR_ATTRIBUTE]
    rgba = np.empty(len(attr.data) * 4)
    attr.data.foreach_get('color', rgba)
    shade = rgba.reshape(-1, 4)[:, 0]
    for value in sorted({round(float(v), 4) for v in shade}):
        group = [p for p in mesh.polygons
                 if round(float(shade[p.vertices[0]]), 4) == value]
        centre = np.mean([positions[v] for p in group for v in p.vertices], axis=0)
        outward = np.mean([np.dot(np.array(p.normal), np.array(p.center) - centre)
                           for p in group])
        assert outward > 0, f'shell {value} winds inwards ({outward:.3f})'

    # a density with only one sign uses the whole ramp for its levels
    fresh_scene()
    positive = import_density_mesh(f'{SCRATCH}/water.cube', 'water.cube',
                                   iso_value=0.05, shells=4, import_atoms=False,
                                   layered=False)
    colors, alphas = shades(positive)
    assert len(colors) == 4 and colors[0] == 0.0 and colors[-1] == 1.0, colors
    assert alphas == colors, (colors, alphas)

step('density_mesh_shells', run_density_mesh_shells)

def run_density_shell_layers():
    """layered=True gives every shell its own render pass and composites
    them by value, with the structure on a pass of its own on top."""
    from importlib import util
    if util.find_spec('skimage') is None:
        print('scikit-image not installed - skipping the actual import')
        return
    from blender_importASE.density_mesh import import_density_mesh
    from blender_importASE.node_networks.compat import (compositor_tree,
                                                        alpha_over_sockets)

    fresh_scene()
    shells = import_density_mesh(f'{SCRATCH}/water.cube', 'water.cube',
                                 iso_value=0.03, shells=4, layered=True,
                                 import_atoms=True, outline=False)
    assert isinstance(shells, list) and len(shells) == 4, shells
    scene = bpy.context.scene
    base = scene.view_layers[0]
    structure_layer = next(vl for vl in scene.view_layers
                           if vl.name.endswith('_structure'))

    def visible(view_layer):
        return {c.collection.name for c in view_layer.layer_collection.children
                if not c.exclude}

    # base + one per shell + the structure's own
    assert len(scene.view_layers) == 6, [vl.name for vl in scene.view_layers]
    shell_collections = {f'{obj.name}_layer' for obj in shells}
    for obj in shells:
        layer = scene.view_layers[f'{obj.name}_layer']
        seen = visible(layer)
        assert f'{obj.name}_layer' in seen, (layer.name, seen)
        assert not (shell_collections - {f'{obj.name}_layer'}) & seen, (layer.name, seen)
        # the structure must not occlude the bands it is composited over
        assert not any(name.startswith('H2O') for name in seen), (layer.name, seen)
    assert any(name.startswith('H2O') for name in visible(structure_layer)), \
        visible(structure_layer)
    assert not shell_collections & visible(structure_layer), visible(structure_layer)
    # the base layer stays the viewport's: it shows everything, and is kept
    # out of the render instead of being emptied (else the structure would
    # vanish from the viewport after a layered import)
    assert shell_collections <= visible(base), visible(base)
    assert any(name.startswith('H2O') for name in visible(base)), visible(base)
    assert not base.use, 'the base layer must not render on its own'

    # composited bottom to top: backdrop, weakest .. strongest, structure
    tree = compositor_tree(scene)
    assert tree is not None, 'no compositor'
    assert scene.render.film_transparent, 'the layers cannot be alpha-overed'
    # 5.x composites through a group output, 4.x through a Composite node
    composite = next(n for n in tree.nodes
                     if n.bl_idname in ('CompositorNodeComposite', 'NodeGroupOutput'))
    order, node = [], composite.inputs[0].links[0].from_node
    while node.bl_idname == 'CompositorNodeAlphaOver':
        background, foreground, _factor = alpha_over_sockets(node)
        order.append(foreground.links[0].from_node.layer)
        node = background.links[0].from_node
    order.append(node.layer)
    # every renderable object here belongs to a pass already, so there is
    # no separate backdrop pass and the weakest shell is the bottom
    expected = ([structure_layer.name]
                + [f'{obj.name}_layer' for obj in reversed(shells)])
    assert order == expected, (order, expected)

step('density_shell_layers', run_density_shell_layers)

def run_global_supercell():
    """One supercell for the whole structure: the atoms through their
    modifier, the density volumes by tiling the grid, and the isosurface
    meshes by recomputing them on a tiled grid."""
    from importlib import util
    if util.find_spec('skimage') is None or (
            util.find_spec('openvdb') is None and util.find_spec('pyopenvdb') is None):
        print('scikit-image or openvdb missing - skipping')
        return
    import numpy as np
    from blender_importASE import controls
    from blender_importASE.density_mesh import import_density_mesh
    from blender_importASE.node_networks.compat import get_mod_input

    fresh_scene()
    import_ase_molecule(f'{SCRATCH}/CHGCAR', 'CHGCAR', representation='nodes',
                        animate=False, read_density=True, outline=False)
    # the mesh importer has to read a CHGCAR too: ase's VaspChargeDensity
    # returns no grids at all for this file, which is what the add-on's own
    # read_vasp_density works around
    isomesh = import_density_mesh(f'{SCRATCH}/CHGCAR', 'CHGCAR', iso_value=1.0,
                                  shells=1, import_atoms=False)
    structure = next(o for o in bpy.data.objects
                     if o.type == 'MESH' and 'atom_color' in o.data.attributes)
    bpy.context.view_layer.objects.active = structure

    modifiers, volumes, meshes = controls._supercell_parts(structure)
    assert modifiers and volumes and meshes, (len(modifiers), len(volumes), len(meshes))

    def span(obj):
        co = np.array([v.co[:] for v in obj.data.vertices])
        return co.max(axis=0) - co.min(axis=0)

    before = span(isomesh)
    assert bpy.ops.ase.global_supercell.poll(), 'not offered on the structure'
    bpy.ops.ase.global_supercell(repeat_x=2, repeat_y=2, repeat_z=1)

    # the atoms repeat through the node group ...
    supercell = next(m for m in structure.modifiers
                     if m.node_group and m.node_group.name.startswith('supercell'))
    names = {i.name: i.identifier for i in supercell.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    assert [get_mod_input(supercell, names[axis])
            for axis in ('repeat_x', 'repeat_y', 'repeat_z')] == [2, 2, 1]
    # ... the volumes by a tiled grid ...
    for volume in volumes:
        assert list(volume['ase_density_repeat']) == [2, 2, 1], volume.name
        assert '2x2x1' in volume.data.filepath, volume.data.filepath
    # ... and the isosurface by being recomputed, not repeated
    after = span(isomesh)
    assert after[0] > 1.8 * before[0] and after[1] > 1.8 * before[1], (before, after)
    assert abs(after[2] - before[2]) < 0.1, (before, after)
    assert list(isomesh['ase_density_repeat']) == [2, 2, 1], isomesh['ase_density_repeat'][:]

step('global_supercell', run_global_supercell)

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

def run_density_supercell():
    """The density volume must be able to span a supercell, and its
    isosurfaces must take their material from the object's slots."""
    from importlib import util
    if util.find_spec('openvdb') is None and util.find_spec('pyopenvdb') is None:
        print('openvdb not installed - skipping')
        return
    try:
        import openvdb as vdb
    except ImportError:
        import pyopenvdb as vdb
    from blender_importASE.import_cubefiles import density_supercell

    def grid_dims(volume_obj):
        grid = vdb.read(bpy.path.abspath(volume_obj.data.filepath), 'density')
        box = grid.evalActiveVoxelBoundingBox()
        import numpy as np
        return tuple(int(n) for n in np.asarray(box[1]) - np.asarray(box[0]) + 1)

    run_import(f'{SCRATCH}/CHGCAR', representation='nodes', animate=False,
               read_density=True)
    vol = next(o for o in bpy.data.objects
               if o.type == 'VOLUME' and 'spin' not in o.name)

    # materials come from the object's own slots now, not from a modifier
    # input: slot 0 is the positive lobe, slot 1 the negative one. The slots
    # must be OBJECT-linked: a material index resolves against the material
    # list the geometry carries, and the isosurface Volume to Mesh builds
    # carries none, so data-linked slots render plain white however right
    # the indices are.
    slots = [(sl.link, sl.material.name if sl.material else None)
             for sl in vol.material_slots]
    assert slots == [('OBJECT', '+ material'), ('OBJECT', '- material')], slots
    group = vol.modifiers[0].node_group
    inputs = [i for i in group.interface.items_tree
              if getattr(i, 'in_out', None) == 'INPUT']
    assert not any(i.socket_type == 'NodeSocketMaterial' for i in inputs), \
        [i.name for i in inputs]
    # both halves are needed: Set Material gives the geometry a material list
    # (Cycles clamps the index to it), Set Material Index makes that index
    # resolve against the object's slots
    assert len([n for n in group.nodes
                if n.bl_idname == 'GeometryNodeSetMaterial']) == 2, 'no Set Material per sign'
    assert any(n.bl_idname == 'GeometryNodeSetMaterialIndex' for n in group.nodes), \
        'no Set Material Index in the density group'

    base = grid_dims(vol)
    assert list(vol['ase_grid_shape']) == list(base), (vol['ase_grid_shape'], base)
    density_supercell(vol, (2, 2, 1))
    assert grid_dims(vol) == (base[0] * 2, base[1] * 2, base[2]), grid_dims(vol)
    assert list(vol['ase_density_repeat']) == [2, 2, 1], vol['ase_density_repeat'][:]
    # always tiled from the single-cell grid, never compounded
    density_supercell(vol, (3, 1, 1))
    assert grid_dims(vol) == (base[0] * 3, base[1], base[2]), grid_dims(vol)
    density_supercell(vol, (1, 1, 1))
    assert grid_dims(vol) == base, grid_dims(vol)
    assert vol.data.filepath == vol['ase_base_vdb'], vol.data.filepath

    # the offset is the supercell group's Offset_x/y/z for densities, and
    # unlike the repeat it stays a live modifier input: a plain translation
    # of the isosurface along the lattice vectors, nothing rewritten
    from blender_importASE.node_networks.compat import set_mod_input, get_mod_input
    mod = vol.modifiers[0]
    names = {i.name: i.identifier for i in mod.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    for axis in 'abc':
        assert f'offset {axis}' in names and f'cell {axis}' in names, sorted(names)
    cell = [list(get_mod_input(mod, names[f'cell {axis}'])) for axis in 'abc']
    assert all(any(abs(v) > 1e-6 for v in vec) for vec in cell), cell
    transform = [n for n in mod.node_group.nodes
                 if n.bl_idname == 'GeometryNodeTransform']
    assert transform and transform[0].inputs['Translation'].is_linked, \
        'the offset does not reach a Transform node'
    before = vol.data.filepath
    set_mod_input(mod, names['offset a'], 1)
    assert vol.data.filepath == before, 'the offset rewrote the grid'
    set_mod_input(mod, names['offset a'], 0)

    # ... and the operator does it for every density of the structure
    structure = next(o for o in bpy.data.objects
                     if o.type == 'MESH' and 'atom_color' in o.data.attributes)
    bpy.context.view_layer.objects.active = structure
    assert bpy.ops.ase.density_supercell.poll(), 'operator not available on the structure'
    bpy.ops.ase.density_supercell(repeat_x=2, repeat_y=1, repeat_z=1, offset_x=-1)
    for volume_obj in [o for o in bpy.data.objects if o.type == 'VOLUME']:
        assert list(volume_obj['ase_density_repeat']) == [2, 1, 1], volume_obj.name
        vmod = volume_obj.modifiers[0]
        vnames = {i.name: i.identifier for i in vmod.node_group.interface.items_tree
                  if getattr(i, 'in_out', None) == 'INPUT'}
        assert get_mod_input(vmod, vnames['offset a']) == -1, volume_obj.name

step('density_supercell', run_density_supercell)

def run_density_node_upgrade():
    """A density imported by an older add-on keeps that version's node
    group - a modifier never swaps group by itself - so the offsets are
    missing until it is upgraded. That is what 'Update density nodes' is
    for, and it must carry the settings and the cell vectors across."""
    from importlib import util
    if util.find_spec('openvdb') is None and util.find_spec('pyopenvdb') is None:
        print('openvdb not installed - skipping')
        return
    from blender_importASE.node_networks.compat import get_mod_input, set_mod_input
    from blender_importASE import controls

    run_import(f'{SCRATCH}/CHGCAR', representation='nodes', animate=False,
               read_density=True)
    vol = next(o for o in bpy.data.objects
               if o.type == 'VOLUME' and 'spin' not in o.name)
    mod = vol.modifiers[0]
    names = {i.name: i.identifier for i in mod.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    set_mod_input(mod, names['isovalue'], 0.125)

    # turn its group back into an older one: no stamp, no offsets
    group = mod.node_group
    group.description = ''
    for item in list(group.interface.items_tree):
        if item.name.startswith(('offset ', 'cell ')):
            group.interface.remove(item)
    assert controls._density_nodes_outdated(vol), 'outdated group not detected'

    bpy.context.view_layer.objects.active = vol
    assert bpy.ops.ase.upgrade_density_nodes.poll(), 'upgrade not offered'
    bpy.ops.ase.upgrade_density_nodes()
    assert not controls._density_nodes_outdated(vol), 'still outdated after the upgrade'

    mod = vol.modifiers[0]
    names = {i.name: i.identifier for i in mod.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    assert abs(get_mod_input(mod, names['isovalue']) - 0.125) < 1e-6, \
        'the isovalue was lost in the upgrade'
    cell = [list(get_mod_input(mod, names[f'cell {axis}'])) for axis in 'abc']
    assert all(any(abs(v) > 1e-6 for v in vec) for vec in cell), cell
    slots = [(sl.link, sl.material.name if sl.material else None)
             for sl in vol.material_slots]
    assert slots == [('OBJECT', '+ material'), ('OBJECT', '- material')], slots

step('density_node_upgrade', run_density_node_upgrade)

def run_density_material_render():
    """Both isosurfaces must actually render in their own color, and the
    object's material slots must drive them - in Cycles, not just EEVEE.

    Renders with Cycles on the CPU on purpose: this is the engine that
    exposed the bug (it clamps a face's material index to the material list
    the geometry itself carries, where EEVEE falls back to the object's
    slots), and it needs no GL context in a headless run.
    """
    from importlib import util
    if util.find_spec('openvdb') is None and util.find_spec('pyopenvdb') is None:
        print('openvdb not installed - skipping')
        return
    import numpy as np
    from blender_importASE.node_networks.compat import set_mod_input

    run_import(f'{SCRATCH}/mo.cube', representation='nodes', animate=False,
               read_density=True, outline=False)
    vol = next(o for o in bpy.data.objects if o.type == 'VOLUME')
    for ob in bpy.data.objects:
        ob.hide_render = (ob is not vol)
    mod = vol.modifiers[0]
    names = {i.name: i.identifier for i in mod.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    set_mod_input(mod, names['isovalue'], 0.02)

    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.samples = 8
    scene.cycles.device = 'CPU'
    scene.render.resolution_x = scene.render.resolution_y = 120
    scene.render.film_transparent = True
    scene.world = bpy.data.worlds.new('render world')
    scene.world.use_nodes = True
    scene.world.node_tree.nodes['Background'].inputs[1].default_value = 1.0
    camera_data = bpy.data.cameras.new('render cam')
    camera_data.type = 'ORTHO'
    camera_data.ortho_scale = 12
    camera = bpy.data.objects.new('render cam', camera_data)
    scene.collection.objects.link(camera)
    scene.camera = camera
    camera.location = (4, 4, 30)

    def channels():
        scene.render.filepath = f'{SCRATCH}/density_material_render.png'
        bpy.ops.render.render(write_still=True)
        img = bpy.data.images.load(scene.render.filepath)
        px = np.array(img.pixels[:]).reshape(img.size[1], img.size[0], 4)
        lit = px[..., :3][px[..., 3] > 0.5]
        bpy.data.images.remove(img)
        assert len(lit), 'the density did not render at all'
        # pixels whose brightest channel is red / green / blue
        return [int((lit.argmax(axis=1) == channel).sum()) for channel in range(3)]

    red, green, blue = channels()
    assert red > 100 and blue > 100, \
        f'the +/- lobes did not both render in color: R{red} G{green} B{blue}'

    # the Material Properties tab has to drive it: slot 0 is the + lobe
    from blender_importASE.node_networks.electron_density_nodes import newShader
    lime = newShader('lime test', 0, 1, 0)
    assert vol.material_slots[0].link == 'OBJECT', vol.material_slots[0].link
    vol.material_slots[0].material = lime
    red, green, blue = channels()
    assert green > 100 and blue < 50, \
        f'swapping slot 0 did not repaint the + lobe: R{red} G{green} B{blue}'

step('density_material_render', run_density_material_render)

def run_density_cutoffs():
    """Every cutoff is a depth in from its own face of the density's
    bounding box: 0 cuts nothing wherever the density sits, and raising
    one eats into that side. Checked by rendering, because only the
    picture says which face actually lost geometry."""
    from importlib import util
    if util.find_spec('openvdb') is None and util.find_spec('pyopenvdb') is None:
        print('openvdb not installed - skipping')
        return
    import numpy as np
    from blender_importASE.node_networks.compat import set_mod_input

    run_import(f'{SCRATCH}/mo.cube', representation='nodes', animate=False,
               read_density=True, outline=False)
    vol = next(o for o in bpy.data.objects if o.type == 'VOLUME')
    for ob in bpy.data.objects:
        ob.hide_render = (ob is not vol)
    mod = vol.modifiers[0]
    cutoffs = [i for i in mod.node_group.interface.items_tree
               if i.name.startswith('cutoff')]
    assert len(cutoffs) == 6, [i.name for i in cutoffs]
    for socket in cutoffs:
        assert (socket.default_value, socket.min_value, socket.max_value) == (0.0, 0.0, 100.0), \
            (socket.name, socket.default_value, socket.min_value, socket.max_value)
        assert socket.subtype == 'DISTANCE', (socket.name, socket.subtype)

    names = {i.name: i.identifier for i in mod.node_group.interface.items_tree
             if getattr(i, 'in_out', None) == 'INPUT'}
    set_mod_input(mod, names['isovalue'], 0.02)
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.samples = 4
    scene.cycles.device = 'CPU'
    scene.render.resolution_x = scene.render.resolution_y = 150
    scene.render.film_transparent = True
    camera_data = bpy.data.cameras.new('cut cam')
    camera_data.type = 'ORTHO'
    camera_data.ortho_scale = 12
    camera = bpy.data.objects.new('cut cam', camera_data)
    scene.collection.objects.link(camera)
    scene.camera = camera
    camera.location = (4, 4, 30)

    def extent():
        scene.render.filepath = f'{SCRATCH}/density_cut.png'
        bpy.ops.render.render(write_still=True)
        img = bpy.data.images.load(scene.render.filepath)
        px = np.array(img.pixels[:]).reshape(img.size[1], img.size[0], 4)
        columns = np.nonzero((px[..., 3] > 0.5).any(axis=0))[0]
        bpy.data.images.remove(img)
        assert len(columns), 'the density did not render'
        scale = camera_data.ortho_scale / px.shape[1]
        return columns.min() * scale, columns.max() * scale

    left, right = extent()
    set_mod_input(mod, names['cut'], True)
    assert extent() == (left, right), 'cut with every cutoff at 0 removed geometry'
    set_mod_input(mod, names['cutoff X'], 2.0)
    cut_left, cut_right = extent()
    assert cut_left - left > 1.5 and abs(cut_right - right) < 0.2, \
        f'cutoff X did not eat in from the low-x face: {left, right} -> {cut_left, cut_right}'
    set_mod_input(mod, names['cutoff X'], 0.0)
    set_mod_input(mod, names['cutoff -X'], 2.0)
    cut_left, cut_right = extent()
    assert right - cut_right > 1.5 and abs(cut_left - left) < 0.2, \
        f'cutoff -X did not eat in from the high-x face: {left, right} -> {cut_left, cut_right}'

step('density_cutoffs', run_density_cutoffs)

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
        if os.path.splitext(fname)[1].lower() not in ('.xyz', '.extxyz', '.cif', '.res', '.ins'):
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
    assert default_element_color('Pb') == (0.0, 0.4, 0.2, 1.0), default_element_color('Pb')

    fresh_scene()
    ase.io.write(f'{SCRATCH}/lead.xyz', ase.Atoms('Pb2', positions=[(0, 0, 0), (3.5, 0, 0)]))
    import_ase_molecule(f'{SCRATCH}/lead.xyz', 'lead.xyz', representation='nodes',
                        animate=False, read_density=False)
    bsdf = next(n for n in bpy.data.materials['Pb'].node_tree.nodes
                if n.type == 'BSDF_PRINCIPLED')
    assert tuple(round(v, 4) for v in bsdf.inputs['Base Color'].default_value) \
        == (0.0, 0.4, 0.2, 1.0), tuple(bsdf.inputs['Base Color'].default_value)
    assert bsdf.inputs['Metallic'].default_value == 1.0, bsdf.inputs['Metallic'].default_value
    assert bsdf.inputs['Roughness'].default_value == 0.5, bsdf.inputs['Roughness'].default_value

step('lead_defaults', run_lead_defaults)
step('element_colors', run_element_colors)
step('operator_via_ops', lambda: (
    fresh_scene(),
    bpy.ops.import_mesh.ase(directory=SCRATCH, files=[{"name": "crystal.cif"}]),
))


def run_render_animations_plan():
    """The 'Render multiple animations' machinery, without actually rendering.

    Covers the two things that are easy to get wrong when a scene holds several
    trajectories: the per-collection frame range must come from that
    collection's own shape keys (the importer rewrites scene.frame_end on every
    import, so the scene range belongs to the last one imported), and restoring
    a collection must not un-hide the importer's hidden helper meshes.
    """
    from blender_importASE import render_vpts as rv

    fresh_scene()
    bpy.ops.object.camera_add(location=(0, 0, 20))
    scene = bpy.context.scene

    # two trajectories of different length in one scene
    traj = ase.io.read(f'{SCRATCH}/traj.xyz', index=':')
    ase.io.write(f'{SCRATCH}/traj_short.xyz', traj[:3])
    for fname in ('traj.xyz', 'traj_short.xyz'):
        import_ase_molecule(f'{SCRATCH}/{fname}', fname, representation='nodes',
                            animate=True, read_density=False)

    by_name = {c.name: c for c in scene.collection.children}
    assert len(by_name) == 2, sorted(by_name)
    ranges = {n: rv.collection_frame_range(c) for n, c in by_name.items()}
    lengths = sorted((hi - lo + 1) for lo, hi in ranges.values())
    assert lengths == [3, 5], f'per-collection ranges wrong: {ranges}'

    cams = rv.scene_cameras(scene)
    assert len(cams) == 1, cams
    jobs = rv.build_jobs(scene, SCRATCH, cams)
    assert len(jobs) == 8, f'expected 3+5 jobs, got {len(jobs)}'
    # single camera -> no <camera>_ prefix, and one subfolder per collection
    assert jobs[0][3].name == '0000.png', jobs[0][3].name
    assert jobs[0][3].parent.name in by_name, jobs[0][3].parent

    # ... and two cameras -> prefixed names, twice the jobs
    bpy.ops.object.camera_add(location=(20, 0, 0))
    cams2 = rv.scene_cameras(scene)
    jobs2 = rv.build_jobs(scene, SCRATCH, cams2)
    assert len(jobs2) == 16, len(jobs2)
    assert jobs2[0][3].name.startswith(cams2[0].name + '_'), jobs2[0][3].name

    # stride and an explicit range
    assert len(rv.build_jobs(scene, SCRATCH, cams, stride=2)) == 2 + 3
    assert len(rv.build_jobs(scene, SCRATCH, cams, start=0, end=1)) == 4

    # hide/restore must reproduce the baseline exactly, helper meshes included
    base = rv.visibility_snapshot(scene)
    hidden = {ob.name for c in scene.collection.children
              for ob in c.objects if ob.hide_render}
    assert hidden, 'expected the importer to leave helper meshes hidden'
    for coll in scene.collection.children:
        rv.set_collection_visible(coll, base[coll.name], False)
    assert all(ob.hide_render for c in scene.collection.children for ob in c.objects)
    for coll in scene.collection.children:
        rv.set_collection_visible(coll, base[coll.name], True)
    after = {ob.name for c in scene.collection.children
             for ob in c.objects if ob.hide_render}
    assert after == hidden, f'visibility not restored: {after} != {hidden}'


step('render_animations_plan', run_render_animations_plan)
step('render_animations_registered', lambda: (
    None if hasattr(bpy.ops.render, 'render_animations')
    else (_ for _ in ()).throw(AssertionError('render.render_animations missing'))
))
step('unregister', blender_importASE.unregister)

print('### SUMMARY')
for k, v in results.items():
    print(f"### {v:4s} {k}")
if any(v == 'FAIL' for v in results.values()):
    sys.exit(1)
