"""Import structures with coordination polyhedra.

For every atom with at least ``min_neighbors`` neighbors, the convex hull
of its neighbor shell is added as mesh faces, so coordination polyhedra
(octahedra, tetrahedra, ...) render as solid faces while the usual
atoms_and_bonds node group draws the atoms and bonds of the same mesh.

Ported from blender_polyhedra.ipynb; the cutoff multipliers it used are
exposed as options.
"""
import bpy
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.data import covalent_radii

from .utils import atomcolors
from .node_networks.nodes_atoms_and_bonds import set_atoms_node_group, atoms_and_bonds, read_structure
from .node_networks.bond_mat import create_bondmat
from .node_networks.electron_density_nodes import newShader
from .node_networks.outline import outline_objects

POLYHEDRA_MATERIAL = 'polyhedra material'


def build_polyhedra_atoms(atoms, expand_cutoff=1.2, trim_cutoff=1.0,
                          poly_cutoff=1.1, min_neighbors=4,
                          include_hydrogen=False):
    """Expand periodic neighbor images and compute convex-hull faces
    around every coordination center.

    Returns (new_atoms, faces): a non-periodic Atoms object whose
    positions are the mesh vertices, and the polyhedra faces as vertex
    index lists.

    expand_cutoff:  covalent-radius multiplier used to pull in periodic
                    neighbor images so boundary polyhedra are closed
    trim_cutoff:    multiplier for removing atoms left without neighbors
                    after the expansion
    poly_cutoff:    multiplier defining the neighbor shell that forms a
                    polyhedron
    min_neighbors:  minimum shell size; smaller shells get no polyhedron
    include_hydrogen: also use H as polyhedra centers/corners
    """
    from scipy.spatial import ConvexHull  # deferred: only this importer needs scipy

    positions = atoms.get_positions()
    cell = atoms.get_cell()

    # pull in periodic images of all neighbors, deduplicated per (atom, image)
    nl = NeighborList([covalent_radii[num] * expand_cutoff for num in atoms.numbers],
                      self_interaction=False, bothways=True)
    nl.update(atoms)
    new_atoms = Atoms()
    new_atoms.pbc = False
    seen = set()

    def append_image(index, offset):
        key = (index, tuple(int(o) for o in offset))
        if key in seen:
            return
        seen.add(key)
        new_atoms.append(atoms[index])
        new_atoms[-1].position = positions[index] + np.dot(offset, cell)

    for i in range(len(atoms)):
        append_image(i, (0, 0, 0))
        indices, offsets = nl.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            if j >= i:
                append_image(j, offset)

    # drop expanded atoms that ended up without any neighbor
    nl = NeighborList([covalent_radii[num] * trim_cutoff for num in new_atoms.numbers],
                      self_interaction=False, bothways=True)
    nl.update(new_atoms)
    for i in reversed(range(len(new_atoms))):
        if len(nl.get_neighbors(i)[0]) == 0:
            del new_atoms[i]

    # convex hull of each coordination shell
    nl = NeighborList([covalent_radii[num] * poly_cutoff for num in new_atoms.numbers],
                      self_interaction=False, bothways=True)
    nl.update(new_atoms)
    positions = new_atoms.get_positions()
    faces = []
    for i in range(len(new_atoms)):
        if not include_hydrogen and new_atoms[i].symbol == 'H':
            continue
        indices, offsets = nl.get_neighbors(i)
        neighbor_indices = [j for j in indices
                            if include_hydrogen or new_atoms[j].symbol != 'H']
        if len(neighbor_indices) < min_neighbors:
            continue
        neighbor_indices = np.array(neighbor_indices)
        neighbor_positions = positions[neighbor_indices]
        try:
            hull = ConvexHull(neighbor_positions - positions[i])
        except Exception:
            # degenerate shells (planar/linear) have no 3d hull
            continue
        for simplex in hull.simplices:
            faces.append([int(v) for v in neighbor_indices[simplex]])
    return new_atoms, faces


def import_polyhedra(filepath, filename, expand_cutoff=1.2, trim_cutoff=1.0,
                     poly_cutoff=1.1, min_neighbors=4, include_hydrogen=False,
                     resolution=16, colorbonds=True, bond_distance=0.66,
                     bond_radius=0.1, outline=False, **kwargs):
    import ase.io
    atoms = ase.io.read(filepath)

    new_atoms, faces = build_polyhedra_atoms(
        atoms, expand_cutoff=expand_cutoff, trim_cutoff=trim_cutoff,
        poly_cutoff=poly_cutoff, min_neighbors=min_neighbors,
        include_hydrogen=include_hydrogen)
    print(f'polyhedra: {len(new_atoms)} atoms, {len(faces)} faces')

    atomcolor = atomcolors()
    atomcolor.setup_materials(atoms, colorbonds=colorbonds)
    my_coll = bpy.data.collections.new(
        name=atoms.get_chemical_formula() + '_polyhedra_' + filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection

    obj, mesh = read_structure(new_atoms,
                               atoms.get_chemical_formula() + '_polyhedra_' + filename.split('.')[0],
                               animate=False, faces=faces)

    set_atoms_node_group()
    elements_name = '_'.join(list(set(new_atoms.get_chemical_symbols())))
    bondmat = create_bondmat(colorbonds=colorbonds, name=elements_name)
    atoms_from_verts = atoms_and_bonds(obj, new_atoms, 'GeometryNodes', bondmat=bondmat)
    obj.modifiers['GeometryNodes'].node_group = atoms_from_verts
    obj.modifiers['GeometryNodes']["Socket_2"] = bond_distance
    obj.modifiers['GeometryNodes']["Socket_3"] = bond_radius
    obj.modifiers['GeometryNodes']["Socket_4"] = resolution

    # the polyhedra faces get their own (semi-transparent) material slot,
    # appended after the element and bond materials; mat_slot is picked up
    # by the Set Material Index node at the end of the atoms_and_bonds tree
    poly_mat = newShader(POLYHEDRA_MATERIAL, 0.4, 0.55, 0.85)
    obj.data.materials.append(poly_mat)
    poly_slot = len(obj.data.materials) - 1
    if 'mat_slot' not in mesh.attributes:
        mesh.attributes.new(name='mat_slot', type='INT', domain='FACE')
    mesh.attributes['mat_slot'].data.foreach_set('value', [poly_slot] * len(mesh.polygons))
    mesh.update()

    if outline:
        outline_objects([obj], modifier='GeometryNodes.001')
    return obj
