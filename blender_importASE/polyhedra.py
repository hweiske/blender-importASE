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
from .node_networks.electron_density_nodes import newMaterial
from .node_networks.outline import outline_objects

POLYHEDRA_MATERIAL = 'polyhedra material'


def _polyhedra_material():
    """Semi-transparent material reading the per-vertex 'atom_color' attribute,
    so every polyhedron is tinted by the element colors of its corner
    atoms (e.g. brown SbBr6 octahedra from the Br corners)."""
    mat = newMaterial(POLYHEDRA_MATERIAL)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    principled = nodes.get('Principled BSDF')
    principled.inputs[1].default_value = 0.0   # metallic
    principled.inputs[2].default_value = 0.6   # roughness
    principled.inputs[3].default_value = 1.45  # IOR
    principled.inputs[4].default_value = 0.6   # alpha
    attr = nodes.get('Attribute')
    if attr is None:
        attr = nodes.new('ShaderNodeAttribute')
        attr.name = 'Attribute'
        attr.location = (-300, 200)
    attr.attribute_name = 'atom_color'
    attr.attribute_type = 'GEOMETRY'
    links.new(attr.outputs['Color'], principled.inputs['Base Color'])
    return mat


def build_polyhedra_atoms(atoms, expand_cutoff=1.2, trim_cutoff=1.0,
                          poly_cutoff=1.1, min_neighbors=4,
                          include_hydrogen=False, single_element_corners=True):
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
    single_element_corners: restrict each polyhedron to corner atoms of a
                    single element (the coordinating counter-ion). Drops
                    same-element neighbors -- e.g. the next-nearest Na-Na
                    contacts that otherwise bloat a Na-centered hull -- so
                    NaCl renders as clean NaCl6 / ClNa6 octahedra. Off by
                    default so same-element clusters (e.g. B6) are unaffected.
    """
    from collections import Counter
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
        if single_element_corners:
            # keep only counter-ion corners: drop neighbors of the center's own
            # element (e.g. next-nearest Na-Na), and if several other elements
            # remain, keep the most common one so every corner is one element.
            center_sym = new_atoms[i].symbol
            corner = [j for j in neighbor_indices
                      if new_atoms[j].symbol != center_sym]
            if corner:
                top = Counter(new_atoms[j].symbol for j in corner).most_common(1)[0][0]
                corner = [j for j in corner if new_atoms[j].symbol == top]
            neighbor_indices = corner
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
                     bond_radius=0.1, outline=False, single_element_corners=True,
                     **kwargs):
    import ase.io
    atoms = ase.io.read(filepath)

    new_atoms, faces = build_polyhedra_atoms(
        atoms, expand_cutoff=expand_cutoff, trim_cutoff=trim_cutoff,
        poly_cutoff=poly_cutoff, min_neighbors=min_neighbors,
        include_hydrogen=include_hydrogen,
        single_element_corners=single_element_corners)
    print(f'polyhedra: {len(new_atoms)} atoms, {len(faces)} faces')

    atomcolor = atomcolors()
    atomcolor.setup_materials(atoms, colorbonds=colorbonds)
    my_coll = bpy.data.collections.new(
        name=atoms.get_chemical_formula() + '_polyhedra_' + filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection

    name = atoms.get_chemical_formula() + '_polyhedra_' + filename.split('.')[0]
    obj, mesh = read_structure(new_atoms, name, animate=False)

    set_atoms_node_group()
    elements_name = '_'.join(list(set(new_atoms.get_chemical_symbols())))
    bondmat = create_bondmat(colorbonds=colorbonds, name=elements_name)
    atoms_from_verts = atoms_and_bonds(obj, new_atoms, 'GeometryNodes', bondmat=bondmat)
    obj.modifiers['GeometryNodes'].node_group = atoms_from_verts
    obj.modifiers['GeometryNodes']["Socket_2"] = bond_distance
    obj.modifiers['GeometryNodes']["Socket_3"] = bond_radius
    obj.modifiers['GeometryNodes']["Socket_4"] = resolution

    # the polyhedra faces live in their own object, so modifiers on the
    # structure (like the outline) never touch them. The material reads
    # the per-vertex 'atom_color' attribute, tinting each face by the
    # element colors of its corner atoms.
    poly_mesh = bpy.data.meshes.new(name + '_faces')
    poly_mesh.from_pydata(new_atoms.get_positions(), [], faces)
    color_attr = poly_mesh.attributes.new(name='atom_color', type='FLOAT_COLOR',
                                          domain='POINT')
    colors = np.empty(len(mesh.vertices) * 4)
    mesh.attributes['atom_color'].data.foreach_get('color', colors)
    color_attr.data.foreach_set('color', colors)
    poly_mesh.materials.append(_polyhedra_material())
    poly_mesh.update()
    poly_obj = bpy.data.objects.new(poly_mesh.name, poly_mesh)
    my_coll.objects.link(poly_obj)

    if outline:
        # atoms and bonds only - the polyhedra object stays outline-free
        outline_objects([obj], modifier='GeometryNodes.001')
    return obj
