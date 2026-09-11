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
from .drawobjects import draw_unit_cell
from .node_networks.nodes_atoms_and_bonds import set_atoms_node_group, atoms_and_bonds, read_structure
from .node_networks.bond_mat import create_bondmat
from .node_networks.electron_density_nodes import newMaterial
from .node_networks.outline import outline_objects
from .node_networks.compat import set_mod_input

POLYHEDRA_MATERIAL = 'polyhedra material'


def _polyhedra_material():
    """Material of the polyhedra faces: a Principled BSDF and a Glass BSDF
    mixed half and half, both tinted by the per-vertex 'atom_color'
    attribute, so every polyhedron takes the element colors of its corner
    atoms (e.g. brown SbBr6 octahedra from the Br corners).

    The principled half carries the color and the transparency, the glass
    half the refraction and the bright edges a solid polyhedron gets where
    it is seen at a glancing angle.

    Only built when the material does not already have that glass half:
    an older 'polyhedra material' is rebuilt, one that is already this
    keeps whatever was tweaked on it.
    """
    mat = newMaterial(POLYHEDRA_MATERIAL)
    tree = mat.node_tree
    nodes, links = tree.nodes, tree.links
    if any(node.bl_idname == 'ShaderNodeBsdfGlass' for node in nodes):
        return mat
    for node in list(nodes):
        nodes.remove(node)

    attribute = nodes.new('ShaderNodeAttribute')
    attribute.name = 'Attribute'
    attribute.attribute_name = 'atom_color'
    attribute.attribute_type = 'GEOMETRY'
    attribute.location = (-320, 60)

    principled = nodes.new('ShaderNodeBsdfPrincipled')
    principled.location = (-60, 300)
    principled.inputs['Metallic'].default_value = 0.0
    principled.inputs['Roughness'].default_value = 0.6
    principled.inputs['IOR'].default_value = 1.45
    principled.inputs['Alpha'].default_value = 0.3

    glass = nodes.new('ShaderNodeBsdfGlass')
    glass.location = (-60, -120)
    glass.distribution = 'MULTI_GGX'
    glass.inputs['Roughness'].default_value = 0.0
    glass.inputs['IOR'].default_value = 1.5

    mix = nodes.new('ShaderNodeMixShader')
    mix.location = (280, 160)
    mix.inputs[0].default_value = 0.5   # index: the socket is 'Fac' on 4.x, 'Factor' on 5.x

    output = nodes.new('ShaderNodeOutputMaterial')
    output.location = (500, 160)

    links.new(attribute.outputs['Color'], principled.inputs['Base Color'])
    links.new(attribute.outputs['Color'], glass.inputs['Color'])
    links.new(principled.outputs['BSDF'], mix.inputs[1])
    links.new(glass.outputs['BSDF'], mix.inputs[2])
    links.new(mix.outputs[0], output.inputs['Surface'])
    return mat


def bond_neighbors(atoms, bond_cutoff, skin=0.0):
    """Periodic neighbor list in which a bond is d < bond_cutoff*(r1+r2).

    skin defaults to 0.0 on purpose. ASE's own default is 0.3, which it
    adds to *each* atom's radius, so the criterion silently becomes
    d < mult*r1 + mult*r2 + 0.6 A. An additive term over-inflates small
    radii - hydrogen goes from 1.3*0.31 = 0.40 to 0.70 A - which counts
    N-H...Br hydrogen bonds (2.30 A in a relaxed hybrid antimony bromide)
    as covalent bonds. That fuses the organic spacers and the Sb2Br10
    anions into one endless network, and then there is no molecule left to
    complete. The expansion, trim and polyhedra-shell searches pass ASE's
    0.3 back in, because their default multipliers were tuned with it.
    """
    return NeighborList([covalent_radii[num] * bond_cutoff for num in atoms.numbers],
                        self_interaction=False, bothways=True, skin=skin)


def _components(neighbors):
    """Connected components of an adjacency list, as lists of indices."""
    components = []
    seen = set()
    for start in range(len(neighbors)):
        if start in seen:
            continue
        component, stack = [], [start]
        seen.add(start)
        while stack:
            i = stack.pop()
            component.append(i)
            for j in neighbors[i]:
                if j not in seen:
                    seen.add(int(j))
                    stack.append(int(j))
        components.append(sorted(component))
    return components


def grow_shells(neighbors, seeds, max_shells=None, stop_on_second_image=True,
                node_cap=None):
    """Grow a cluster shell by shell over (atom index, cell image) nodes.

    Each step adds the bonded neighbors of the shell before it, carrying
    the accumulated image offset along every bond, so the cluster grows in
    unwrapped space: a molecule the cell cuts in half grows straight across
    the boundary and comes out in one piece.

    Returns (nodes, shells, closed). `closed` is True when a shell added
    nothing new - one more shell would not grow the cluster, so the
    molecule is complete. Reaching a *second image* of an atom the cluster
    already holds means a bond path came back to that atom with a nonzero
    net lattice translation: the component is an extended framework (chain,
    layer, 3d network) with no finite molecule to import, and the growth
    stops with closed=False. That test is off when growing shells around a
    framework, where several images of one atom are exactly the point.

    `max_shells` stops after that many shells and `node_cap` is a backstop
    against a runaway cluster; both report closed=False.
    """
    nodes = set(seeds)
    images_of = {}
    for index, offset in nodes:
        images_of.setdefault(index, set()).add(offset)
    frontier = list(nodes)
    shells = 0
    while frontier:
        if max_shells is not None and shells >= max_shells:
            return nodes, shells, False
        new = []
        for index, offset in frontier:
            for other, image in zip(*neighbors[index]):
                other = int(other)
                shifted = (offset[0] + int(image[0]),
                           offset[1] + int(image[1]),
                           offset[2] + int(image[2]))
                if shifted in images_of.get(other, ()):
                    continue
                if stop_on_second_image and other in images_of:
                    return nodes, shells, False
                nodes.add((other, shifted))
                images_of.setdefault(other, set()).add(shifted)
                new.append((other, shifted))
        if node_cap is not None and len(nodes) > node_cap:
            print(f'polyhedra: a molecule grew past {node_cap} atoms - '
                  'treating it as an extended framework')
            return nodes, shells, False
        frontier = new
        shells += 1
    return nodes, shells, True


def select_complete_molecules(atoms, bond_cutoff=1.3, cell_margin=0.0,
                              all_images=True, framework_shells=1):
    """Whole molecules out of a periodic cell, grown shell by shell.

    Every molecule of the cell is grown from one of its atoms over
    (atom, cell image) nodes until one more shell would not add anything
    (see `grow_shells`). The result is one template per molecule -
    unwrapped, so a molecule cut by a cell face is whole, and of any size,
    including molecules longer than the cell itself.

    Each atom already lies in some lattice cell, `home` below; an atom of a
    template moved by the lattice translation t therefore sits inside the
    central cell exactly when `offset + t == -home`. The translations that
    place at least one of a molecule's atoms in the cell are just the
    negated (offset + home) values, so the copies fall out of the growth
    itself - no supercell is ever built:

    - `all_images=True` imports the molecule at every one of them, so a
      molecule straddling a cell face arrives twice, one straddling an edge
      four times, and the cell stays fully populated instead of keeping the
      hole a cut molecule came from.
    - `all_images=False` imports one copy only, the image holding the most
      of the cell's own atoms (ties go to the shortest translation, so the
      choice is reproducible: the Sb2Br10 dimers of a hybrid bromide split
      6/6 over two images and 3/3/3/3 over four).

    `cell_margin` grows the region a molecule has to reach, in angstrom,
    pulling in the surrounding molecules that come that close to the cell.

    A component that never closes is an extended framework, which has no
    molecule to complete: its atoms in the cell are grown by
    `framework_shells` bonded shells instead - enough to close their
    coordination polyhedra at the cell boundary at the default of 1.

    Returns (new_atoms, trimmable): the non-periodic selection, and a
    per-atom flag for the caller's trim pass. A molecule is complete by
    construction and never trimmable; a framework's added shells are, just
    as they were before.
    """
    import itertools

    cell = np.asarray(atoms.get_cell())
    positions = atoms.get_positions()
    nl = bond_neighbors(atoms, bond_cutoff)
    nl.update(atoms)
    neighbors = [nl.get_neighbors(i) for i in range(len(atoms))]

    # the lattice cell each atom already lies in. Deliberately not
    # atoms.wrap(): several fixtures store atoms at fractional -0.5 or a
    # hair outside [0, 1), and moving those would shift the structure the
    # importer has always drawn.
    fractional = np.linalg.solve(cell.T, positions.T).T
    home = np.floor(fractional + 1e-8).astype(int)
    margin = cell_margin * np.linalg.norm(np.asarray(atoms.get_cell().reciprocal()), axis=1)

    def translations(template):
        """The lattice translations of one molecule that get imported."""
        candidates = {tuple(-(np.asarray(offset) + home[index]))
                      for index, offset in template}
        if margin.any():
            # also the copies that merely come within cell_margin of the
            # cell: widen by one cell each way and test geometrically
            candidates = {tuple(np.asarray(t) + d)
                          for t in candidates
                          for d in itertools.product((-1, 0, 1), repeat=3)}
            candidates = {t for t in candidates
                          if any(np.all((fractional[index] + offset + t > -margin)
                                        & (fractional[index] + offset + t < 1 + margin))
                                 for index, offset in template)}
        if not all_images and candidates:
            inside = lambda t: sum(  # noqa: E731 - local scoring helper
                1 for index, offset in template
                if not (np.asarray(offset) + t + home[index]).any())
            return [max(sorted(candidates),
                        key=lambda t: (inside(t), -sum(abs(x) for x in t)))]
        return sorted(candidates)

    molecules, framework = [], []
    for component in _components([indices for indices, _ in neighbors]):
        template, _shells, closed = grow_shells(
            neighbors, [(component[0], (0, 0, 0))],
            node_cap=max(50 * len(atoms), 10000))
        if closed:
            molecules.append(sorted(template))
        else:
            # retire the whole component: re-seeding the atoms the aborted
            # growth never reached would grow bogus partial molecules out
            # of the same endless network
            framework.append(component)

    nodes, trimmable = [], []
    for template in molecules:
        for translation in translations(template):
            for index, offset in template:
                nodes.append((index, tuple(np.asarray(offset) + translation)))
                trimmable.append(False)

    if framework:
        seeds = [(index, (0, 0, 0)) for component in framework for index in component]
        grown, _shells, _closed = grow_shells(
            neighbors, seeds, max_shells=framework_shells,
            stop_on_second_image=False)
        placed = set(nodes)
        nodes.extend(seeds)
        trimmable.extend([False] * len(seeds))
        shell = sorted(set(grown) - set(seeds) - placed)
        nodes.extend(shell)
        trimmable.extend([True] * len(shell))

    # deduplicate, first occurrence wins, and an atom that any cluster
    # needs for its own sake is never trimmable
    unique = {}
    for node, trim in zip(nodes, trimmable):
        unique[node] = trim if node not in unique else (unique[node] and trim)
    keep = list(unique)
    trimmable = [unique[node] for node in keep]
    new_atoms = Atoms(numbers=[atoms.numbers[index] for index, _ in keep],
                      positions=[positions[index] + np.asarray(offset) @ cell
                                 for index, offset in keep])
    new_atoms.pbc = False
    new_atoms.cell = None
    return new_atoms, trimmable


def build_polyhedra_atoms(atoms, expand_cutoff=1.2, trim_cutoff=1.0,
                          poly_cutoff=1.1, min_neighbors=4,
                          include_hydrogen=False, single_element_corners=True,
                          complete_molecules=True, bond_cutoff=1.3,
                          cell_margin=0.0, all_images=True,
                          framework_shells=1):
    """Extend the structure past the cell and compute convex-hull faces
    around every coordination center.

    Returns (new_atoms, faces): a non-periodic Atoms object whose
    positions are the mesh vertices, and the polyhedra faces as vertex
    index lists.

    complete_molecules: import whole molecules -- grow each one shell by
                    shell out of the periodic cell, so a molecule the cell
                    cuts in half arrives in one piece (see
                    select_complete_molecules). Off falls back to the plain
                    image expansion below, as does input without a 3d cell.
    bond_cutoff:    what counts as a bond while growing molecules: two
                    atoms belong to the same molecule below this multiple
                    of the sum of their covalent radii. Skin-free, unlike
                    the multipliers below - see bond_neighbors
    cell_margin:    angstrom to grow the selection region around the cell
                    by, pulling in the surrounding whole molecules that
                    come within that distance of it
    all_images:     import every image of a molecule that reaches into the
                    cell (off: one whole copy per molecule)
    framework_shells: bonded shells to grow around the cell's atoms of an
                    extended framework, which has no molecule to complete
    expand_cutoff:  covalent-radius multiplier pulling in periodic neighbor
                    images in the plain expansion, so boundary polyhedra
                    close there too
    trim_cutoff:    multiplier for removing atoms left without neighbors
                    after the expansion (never applied to a complete
                    molecule, which is whole by construction)
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

    if complete_molecules and atoms.cell.rank == 3:
        new_atoms, trimmable = select_complete_molecules(
            atoms, bond_cutoff=bond_cutoff, cell_margin=cell_margin,
            all_images=all_images, framework_shells=framework_shells)
    else:
        # the plain expansion: pull in the periodic images of every
        # neighbor, deduplicated per (atom, image)
        positions = atoms.get_positions()
        cell = atoms.get_cell()
        # ASE's 0.3 skin kept here on purpose: the expansion defaults were
        # tuned with it (see bond_neighbors for what it does to the criterion)
        nl = bond_neighbors(atoms, expand_cutoff, skin=0.3)
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
        trimmable = [True] * len(new_atoms)

    # drop expanded atoms that ended up without any neighbor
    # skin kept: trim_cutoff=1.0 is only generous enough to keep a bonded
    # image because ASE's 0.3 is added to each radius
    nl_trim = bond_neighbors(new_atoms, trim_cutoff, skin=0.3)
    nl_trim.update(new_atoms)
    for i in reversed(range(len(new_atoms))):
        if trimmable[i] and len(nl_trim.get_neighbors(i)[0]) == 0:
            del new_atoms[i]

    # convex hull of each coordination shell
    # skin kept: the hull search needs ~3.45 A to catch a bridging Br of an
    # SbBr6 octahedron, which poly_cutoff=1.1 only reaches with it
    nl = bond_neighbors(new_atoms, poly_cutoff, skin=0.3)
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
                     complete_molecules=True, bond_cutoff=1.3, cell_margin=0.0,
                     all_images=True, framework_shells=1, unit_cell=False,
                     **kwargs):
    import ase.io
    atoms = ase.io.read(filepath)

    new_atoms, faces = build_polyhedra_atoms(
        atoms, expand_cutoff=expand_cutoff, trim_cutoff=trim_cutoff,
        poly_cutoff=poly_cutoff, min_neighbors=min_neighbors,
        include_hydrogen=include_hydrogen,
        single_element_corners=single_element_corners,
        complete_molecules=complete_molecules, bond_cutoff=bond_cutoff,
        cell_margin=cell_margin, all_images=all_images,
        framework_shells=framework_shells)
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
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_2", bond_distance)
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_3", bond_radius)
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_4", resolution)

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
    # carry the atomic numbers over as well: 'element' next to 'atom_color'
    # is what lets a later color change (ASE sidebar swatch) find the
    # points of one element and repaint the faces with them
    element_attr = poly_mesh.attributes.new(name='element', type='FLOAT',
                                            domain='POINT')
    numbers = np.empty(len(mesh.vertices))
    mesh.attributes['element'].data.foreach_get('value', numbers)
    element_attr.data.foreach_set('value', numbers)
    poly_mesh.materials.append(_polyhedra_material())
    poly_mesh.update()
    poly_obj = bpy.data.objects.new(poly_mesh.name, poly_mesh)
    my_coll.objects.link(poly_obj)

    if outline:
        # atoms and bonds only - the polyhedra object stays outline-free
        outline_objects([obj], modifier='GeometryNodes.001')

    if unit_cell and atoms.cell.rank == 3:
        # into this structure's collection, and after the outline pass so
        # the flat black cell edges never get an outline of their own
        draw_unit_cell(atoms)
    return obj
