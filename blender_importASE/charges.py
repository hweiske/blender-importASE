"""Import a structure with per-atom partial charges from a csv file.

The csv contains one charge per atom, in the same order as the atoms in
the structure file (either a bare column of numbers or rows whose last
numeric field is the charge; header lines are skipped). The values are
stored as a 'charge' float point attribute, transferred onto the atom
spheres and bond endpoints by the atoms_and_bonds node group
(start_charge / end_charge / CHARGE_CURVE attributes), and visualized by
two generated materials:

- 'charge_atoms': atoms colored by charge (red negative, white neutral,
  blue positive)
- 'color_curve_charge': bonds with a charge gradient between their two
  atoms, analogous to the element color curve

A 'charge_colors' switch on the modifier (visible in the ASE panel)
toggles between the charge materials and the normal element colors.
"""
import bpy
import numpy as np

from .utils import atomcolors
from .node_networks.nodes_atoms_and_bonds import set_atoms_node_group, atoms_and_bonds, read_structure
from .node_networks.bond_mat import create_bondmat
from .node_networks.electron_density_nodes import newMaterial
from .node_networks.outline import outline_objects
from .node_networks.compat import set_mod_input

CHARGE_ATOMS_MATERIAL = 'charge_atoms'
CHARGE_BONDS_MATERIAL = 'color_curve_charge'


def read_charges_csv(filepath):
    """One charge per row, same order as the atoms. Accepts a bare number
    per line or csv/semicolon rows whose last numeric field is the charge
    (e.g. 'C, 0.23'); non-numeric rows (headers) are skipped."""
    values = []
    with open(filepath) as fh:
        for line in fh:
            fields = line.replace(';', ',').split(',')
            value = None
            for field in reversed(fields):
                try:
                    value = float(field)
                    break
                except ValueError:
                    continue
            if value is not None:
                values.append(value)
    return np.array(values)


def _charge_material(name, attribute):
    """Material reading a float attribute through a symmetric
    red-white-blue ramp (red negative, blue positive)."""
    mat = newMaterial(name)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    principled = nodes.get('Principled BSDF')

    attr = nodes.get('Attribute')
    if attr is None:
        attr = nodes.new('ShaderNodeAttribute')
        attr.name = 'Attribute'
        attr.location = (-700, 200)
    attr.attribute_name = attribute
    attr.attribute_type = 'GEOMETRY'

    map_range = nodes.get('Map Range')
    if map_range is None:
        map_range = nodes.new('ShaderNodeMapRange')
        map_range.name = 'Map Range'
        map_range.location = (-500, 200)
        map_range.clamp = True

    ramp = nodes.get('Color Ramp')
    if ramp is None:
        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.name = 'Color Ramp'
        ramp.location = (-300, 200)
        ramp.color_ramp.elements[0].color = (0.85, 0.10, 0.10, 1)  # negative
        ramp.color_ramp.elements[1].color = (0.05, 0.15, 0.85, 1)  # positive
        mid = ramp.color_ramp.elements.new(0.5)
        mid.color = (0.95, 0.95, 0.95, 1)

    links.new(attr.outputs['Fac'], map_range.inputs['Value'])
    links.new(map_range.outputs['Result'], ramp.inputs['Fac'])
    links.new(ramp.outputs['Color'], principled.inputs['Base Color'])
    return mat


def make_charge_materials(charges):
    """Create/update the charge materials with a symmetric value range
    around zero, so zero charge is always white."""
    limit = float(np.abs(charges).max()) or 1.0
    materials = []
    for name, attribute in ((CHARGE_ATOMS_MATERIAL, 'charge'),
                            (CHARGE_BONDS_MATERIAL, 'CHARGE_CURVE')):
        mat = _charge_material(name, attribute)
        map_range = mat.node_tree.nodes['Map Range']
        map_range.inputs['From Min'].default_value = -limit
        map_range.inputs['From Max'].default_value = limit
        materials.append(mat)
    return materials


def import_charges(filepath, filename, charge_filepath, resolution=16,
                   colorbonds=True, bond_distance=0.66, bond_radius=0.1,
                   outline=False, **kwargs):
    import ase.io
    atoms = ase.io.read(filepath)

    charges = read_charges_csv(charge_filepath)
    if len(charges) != len(atoms):
        raise ValueError(
            f'{charge_filepath} contains {len(charges)} charges but '
            f'{filepath} has {len(atoms)} atoms')

    atomcolor = atomcolors()
    atomcolor.setup_materials(atoms, colorbonds=colorbonds)
    my_coll = bpy.data.collections.new(
        name=atoms.get_chemical_formula() + '_charges_' + filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection

    obj, mesh = read_structure(atoms,
                               atoms.get_chemical_formula() + '_charges_' + filename.split('.')[0],
                               animate=False)

    # per-atom charges as a float point attribute
    mesh.attributes.new(name='charge', type='FLOAT', domain='POINT')
    mesh.attributes['charge'].data.foreach_set('value', charges)
    mesh.update()

    set_atoms_node_group()
    elements_name = '_'.join(list(set(atoms.get_chemical_symbols())))
    bondmat = create_bondmat(colorbonds=colorbonds, name=elements_name)
    atoms_from_verts = atoms_and_bonds(obj, atoms, 'GeometryNodes',
                                       bondmat=bondmat, with_charges=True)
    obj.modifiers['GeometryNodes'].node_group = atoms_from_verts
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_2", bond_distance)
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_3", bond_radius)
    set_mod_input(obj.modifiers['GeometryNodes'], "Socket_4", resolution)

    # the charge materials occupy the two slots after the bond material;
    # the node tree's charge_colors switch points the faces at them
    for mat in make_charge_materials(charges):
        obj.data.materials.append(mat)

    if outline:
        outline_objects([obj], modifier='GeometryNodes.001')
    return obj
