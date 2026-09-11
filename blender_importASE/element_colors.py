"""Editing the per-element colors of an imported structure.

An element's color lives in two places once a structure is imported:

- the per-element materials ('Sb', 'Sb-bond', 'Sb-C-bond', ...), which
  shade the atom spheres and the ball'n'stick bonds, and
- the 'atom_color' point attribute of the structure mesh, which the
  geometry-node bonds sample at both ends and blend along the curve (the
  colored-bonds option), and which tints polyhedra faces.

Changing a material color therefore only ever repainted the atoms - the
node bonds kept the color the structure was imported with, and the
attribute itself is not something the Blender UI lets you edit per
element. This module keeps the two in sync:

- the ASE sidebar draws one color swatch per element of the active
  structure (`draw_element_colors`), bound straight to the element
  material, and
- a depsgraph handler notices element-material color changes - from that
  swatch, the Material Properties tab, a driver or a script - and rewrites
  the matching 'atom_color' entries of every structure in the file.

`ASE_OT_sync_element_colors` does the same pass on demand, and
`ASE_OT_reset_element_colors` puts an element back to the add-on's
default color.
"""
import re

import bpy
import numpy as np
from ase.data import atomic_numbers, chemical_symbols
from bpy.app.handlers import persistent

from .utils import default_element_color

# 'Sb-C-bond' / 'Sb-C-bond-smooth': the two-sided gradient material of a
# bond between two different elements (see atomcolors.create_bondmat)
PAIR_BOND_RE = re.compile(r'^([A-Z][a-z]?)-([A-Z][a-z]?)-bond(-smooth)?$')
# node names create_bondmat gives the two ends of such a gradient
PAIR_BOND_NODES = ('Atom 1 Principled BSDF', 'Atom 2 Principled BSDF')

BASE_COLOR = 0  # 'Base Color' input of a Principled BSDF, every version

_SYMBOLS = frozenset(chemical_symbols[1:])

# element materials whose color the handler watches, and the colors it last
# saw them at. Rebuilt when the material count changes (a new import), so
# the handler itself never has to walk bpy.data.
_watched = []
_watched_key = None
_seen_colors = {}
_syncing = False


def _bsdf(mat, name=None):
    """The Principled BSDF of a material, by node name where the material
    builder gave it one, else the first principled node in the tree."""
    if not mat.use_nodes or mat.node_tree is None:
        return None
    nodes = mat.node_tree.nodes
    if name is not None:
        return nodes.get(name)
    node = nodes.get('Principled BSDF')
    if node is not None:
        return node
    return next((n for n in nodes if n.type == 'BSDF_PRINCIPLED'), None)


def _color_sockets(symbol):
    """Every base-color socket that paints `symbol`: its atom material, its
    same-element bond material, and its side of each mixed-pair bond
    material."""
    for name in (symbol, f'{symbol}-bond'):
        mat = bpy.data.materials.get(name)
        node = _bsdf(mat) if mat is not None else None
        if node is not None:
            yield node.inputs[BASE_COLOR]
    for mat in bpy.data.materials:
        match = PAIR_BOND_RE.match(mat.name)
        if match is None:
            continue
        for sym, node_name in zip(match.group(1, 2), PAIR_BOND_NODES):
            if sym != symbol:
                continue
            node = _bsdf(mat, node_name)
            if node is not None:
                yield node.inputs[BASE_COLOR]


def get_element_color(symbol):
    """The color `symbol` is currently drawn in: its atom material's base
    color, or the add-on default when no material exists (yet)."""
    mat = bpy.data.materials.get(symbol)
    node = _bsdf(mat) if mat is not None else None
    if node is None:
        return default_element_color(symbol)
    return tuple(node.inputs[BASE_COLOR].default_value)


def set_element_color(symbol, color):
    """Paint `symbol` in `color` (RGB or RGBA): every material that shades
    it, plus the 'atom_color' attribute the node bonds read."""
    rgba = tuple(color[:3]) + (color[3] if len(color) > 3 else 1.0,)
    for socket in _color_sockets(symbol):
        socket.default_value = rgba
    _seen_colors[symbol] = rgba
    sync_atom_color_attributes([symbol])


def element_color_socket(symbol):
    """The socket the sidebar swatch edits - the atom material's base
    color - or None when this element has no material."""
    mat = bpy.data.materials.get(symbol)
    node = _bsdf(mat) if mat is not None else None
    return None if node is None else node.inputs[BASE_COLOR]


def structure_symbols(obj):
    """Elements of the structure `obj` belongs to, as symbols.

    'ase_elements' is written by the node importer; fall back to the mesh's
    own 'element' attribute and finally to the element-named atom objects of
    the collection (the ball'n'stick and 3D-print representations).
    """
    numbers = list(obj.get('ase_elements', []))
    if not numbers:
        mesh = obj.data if obj.type == 'MESH' else None
        if mesh is not None and 'element' in mesh.attributes:
            data = mesh.attributes['element'].data
            values = np.empty(len(data), dtype=np.float32)
            data.foreach_get('value', values)
            numbers = np.unique(np.rint(values).astype(int)).tolist()
    symbols = [chemical_symbols[int(z)] for z in sorted(numbers)
               if 0 < int(z) < len(chemical_symbols)]
    if not symbols and obj.users_collection:
        symbols = sorted({o.name.split('.')[0]
                          for o in obj.users_collection[0].all_objects
                          if o.name.split('.')[0] in _SYMBOLS})
    return symbols


def sync_atom_color_attributes(symbols):
    """Rewrite the 'atom_color' entries of the given elements from their
    material colors, in every mesh of the file that carries one.

    Structure meshes (and the polyhedra faces mesh) pair 'atom_color' with
    an 'element' attribute on the same point domain, which is what says
    which points belong to which element.
    """
    wanted = {atomic_numbers[sym]: np.array(get_element_color(sym), dtype=np.float32)
              for sym in symbols if sym in atomic_numbers}
    if not wanted:
        return []
    touched = []
    for mesh in bpy.data.meshes:
        if 'atom_color' not in mesh.attributes or 'element' not in mesh.attributes:
            continue
        color_data = mesh.attributes['atom_color'].data
        element_data = mesh.attributes['element'].data
        count = len(color_data)
        if count == 0 or len(element_data) != count:
            continue  # not a structure mesh: different domains
        numbers = np.empty(count, dtype=np.float32)
        element_data.foreach_get('value', numbers)
        numbers = np.rint(numbers).astype(int)
        colors = np.empty(count * 4, dtype=np.float32)
        color_data.foreach_get('color', colors)
        colors = colors.reshape(count, 4)
        changed = False
        for number, rgba in wanted.items():
            mask = numbers == number
            if mask.any() and not np.allclose(colors[mask], rgba):
                colors[mask] = rgba
                changed = True
        if changed:
            color_data.foreach_set('color', colors.ravel())
            mesh.update()
            touched.append(mesh)
    if touched:
        # make the geometry-node bonds re-read the attribute
        for obj in bpy.data.objects:
            if obj.data in touched:
                obj.update_tag()
    return touched


def _watched_materials():
    """Element materials to watch for color edits. Cached: the handler runs
    on every depsgraph update, so it must not walk bpy.data each time."""
    global _watched, _watched_key
    key = len(bpy.data.materials)
    if key != _watched_key:
        _watched = [mat.name for mat in bpy.data.materials if mat.name in _SYMBOLS]
        _watched_key = key
    return _watched


@persistent
def _track_element_colors(scene, depsgraph=None):
    """Push element-material color edits into the 'atom_color' attributes.

    Polls the watched materials rather than filtering depsgraph.updates, so
    it catches the color changing from anywhere (sidebar swatch, Material
    Properties, driver, script) on every Blender version. Writing the
    attributes tags the structures, which fires this handler again - the
    colors then match what was last seen, so it stops there.
    """
    global _syncing
    if _syncing:
        return
    changed = []
    for symbol in _watched_materials():
        color = get_element_color(symbol)
        if _seen_colors.get(symbol) != color:
            _seen_colors[symbol] = color
            changed.append(symbol)
    if not changed:
        return
    _syncing = True
    try:
        sync_atom_color_attributes(changed)
    finally:
        _syncing = False


@persistent
def _reset_cache(*args):
    """Another file's materials are unrelated to the ones cached here."""
    global _watched_key
    _watched_key = None
    _seen_colors.clear()


class ASE_OT_sync_element_colors(bpy.types.Operator):
    """Repaint the bonds from the current element colors"""
    bl_idname = 'ase.sync_element_colors'
    bl_label = 'Apply colors to bonds'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        _reset_cache()  # also picks up element materials added since
        symbols = structure_symbols(context.active_object)
        touched = sync_atom_color_attributes(symbols)
        if touched:
            self.report({'INFO'}, f'repainted {len(touched)} structure mesh(es)')
        else:
            self.report({'INFO'}, 'bonds already match the element colors')
        return {'FINISHED'}


class ASE_OT_reset_element_colors(bpy.types.Operator):
    """Put this structure's elements back to their default colors"""
    bl_idname = 'ase.reset_element_colors'
    bl_label = 'Default colors'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for symbol in structure_symbols(context.active_object):
            set_element_color(symbol, default_element_color(symbol))
        return {'FINISHED'}


def draw_element_colors(layout, obj):
    """Sidebar section: one color swatch per element of this structure.

    The swatch edits the element's atom material directly, so it also
    covers structures imported before this panel existed; the handler above
    carries the change over to the bonds.
    """
    symbols = structure_symbols(obj)
    sockets = [(sym, element_color_socket(sym)) for sym in symbols]
    sockets = [(sym, socket) for sym, socket in sockets if socket is not None]
    if not sockets:
        return
    box = layout.box()
    box.label(text='Element colors')
    flow = box.grid_flow(row_major=True, columns=2, align=True)
    for symbol, socket in sockets:
        row = flow.row(align=True)
        row.label(text=symbol)
        row.prop(socket, 'default_value', text='')
    row = box.row(align=True)
    row.operator('ase.sync_element_colors', icon='FILE_REFRESH', text='Apply to bonds')
    row.operator('ase.reset_element_colors', icon='LOOP_BACK', text='Defaults')


classes = (ASE_OT_sync_element_colors, ASE_OT_reset_element_colors)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    if _track_element_colors not in bpy.app.handlers.depsgraph_update_post:
        bpy.app.handlers.depsgraph_update_post.append(_track_element_colors)
    if _reset_cache not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(_reset_cache)


def unregister():
    for handlers, handler in ((bpy.app.handlers.depsgraph_update_post, _track_element_colors),
                              (bpy.app.handlers.load_post, _reset_cache)):
        if handler in handlers:
            handlers.remove(handler)
    _reset_cache()
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
