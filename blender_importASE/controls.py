"""Control tables and UI for per-element/per-pair settings.

Per-element and per-element-pair options used to be represented as one
modifier socket plus a dedicated node sub-network each, which made node
group creation scale quadratically with the number of distinct elements.
They are now stored as data: small hidden table meshes whose point
attributes are looked up from a constant-size node chain.

- pair table:    one point per element pair id (min(Z1,Z2)*119 + max),
                 boolean attribute 'cut' - True hides bonds of that pair.
- element table: one point per atomic number, integer attribute
                 'radius_mode' - 0 = covalent radius, 1 = vdW radius.

The tables are exposed in the sidebar (N panel) of the 3D viewport under
the 'ASE' tab.
"""
import bpy
from ase.data import chemical_symbols

PAIR_STRIDE = 119  # > max atomic number, so pair ids are unique


def pair_id(z1, z2):
    return min(z1, z2) * PAIR_STRIDE + max(z1, z2)


def make_control_tables(name, numbers):
    """Create the hidden pair/element table objects for one structure."""
    zmax = max(numbers)

    pair_mesh = bpy.data.meshes.new(f'{name}_pair_table')
    pair_mesh.from_pydata([(0, 0, 0)] * (pair_id(zmax, zmax) + 1), [], [])
    pair_mesh.attributes.new(name='cut', type='BOOLEAN', domain='POINT')
    pair_obj = bpy.data.objects.new(pair_mesh.name, pair_mesh)

    element_mesh = bpy.data.meshes.new(f'{name}_element_table')
    element_mesh.from_pydata([(0, 0, 0)] * (zmax + 1), [], [])
    element_mesh.attributes.new(name='radius_mode', type='INT', domain='POINT')
    element_obj = bpy.data.objects.new(element_mesh.name, element_mesh)

    for ob in (pair_obj, element_obj):
        bpy.context.collection.objects.link(ob)
        ob.hide_viewport = True
        ob.hide_render = True
    return pair_obj, element_obj


def find_ase_modifier(obj):
    """Return (modifier, {socket name: identifier}) of the atoms_and_bonds
    modifier on obj, or (None, None)."""
    if obj is None:
        return None, None
    for m in obj.modifiers:
        if m.type != 'NODES' or m.node_group is None:
            continue
        idents = {}
        for item in m.node_group.interface.items_tree:
            if getattr(item, 'in_out', None) == 'INPUT':
                idents[item.name] = item.identifier
        if 'pair_table' in idents and 'element_table' in idents:
            return m, idents
    return None, None


def _table_object(context, which):
    mod, idents = find_ase_modifier(context.active_object)
    if mod is None:
        return None
    return mod.get(idents[which])


def _touch(table_obj):
    # make dependent geometry-nodes modifiers re-evaluate
    table_obj.data.update()
    table_obj.update_tag()
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            area.tag_redraw()


class ASE_OT_toggle_pair_cut(bpy.types.Operator):
    """Show or hide bonds between this element pair"""
    bl_idname = 'ase.toggle_pair_cut'
    bl_label = 'Toggle bonds of element pair'
    bl_options = {'REGISTER', 'UNDO'}

    pair_id: bpy.props.IntProperty()

    def execute(self, context):
        table = _table_object(context, 'pair_table')
        if table is None:
            self.report({'WARNING'}, 'no ASE pair table on active object')
            return {'CANCELLED'}
        data = table.data.attributes['cut'].data
        data[self.pair_id].value = not data[self.pair_id].value
        _touch(table)
        return {'FINISHED'}


class ASE_OT_set_radius_mode(bpy.types.Operator):
    """Switch this element between covalent and vdW radius"""
    bl_idname = 'ase.set_radius_mode'
    bl_label = 'Toggle element radius mode'
    bl_options = {'REGISTER', 'UNDO'}

    number: bpy.props.IntProperty()

    def execute(self, context):
        table = _table_object(context, 'element_table')
        if table is None:
            self.report({'WARNING'}, 'no ASE element table on active object')
            return {'CANCELLED'}
        data = table.data.attributes['radius_mode'].data
        data[self.number].value = 0 if data[self.number].value else 1
        _touch(table)
        return {'FINISHED'}


# node group name prefix -> panel section title
GROUP_TITLES = (
    ('atoms_and_bonds', 'Atoms & bonds'),
    ('hide atoms', 'Hide atoms'),
    ('supercell', 'Supercell'),
    ('outline', 'Outline'),
    ('BONDS_GEOMETRY', 'Bonds'),
    ('visualize_edensity', 'Electron density'),
)

# socket types that make no sense as a panel slider
SKIP_SOCKET_TYPES = {'NodeSocketGeometry', 'NodeSocketObject',
                     'NodeSocketCollection', 'NodeSocketMaterial',
                     'NodeSocketImage', 'NodeSocketMenu'}


def _group_title(node_group_name):
    for prefix, title in GROUP_TITLES:
        if node_group_name.startswith(prefix):
            return title
    return None


def _draw_modifier_inputs(layout, mod):
    """Draw all scalar inputs of a geometry-nodes modifier as regular,
    keyframeable modifier properties."""
    keys = mod.keys()
    col = layout.column(align=True)
    for item in mod.node_group.interface.items_tree:
        if getattr(item, 'in_out', None) != 'INPUT':
            continue
        if item.socket_type in SKIP_SOCKET_TYPES:
            continue
        if item.identifier not in keys:
            continue
        col.prop(mod, f'["{item.identifier}"]', text=item.name)


def _iter_gn_modifiers(obj):
    for mod in obj.modifiers:
        if mod.type == 'NODES' and mod.node_group is not None:
            title = _group_title(mod.node_group.name)
            if title:
                yield mod, title


def _find_density_modifiers(obj):
    """Find all visualize_edensity modifiers on sibling objects (the density
    volumes imported alongside this structure, e.g. total charge + spin)."""
    found = []
    seen = set()
    for coll in obj.users_collection:
        for ob in coll.objects:
            if ob == obj or ob.name in seen:
                continue
            for mod in ob.modifiers:
                if (mod.type == 'NODES' and mod.node_group is not None
                        and mod.node_group.name.startswith('visualize_edensity')):
                    found.append((ob, mod))
                    seen.add(ob.name)
    return found


class ASE_PT_controls(bpy.types.Panel):
    bl_label = 'ASE structure'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ASE'

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and next(_iter_gn_modifiers(obj), None) is not None

    def draw(self, context):
        obj = context.active_object
        for mod, title in _iter_gn_modifiers(obj):
            box = self.layout.box()
            box.label(text=title)
            _draw_modifier_inputs(box, mod)
            if mod.node_group.name.startswith('atoms_and_bonds'):
                self.draw_tables(context, box)

        # electron densities live on sibling objects in the same collection
        # (a spin-polarized CHGCAR yields a total and a spin volume)
        if not any(t == 'Electron density' for _, t in _iter_gn_modifiers(obj)):
            for density_obj, density_mod in _find_density_modifiers(obj):
                box = self.layout.box()
                title = 'Spin difference' if '_spin' in density_obj.name else 'Electron density'
                box.label(text=f'{title} ({density_obj.name})')
                _draw_modifier_inputs(box, density_mod)

    def draw_tables(self, context, layout):
        obj = context.active_object
        numbers = sorted(obj.get('ase_elements', []))
        if not numbers:
            return

        pair_table = _table_object(context, 'pair_table')
        element_table = _table_object(context, 'element_table')

        if element_table is not None:
            modes = element_table.data.attributes['radius_mode'].data
            layout.label(text='Atom radius (covalent / vdW)')
            flow = layout.grid_flow(row_major=True, columns=3, align=True)
            for z in numbers:
                mode = bool(modes[z].value)
                op = flow.operator('ase.set_radius_mode',
                                   text=f'{chemical_symbols[z]}: {"vdW" if mode else "cov"}',
                                   depress=mode)
                op.number = z

        if pair_table is not None:
            cuts = pair_table.data.attributes['cut'].data
            layout.label(text='Hide bonds between')
            flow = layout.grid_flow(row_major=True, columns=3, align=True)
            for i, zi in enumerate(numbers):
                for zj in numbers[i:]:
                    pid = pair_id(zi, zj)
                    cut = bool(cuts[pid].value)
                    op = flow.operator('ase.toggle_pair_cut',
                                       text=f'{chemical_symbols[zi]}-{chemical_symbols[zj]}',
                                       depress=cut)
                    op.pair_id = pid


classes = (ASE_OT_toggle_pair_cut, ASE_OT_set_radius_mode, ASE_PT_controls)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
