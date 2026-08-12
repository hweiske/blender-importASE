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

from .node_networks.compat import get_mod_input, mod_input_keys, mod_input_ui

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
    return get_mod_input(mod, idents[which])


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
        if self.pair_id < 0 or self.pair_id >= len(data):
            self.report({'WARNING'}, f'pair_id {self.pair_id} out of range')
            return {'CANCELLED'}
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
        if self.number < 0 or self.number >= len(data):
            self.report({'WARNING'}, f'element number {self.number} out of range')
            return {'CANCELLED'}
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
    keys = mod_input_keys(mod)
    col = layout.column(align=True)
    for item in mod.node_group.interface.items_tree:
        if getattr(item, 'in_out', None) != 'INPUT':
            continue
        if item.socket_type in SKIP_SOCKET_TYPES:
            continue
        if item.identifier not in keys:
            continue
        data, prop_path = mod_input_ui(mod, item.identifier)
        col.prop(data, prop_path, text=item.name)


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


class ASE_OT_rebuild_supports(bpy.types.Operator):
    """(Re)generate the 3D-print supports of the active structure; adjust
    the parameters in the redo panel (F9) and the supports rebuild live"""
    bl_idname = 'ase.rebuild_supports'
    bl_label = 'Rebuild 3D-print supports'
    bl_options = {'REGISTER', 'UNDO'}

    base_radius: bpy.props.FloatProperty(
        name="base radius", default=0.25, min=0.01, soft_max=1.0,
        description="pillar radius at the plate")
    tip_radius: bpy.props.FloatProperty(
        name="contact radius", default=0.1, min=0.01, soft_max=1.0,
        description="pillar radius at the atom contact point")
    support_layer: bpy.props.FloatProperty(
        name="support drop", default=0.3, min=0.0, soft_max=5.0,
        description="minimum height a bonded/touching neighbor must sit below an atom to hold it up; larger adds more pillars")
    plate_thickness: bpy.props.FloatProperty(
        name="plate thickness", default=0.6, min=0.1, soft_max=3.0)
    plate_holes: bpy.props.BoolProperty(
        name="plate holes", default=True,
        description="regular grid of square holes in the plate to save material")
    plate_gap: bpy.props.FloatProperty(
        name="plate gap", default=2.0, min=0.0, soft_max=10.0,
        description="distance from the base plate up to the lowest atom")

    def execute(self, context):
        from .exports import rebuild_supports
        try:
            rebuild_supports(context, base_radius=self.base_radius,
                             tip_radius=self.tip_radius,
                             support_layer=self.support_layer,
                             plate_thickness=self.plate_thickness,
                             plate_holes=self.plate_holes,
                             pillar_length=self.plate_gap)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


def _selected_atom_indices(obj):
    """Indices of the selected atoms (= selected vertices) of a structure,
    read from the edit-mode cage when the object is in edit mode."""
    if obj.mode == 'EDIT':
        import bmesh
        bm = bmesh.from_edit_mesh(obj.data)
        return [v.index for v in bm.verts if v.select]
    return [v.index for v in obj.data.vertices if v.select]


class ASE_OT_add_dotted_bond(bpy.types.Operator):
    """Draw a dotted bond between two atoms: select exactly two atoms
    (vertices) of the structure, then click. With 'replace solid bond' the
    normal bond between them is hidden, otherwise the dots are simply added
    - useful for partial bonds in a transition state or hydrogen bonds"""
    bl_idname = 'ase.add_dotted_bond'
    bl_label = 'Add dotted bond'
    bl_options = {'REGISTER', 'UNDO'}

    bond_type: bpy.props.EnumProperty(
        name="bond type",
        description="how the bond between the two atoms is drawn",
        items=[
            ('DOTTED', 'Dotted', 'A row of spheres between the two atoms'),
            ('SCALED', 'Scaled',
             'A solid bond that gets thinner the longer it is, up to the chosen radius'),
            ('DASHED', 'Dashed', 'Alternating cylinder segments between the two atoms'),
        ],
        default='DOTTED')
    atom_a: bpy.props.IntProperty(
        name="atom A", default=-1, min=-1,
        description="first atom index; -1 uses the current selection")
    atom_b: bpy.props.IntProperty(
        name="atom B", default=-1, min=-1,
        description="second atom index; -1 uses the current selection")
    segments: bpy.props.IntProperty(
        name="dots / dashes", default=10, min=1, soft_max=60,
        description="number of dots (dotted) or dashes (dashed); unused for scaled")
    radius: bpy.props.FloatProperty(
        name="radius", default=0.08, min=0.0, soft_max=1.0,
        description="dot radius, or the maximum bond radius for the scaled style")
    reference_length: bpy.props.FloatProperty(
        name="reference length", default=1.5, min=0.0001, soft_max=10.0,
        description="scaled style only: bonds at or below this length get the full "
                    "radius, longer ones get proportionally thinner")
    resolution: bpy.props.IntProperty(
        name="resolution", default=2, min=1, soft_max=32,
        description="icosphere subdivisions (dotted) or profile vertices "
                    "(scaled/dashed)")
    replace: bpy.props.BoolProperty(
        name="replace solid bond", default=False,
        description="also hide the normal bond between the two atoms, so the "
                    "dotted bond takes its place")
    outline: bpy.props.BoolProperty(name="outline", default=True)

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and find_ase_modifier(obj)[0] is not None

    def execute(self, context):
        obj = context.active_object
        if self.atom_a >= 0 and self.atom_b >= 0:
            index_a, index_b = self.atom_a, self.atom_b
        else:
            selected = _selected_atom_indices(obj)
            if len(selected) != 2:
                self.report({'ERROR'},
                            f'select exactly two atoms of "{obj.name}" '
                            f'(currently {len(selected)}) - enter edit mode, '
                            'pick the two vertices, then run this again')
                return {'CANCELLED'}
            index_a, index_b = selected
            # write them back so the redo panel (F9) shows and can tweak them
            self.atom_a, self.atom_b = index_a, index_b

        from .dotted_bond import add_bond
        # the dotted style's default resolution (icosphere subdivisions) is
        # far lower than a profile resolution, so only pass it on when the
        # user actually changed it away from the property default
        resolution = self.resolution if self.resolution != 2 else None
        try:
            bond = add_bond(obj, index_a, index_b, style=self.bond_type,
                            segments=self.segments, radius=self.radius,
                            resolution=resolution,
                            reference_length=self.reference_length,
                            replace=self.replace, outline=self.outline)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'},
                    f'{self.bond_type.lower()} bond {index_a}-{index_b}: {bond.name}')
        return {'FINISHED'}

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'bond_type')
        layout.prop(self, 'radius')
        if self.bond_type in {'DOTTED', 'DASHED'}:
            layout.prop(self, 'segments')
        if self.bond_type == 'SCALED':
            layout.prop(self, 'reference_length')
        layout.prop(self, 'resolution')
        layout.prop(self, 'replace')
        layout.prop(self, 'outline')
        row = layout.row(align=True)
        row.prop(self, 'atom_a')
        row.prop(self, 'atom_b')


class ASE_OT_reset_custom_bonds(bpy.types.Operator):
    """Restore every solid bond that a dotted/scaled/dashed bond replaced,
    and delete those custom bonds again"""
    bl_idname = 'ase.reset_custom_bonds'
    bl_label = 'Reset custom bonds'
    bl_options = {'REGISTER', 'UNDO'}

    remove_objects: bpy.props.BoolProperty(
        name="delete bond objects", default=True,
        description="also delete the dotted/scaled/dashed bond objects of this "
                    "structure; untick to only bring the solid bonds back")

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and find_ase_modifier(obj)[0] is not None

    def execute(self, context):
        from .dotted_bond import reset_custom_bonds
        restored, removed = reset_custom_bonds(context.active_object,
                                               remove_objects=self.remove_objects)
        self.report({'INFO'},
                    f'restored solid bonds on {restored} atom(s), '
                    f'removed {removed} bond object(s)')
        return {'FINISHED'}


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
                selected = len(_selected_atom_indices(obj))
                col = box.column(align=True)
                col.operator('ase.add_dotted_bond', icon='PARTICLES')
                col.operator('ase.reset_custom_bonds', icon='LOOP_BACK')
                if selected != 2:
                    box.label(text=f'select 2 atoms ({selected} selected)',
                              icon='INFO')

        # 3D-print supports: offered when the collection has real atom
        # meshes (the 3D print representation)
        from ase.data import chemical_symbols as _symbols
        if any(o.type == 'MESH' and o.name.split('.')[0] in _symbols
               for o in obj.users_collection[0].all_objects):
            box = self.layout.box()
            box.label(text='3D printing')
            box.operator('ase.rebuild_supports', icon='MOD_LATTICE')

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


classes = (ASE_OT_toggle_pair_cut, ASE_OT_set_radius_mode,
           ASE_OT_rebuild_supports, ASE_OT_add_dotted_bond,
           ASE_OT_reset_custom_bonds, ASE_PT_controls)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
