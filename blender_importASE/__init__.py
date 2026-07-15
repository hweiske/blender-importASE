#!/usr/bin/env python
import os.path
import importlib
import bpy
import sys
import subprocess
from bpy_extras.io_utils import ImportHelper, ExportHelper
from os.path import join
from importlib import util


__author__ = "Hendrik Weiske"
__credits__ = ["Franz Thiemann"]
__version__ = "2.2"
__maintainer__ = "Hendrik Weiske"
__email__ = "hendrik.weiske@uni-leipzig.de"

bl_info = {
    "name": "ASE Importer",
    "description": "Import molecules using ASE",
    "author": "Hendrik Weiske",
    "version": (2, 2),
    "blender": (4, 4, 0),
    "location": "File > Import",
    "category": "Import-Export",
}

class ImportASEMolecule(bpy.types.Operator, ImportHelper):
    bl_idname = "import_mesh.ase"
    bl_label = "Import ASE Molecule"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".*"

    scale: bpy.props.FloatProperty(
        name="Scale",
        description="scaling the atoms",
        default=0.5,
        min=0.0,
        soft_max=10,
    )

    resolution: bpy.props.IntProperty(
        name='resolution',
        description='resolution of bonds and atoms. Breaks for Balls\'n\'Sticks and Licorice with Longbonds',
        default=32,
    )

    colorbonds: bpy.props.BoolProperty(
        name='colorbonds',
        description="Color the bonds according to the connecting atoms",
        default=False,
    )
    long_bonds: bpy.props.BoolProperty(
        name='Use Longbonds',
        description="Mitigates lines in the middle of bonds where individual bonds meet,  relevant for Balls'n'Sticks and Licorice",
        default=True,
    )
    color: bpy.props.FloatProperty(
        name="color",
        description="color for gray bonds in BW-scale, relevant for Balls'n'Sticks and Licorice",
        default=0.6,
        min=0.0,
        max=1.0,
    )
    unit_cell: bpy.props.BoolProperty(
        name='unit_cell',
        description="Draw unit cell",
        default=False,
    )

    representation: bpy.props.EnumProperty(
        name="representation",
        description="select the representation for your structure",
        items=[
            ("Balls'n'Sticks", "Balls'n'Sticks", "Balls and sticks representaiton"),
            ("Licorice", "Licorice", "Licorice representation"),
            ('VDW', 'VDW', 'VDW Radii, no bonds'),
            ('3D_print', '3D print (spheres + bonds)', 'Real sphere meshes plus geometry-node bonds with icosphere joints - suited for 3D printing (was: bonds_fromnodes)'),
            ('nodes', 'nodes', 'Everything from geometrynodes. Fastest'),
        ],
        default="nodes"
    )
    read_density: bpy.props.BoolProperty(
        name='load e-density',
        description="load electron-density as volume and use a node-tree for the creation of isosurfaces (.cube files and VASP CHGCAR/PARCHG/AECCAR)",
        default=True,
    )
    zero_cell: bpy.props.BoolProperty(
        name='Assume zero centered cell',
        description="Calculations in AMS result in zero-centered cells. This will confuse the longbond algorithm and cause" \
                    "the atoms to lie outside the unit cell (if drawn). This option compensates for that. relevant for Balls'n'Sticks and Licorice",
        default=False,
    )
    animate: bpy.props.BoolProperty(
        name='animate',
        description="animate trajectory (all traj-files readable by ASE)",
        default=True,
    )
    imageslice: bpy.props.IntProperty(
        name='nth-image',
        description='when loading long trajectories it is recommended not to use all images, since that will scale poorly depending on the number of bonds in the molecule and drastically influence performance',
        default=1
    )
    overwrite: bpy.props.BoolProperty(
        name='overwrite',
        description='overwrite representation to "nodes" for animations. Recommended',
        default=True
    )
    outline: bpy.props.BoolProperty(
        name='outline',
        description='add outline modifier for atoms and bonds',
        default=True
    )
    add_supercell: bpy.props.BoolProperty(
        name='add_supercell',
        description='add supercell modifier when PBC',
        default=True
    )
    files: bpy.props.CollectionProperty(
        type=bpy.types.OperatorFileListElement,
        options={'HIDDEN', 'SKIP_SAVE'},
        description='List of files to be imported'
    )
    directory: bpy.props.StringProperty(
        name='folder',
        description='directory of file',
        subtype='DIR_PATH'
    )
    

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "resolution")
        layout.prop(self, "scale")
        layout.prop(self, 'outline')
        layout.prop(self, 'add_supercell')
        layout.prop(self, 'colorbonds')
        layout.prop(self, 'long_bonds')
        layout.prop(self, 'representation')
        layout.prop(self, 'color')
        layout.prop(self, 'unit_cell')
        layout.prop(self, 'read_density')
        layout.prop(self, 'zero_cell')
        layout.prop(self, 'animate')
        layout.prop(self, 'overwrite')
        layout.prop(self,'imageslice')

    def execute(self, context):
        # When invoked from the GUI file dialog, ImportHelper populates
        # self.files + self.directory.  When called programmatically with
        # bpy.ops.import_mesh.ase(filepath=...) — the documented single-file
        # form — only self.filepath is populated and self.files stays empty.
        # Normalise both into a (directory, [name, ...]) pair before the loop.
        if self.files:
            directory = self.directory
            names = [f.name for f in self.files]
        elif self.filepath:
            directory, name = os.path.split(self.filepath)
            names = [name]
        else:
            self.report({'ERROR'}, "No filepath or files provided")
            return {'CANCELLED'}

        from .ui import import_ase_molecule
        for name in names:
            filepath = join(directory, name)
            try:
                import_ase_molecule(
                    filepath, name,
                    resolution=self.resolution,
                    color=self.color, colorbonds=self.colorbonds,
                    long_bonds=self.long_bonds, scale=self.scale,
                    unit_cell=self.unit_cell, representation=self.representation,
                    read_density=self.read_density,
                    shift_cell=self.zero_cell, imageslice=self.imageslice,
                    animate=self.animate, outline=self.outline,
                    overwrite=self.overwrite, add_supercell=self.add_supercell,
                )
            except ValueError as exc:
                self.report({'ERROR'}, str(exc))
                return {'CANCELLED'}
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class ImportASEPolyhedra(bpy.types.Operator, ImportHelper):
    """Import a structure with coordination polyhedra: the convex hull of
    every coordination shell is drawn as solid faces"""
    bl_idname = "import_mesh.ase_polyhedra"
    bl_label = "Import ASE Polyhedra"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".*"

    expand_cutoff: bpy.props.FloatProperty(
        name="expansion cutoff",
        description="covalent-radius multiplier used to pull in periodic neighbor images so polyhedra at the cell boundary are closed",
        default=1.2,
        min=0.5,
        soft_max=2.0,
    )
    trim_cutoff: bpy.props.FloatProperty(
        name="trim cutoff",
        description="covalent-radius multiplier below which expanded atoms without neighbors are removed again",
        default=1.0,
        min=0.5,
        soft_max=2.0,
    )
    poly_cutoff: bpy.props.FloatProperty(
        name="polyhedra cutoff",
        description="covalent-radius multiplier defining the neighbor shell that forms a polyhedron",
        default=1.1,
        min=0.5,
        soft_max=2.0,
    )
    min_neighbors: bpy.props.IntProperty(
        name="min neighbors",
        description="minimum number of neighbors an atom needs to get a coordination polyhedron",
        default=4,
        min=4,
    )
    include_hydrogen: bpy.props.BoolProperty(
        name="include hydrogen",
        description="also use hydrogen atoms as polyhedra centers and corners",
        default=False,
    )
    single_element_corners: bpy.props.BoolProperty(
        name="single-element corners",
        description="restrict each polyhedron to corner atoms of a single element "
                    "(the coordinating counter-ion), so ionic solids like NaCl "
                    "render as clean NaCl6 / ClNa6 octahedra instead of hulls that "
                    "swallow next-nearest same-element atoms. Turn off for "
                    "same-element clusters such as B6",
        default=True,
    )
    resolution: bpy.props.IntProperty(
        name='resolution',
        description='resolution of bonds and atoms',
        default=16,
    )
    colorbonds: bpy.props.BoolProperty(
        name='colorbonds',
        description="Color the bonds according to the connecting atoms",
        default=True,
    )
    bond_distance: bpy.props.FloatProperty(
        name="bond distance",
        description="bond distance criterion passed to the atoms_and_bonds node group",
        default=0.66,
        min=0.0,
        soft_max=2.0,
    )
    bond_radius: bpy.props.FloatProperty(
        name="bond radius",
        description="bond cylinder radius",
        default=0.1,
        min=0.0,
        soft_max=1.0,
    )
    outline: bpy.props.BoolProperty(
        name='outline',
        description='add outline modifier to the atoms and bonds (the polyhedra faces stay outline-free)',
        default=True,
    )
    files: bpy.props.CollectionProperty(
        type=bpy.types.OperatorFileListElement,
        options={'HIDDEN', 'SKIP_SAVE'},
        description='List of files to be imported'
    )
    directory: bpy.props.StringProperty(
        name='folder',
        description='directory of file',
        subtype='DIR_PATH'
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'expand_cutoff')
        layout.prop(self, 'trim_cutoff')
        layout.prop(self, 'poly_cutoff')
        layout.prop(self, 'min_neighbors')
        layout.prop(self, 'include_hydrogen')
        layout.prop(self, 'single_element_corners')
        layout.prop(self, 'resolution')
        layout.prop(self, 'colorbonds')
        layout.prop(self, 'bond_distance')
        layout.prop(self, 'bond_radius')
        layout.prop(self, 'outline')

    def execute(self, context):
        if self.files:
            directory = self.directory
            names = [f.name for f in self.files]
        elif self.filepath:
            directory, name = os.path.split(self.filepath)
            names = [name]
        else:
            self.report({'ERROR'}, "No filepath or files provided")
            return {'CANCELLED'}

        from .polyhedra import import_polyhedra
        for name in names:
            import_polyhedra(
                join(directory, name), name,
                expand_cutoff=self.expand_cutoff,
                trim_cutoff=self.trim_cutoff,
                poly_cutoff=self.poly_cutoff,
                min_neighbors=self.min_neighbors,
                include_hydrogen=self.include_hydrogen,
                single_element_corners=self.single_element_corners,
                resolution=self.resolution,
                colorbonds=self.colorbonds,
                bond_distance=self.bond_distance,
                bond_radius=self.bond_radius,
                outline=self.outline,
            )
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


# dynamic EnumProperty items must stay referenced from python or Blender
# shows garbage strings - module-level caches hold them alive
_density_file_items_cache = []
_csv_file_items_cache = []


def _sibling_file_items(operator, cache, match):
    """Enum items: files next to the one selected in the import browser."""
    directory = operator.directory or os.path.dirname(operator.filepath)
    items = [('NONE', '(none)', 'no file selected')]
    try:
        current = {f.name for f in operator.files} | {os.path.basename(operator.filepath)}
        for fname in sorted(os.listdir(directory)):
            if fname in current:
                continue
            if match(fname):
                items.append((fname, fname, os.path.join(directory, fname)))
    except OSError:
        pass
    cache[:] = items
    return cache


class ImportASEDensityMesh(bpy.types.Operator, ImportHelper):
    """Import the +/- isosurfaces of a density file as a real mesh
    (marching cubes), optionally colored by a second density file"""
    bl_idname = "import_mesh.ase_density_mesh"
    bl_label = "Import ASE Density as Mesh"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".*"

    iso_value: bpy.props.FloatProperty(
        name="isovalue",
        description="isosurface level; both +isovalue and -isovalue surfaces are generated when present in the data",
        default=0.03,
        precision=4,
        soft_min=0.0001,
        soft_max=10.0,
    )
    color_choice: bpy.props.EnumProperty(
        name="color density",
        description="second density file from the same folder; its values are sampled on the isosurface and drive the color ramp",
        items=lambda self, context: _sibling_file_items(
            self, _density_file_items_cache,
            lambda f: f.lower().endswith('.cube') or f.upper().startswith(('CHGCAR', 'CHG', 'PARCHG', 'AECCAR'))),
    )
    color_min: bpy.props.FloatProperty(
        name="color min",
        description="value of the color density mapped to the low end of the color ramp; leave min = max for automatic normalization to the sampled range",
        default=0.0,
        precision=4,
    )
    color_max: bpy.props.FloatProperty(
        name="color max",
        description="value of the color density mapped to the high end of the color ramp; values outside the range are clamped",
        default=0.0,
        precision=4,
    )
    sample_interior: bpy.props.BoolProperty(
        name="sample interior",
        description="color each surface point with the strongest (largest magnitude) color-density value found along the surface normal through the whole volume, instead of the value directly on the surface - projects buried features onto the isosurface",
        default=False,
    )
    shade_smooth: bpy.props.BoolProperty(
        name="shade smooth",
        description="smooth-shade the isosurface",
        default=True,
    )
    import_atoms: bpy.props.BoolProperty(
        name="import atoms",
        description="also import the structure from the density file as the nodes representation",
        default=True,
    )
    preset: bpy.props.EnumProperty(
        name="shader preset",
        description="initial color ramp of the generated isosurface material (one material per preset; edits survive re-imports)",
        items=[
            ('DEFAULT', 'red-white-blue', 'soft red to white to blue ramp'),
            ('ELSTAT', 'elstat. potential', 'pure blue to white to red (electrostatic potential map)'),
            ('LED', 'LED', 'red, green, blue at ramp positions 0.8, 0.9, 1.0'),
        ],
        default='DEFAULT',
    )
    files: bpy.props.CollectionProperty(
        type=bpy.types.OperatorFileListElement,
        options={'HIDDEN', 'SKIP_SAVE'},
        description='List of files to be imported'
    )
    directory: bpy.props.StringProperty(
        name='folder',
        description='directory of file',
        subtype='DIR_PATH'
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'iso_value')
        layout.prop(self, 'color_choice')
        row = layout.row(align=True)
        row.prop(self, 'color_min')
        row.prop(self, 'color_max')
        layout.prop(self, 'sample_interior')
        layout.prop(self, 'preset')
        layout.prop(self, 'import_atoms')
        layout.prop(self, 'shade_smooth')

    def execute(self, context):
        if self.files:
            directory = self.directory
            names = [f.name for f in self.files]
        elif self.filepath:
            directory, name = os.path.split(self.filepath)
            names = [name]
        else:
            self.report({'ERROR'}, "No filepath or files provided")
            return {'CANCELLED'}

        from .density_mesh import import_density_mesh
        if self.color_choice and self.color_choice != 'NONE':
            color_filepath = join(directory, self.color_choice)
        else:
            color_filepath = None
        for name in names:
            try:
                import_density_mesh(
                    join(directory, name), name,
                    color_filepath=color_filepath,
                    iso_value=self.iso_value,
                    shade_smooth=self.shade_smooth,
                    preset=self.preset,
                    import_atoms=self.import_atoms,
                    color_min=self.color_min,
                    color_max=self.color_max,
                    sample_interior=self.sample_interior,
                )
            except ValueError as exc:
                self.report({'ERROR'}, str(exc))
                return {'CANCELLED'}
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class ImportASECharges(bpy.types.Operator, ImportHelper):
    """Import a structure with per-atom partial charges from a csv file
    (one charge per atom, same order as the structure file); atoms and
    bonds can be colored by charge via the 'charge_colors' switch"""
    bl_idname = "import_mesh.ase_charges"
    bl_label = "Import ASE Charges"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".*"

    charge_choice: bpy.props.EnumProperty(
        name="charges csv",
        description="csv file from the same folder with one partial charge per atom, in the same order as the atoms in the structure file",
        items=lambda self, context: _sibling_file_items(
            self, _csv_file_items_cache,
            lambda f: f.lower().endswith(('.csv', '.txt', '.dat'))),
    )
    resolution: bpy.props.IntProperty(
        name='resolution',
        description='resolution of bonds and atoms',
        default=16,
    )
    colorbonds: bpy.props.BoolProperty(
        name='colorbonds',
        description="Color the bonds according to the connecting atoms (used when charge colors are switched off)",
        default=True,
    )
    bond_distance: bpy.props.FloatProperty(
        name="bond distance",
        description="bond distance criterion passed to the atoms_and_bonds node group",
        default=0.66,
        min=0.0,
        soft_max=2.0,
    )
    bond_radius: bpy.props.FloatProperty(
        name="bond radius",
        description="bond cylinder radius",
        default=0.1,
        min=0.0,
        soft_max=1.0,
    )
    outline: bpy.props.BoolProperty(
        name='outline',
        description='add outline modifier to the atoms and bonds',
        default=True,
    )
    files: bpy.props.CollectionProperty(
        type=bpy.types.OperatorFileListElement,
        options={'HIDDEN', 'SKIP_SAVE'},
        description='List of files to be imported'
    )
    directory: bpy.props.StringProperty(
        name='folder',
        description='directory of file',
        subtype='DIR_PATH'
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'charge_choice')
        layout.prop(self, 'resolution')
        layout.prop(self, 'colorbonds')
        layout.prop(self, 'bond_distance')
        layout.prop(self, 'bond_radius')
        layout.prop(self, 'outline')

    def execute(self, context):
        if self.files:
            directory = self.directory
            names = [f.name for f in self.files]
        elif self.filepath:
            directory, name = os.path.split(self.filepath)
            names = [name]
        else:
            self.report({'ERROR'}, "No filepath or files provided")
            return {'CANCELLED'}
        if self.charge_choice and self.charge_choice != 'NONE':
            charge_filepath = join(directory, self.charge_choice)
        else:
            self.report({'ERROR'}, "No charges csv file selected")
            return {'CANCELLED'}

        from .charges import import_charges
        for name in names:
            try:
                import_charges(
                    join(directory, name), name,
                    charge_filepath=charge_filepath,
                    resolution=self.resolution,
                    colorbonds=self.colorbonds,
                    bond_distance=self.bond_distance,
                    bond_radius=self.bond_radius,
                    outline=self.outline,
                )
            except ValueError as exc:
                self.report({'ERROR'}, str(exc))
                return {'CANCELLED'}
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class ExportASEXyz(bpy.types.Operator, ExportHelper):
    """Export the active nodes-representation structure to a .xyz file
    (vertex positions with their stored element numbers as symbols)"""
    bl_idname = "export_mesh.ase_xyz"
    bl_label = "Export ASE xyz"

    filename_ext = ".xyz"
    filter_glob: bpy.props.StringProperty(default="*.xyz", options={'HIDDEN'})

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj is not None and obj.type == 'MESH'
                and 'element' in obj.data.attributes)

    def execute(self, context):
        from .exports import export_xyz
        try:
            n = export_xyz(context.active_object, self.filepath)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, f'wrote {n} atoms to {self.filepath}')
        return {'FINISHED'}


class ExportASE3DPrint(bpy.types.Operator, ExportHelper):
    """Export the active structure's collection for 3D printing: one STL
    per element (atoms joined), the bonds, and resin supports, zipped
    into a single archive"""
    bl_idname = "export_mesh.ase_3dprint"
    bl_label = "Export ASE 3D print"

    filename_ext = ".zip"
    filter_glob: bpy.props.StringProperty(default="*.zip", options={'HIDDEN'})

    generate_supports: bpy.props.BoolProperty(
        name="generate supports",
        description="generate supports (with the parameters below) only if none exist yet; existing supports - e.g. from 'Rebuild 3D-print supports' in the ASE panel - are always exported as-is. Turn off to export without supports when none exist",
        default=True,
    )
    base_radius: bpy.props.FloatProperty(
        name="base radius",
        description="pillar radius at the plate (structure units, i.e. Angstrom)",
        default=0.25,
        min=0.01,
        soft_max=1.0,
    )
    tip_radius: bpy.props.FloatProperty(
        name="contact radius",
        description="pillar radius at the atom contact point",
        default=0.1,
        min=0.01,
        soft_max=1.0,
    )
    support_layer: bpy.props.FloatProperty(
        name="support drop",
        description="minimum height a bonded/touching neighbor must sit below an atom to hold it up; atoms without such a lower neighbor (islands, horizontal or upward-only branches) get their own pillar. Larger values add more pillars",
        default=0.3,
        min=0.0,
        soft_max=5.0,
    )
    plate_thickness: bpy.props.FloatProperty(
        name="plate thickness",
        description="thickness of the base plate",
        default=0.6,
        min=0.1,
        soft_max=3.0,
    )
    plate_holes: bpy.props.BoolProperty(
        name="plate holes",
        description="punch a regular grid of square holes into the base plate to save material (pillars landing on a hole are moved onto material)",
        default=True,
    )
    plate_gap: bpy.props.FloatProperty(
        name="plate gap",
        description="distance from the base plate up to the lowest atom (pillar length below the model)",
        default=2.0,
        min=0.0,
        soft_max=10.0,
    )

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'generate_supports')
        layout.prop(self, 'base_radius')
        layout.prop(self, 'tip_radius')
        layout.prop(self, 'support_layer')
        layout.prop(self, 'plate_thickness')
        layout.prop(self, 'plate_holes')
        layout.prop(self, 'plate_gap')

    def execute(self, context):
        from .exports import export_3dprint
        try:
            files = export_3dprint(context, self.filepath,
                                   generate_supports=self.generate_supports,
                                   base_radius=self.base_radius,
                                   tip_radius=self.tip_radius,
                                   support_layer=self.support_layer,
                                   plate_thickness=self.plate_thickness,
                                   plate_holes=self.plate_holes,
                                   plate_gap=self.plate_gap)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, f'wrote {", ".join(files)} to {self.filepath}')
        return {'FINISHED'}


class ASEAddonPreferences(bpy.types.AddonPreferences):
    bl_idname = __name__

    install_failed: bpy.props.BoolProperty(default=False)

    def draw(self, context):
        layout = self.layout
        if self.install_failed:
            layout.label(text="ASE installation failed. Please check your internet connection.", icon='ERROR')
        else:
            layout.label(text="ASE installation successful.", icon='CHECKMARK')


# (import name, pip name, required). ase is required - the operators are only
# registered when it is available. scipy (polyhedra) and scikit-image (density
# mesh) back individual features; they are installed up front like ase so the
# first use of those importers does not stall, but a failure to install them
# only disables their feature rather than the whole add-on.
DEPENDENCIES = [
    ("ase", "ase", True),
    ("scipy", "scipy", False),
    ("skimage", "scikit-image", False),
]


def _install_package(import_name, pip_name):
    """pip install pip_name into Blender's user modules path and return whether
    import_name is importable afterwards."""
    install_path = os.path.join(bpy.utils.script_path_user(), "modules")
    subprocess.check_call([sys.executable, "-m", "pip", "install",
                           "--target", install_path, pip_name])
    if install_path not in sys.path:
        sys.path.append(install_path)
    importlib.invalidate_caches()
    return util.find_spec(import_name) is not None


def check_dependency():
    """Make sure ase and the optional feature dependencies (scipy,
    scikit-image) are importable, installing any that are missing into
    Blender's user modules path. Returns True when the required ase is
    available; optional dependencies only print a warning on failure."""
    ase_available = True
    for import_name, pip_name, required in DEPENDENCIES:
        if util.find_spec(import_name) is not None:
            continue
        print(f"{pip_name} not present in Blender python. Attempting install. "
              "This could take a moment...")
        try:
            installed = _install_package(import_name, pip_name)
        except subprocess.CalledProcessError:
            installed = False
        if installed:
            print(f"Installed {pip_name}")
        elif required:
            print(f"Failed to install {pip_name}. Please check your internet "
                  "connection and try again or install manually.")
            ase_available = False
        else:
            print(f"Could not install {pip_name}; the feature that needs it "
                  "will be unavailable until it is installed.")
    return ase_available

def menu_func_import(self, context):
    self.layout.operator(ImportASEMolecule.bl_idname, text="ASE Molecule (.*)")
    self.layout.operator(ImportASEPolyhedra.bl_idname, text="ASE Polyhedra (.*)")
    self.layout.operator(ImportASEDensityMesh.bl_idname, text="ASE Density as Mesh (.*)")
    self.layout.operator(ImportASECharges.bl_idname, text="ASE Charges (.*)")

def menu_func_export(self, context):
    self.layout.operator(ExportASEXyz.bl_idname, text="ASE xyz (.xyz)")
    self.layout.operator(ExportASE3DPrint.bl_idname, text="ASE 3D print (.zip)")

def register():
    bpy.utils.register_class(ASEAddonPreferences)
    dependency = check_dependency()
    if dependency:
        bpy.utils.register_class(ImportASEMolecule)
        bpy.utils.register_class(ImportASEPolyhedra)
        bpy.utils.register_class(ImportASEDensityMesh)
        bpy.utils.register_class(ImportASECharges)
        bpy.utils.register_class(ExportASEXyz)
        bpy.utils.register_class(ExportASE3DPrint)
        bpy.types.TOPBAR_MT_file_import.append(menu_func_import)
        bpy.types.TOPBAR_MT_file_export.append(menu_func_export)
        # deferred so the addon can load (and show its preferences) when
        # ase is not installed yet - controls imports ase at module level
        from . import controls
        controls.register()
    else:
        prefs = bpy.context.preferences.addons[__name__].preferences
        prefs.install_failed = True

def unregister():
    try:
        from . import controls
        controls.unregister()
    except Exception:
        print("ASE controls were not registered, skipping.")
    for cls in (ImportASEMolecule, ImportASEPolyhedra, ImportASEDensityMesh,
                ImportASECharges, ExportASEXyz, ExportASE3DPrint):
        try:
            bpy.utils.unregister_class(cls)
        except RuntimeError:
            print(f"{cls.__name__} was not registered, skipping.")
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)
    bpy.utils.unregister_class(ASEAddonPreferences)


if __name__ == "__main__":
    register()
