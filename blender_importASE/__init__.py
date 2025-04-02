#!/usr/bin/env python

import bpy
from bpy_extras.io_utils import ImportHelper
from os.path import join
from .ui import import_ase_molecule


__author__ = "Hendrik Weiske"
__credits__ = ["Franz Thiemann"]
__version__ = "1.3" 
__maintainer__ = "Hendrik Weiske"
__email__ = "hendrik.weiske@uni-leipzig.de"

bl_info = {
    "name": "ASE Importer",
    "description": "Import molecules using ASE",
    "author": "Hendrik Weiske",
    "version": (1, 3),
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
        description='resolution of bonds and atoms',
        default=16,
    )

    colorbonds: bpy.props.BoolProperty(
        name='colorbonds',
        description="Color the bonds according to the surrounding atoms",
        default=False,
    )
    fix_bonds: bpy.props.BoolProperty(
        name='Use Longbonds',
        description="Mitigates lines in the middle of bonds where individual bonds meet",
        default=True,
    )
    color: bpy.props.FloatProperty(
        name="color",
        description="color for gray bonds in BW-scale",
        default=0.6,
        min=0.0,
        max=1.0,
    )
    unit_cell: bpy.props.BoolProperty(
        name='unit_cell',
        description="Draw unit cell",
        default=False,
    )
    separate_collections: bpy.props.BoolProperty(
        name='separate_collections',
        description="separate collections by unit cell, atoms and bonds",
        default=False,
    )
    representation: bpy.props.EnumProperty(
        name="representation",
        description="select the representation for your structure",
        items=[
            ("Balls'n'Sticks", "Balls'n'Sticks", "Balls and sticks representaiton"),
            ("Licorice", "Licorice", "Licorice representation"),
            ('VDW', 'VDW', 'VDW Radii, no bonds'),
            ('bonds_fromnodes', 'bonds_fromnodes', 'bonds from geometrynodes'),
            ('nodes', 'nodes', 'Everything from geometrynodes. Fastest'),
        ],
        default="nodes"
    )
    read_density: bpy.props.BoolProperty(
        name='load e-density',
        description="load electron-density as volume and use a node-tree for the creation of isosurfaces (only .cube-files)",
        default=True,
    )
    zero_cell: bpy.props.BoolProperty(
        name='Assume zero centered cell',
        description="Calculations in AMS result in zero-centered cells. This will confuse the longbond algorithm and cause" \
                    "the atoms to lie outside the unit cell (if drawn). This option compensates for that.",
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
        description='overwrite representation to "nodes" for animations',
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
        layout.prop(self, 'fix_bonds')
        layout.prop(self, 'representation')
        layout.prop(self, 'color')
        layout.prop(self, 'unit_cell')
        layout.prop(self, 'separate_collections')
        layout.prop(self, 'read_density')
        layout.prop(self, 'zero_cell')
        layout.prop(self, 'animate')
        layout.prop(self, 'overwrite')
        layout.prop(self,'imageslice')

    def execute(self, context):
        
        for file in self.files:
            filepath = join(self.directory, file.name)
            # this section causes the representation to be ignored when overwrite is checked
            # we should come up with something else for now set default to false

            
            import_ase_molecule(filepath, file.name,
                            
                                resolution=self.resolution,
                                color=self.color, colorbonds=self.colorbonds, fix_bonds=self.fix_bonds, scale=self.scale,
                                unit_cell=self.unit_cell, representation=self.representation,
                                separate_collections=self.separate_collections,
                                read_density=self.read_density, 
                                shift_cell=self.zero_cell,imageslice=self.imageslice,
                                animate=self.animate, outline=self.outline,
                                overwrite=self.overwrite, add_supercell=self.add_supercell
                                )
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


def menu_func_import(self, context):
    self.layout.operator(ImportASEMolecule.bl_idname, text="ASE Molecule (.*)")


def register():
    bpy.utils.register_class(ImportASEMolecule)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportASEMolecule)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()
