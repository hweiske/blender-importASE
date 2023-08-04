#!/usr/bin/env python

__author__ = "Hendrik Weiske"
__credits__ = ["Tassem El-Sayed"]
__version__ = "1.0.0"
__maintainer__ = "Hendrik Weiske"
__email__ = "hendrik.weiske@uni-leipzig.de"

bl_info = {
    "name": "ASE Importer",
    "description": "Import molecules using ASE",
    "author": "Your Name",
    "version": (1, 0),
    "blender": (2, 92, 0),
    "location": "File > Import",
    "category": "Import-Export",
}

import bpy
from bpy_extras.io_utils import ImportHelper
from ase import io
import ase
from ase.data import covalent_radii,colors
from ase.build import make_supercell
import numpy as np
from os.path import join
class ImportASEMolecule(bpy.types.Operator, ImportHelper):
    bl_idname = "import_mesh.ase"
    bl_label = "Import ASE Molecule"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".*"
   
    supercell1: bpy.props.IntVectorProperty(
        name="Supercell-vector",
        size=9,
        default=[1, 0, 0, 0,1,0, 0,0,1],
    )
    scale: bpy.props.FloatProperty(
        name="Scale",
        description="scaling th atoms",
        default=1.0,
        min=0.0,
        soft_max=10,
    )
    colorbonds: bpy.props.BoolProperty(
        name='colorbonds',
        description="Color the bonds according to the surrounding atoms",
        default=False,
            )
    color: bpy.props.FloatProperty(
        name="color",
        description="color for gray bonds in BW-scale",
        default=0.2,
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
            ("Balls'n'Sticks","Balls'n'Sticks", "Balls and sticks representaiton"),
            ("Licorice","Licorice", "Licorice representation"),
            ('VDW','VDW','VDW Radii, no bonds'),
        ],
        default="Balls'n'Sticks"
    )
    files: bpy.props.CollectionProperty(
        type=bpy.types.OperatorFileListElement,
        options={'HIDDEN','SKIP_SAVE'},
        description='List of files to be imported'
    )
    directory:bpy.props.StringProperty(
        name='folder',
        description='where to put the images',
        subtype='DIR_PATH'
    )
    def draw(self, context):
        layout = self.layout
        box= layout.box()
        box.label(text="")
        box.separator()
        box.operator("import_scene.my_format", text="Import")
        for i in range(3):
            row = box.row(align=True)
            for j in range(3):
                row.prop(self, "supercell1", index=i * 3 + j, emboss=False, slider=True)
        layout.prop(self,"scale")
        layout.prop(self,'colorbonds')
        layout.prop(self,'representation')
        layout.prop(self,'color')
        layout.prop(self,'unit_cell')
        layout.prop(self,'separate_collections')
    def execute(self, context):
        for file in self.files:
            filepath = join(self.directory,file.name)
            matrix = np.array(self.supercell1).reshape((3, 3))
            print(matrix)
            import_ase_molecule(filepath,file.name,matrix,
                                color=self.color,colorbonds=self.colorbonds,scale=self.scale,
                                unit_cell=self.unit_cell,representation=self.representation,
                                separate_collections=self.separate_collections)
        return {"FINISHED"}
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self) 
        return {'RUNNING_MODAL'} 

import ase.neighborlist

def draw_atoms(atoms,scale=1,representation="Balls'n'Sticks"):
    cnt = 0
    #bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0),segments = 16 ,ring_count = 16)
    bpy.ops.object.shade_smooth()
    sphere = bpy.context.object
    sphere.name = 'ref_sphere'
    for atom in atoms:
        ob = sphere.copy()
        ob.data = sphere.data.copy()
        ob.location = atom.position
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob) 
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].name = atom.symbol
        if representation == "Balls'n'Sticks":
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [covalent_radii[atom.number]*0.5*scale,]*3
        elif representation == 'Licorice':
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [0.1]*3
        else:
            bpy.context.view_layer.active_layer_collection.collection.objects[-1].scale = [covalent_radii[atom.number]*scale,]*3
        print(bpy.data.node_groups)
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials[atom.symbol])
        cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_sphere'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None

def draw_bonds(atoms):
    nl = ase.neighborlist.NeighborList([covalent_radii[atomic_number]*0.9  for atomic_number in atoms.numbers], self_interaction = False,bothways=True)
    nl.update(atoms)
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='bonds')
    except:
        None
#    bpy.ops.surface.primitive_nurbs_surface_cylinder_add(radius=1.0, enter_editmode=False, align='WORLD', location=(0.0, 0.0, 0.0), rotation=(0.0, 0.0, 0.0), scale=(0.0, 0.0, 0.0))
    bpy.ops.mesh.primitive_cylinder_add(vertices=16)
    bpy.ops.object.shade_smooth()
    bond = bpy.context.object
    bond.name = 'ref_bond'
    cnt = 0
    for atom in atoms:
        if nl.get_neighbors(atom.index)[0].size > 0:
            neighbors,offsets=nl.get_neighbors(atom.index)
            for neighbor,offset in zip(neighbors, offsets):
                displacements =[ 0.5*(atoms.positions[neighbor] - atom.position + np.dot(offset, atoms.cell)), 0.5*( atom.position - np.dot(offset, atoms.cell) - atoms.positions[neighbor] )]
                for n,displacement in enumerate(displacements):
                    if n == 0:
                        location= atom.position + (displacement/2)
    #                else:
     #                   location=atoms[neighbor].position + (displacement/2)
                    distance = atoms.get_distance(atom.index,neighbor,mic=True)/2
                    ob = bond.copy()
                    ob.data = bond.data.copy()
                    bpy.context.view_layer.active_layer_collection.collection.objects.link(ob) 
                    ob.name = f'{atom.symbol}{atom.index}-{atoms[neighbor].symbol}{neighbor}'
                    bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials[f'{atom.symbol}-bond'])
                    ob.location = location
                    ob.scale = (1,1,1)
                    #if n == 0:
                    ob.dimensions = (0.2,0.2,distance)
                    phi = np.arctan2(displacement[1], displacement[0]) 
                    theta = np.arccos(displacement[2] / distance) 
                    ob.rotation_euler[1] = theta
                    ob.rotation_euler[2] = phi
                    break
                cnt += 1
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_bond'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None
def draw_unit_cell(atoms):
    bpy.ops.object.select_all(action='DESELECT')
    try:
        bpy.ops.group.create(name='cell')
    except:
        None
    
    #SETUP MATERIAL
    matu=bpy.data.materials.new(name = 'unit_cell')
    matu.use_nodes=True
    tu=matu.node_tree
    su=tu.nodes['Principled BSDF']
    COL=(0.1,0.1,0.1,1)
    su.inputs[0].default_value=COL
    bpy.ops.mesh.primitive_cylinder_add(vertices=16)
    bpy.ops.object.shade_smooth()
    cell = bpy.context.object
    cell.name = 'ref_cell'
    cnt = 0
    X=[0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0]
    Y=[0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,0]
    Z=[0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1]
    for n in range(1,len(X)):
        pos1=np.array([X[n-1],Y[n-1],Z[n-1]])
        pos2=np.array([X[n],Y[n],Z[n]])
        location1=np.dot(pos1,atoms.cell)
        location2=np.dot(pos2,atoms.cell)
        print(n,location1,location2,pos1,pos2)
        displacement=location2-location1
        distance=np.linalg.norm(location1-location2)
        ob=cell.copy()
        ob.data=cell.data.copy()
        bpy.context.view_layer.active_layer_collection.collection.objects.link(ob)
        ob.name = f'unitcell-cylinder'
        bpy.context.view_layer.active_layer_collection.collection.objects[-1].data.materials.append(bpy.data.materials['unit_cell'])
        ob.location=location1+(displacement/2)
        ob.scale=(1,1,1)
        ob.dimensions=(0.1,0.1,distance)
        phi=np.arctan2(displacement[1],displacement[0])
        print(displacement,distance)
        theta=np.arccos(displacement[2]/distance)
        ob.rotation_euler[1]=theta
        ob.rotation_euler[2]=phi
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects['ref_cell'].select_set(True)
    bpy.ops.object.delete()
    bpy.ops.object.select_all(action='DESELECT')
    return None
def group_atoms(atoms):
    atom_types = set(atoms.get_chemical_symbols())
    bpy.ops.object.select_all(action='DESELECT')
    for atom_type in atom_types:
        try:
            bpy.ops.group.create(name=atom_type)
        except:
            None
    return None
    
def setup_materials(atoms,colorbonds=False,color=0.2):
    color_dict={  'H'     :(      1.00, 1.00, 1.00        ),
        'C'     :(      0.05, 0.05, 0.05        ),
        'Si'    :(      0.02, 0.38, 0.67        ),
        'Ge'    :(      0.05,0.45,0.45  ),
        'Ga'    :(      0.33, 0.71, 0.09        ),
        'In'    :(     0,0,0            ),
        'N'     :(      0.00, 0.00, 1.00        ),
        'P'     :(      1.00, 0.50, 0.00        ),
        'As'    :(      0.75, 0.54, 0.00        ),
        'Sb'    :(      0.74, 0.46, 0.17        ),
        'Bi'    :(      0.82, 0.71, 0.55        ),
        'O'     :(      1.00, 0.00, 0.00        ),
        'S'     :(      1.00, 1.00, 0.00        ),
        'F'     :(      0.00, 1.00, 0.00        ),
        'Cl'    :(      0.50, 1.00, 0.00        ),
        'Br'    :(      0.39, 0.15, 0.03        ),
        'I'     :(      1.00, 0.00, 1.00        ),
        'Ti'    :(      0.25, 1.75, 0.75        )}
    atom_types = set(atoms.get_chemical_symbols())
    atom_n=list(set(atoms.numbers))
    for n,atom_type in enumerate(atom_types):
            matat=bpy.data.materials.new(name = str(atom_type))
            matb=bpy.data.materials.new(name = f'{atom_type}-bond')
            matb.use_nodes=True
            matat.use_nodes=True
            ta=matat.node_tree
            tb=matb.node_tree
            sa=ta.nodes['Principled BSDF']
            sb=tb.nodes['Principled BSDF']
            if atom_type in color_dict:
                COL=list(color_dict[atom_type]) + [1]
            else:
                COL=list(colors.jmol_colors[atom_n[n]]) + [1]
#            bpy.data.materials[str(atom_type)].diffuse_color = COL
            sa.inputs[0].default_value=COL
            if colorbonds == False:
 #               bpy.data.materials[f'{atom_type}-bond'].diffuse_color = [color,color,color,1]
                sb.inputs[0].default_value=[color,color,color,1]
            else:
                sb.inputs[0].default_value=COL
  #              bpy.data.materials[f'{atom_type}-bond'].diffuse_color = COL
    #        bpy.data.materials[str(atom_type)].metallic = 0.2
   #         bpy.data.materials[f'{atom_type}-bond'].metallic = 0.2

            
    return None 
def import_ase_molecule(filepath,filename,matrix,colorbonds=False,color=0.2,scale=1,unit_cell=False, representation="Balls'n'Sticks",separate_collections=False):
    atoms = ase.io.read(filepath)
    atoms=make_supercell(atoms,matrix)
    if separate_collections == True:
        my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_atoms')
    else:
        my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0])
    bpy.context.scene.collection.children.link(my_coll)
    layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
    bpy.context.view_layer.active_layer_collection = layer_collection
    setup_materials(atoms,colorbonds=colorbonds,color=color)
    group_atoms(atoms)
    draw_atoms(atoms,scale=scale,representation=representation)
    if representation != 'VDW':
        if separate_collections == True:
            my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_bonds')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        draw_bonds(atoms)
    if unit_cell == True and atoms.pbc.all() != False:
        if separate_collections == True:
            my_coll=bpy.data.collections.new(name=atoms.get_chemical_formula()+'_'+filename.split('.')[0]+'_cell')
            bpy.context.scene.collection.children.link(my_coll)
            layer_collection = bpy.context.view_layer.layer_collection.children[my_coll.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
        draw_unit_cell(atoms)

    
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

