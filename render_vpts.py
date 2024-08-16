import bpy
from os.path import join
from pathlib import Path
from bpy_extras.io_utils import ExportHelper

bl_info = {
    "name": "Render vpts",
    "blender": (2, 80, 0),
    "category": "Render",
}

class RenderImageOperator(bpy.types.Operator,ExportHelper):
    bl_idname = "render.render_vpts"
    bl_label = "Render structure vpts"
    directory = bpy.props.StringProperty(
        name='folder',
        description='where to put the images',
        subtype='DIR_PATH'
    )
    #imagename: bpy.props.StringProperty(
    #    name='imagename',
    #    description="basename of the image",
    #        default='image'
    #        )
    #imagepath:bpy.props.StringProperty(
     #   name='imagepath',
     #   description="path for the image",
     #       default='',
     #       subtype='FILE_PATH'
    #        )
    #imagepath = bpy.props.StringProperty(subtype="FILE_PATH") 
    def execute(self, context):
        # Set up scene
        scene = bpy.context.scene
        cameras=self.get_camera_list()
        collections=scene.collection.children
        #print(collections)
        for camera in cameras:
            self.toggle(camera,SET=True)
        for collection in collections:
            self.toggle_collection(collection,SET=True)
        for collection in collections:
            self.toggle_collection(collection,SET=False)
            #print(collection.name)
            for camera in cameras:
                self.toggle(camera,SET=False)
                bpy.context.scene.camera = camera   
                self.RENDER(FILEPATH=str(Path(f'{join(self.directory,collection.name)}_{camera.name}.png')))
                #print(camera.name)
                self.toggle(camera,SET=True)
            self.toggle_collection(collection,SET=True)
        for camera in cameras:
            self.toggle(camera,SET=False)
        for collection in collections:
            self.toggle_collection(collection,SET=False)
        return {'FINISHED'}
    #def draw_func(self, context):
    #    layout = self.layout
    #    layout.operator("render.render_image", text="Render Image")
    def toggle(self,obj,SET=True):
        obj.hide_render = SET
        obj.hide_viewport = SET  # Optional: hide in the viewport as well
        for child in obj.children:
                child.hide_render = SET
                child.hide_viewport = SET  # Optional: hide in the viewport as well
        return(None)
    def toggle_collection(self,coll_obj,SET=True):
        coll_obj.hide_render = SET
        #coll_obj.hide_viewport = SET
        #bpy.context.view_layer.active_layer_collection = 
        for obj in coll_obj.objects:
            obj.hide_render = SET
        #     obj.hide_viewport = SET  # Optional: hide in the viewport as well
    def RENDER(self,FILEPATH='./'):
        render_settings = bpy.context.scene.render
        render_settings.filepath = FILEPATH  # Replace with desired output file path
#        render_settings.engine = 'CYCLES'  # Replace with desired render engine
        #render_settings.resolution_x = 2000  # Replace with desired resolution
        #render_settings.resolution_y = 2000
        bpy.ops.render.render(write_still=True)
    def get_camera_list(self):
        camera_list = []
        scene = bpy.context.scene
        for obj in scene.objects:
            if obj.type == 'CAMERA':
                camera_list.append(obj)
        return camera_list
    @classmethod
    def poll(cls, context):
        return context.scene is not None
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self) 
        return {'RUNNING_MODAL'} 

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "imagepath")
    
    def check(self, context):
        return self.filepath != ""

def register():
    bpy.utils.register_class(RenderImageOperator)
    bpy.types.TOPBAR_MT_render.append(menu_func)
   # bpy.ops.render.render_structure('INVOKE_DEFAULT')


def unregister():
    bpy.utils.unregister_class(RenderImageOperator)
    bpy.types.TOPBAR_MT_render.remove(menu_func)


def menu_func(self, context):
    self.layout.operator(RenderImageOperator.bl_idname)
if __name__ == "__main__":
    register()
    
