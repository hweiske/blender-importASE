# Collection of Blender Addons for Molecular Structures

## Dependencies
* ASE: You just need to find the location of your blender installation and use pip to install ase for blender to find it. For example when installed with snap:

  `/snap/blender/xxxx/3.x/python/bin/python3.x -m pip install ase`

## Installation
To use the addons in Blender simply download the zip file for yor version `blender_importASE_x_x.zip` go to edit -> preferences -> addons; click install; find the zip file and install it install them. Then activate the new addon in the list. If you want to use the automatic rendering of viewpoints, also download the file `render_vpts.py` and install and activate the same way.

### Developement Install
Symlink the `render_vpts.py` file and the `new_importASE` folder into your addon directory (by default under linux `~/.config/blender/x.x/scripts/addons`).

## Usage

You can now import molecules from the File -> import tab and use render -> render vpts to render all collections seperately for your list of cameras.

Images will be put in the folder with the collection name and the name of the camera (name them top, side, front. camera.001 and camera.002 won't help you understand it).
 
