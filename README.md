# Collection of Blender Addons for Molecular Structures

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10776696.svg)](https://doi.org/10.5281/zenodo.10776696)

## Dependencies

* ASE: You just need to find the location of your blender installation and use pip to install ase for blender to find it. For example when installed with snap:

  `/snap/blender/xxxx/3.x/python/bin/python3.x -m pip install ase`

  * to find your blender-python for installing ase: open blender; open the python console, type "import sys;sys.executable" and hit enter. now the executable is printed - use as above

## Installation

To use the addons in Blender simply download the zip file for yor version `blender_importASE.zip` from the latest release. In Blender go to edit -> preferences -> addons; click install; find the zip file and install it. Then activate the new addon in the list. If you want to use the automatic rendering of viewpoints, also download the file `render_vpts.py` and install and activate the same way.

### Developement Install

Symlink the `render_vpts.py` file and the `blender_importASE` folder into your addon directory (by default under linux `~/.config/blender/x.x/scripts/addons`).

## Usage

You can now import molecules from the File -> import tab and use render -> render vpts to render all collections seperately for your list of cameras.

Images will be put in the folder with the collection name and the name of the camera (name them top, side, front. camera.001 and camera.002 won't help you understand it).
