"""
Regression test for programmatic operator call: bpy.ops.import_mesh.ase(filepath=...).

Run as:
    blender --background --python test/run_in_blender.py

Exits 0 and prints "OK" on success.
"""
import sys
import os
import bpy

HERE = os.path.dirname(os.path.abspath(__file__))
XYZ = os.path.join(HERE, "fixtures", "water.xyz")

bpy.ops.preferences.addon_enable(module='blender_importASE')

# Case 1 — programmatic single-file call (the regression)
n0 = len(bpy.data.objects)
r = bpy.ops.import_mesh.ase(filepath=XYZ, representation="Balls'n'Sticks")
n1 = len(bpy.data.objects)
assert r == {'FINISHED'}, f"filepath= path: result={r}"
assert n1 > n0, f"filepath= path imported zero objects (n0={n0}, n1={n1})"

# Case 2 — directory + files (the GUI path) still works
for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)
n0 = len(bpy.data.objects)
r = bpy.ops.import_mesh.ase(
    directory=os.path.dirname(XYZ),
    files=[{"name": os.path.basename(XYZ)}],
    representation="Balls'n'Sticks",
)
n1 = len(bpy.data.objects)
assert r == {'FINISHED'}, f"files= path: result={r}"
assert n1 > n0, f"files= path imported zero objects (n0={n0}, n1={n1})"

print("OK")
