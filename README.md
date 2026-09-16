# Collection of examples

[![CI](https://github.com/hweiske/blender-importASE/actions/workflows/python-app.yml/badge.svg)](https://github.com/hweiske/blender-importASE/actions/workflows/python-app.yml)
[![Latest release](https://img.shields.io/github/v/release/hweiske/blender-importASE)](https://github.com/hweiske/blender-importASE/releases/latest)
[![Blender](https://img.shields.io/badge/blender-4.4%2B-orange?logo=blender&logoColor=white)](https://www.blender.org/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10776696.svg)](https://doi.org/10.5281/zenodo.10776696)
[![Downloads](https://img.shields.io/github/downloads/hweiske/blender-importASE/total)](https://github.com/hweiske/blender-importASE/releases)

Import molecules, crystals, trajectories, and volumetric data (electron densities, molecular orbitals) into Blender through [ASE](https://gitlab.com/ase/ase) — with geometry-nodes representations, coordination polyhedra, and isosurfaces.

[![Download Add-on](https://img.shields.io/badge/Download-blender__importASE.zip-blue?style=for-the-badge&logo=blender&logoColor=white)](https://github.com/hweiske/blender-importASE/releases/download/latest-build/blender_importASE.zip)

<table>
  <tr>
    <td align="center" width="50%">
      <img src="docs/images/molecule.jpg" alt="Molecules adsorbed on a germanium surface"/><br/>
      <b>Molecules and surfaces</b> — geometry-nodes atoms and bonds, periodic slabs grown shell by shell past the cell
    </td>
    <td align="center" width="50%">
      <img src="docs/images/polyhedra.jpg" alt="Crystal with coordination polyhedra"/><br/>
      <b>Coordination polyhedra</b> — convex hulls of coordination shells as solid faces
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/images/orbital_volume.jpg" alt="Molecular orbital isosurface"/><br/>
      <b>Molecular orbitals &amp; densities</b> — .cube / VASP volumes with node-based isosurfaces
    </td>
    <td align="center" width="50%">
      <img src="docs/images/density_mesh.jpg" alt="Density isosurface as mesh, colored by a second density"/><br/>
      <b>Density as mesh</b> — marching-cubes isosurfaces, optionally colored by a second density file
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/images/charges.jpg" alt="Molecule colored by partial charges"/><br/>
      <b>Partial charges</b> — per-atom charges from a csv file, red-white-blue on atoms and bonds
    </td>
    <td align="center" width="50%">
      <img src="docs/images/trajectory.gif" alt="Animated trajectory"/><br/>
      <b>Trajectories</b> — any ASE-readable trajectory, animated frame by frame (including varying atom counts)
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/images/adps.jpg" alt="Urea from neutron data with thermal ellipsoids"/><br/>
      <b>Thermal ellipsoids</b> — anisotropic displacement parameters from CIF or SHELX .res/.ins, with principal-axis rings
    </td>
    <td align="center" width="50%">
      <img src="docs/images/print_supports.jpg" alt="Molecule with generated resin supports"/><br/>
      <b>3D printing</b> — atoms and bonds with generated resin supports, exported as per-element STLs in one zip
    </td>
  </tr>
</table>

## Dependencies

Dependencies are automatically installed upon activation of the addon using `pip` if an internet connection is present.
In case no internet connection is available. [ASE](https://gitlab.com/ase/ase) needs to be installed manually.

### Manual dependency installation
* Use the blender scripting view to get the module directory: `bpy.utils.script_path_user() + "/modules"`
* Install ASE to the path using pip: `pip install ase --target <install_dir>
* Restart Blender
* 
## Installation

Click the **Download Add-on** button above for `blender_importASE.zip` built straight from the current `main` — rebuilt and republished on every push (see `.github/workflows/build-latest.yml`). For a stable, versioned copy instead, grab it from a [tagged release](https://github.com/hweiske/blender-importASE/releases/latest). Either way: in Blender go to edit -> preferences -> addons; click install; find the zip file and install it. Then activate the new addon in the list. Viewpoint rendering (render -> render vpts, render -> render multiple animations) is part of the addon, so there is nothing else to install.

Every rolling build carries a `BUILD_INFO.txt` naming the commit it came from, so if you ever wonder whether a download is current, open the zip and check it against the [latest commit on `main`](https://github.com/hweiske/blender-importASE/commits/main). If it is behind, your browser served a cached copy — hard-reload the download (Ctrl+Shift+R) or grab it from the [`latest-build` pre-release](https://github.com/hweiske/blender-importASE/releases/tag/latest-build) page directly.

### Developement Install

Symlink the `blender_importASE` folder into your addon directory (by default under linux `~/.config/blender/x.x/scripts/addons`).

## Usage

You can now import molecules from the File -> import tab and use render -> render vpts to render all collections seperately for your list of cameras.

Images will be put in the folder with the collection name and the name of the camera (name them top, side, front. camera.001 and camera.002 won't help you understand it).

### Electron densities

With "load e-density" enabled, volumetric data is imported as a Blender volume
with a node-based isosurface (adjustable isovalue and directional cutoffs):

* `.cube` files (Gaussian cube format)
* VASP files: `CHGCAR`, `CHG`, `PARCHG`, `AECCAR*`. For spin-polarized
  calculations a second volume with the spin difference is created, with green
  (spin-up excess) and pink (spin-down excess) isosurfaces.

Note that densities are shown in e/A^3 (ASE convention), so isovalues from
tools that use the raw CHGCAR values (e.g. VESTA) do not transfer directly.

### Thermal ellipsoids (ADPs) and SHELX files

Structures from single-crystal refinements can be drawn with their anisotropic
displacement parameters as thermal ellipsoids, ORTEP style: each atom becomes
the ellipsoid holding it with a chosen probability (50 % by default), with
black rings along its three principal sections.

* **Automatic.** The "ADPs" option is on by default in both the regular import
  (nodes representation) and the polyhedra import. It only takes effect when
  the file actually carries an anisotropic displacement table; every other
  file imports exactly as before.
* **Sources.** CIF files (`_atom_site_aniso_U_ij`, `B_ij` or `beta_ij`) and
  SHELX `.res` / `.ins` files. Atoms with only an isotropic value become
  spheres of that size; riding hydrogens in a SHELX file get their multiple of
  the parent atom's U_eq. Symmetry-generated atoms get the tensor rotated by
  the operation that generated them.
* **Hydrogens.** Riding hydrogens only carry an isotropic value, which at 50 %
  draws them larger than the atoms they sit on - so they stay the usual small
  spheres unless "hydrogen ADPs" is ticked (e.g. for neutron data with
  anisotropic hydrogens, as in the urea above).
* **Afterwards.** The modifier inputs `adps` (ellipsoids / normal spheres),
  `hydrogen_adps` and `adp_scale` (the probability level, 1.538 = 50 %) switch
  things live; the last modifier, `adp_rings`, sets the ring thickness and
  their material (a black emission). Supercell and per-element hiding work on
  the ellipsoids as on normal atoms.

SHELX `.res` / `.ins` files can be imported in general - with or without
ADPs - since ASE itself has no reader for them: cell, symmetry (LATT/SYMM),
SFAC, free variables and riding hydrogens are read, the difference-map Q peaks
of a `.res` are skipped.

### ASE panel

Imported structures get a panel in the 3D viewport sidebar (N key -> ASE tab)
with the most important settings in one place: per-element switching between
covalent and vdW radii, hiding bonds per element pair, bond distance/radius and
resolution, supercell repeats, outline thickness, per-element visibility, and
the isovalues of any imported densities.

### Custom bonds (dotted / scaled / dashed)

Select two atoms of an imported structure (edit mode, pick the two vertices)
and press "Add custom bond" in the ASE sidebar. The *bond type* dropdown in
the redo panel (F9) picks the style:

* **Dotted** - a row of spheres between the two atoms
* **Scaled** - a solid bond that gets thinner the longer it is, capped at the
  chosen radius. It is measured against the bond's natural length (the two
  atoms' covalent radii added), so a normal-length bond is full thickness and
  a stretched or partial one thins in proportion
* **Dashed** - alternating cylinder segments

Use them for partial bonds in a transition state, hydrogen bonds, or any
interaction the distance-based bond search does not draw. All three take the
bond colour (blended between the two atoms), match the structure's own bond
radius unless you set one, and get the outline. Each bond is listed in the ASE
panel of the structure, so the two atoms can be changed there afterwards. The bond
samples the atom positions live, so it follows the structure and its
trajectory.

"replace solid bond" additionally hides the normal bond between the two atoms,
so the custom one takes its place. Note this is one replacement per atom: if
you replace 0-1 and then 0-2, the 0-1 bond reappears. "Reset custom bonds"
brings every replaced solid bond back and deletes the custom bond objects
again.

### Materials

The materials used by the geometry-node representations are the ones in the
object's Material Properties tab (sorted by element, bond material last), so
you can swap or edit them there and the viewport/render follows.

### 3D printing

The "3D print" representation (formerly `bonds_fromnodes`) imports real
sphere meshes plus geometry-node bond tubes with icospheres at every atom
position ("joint radius" on the modifier), so bonds fuse into a printable
solid. File -> Export -> "ASE 3D print (.zip)" then writes one STL per
element (atoms joined), the bonds, and simple resin supports (base plate +
tapered pillars under the lowest atoms; skipped if the collection already
contains your own "supports" object) into a single zip for the slicer.

### Export to xyz

File -> Export -> "ASE xyz (.xyz)" writes the active nodes-representation
structure back to a plain xyz file, using the vertex positions (in world
coordinates, i.e. including any edits) and the stored element numbers.
