# blender-importASE — guide for Claude instances

This is a Blender add-on for importing atomistic structures (via [ASE](https://wiki.fysik.dtu.dk/ase/)) and turning them into publication-quality renders: molecules, crystals, coordination polyhedra, electron-density isosurfaces (volume or mesh), partial-charge colorings, and 3D-printable models. This document is the reference for driving it — both from the Blender GUI and from Python scripts. Everything the GUI does calls the same functions you can call directly, so scripting and clicking are interchangeable.

- **Package:** `blender_importASE/` (add-on version 2.5.0, min Blender 4.4; tested on 4.4, 5.1 and 5.2).
- **Dependencies:** `ase`, plus `scipy` (polyhedra), `scikit-image` (density-as-mesh), `scm.plams` (AMS TAPE41 volumes), `openvdb`/`pyopenvdb` (volumetric density) — none installed automatically; each has its own "Install" button in the add-on preferences. See [§8](#8-dependencies).

---

## 1. Two ways to drive it

**From the GUI.** *File ▸ Import* gains four entries and *File ▸ Export* gains two (see [§2](#2-operators-file--importexport)). After importing, the **N-panel ▸ ASE tab** (`ASE_PT_controls`) exposes live controls for the active structure: per-element radius mode, per-pair bond hiding, per-element colors, and the 3D-print support rebuilder ([§6](#6-live-controls-the-ase-n-panel)).

**From Python.** Every operator is a thin wrapper over a module-level function. In a headless render or a notebook:

```python
import sys; sys.path.insert(0, '/path/to/blender-importASE')
import bpy, blender_importASE
blender_importASE.register()                       # checks only - see §8, install ase first
from blender_importASE.ui import import_ase_molecule
import_ase_molecule('/data/mol.cube', 'mol.cube', representation='nodes',
                    read_density=True, outline=True, add_supercell=False)
```

Run headless with:
```
blender -b --factory-startup --python your_script.py
```

Note the calling convention shared by all importers: `(filepath, filename, ...)` — `filepath` is the full path, `filename` is the basename (used to detect file type and name the collection). Sibling files (a second density for coloring, a charge CSV) are found in the same directory as `filepath`.

---

## 2. Operators (File ▸ Import/Export)

| Menu entry | idname | Function | What it does |
|---|---|---|---|
| ASE Molecule (.*) | `import_mesh.ase` | `ui.import_ase_molecule` | Any ASE-readable structure/trajectory → atoms + bonds |
| ASE Polyhedra (.*) | `import_mesh.ase_polyhedra` | `polyhedra.import_polyhedra` | Coordination polyhedra (convex hulls of neighbor shells) |
| ASE Density as Mesh (.*) | `import_mesh.ase_density_mesh` | `density_mesh.import_density_mesh` | Marching-cubes isosurface mesh, optionally color-sampled |
| ASE Charges (.*) | `import_mesh.ase_charges` | `charges.import_charges` | Structure with per-atom charge attribute + charge coloring |
| ASE xyz (.xyz) | `export_mesh.ase_xyz` | `exports.export_xyz` | Active structure → .xyz (world-space coords) |
| ASE 3D print (.zip) | `export_mesh.ase_3dprint` | `exports.export_3dprint` | Per-element STLs + bonds + supports, zipped |

**Render ▸ Render structure vpts** (`render.render_vpts`, `render_vpts.RenderImageOperator`) renders every collection separately for every camera, writing `<collection>_<camera>.png` into a chosen folder. It is part of the add-on (it used to be a separate `render_vpts.py` you installed alongside), and needs no dependencies, so it registers even if the ASE install fails.

The importers accept multi-file selection (`files`/`directory` props). Exporters use the active object.

---

## 3. Importing molecules — `import_ase_molecule`

```python
import_ase_molecule(filepath, filename, overwrite=True, add_supercell=True,
    resolution=32, colorbonds=False, long_bonds=False, color=0.2, scale=1,
    unit_cell=False, representation="Balls'n'Sticks", read_density=True,
    shift_cell=False, imageslice=1, frame_interpolation=1, animate=True,
    outline=True, **kwargs)
```
Reads via `ase.io.read(index=':')` (VASP CHGCAR-family via `read_vasp_density`), builds a collection named `<formula>_<stem>`, dispatches by `representation`, optionally draws the unit cell and reads the density volume.

**`representation` values:**
- `"nodes"` *(default in the GUI, recommended)* — the whole structure lives in one geometry-node modifier stack on a single-vertex mesh. Fastest, supports animated trajectories and the live control tables ([§6](#6-live-controls-the-ase-n-panel)). **Use this unless you specifically need real mesh geometry.**
- `"Balls'n'Sticks"` — real sphere + cylinder meshes.
- `"Licorice"` — same path, licorice styling.
- `"VDW"` — van-der-Waals spheres only, no bonds.
- `"3D_print"` — real spheres in an `atoms` sub-collection + node bonds with icosphere joints, so everything fuses watertight for printing/export ([§7](#7-3d-printing)).

**Key options:**
- `outline=True` — adds the dark-rim outline modifier. **House style for this project: always render molecules with `outline=True`.**
- `read_density=True` — if the file carries a volume (`.cube`, CHGCAR/PARCHG/AECCAR, AMS TAPE41), builds a volumetric density object (needs openvdb). Isosurface materials `'+ material'`/`'- material'`.
- `add_supercell=True` — adds the supercell modifier when the cell is periodic; repeat counts live on `Socket_2/3/4` of that modifier ([§5](#5-the-nodes-modifier-stack-scripting-internals)).
- `animate=True`, `imageslice=n` — for trajectories, import only every *n*th image (`overwrite=True` forces `nodes`). Use this to thin out long trajectories.
- `frame_interpolation=n` — spacing of the imported images on the timeline. `1` (default) puts each image on its own frame; `10` leaves 9 empty frames between images for Blender to interpolate, turning a short path (e.g. a 6-image NEB) into a smooth animation. Note this is the opposite of `imageslice`: that one *removes* images, this one *adds* in-between frames. Caveat: on a trajectory whose atom count changes, atoms that appear/disappear slide in from their parked position across the interpolated frames (the images themselves stay exact); the importer prints a warning in that case.
- element colors, roughness and metallic come from `utils.atomcolors` (lead, say, is the violet
  `(0.0, 0.4, 0.2)` — `#00AA7C` in the picker — fully metallic at roughness 0.5); see
  [§6](#6-live-controls-the-ase-n-panel) for editing them per structure.
- `colorbonds=True` — color bond halves by their atoms; `unit_cell=True` draws the cell box as cylinders joined into one object, shaded by the `'unit_cell'` material: flat black (0,0,0) wired straight into the Surface output, so the edges read like the outline instead of catching lights. Editing that RGB node re-colors every later import; the old shaded Principled version is rebuilt on the next import.

The GUI operator passes different defaults (`scale=0.5`, `color=0.6`, `representation="nodes"`; its `zero_cell` maps to `shift_cell`).

---

## 4. The other importers

### Polyhedra — `polyhedra.import_polyhedra`
```python
import_polyhedra(filepath, filename, expand_cutoff=1.2, trim_cutoff=1.0,
    poly_cutoff=1.1, min_neighbors=4, include_hydrogen=False, resolution=16,
    colorbonds=True, bond_distance=0.66, bond_radius=0.1, outline=False,
    single_element_corners=True, complete_molecules=True, bond_cutoff=1.3,
    cell_margin=0.0, all_images=True, framework_shells=1, unit_cell=False, **kwargs)
```
Builds a coordination polyhedron (convex hull, `scipy.spatial.ConvexHull`) around every atom with ≥`min_neighbors` neighbors within `poly_cutoff`×covalent radius. Produces the atoms/bonds structure object **plus** a separate `<name>_faces` mesh carrying `atom_color`/`element` attributes and the `'polyhedra material'` — a Principled BSDF (roughness 0.6, IOR 1.45, alpha 0.3) and a Glass BSDF (multiscatter GGX, roughness 0, IOR 1.5) mixed half and half, both tinted by `atom_color`: the principled half carries the color and the transparency, the glass half the refraction and the bright edges. It is only built when the material does not already have that glass half, so an older one is rebuilt and tweaks to a current one survive a re-import. Outline (when on) goes on the atoms/bonds only, never the faces. `unit_cell=True` draws the cell (skipped when the file has no 3d cell).

**How the structure is extended past the cell** (`build_polyhedra_atoms`). With
`complete_molecules=True` (default), `select_complete_molecules` grows every molecule **shell by
shell** over `(atom index, cell image)` nodes, using only the neighbor list of the N-atom cell —
each step carries the accumulated image offset along every bond, so growth crosses the cell
boundary in unwrapped space. A molecule is finished when one more shell would not add anything;
if a bond path returns to an atom already in the cluster at a *different* image, the component is
an extended framework (chain, layer, 3d network) and has no finite molecule at all. No supercell
is built, and molecules **larger than one cell** work.

Each atom already lies in some lattice cell (`home = floor(fractional)`), so a molecule moved by
the lattice translation `t` puts an atom inside the central cell exactly when `offset + t ==
-home`. The translations that place at least one atom in the cell are therefore just the negated
`offset + home` values — the copy set falls out of the growth itself.

- `all_images=True` (default, "all molecule images") imports the molecule at every one of them:
  a molecule on a cell face arrives twice, on an edge four times, on a corner eight. The cell
  stays fully populated instead of keeping the hole a cut molecule came from.
- `all_images=False` imports one copy per molecule — the image holding most of the cell's own
  atoms, ties going to the shortest translation so the result is reproducible.
- `cell_margin` (Å) grows the region a molecule must reach, pulling in surrounding molecules.
- `framework_shells` (default 1) grows an extended framework's cell atoms by that many bonded
  shells; 1 closes the coordination polyhedra at the boundary, 0 cuts at the cell.

**The bond criterion, and ASE's skin.** `bond_neighbors()` builds every neighbor list in the
module and takes `skin` explicitly, because ASE's `NeighborList` defaults to `skin=0.3` and adds
it to *each* atom's radius: the criterion silently becomes `d < mult·r₁ + mult·r₂ + 0.6 Å`. An
additive term over-inflates small radii — H goes from `1.3·0.31 = 0.40` to `0.70 Å` — so
N–H···Br hydrogen bonds count as covalent bonds. In a relaxed hybrid antimony bromide that fuses
the organic spacers and the Sb₂Br₁₀ anions into one endless network, and nothing can be completed.
Molecule growth therefore uses **`bond_cutoff` with skin 0** (default 1.3, i.e. `d < 1.3·(r₁+r₂)`);
the expansion, trim and hull searches pass ASE's 0.3 back in, since `expand_cutoff=1.2`,
`trim_cutoff=1.0` and `poly_cutoff=1.1` were tuned with it (the hull needs ≈3.45 Å to catch the
bridging Br of an SbBr₆ octahedron). 1.3 also matches what the geometry nodes actually draw —
3.37 Å for Sb–Br versus the node tree's 3.42 Å at `bond_distance=0.66`.

| structure (cell content) | complete molecules | one copy (`all_images=False`) | expansion only |
|---|---|---|---|
| relaxed hybrid Sb₂Br₁₀ + spacers (120 atoms) | 360 atoms, 96 faces, **every C/N 4-coordinate, every Sb with 6 Br** | 120 atoms, 32 faces, also clean | 136 atoms, 18 of 34 C miscoordinated |
| `crystal.cif` (120 atoms) | 360 atoms, 96 faces | 120 atoms | 136 atoms, cut cations |
| `crystal.cif`, `cell_margin=3.0` | 840 atoms, 352 faces | — | — |
| benzene on the cell corner (12 atoms) | 48 atoms, 4 whole rings | 12 atoms | 19 atoms, rings cut |
| rocksalt 2×2×2 (framework, 16 atoms) | 71 atoms / 244 faces at `framework_shells=1` (16/16 at 0, 229/1124 at 2) | same (framework path) | 44 atoms, 114 faces |

`complete_molecules=False` restores the older behavior for everything: images of higher-indexed
neighbors only, every atom trimmable. Input without a 3d cell always takes that path — there is
nothing to complete.

### Density as mesh — `density_mesh.import_density_mesh`
```python
import_density_mesh(filepath, filename, color_filepath=None, iso_value=0.03,
    shade_smooth=True, preset='DEFAULT', import_atoms=True, color_min=None,
    color_max=None, sample_interior=False, outline=True, **kwargs)
```
Runs marching cubes on the ±`iso_value` levels (`ValueError` if the iso is outside the data range). If `color_filepath` is given, samples that second density onto the surface into a `density_color` vertex attribute; `sample_interior=True` takes the strongest value along the surface normal through the whole volume rather than at the surface point. `color_min==color_max` (default 0) auto-normalizes.

**`preset`** picks the color-ramp material:
- `'DEFAULT'` → `'density_mesh material'`, red (0.0) → white (0.5) → blue (1.0)
- `'ELSTAT'` → `'elstat_potential material'`, blue → white → red
- `'LED'` → `'LED material'`, red (0.8) → green (0.9) → blue (1.0)

To make the isosurface semi-transparent, set the material's Principled BSDF `Alpha` after import:
```python
mat = bpy.data.materials['LED material']
mat.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = 0.4
```
`import_atoms=True` also imports the structure as `nodes`, outlined by default (`outline=True`).

### AMS TAPE41 volumes — `import_cubefiles.tape41_import`

A TAPE41 (the binary KF file AMS/BAND writes for grid-based quantities -
e.g. `NOCVdRhoPlot` off a PEDANOCV restart) is not an ase.io format at all,
so it takes a separate path through `import_ase_molecule`, detected by
`import_cubefiles.is_ams_tape41` (matches a bare `TAPE41`, `*.t41`, or any
`*TAPE41`-suffixed name - the common case of copying it out of a results
directory with the original name kept, e.g. `restart.TAPE41`). Both atoms
and every named volume come from one `import_cubefiles.read_tape41`
pass (`scm.plams.tools.kftools.KFFile`, reading the `Geometry`/`Grid`/`FOO`
sections directly - no external `densf`/`amsvol2cube` conversion step, and
no intermediate `.cube` file):

```python
import_ase_molecule('/path/to/restart.TAPE41', 'restart.TAPE41',
                    representation='nodes', read_density=True)
```

A TAPE41 can carry several named volumes at once (e.g. two NOCV pairs,
`dRhoNOCV=1,k=1` and `dRhoNOCV=2,k=1`, from one `NOCVdRhoPlot: 1 Band 1 2`
restart request) - every one becomes its own volume object, named after
its KF variable, each independently toggleable and using the same
`'+ material'`/`'- material'` convention as `cube2vol`. Call
`read_tape41(filepath, volumes=[...])` directly to select specific names
instead of importing every volume in the file, or `is_ams_tape41(filename)`
to test a path before deciding how to import it.

### Electron density volumes — supercell and materials

A density volume is a VDB grid with the `visualize_edensity` node group on it (Volume to Mesh at
`isovalue`, the `- material` lobe at the mirrored threshold, optional cutoff planes behind the
`cut` switch).

**Cutting.** All six cutoffs are a **depth in angstrom measured in from their own face** of the
density's bounding box: `cutoff X/Y/Z` add to the box minimum, `cutoff -X/-Y/-Z` subtract from the
maximum, and everything past the plane is deleted. 0 therefore cuts nothing *wherever* the density
sits — which an absolute coordinate could not do, since a cube file centred on the origin runs into
negative x, y and z. They share one range (0–100 Å, `DISTANCE` subtype so they read in the scene's
unit) instead of the 100/1000/10000 mixture they had. Measured on `mo.cube`: `cut` on with every
cutoff at 0 renders pixel-identical to `cut` off, `cutoff X = 2` moves the low-x edge in by 2.0 Å
and `cutoff -X = 2` moves the high-x edge in by 2.1 Å.

**Supercell.** Repeating the *isosurface mesh* of one cell cannot extend a density: every copy is
still capped at the cell face it was generated in, so a lobe crossing the boundary keeps its flat
cut. The grid of a periodic calculation is itself periodic — one cell's worth of samples, far
plane not repeated — so tiling the grid *is* the supercell density, and the marching cubes then
runs across the interior boundaries. `import_cubefiles.density_supercell(volume_obj, (nx, ny, nz))`
does that, writing `<name>_<nx>x<ny>x<nz>.vdb` next to the original and repointing the volume at
it. The ASE sidebar's **Density supercell** button (`ase.density_supercell`) applies it to every
density of the structure's collection and defaults its repeats to the structure's own supercell
modifier, so one click matches the two. It always tiles from the single-cell grid recorded at
import (`ase_base_vdb`, with `ase_grid_shape` for the true grid extent), so repeats never
compound; `(1, 1, 1)` restores the original file. Measured on the CHGCAR fixture: a 2.88 Å
footprint becomes 5.72 × 5.72 Å at 2×2×1 and 8.56 × 2.88 Å at 3×1×1, with the quarter-lobes at the
old cell corners joining into whole ones.

**Offset.** `offset a/b/c` (int, one per lattice vector) matches the supercell group's
`Offset_x/y/z`: it translates the finished isosurface by whole cells along `cell a/b/c` (vector
sockets the import fills in from the grid spacing x sample count, since one group is shared by
every density in the file). Being a plain transform it stays **live in the modifier** — nothing is
rewritten, and setting it back to 0 restores the original position. Clicking **Density supercell**
copies both the repeats and the offsets from the structure's own supercell modifier.

**Why the repeat is a rebuild and the offset is not.** A volume cannot be tiled inside geometry
nodes on either supported version: the nodes that could resample a grid live (Sample Grid feeding
a Volume Cube) do not evaluate in 4.4 or in 5.2.1, and tiling the *isosurface mesh* instead gives
what a plain supercell node network gives — every copy still capped at the cell face it was
generated in. Deleting those cap faces after realizing the instances was measured to give the
right extent (5.71 A vs the grid tiling's 5.72 A at 2x2x1) but leaves visible seams where the
half-lobes of neighbouring copies meet, instead of one continuous surface. Rewriting the grid is
therefore the only way to a genuinely cut-free density supercell; it is kept reversible (the base
grid is never touched, `(1,1,1)` restores it) rather than live.

**Upgrading an older density.** A geometry-nodes modifier keeps the group it was created with, so
a .blend saved by an earlier version never gains what a new group revision adds — its densities
have no `offset a/b/c` at all, and still take their materials from modifier sockets. The sidebar
detects that (`_density_nodes_outdated`) and offers **Update density nodes**
(`ase.upgrade_density_nodes`, also run automatically by the supercell button): it repoints the
modifier at the current group, carries every setting over *by socket name*, recovers `cell a/b/c`
from the grid's own transform (`density_cell_vectors`, one `indexToWorld` step per axis × the
sample count, so a triclinic cell works too) and leaves the material slots as they are. Verified
against a scene saved by the previous version: isovalue preserved, cell recovered as 2.831 Å, and
the offset then moves the isosurface by exactly one cell.

**Materials.** The two isosurfaces are tagged with a `mat_slot` face attribute (0 = positive,
1 = negative) that a Set Material Index node at the end of the group turns into the real material
index, exactly as `atoms_and_bonds` does — so the Material Properties tab of the volume object is
in control (slot 0 `+ material`, slot 1 `- material`, or the spin pair for a CHGCAR's second
volume) instead of a material picked in the modifier.

Making that work in **both render engines** needs three pieces, and dropping any one of them
breaks it in a way only a render shows (measured on `mo.cube`, pixels of the two lobes):

| setup | EEVEE | Cycles |
|---|---|---|
| data-linked slots, index only | 12648 grey | 12648 grey |
| object-linked slots, index only | 2792 blue + 2783 red | **5601 blue, 0 red** |
| object-linked + Set Material per sign + Set Material Index | 2792 blue + 2783 red | 2808 blue + 2786 red |

- **A Set Material per sign** inside the group (the `+`/`- material` defaults) gives the geometry a
  two-entry material list. A material index only means something within the list the *geometry*
  carries, and a mesh built by Volume to Mesh starts with an empty one: EEVEE quietly falls back to
  the object's slots, Cycles clamps every face to the first entry — which is why Cycles painted
  both lobes with slot 0.
- **Set Material Index** from `mat_slot` at the end is what makes that index resolve against the
  *object's* slots. Without it the geometry keeps the materials the group assigned and the
  Properties tab has no effect (measured: swapping slot 0 changed nothing).
- **`slot.link = 'OBJECT'`** (`density_materials`), because the isosurface carries no material list
  to link against; a mesh structure gets away with data links since the geometry flowing through
  its tree is the object's own mesh.

With all three, swapping the material in slot 0 repaints the positive lobe in EEVEE and in Cycles —
that is the point of putting the Properties tab in charge, and `density_material_render` in the
smoke test renders it with Cycles to keep it that way. The group carries a revision stamp in its
description; a `visualize_edensity` from an older version is renamed `visualize_edensity_old` and
rebuilt on the next import rather than silently reused.

### Charges — `charges.import_charges`
```python
import_charges(filepath, filename, charge_filepath, resolution=16,
    colorbonds=True, bond_distance=0.66, bond_radius=0.1, outline=False, **kwargs)
```
Reads one charge per atom from a CSV/txt/dat (`read_charges_csv` takes the last numeric field per row; `ValueError` if the count ≠ atom count), stores a per-atom `charge` float attribute, and builds the structure with a `charge_colors` switch that toggles between element colors and a symmetric ±limit charge ramp (materials `'charge_atoms'` and `'color_curve_charge'`). Pass `outline=True` for the house style.

---

## 5. The `nodes` modifier stack (scripting internals)

A `representation='nodes'` import stacks these modifiers on the single-vertex mesh, named `GeometryNodes`, `GeometryNodes.001`, … in order:

1. `GeometryNodes` → `hide atoms`
2. `GeometryNodes.001` → `supercell` — **only if periodic and `add_supercell`**; when skipped, all later indices shift down by one.
3. next → `atoms_and_bonds_<formula>` — the main modifier
4. next → `outline` (if `outline=True`)

**Don't hardcode modifier indices** — a non-periodic import has no supercell modifier. Use `controls.find_ase_modifier(obj)` which returns `(modifier, {socket_name: identifier})`.

**Setting/reading modifier inputs:** always go through the compat helpers — Blender 5.2 moved geometry-nodes modifier inputs off the modifier (`mod["Socket_N"]` no longer works there):
```python
from blender_importASE.node_networks.compat import set_mod_input, get_mod_input
set_mod_input(mod, 'Socket_2', 0.8)      # works on 4.4 / 5.1 / 5.2, tags the depsgraph
value = get_mod_input(mod, 'Socket_2')
```

**atoms_and_bonds sockets:** `Socket_2` = bond distance, `Socket_3` = bond radius, `Socket_4` = resolution, `atom_scale` = overall atom-size multiplier (default 1.0). `pair_table`/`element_table` are Object sockets — reach them by identifier from `find_ase_modifier`, never by hardcoded name.

**supercell sockets:** `Socket_2/3/4` = repeat x/y/z (int, default 1).

**Control tables** (per structure, edited live or by script):
- `pair_table` mesh: one point per element pair, BOOLEAN `cut` attribute — True hides that pair's bonds. Index with `controls.pair_id(z1, z2)` (`= min*119 + max`).
- `element_table` mesh: one point per atomic number, INT `radius_mode` — 0 = covalent, 1 = vdW.

Script example (hide Cu–Cu bonds, show Cu as vdW):
```python
from blender_importASE import controls
from blender_importASE.controls import pair_id
from blender_importASE.node_networks.compat import get_mod_input
mod, idents = controls.find_ase_modifier(obj)
pair_table = get_mod_input(mod, idents['pair_table'])
element_table = get_mod_input(mod, idents['element_table'])
pair_table.data.attributes['cut'].data[pair_id(29, 29)].value = True
element_table.data.attributes['radius_mode'].data[29].value = 1
pair_table.data.update(); element_table.data.update()
obj.update_tag()
```

---

## 6. Live controls (the ASE N-panel)

`ASE_PT_controls` (VIEW_3D ▸ N-panel ▸ **ASE** tab) appears when the active object has a recognized geometry-node modifier. It draws:
- a per-element radius-mode grid (`ase.set_radius_mode`, prop `number`) — covalent ↔ vdW,
- a per-pair bond-cut grid (`ase.toggle_pair_cut`, prop `pair_id`),
- an **Element colors** box: one swatch per element of the structure (see below),
- a **3D printing** box with **Rebuild 3D-print supports** (`ase.rebuild_supports`) when the collection holds real element meshes,
- one box per sibling density/geometry modifier.

**Custom bonds** (`ase.add_dotted_bond`, `dotted_bond.add_bond`) draw a bond between two atoms that the distance-based search does not - a partial bond in a transition state, a hydrogen bond, and so on. Select exactly two atoms (vertices) and click *Add custom bond*; the `bond type` dropdown in the redo panel picks the style:

- `DOTTED` - a row of spheres (`dotted_bond` group)
- `SCALED` - a solid bond thinning with length, capped at `radius`. The reference is the bond's *natural* length (the two atoms' covalent radii added), so a bond at its normal length is full thickness and a stretched/partial one thins in proportion (`scaled_bond` group)
- `DASHED` - alternating cylinder segments (`dashed_bond` group)

```python
from blender_importASE.custom_bonds import add_bond, reset_custom_bonds
add_bond(structure_obj, 3, 17, style='DASHED', segments=8, replace=True)
# radius=None (the default) matches the structure's own bond radius
reset_custom_bonds(structure_obj)          # restore solid bonds, delete the custom ones
```

Each style is its own node group sharing one front-end: both atom positions are sampled from the structure, so the bond follows it (including trajectory animation). All write the `COLOR_CURVE` attribute the bond material reads, blended between the two atoms' colors, and set the material index of the structure's bond slot; the outline modifier is added as usual. Dashes are made by resampling the line, dropping every other edge and running the rest through Curve to Mesh with a circle profile - the segments come out aligned with the bond with no rotation maths.

`replace` stores a per-atom `dotted_partner` int on the structure mesh (the other atom's index + 1; 0/missing = none) which atoms_and_bonds ORs into its bond delete selection. That is **one partner per atom**: replacing 0-1 then 0-2 leaves the 0-1 bond visible again. `ase.reset_custom_bonds` / `reset_custom_bonds()` clears them all and removes the bond objects (`remove_objects=False` keeps the objects). Because the delete-selection wiring lives in atoms_and_bonds, structures imported before this feature need a re-import for `replace` (adding the bond itself works either way).

**Element colors** (`element_colors.py`). An element's color lives in two places: its materials (`'Sb'`, `'Sb-bond'`, and the two-sided `'Sb-C-bond'` gradients) shade the atoms and the ball'n'stick bonds, while the structure mesh's `atom_color` point attribute is what the geometry-node bonds sample at both ends and blend along the curve (`colorbonds`) - and what tints polyhedra faces. The panel's swatch edits the element's atom material; a `depsgraph_update_post` handler notices the change - from the swatch, the Material Properties tab, a driver or a script - and rewrites the `atom_color` entries of that element in every mesh of the file, so the bonds follow. Any mesh pairing `atom_color` with an `element` attribute on the same point domain is covered.

```python
from blender_importASE.element_colors import (set_element_color, get_element_color,
                                              sync_atom_color_attributes)
set_element_color('Sb', (0.571125, 0.109462, 0.274677))   # materials + attributes
get_element_color('Sb')                                    # what it is drawn in now
sync_atom_color_attributes(['Sb', 'C'])                    # attributes from the materials
```

Colors are **linear**, the way Blender reads a base color or a FLOAT_COLOR attribute; the hex codes noted in `utils.atomcolors.color_dict` are their sRGB equivalents, i.e. what the color picker shows. `utils.default_element_color(symbol)` gives the add-on's default (jmol colors for elements the scheme doesn't cover), `ase.reset_element_colors` puts the structure's elements back to it, and `ase.sync_element_colors` runs the sync pass on demand. An element material that already exists in the file keeps its color when another structure is imported, so a re-import never resets a picked color.

`ase.rebuild_supports` is live-adjustable in the F9 redo panel: `base_radius` (0.25), `tip_radius`/"contact radius" (0.1), `support_layer`/"support drop" (0.3), `plate_thickness` (0.6), `plate_holes` (True), `plate_gap` (2.0). It removes existing auto-supports and rebuilds via `exports.build_supports`.

---

## 7. 3D printing

Import with `representation='3D_print'` (real spheres + fused node bonds with icosphere joints). Then either use the N-panel **Rebuild supports** button, or export directly:

```python
exports.export_3dprint(context, filepath, generate_supports=True,
    base_radius=0.25, tip_radius=0.1, support_layer=0.8, plate_thickness=0.6,
    plate_holes=True, plate_gap=2.0)
```
Writes one STL per element + `bonds.stl` + `supports.stl` and zips them.

**Support behavior to know:** supports are generated **only if none already exist**. If you (or the Rebuild button) already made supports, `export_3dprint` exports them **as-is** and ignores its own support parameters. So the workflow is: rebuild supports in the N-panel until they look right, *then* export — the export reuses exactly what you see. Support generation uses a BVH of the real drawn bond geometry as the single source of truth for connectivity and obstacle avoidance (never skewers an atom; may graze bonds), with bottom-up island grounding and an optional holed base plate.

`exports.build_supports(atom_objects, collection, ...)` is the underlying builder (defaults differ slightly: `support_layer=0.8`); pass `bond_objects=` the bond meshes so it can avoid them.

---

## 8. Dependencies

**Nothing installs automatically.** `register()` only *checks* what's importable
(`check_dependency()` — no pip call in it at all) and registers whatever operators that
allows; installing is a deliberate, per-package action from the add-on's preferences
panel (Edit ▸ Preferences ▸ Add-ons ▸ ASE Importer), which lists every `DEPENDENCIES`
entry with a ✓/✗ status and an **Install `<name>`** button when missing
(`ASEInstallDependency`, `bl_idname = "ase.install_dependency"`). Scripts can still call
it directly: `bpy.ops.ase.install_dependency(import_name='ase', pip_name='ase')`.

- **ase** — required; without it the import/export operators, and the whole *File ▸
  Import/Export* menu entries, stay unregistered — clicking **Install ase** registers
  them immediately in the same session (`_register_feature_operators()`, no restart or
  add-on re-toggle needed).
- **scipy** — polyhedra only.
- **scikit-image** — density-as-mesh.
- **scm.plams** (pip name `plams`) — AMS TAPE41 volumes (`import_cubefiles.read_tape41`/`tape41_import`). Pure Python (no compiled extensions), so it never triggers the native-package mismatch cleanup below the way scipy/scikit-image can.
- **openvdb/pyopenvdb** — volumetric density (`import_cubefiles.data2vol`); has no pip wheel worth auto-offering either, so importing it is always a manual, `ImportError`-on-failure story regardless of this section.

`check_dependency()` (and the install button) check with a real `importlib.import_module()` (not `find_spec`), so an installed-but-broken package is correctly treated as missing and reinstallable — this matters because `find_spec` only checks that a module is *findable*, not that it actually imports.

Before checking anything, `check_dependency()` also scans every site-packages-like directory that could hold a stale, wrong-Python compiled package — `_native_package_roots()` returns both the user `modules` folder (which our own installs write into) *and* Blender's own bundled interpreter's `site-packages` (via `sysconfig.get_paths()`, since numpy ships as part of Blender itself) — and deletes any mismatched files it finds in each. This guards against two related failure modes: our install path is keyed on Blender's *version* folder, not its bundled Python (a build bump or a "copy previous settings" migration can leave an old-interpreter numpy/scipy/scikit-image build sitting in the new version's `modules` folder); and Blender's *own* bundled numpy can end up half-updated the same way (e.g. a partial Steam update that bumps the embedded Python but doesn't cleanly replace every compiled file). Either one raises `ImportError: ... Importing the numpy C-extensions failed ... incompatible with python 'cpython-31X'` deep inside the package — `find_spec` doesn't catch this because the files are still *there*, just unloadable.

The scan reads each package's `RECORD` file (written by pip) for `cpXY`-tagged extension filenames; on a mismatch it deletes only the files that record lists which no other, correctly-tagged `RECORD` in the same directory also claims (pip's `--target` mode doesn't clean up a prior conflicting install of the same package the way a normal `--upgrade` would, so two dist-infos can end up describing overlapping files in one shared package directory - a later install overwrites same-path files in place, so an old record can still list a path whose current content actually belongs to the good install). If a file can't be deleted (no write permission - typical for Blender's own install directory), or if something got purged from Blender's own bundled Python rather than our own `modules` folder (where nothing will reinstall a working copy), `ASEAddonPreferences.stale_native_packages` reports it so it's visible in the add-on's preferences panel instead of a bare console print.

Heads-up: if the user has a dev ASE checkout on `sys.path`, it can shadow the pip ASE and cause version-skew surprises (e.g. CHGCAR species-line parsing). `read_vasp_density` is gzip-aware and retries stripping `/` from POTCAR-style species lines.

---

## 9. Rendering conventions for this repo

Gallery panels and GIFs are rendered from `blender_startup.blend` (Cycles + packed HDRI world) with an orthographic camera framed to fit, 1080×1080 for panels. **Always render molecules with `outline=True`.** Reproducible render scripts live in `docs/` (e.g. `docs/render_trajectory_gif.py`). The scratchpad `render_panels.py` regenerates the whole `docs/images/` gallery.

Each custom bond is also listed in the ASE panel of its **structure** (a *Custom bonds* box), so `atom A` / `atom B` and the style's settings can be changed there without selecting the bond object.
