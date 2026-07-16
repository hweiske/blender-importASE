# blender-importASE — guide for Claude instances

This is a Blender add-on for importing atomistic structures (via [ASE](https://wiki.fysik.dtu.dk/ase/)) and turning them into publication-quality renders: molecules, crystals, coordination polyhedra, electron-density isosurfaces (volume or mesh), partial-charge colorings, and 3D-printable models. This document is the reference for driving it — both from the Blender GUI and from Python scripts. Everything the GUI does calls the same functions you can call directly, so scripting and clicking are interchangeable.

- **Package:** `blender_importASE/` (add-on version 2.2, min Blender 4.4; tested on 4.4 and 5.1).
- **Dependencies:** `ase` (auto-installed on first `register()`), plus `scipy` (polyhedra), `scikit-image` (density-as-mesh, auto-installed), `openvdb`/`pyopenvdb` (volumetric density). See [§8](#8-dependencies).

---

## 1. Two ways to drive it

**From the GUI.** *File ▸ Import* gains four entries and *File ▸ Export* gains two (see [§2](#2-operators-file--importexport)). After importing, the **N-panel ▸ ASE tab** (`ASE_PT_controls`) exposes live controls for the active structure: per-element radius mode, per-pair bond hiding, and the 3D-print support rebuilder ([§6](#6-live-controls-the-ase-n-panel)).

**From Python.** Every operator is a thin wrapper over a module-level function. In a headless render or a notebook:

```python
import sys; sys.path.insert(0, '/path/to/blender-importASE')
import bpy, blender_importASE
blender_importASE.register()                       # installs ase if missing
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

The importers accept multi-file selection (`files`/`directory` props). Exporters use the active object.

---

## 3. Importing molecules — `import_ase_molecule`

```python
import_ase_molecule(filepath, filename, overwrite=True, add_supercell=True,
    resolution=32, colorbonds=False, long_bonds=False, color=0.2, scale=1,
    unit_cell=False, representation="Balls'n'Sticks", read_density=True,
    shift_cell=False, imageslice=1, animate=True, outline=True, **kwargs)
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
- `read_density=True` — if the file carries a volume (`.cube`, CHGCAR/PARCHG/AECCAR), builds a volumetric density object (needs openvdb). Isosurface materials `'+ material'`/`'- material'`.
- `add_supercell=True` — adds the supercell modifier when the cell is periodic; repeat counts live on `Socket_2/3/4` of that modifier ([§5](#5-the-nodes-modifier-stack-scripting-internals)).
- `animate=True`, `imageslice=n` — for trajectories, animate every *n*th frame (`overwrite=True` forces `nodes`).
- `colorbonds=True` — color bond halves by their atoms; `unit_cell=True` draws the cell box.

The GUI operator passes different defaults (`scale=0.5`, `color=0.6`, `representation="nodes"`; its `zero_cell` maps to `shift_cell`).

---

## 4. The other importers

### Polyhedra — `polyhedra.import_polyhedra`
```python
import_polyhedra(filepath, filename, expand_cutoff=1.2, trim_cutoff=1.0,
    poly_cutoff=1.1, min_neighbors=4, include_hydrogen=False, resolution=16,
    colorbonds=True, bond_distance=0.66, bond_radius=0.1, outline=False, **kwargs)
```
Builds a coordination polyhedron (convex hull, `scipy.spatial.ConvexHull`) around every atom with ≥`min_neighbors` neighbors within `poly_cutoff`×covalent radius. `expand_cutoff` pulls in periodic images so boundary polyhedra close; `trim_cutoff` removes stragglers. Produces the atoms/bonds structure object **plus** a separate `<name>_faces` mesh carrying an `atom_color` attribute and the semi-transparent `'polyhedra material'`. Outline (when on) goes on the atoms/bonds only, never the faces.

### Density as mesh — `density_mesh.import_density_mesh`
```python
import_density_mesh(filepath, filename, color_filepath=None, iso_value=0.03,
    shade_smooth=True, preset='DEFAULT', import_atoms=True, color_min=None,
    color_max=None, sample_interior=False, **kwargs)
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
`import_atoms=True` also imports the structure as `nodes`. To get outlined atoms alongside the mesh, import them separately with `import_ase_molecule(..., outline=True)` and pass `import_atoms=False`.

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
- a **3D printing** box with **Rebuild 3D-print supports** (`ase.rebuild_supports`) when the collection holds real element meshes,
- one box per sibling density/geometry modifier.

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

- **ase** — required; auto-pip-installed into Blender's user `modules` path on `register()`. If install fails, no operators register and `ASEAddonPreferences.install_failed` is set.
- **scipy** — polyhedra only; lazy import, no auto-install.
- **scikit-image** — density-as-mesh; auto-installed on demand.
- **openvdb/pyopenvdb** — volumetric density (`import_cubefiles.data2vol`); no auto-install, clear `ImportError` if missing.

Heads-up: if the user has a dev ASE checkout on `sys.path`, it can shadow the pip ASE and cause version-skew surprises (e.g. CHGCAR species-line parsing). `read_vasp_density` is gzip-aware and retries stripping `/` from POTCAR-style species lines.

---

## 9. Rendering conventions for this repo

Gallery panels and GIFs are rendered from `blender_startup.blend` (Cycles + packed HDRI world) with an orthographic camera framed to fit, 1080×1080 for panels. **Always render molecules with `outline=True`.** Reproducible render scripts live in `docs/` (e.g. `docs/render_trajectory_gif.py`). The scratchpad `render_panels.py` regenerates the whole `docs/images/` gallery.
