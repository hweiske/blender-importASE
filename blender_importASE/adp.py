"""Anisotropic displacement parameters (ADPs) from CIF and SHELX files.

The CIF lists one displacement tensor per site of the asymmetric unit,
either anisotropic (_atom_site_aniso_U_ij, B_ij or beta_ij, matched to the
sites by label) or isotropic (_atom_site_U_iso_or_equiv / B_iso_or_equiv).
Every atom ASE generates from a site by a space-group operation gets that
site's tensor rotated by the operation, in Cartesian coordinates:

    U* = N U_cif N              N = diag(|a*|, |b*|, |c*|)
    U_site = A U* A^T           A: lattice vectors as columns
    U_atom = R U_site R^T       R = A W A^-1, W the fractional rotation

The principal axes of U_atom (its eigenvectors) and the RMS displacements
along them (square roots of the eigenvalues) give the thermal ellipsoid;
the probability level only scales it (see probability_scale).
"""
import re

import numpy as np

ADP_RINGS_MATERIAL = 'adp_rings'
EIGHT_PI2 = 8 * np.pi ** 2
TWO_PI2 = 2 * np.pi ** 2
_TENSOR_KEYS = ('11', '22', '33', '12', '13', '23')


def probability_scale(probability):
    """Factor from the RMS displacements to the semi-axes of the ellipsoid
    holding the atom with the given probability. The squared Mahalanobis
    distance of a 3d Gaussian is chi-square distributed with 3 degrees of
    freedom, so 0.5 (the ORTEP default) gives 1.5382. Probabilities outside
    (0, 1) have no ellipsoid and are clamped just inside."""
    from scipy.stats import chi2
    probability = min(max(float(probability), 1e-6), 1 - 1e-6)
    return float(np.sqrt(chi2.ppf(probability, 3)))


def cif_to_cartesian(tensor, cell):
    """U_ij as a CIF or SHELX file lists them (in the basis normalized to the
    reciprocal axis lengths) as a Cartesian tensor."""
    lattice = np.asarray(cell).T
    N = np.diag(np.linalg.norm(np.asarray(cell.reciprocal()), axis=1))
    return lattice @ N @ tensor @ N @ lattice.T


def _floats(values):
    """CIF values as floats: standard uncertainties like '0.0234(5)' are
    dropped, unknown values ('?', '.') become NaN."""
    out = []
    for value in values:
        if isinstance(value, str):
            value = re.sub(r'\(\d+\)$', '', value.strip())
        try:
            out.append(float(value))
        except (TypeError, ValueError):
            out.append(np.nan)
    return np.array(out, dtype=float)


def _n_sites(info):
    for key in ('_atom_site_label', '_atom_site_type_symbol',
                '_atom_site_fract_x', '_atom_site_cartn_x'):
        if key in info:
            return len(info[key])
    raise ValueError('the CIF has no _atom_site loop')


def site_adps(info, cell):
    """Cartesian displacement tensor per row of the _atom_site loop, as an
    (n_sites, 3, 3) array; NaN for sites without any displacement data.

    `info` is the tag dictionary ASE stores with read(..., store_tags=True).
    """
    n_sites = _n_sites(info)
    lattice = np.asarray(cell).T
    U = np.full((n_sites, 3, 3), np.nan)

    if '_atom_site_u_iso_or_equiv' in info:
        uiso = _floats(info['_atom_site_u_iso_or_equiv'])
    elif '_atom_site_b_iso_or_equiv' in info:
        uiso = _floats(info['_atom_site_b_iso_or_equiv']) / EIGHT_PI2
    else:
        uiso = None
    if uiso is not None:
        # isotropic in any basis; NaN * 0 keeps a missing value all-NaN
        U[:] = uiso[:, None, None] * np.eye(3)

    labels = [str(label) for label in info.get('_atom_site_label', [])]
    aniso_labels = info.get('_atom_site_aniso_label')
    if not labels or aniso_labels is None:
        return U
    for kind in ('u', 'b', 'beta'):
        keys = [f'_atom_site_aniso_{kind}_{ij}' for ij in _TENSOR_KEYS]
        if all(key in info for key in keys):
            break
    else:
        return U

    u11, u22, u33, u12, u13, u23 = (_floats(info[key]) for key in keys)
    row_of = {label: row for row, label in enumerate(labels)}
    for n, label in enumerate(aniso_labels):
        row = row_of.get(str(label))
        if row is None:
            continue
        tensor = np.array([[u11[n], u12[n], u13[n]],
                           [u12[n], u22[n], u23[n]],
                           [u13[n], u23[n], u33[n]]])
        if kind == 'beta':
            # beta_ij = 2 pi^2 a*_i a*_j U_ij, i.e. already U* up to 2 pi^2
            U[row] = lattice @ (tensor / TWO_PI2) @ lattice.T
        else:
            U[row] = cif_to_cartesian(tensor / EIGHT_PI2 if kind == 'b' else tensor, cell)
    return U


def _site_fractional(info, cell):
    if '_atom_site_fract_x' in info:
        return np.column_stack([_floats(info[f'_atom_site_fract_{x}'])
                                for x in 'xyz'])
    cartesian = np.column_stack([_floats(info[f'_atom_site_cartn_{x}'])
                                 for x in 'xyz'])
    return np.linalg.solve(np.asarray(cell).T, cartesian.T).T


def atom_adps(atoms, symprec=2e-3):
    """Cartesian displacement tensor of every atom of a CIF read with
    ase.io.read(..., store_tags=True): the tensor of the site it was
    generated from, rotated by the space-group operation that generated it.
    Returns an (n_atoms, 3, 3) array, NaN where the CIF has no data."""
    info = atoms.info
    site_U = site_adps(info, atoms.cell)
    kinds = atoms.arrays.get('spacegroup_kinds', np.arange(len(atoms)))
    U = site_U[kinds]
    spacegroup = info.get('spacegroup')
    if spacegroup is None or atoms.cell.rank != 3:
        return U

    lattice = np.asarray(atoms.cell).T
    to_fractional = np.linalg.inv(lattice)
    sites = _site_fractional(info, atoms.cell)
    rotations, translations = spacegroup.get_op()
    fractional = atoms.get_scaled_positions(wrap=False)
    unmatched = 0
    for j, kind in enumerate(kinds):
        if not np.isfinite(U[j]).all():
            continue
        diff = rotations @ sites[kind] + translations - fractional[j]
        diff -= np.rint(diff)
        match = np.flatnonzero(np.all(np.abs(diff) < symprec, axis=1))
        if not match.size:
            unmatched += 1
            continue
        # on a special position several operations match; the site
        # symmetry leaves the tensor the same under any of them
        R = lattice @ rotations[match[0]] @ to_fractional
        U[j] = R @ U[j] @ R.T
    if unmatched:
        print(f'adp: no symmetry operation found for {unmatched} atoms - '
              'their ellipsoids are left unrotated')
    return U


def read_atoms(filepath, index=':'):
    """ase.io.read, extended by SHELX .res / .ins files, keeping a CIF's
    atom-site tags (store_tags) so its displacement parameters can be read.
    Returns a list of Atoms for index=':', else one Atoms."""
    import ase.io
    from .shelx import SHELX_EXTENSIONS, read_shelx

    lower = filepath.lower()
    if lower.endswith(SHELX_EXTENSIONS):
        atoms = read_shelx(filepath)
        return [atoms] if index == ':' else atoms
    if lower.endswith('.cif'):
        return ase.io.read(filepath, index=index, format='cif', store_tags=True)
    return ase.io.read(filepath, index=index)


def anisotropic_adps(atoms):
    """atom_adps for a structure that carries an anisotropic displacement
    table, None for one without (no CIF / SHELX data, or isotropic values
    only) - the importers draw ellipsoids only when there are some."""
    if '_atom_site_aniso_label' not in atoms.info:
        return None
    try:
        U = atom_adps(atoms)
    except (KeyError, ValueError) as exc:
        print(f'adp: displacement parameters not readable ({exc}) - ignored')
        return None
    return U if np.isfinite(U).all(axis=(1, 2)).any() else None


def ellipsoids(U):
    """Principal axes of displacement tensors.

    Returns (rotations, axes, valid): unit quaternions (w, x, y, z) turning
    the x, y and z axes onto the principal axes, the RMS displacements along
    them in angstrom, and whether the tensor describes an ellipsoid at all -
    False where data is missing or the tensor is not positive definite (a
    refinement gone wrong), which the node tree draws as a plain sphere.
    """
    from scipy.spatial.transform import Rotation

    n = len(U)
    rotations = np.tile([1.0, 0.0, 0.0, 0.0], (n, 1))
    axes = np.zeros((n, 3))
    valid = np.zeros(n, dtype=bool)
    finite = np.flatnonzero(np.isfinite(U).all(axis=(1, 2)))
    if not finite.size:
        return rotations, axes, valid
    values, vectors = np.linalg.eigh(U[finite])
    # eigh may hand back a reflection; flip one axis to make it a rotation
    vectors[np.linalg.det(vectors) < 0, :, 2] *= -1
    rotations[finite] = Rotation.from_matrix(vectors).as_quat()[:, [3, 0, 1, 2]]
    axes[finite] = np.sqrt(np.clip(values, 0.0, None))
    valid[finite] = values.min(axis=1) > 0
    return rotations, axes, valid
