"""Structures and displacement parameters from SHELX .res / .ins files.

ASE has no reader for the SHELXL instruction format, so this one builds the
structure the way ASE's CIF reader does - the sites of the file expanded by
the symmetry operations through ase.spacegroup.crystal - and stores the
atom-site data as the CIF tags adp.py reads (as with read(..., store_tags=
True) on a CIF): the displacement parameters of both formats go through
the same code.

What is read:

- CELL (the wavelength in front is skipped), LATT (sign: centrosymmetric or
  not; |n|: the lattice centering), SYMM (the identity is implied), SFAC
  (short form 'SFAC C H O' or one long-form element per line), FVAR
- atoms: 'name sfac x y z sof U11 U22 U33 U23 U13 U12' or '... sof Uiso',
  continued over lines ending in '='. Parameters are decoded like SHELXL
  does: 10 + p is p held fixed, 10m + p (m > 1) is p times free variable m,
  -(10m + p) is p times (free variable m - 1). A Uiso between -5 and -0.5
  rides on the previous non-hydrogen atom: -Uiso times its U_eq.

Everything after END (the difference-map Q peaks of a .res) is ignored, as
are the other instructions.
"""
import math
import re

import numpy as np

from .adp import cif_to_cartesian

SHELX_EXTENSIONS = ('.res', '.ins')

# every SHELXL instruction: an atom name may not be one of them, which is
# how atom lines are told apart
_INSTRUCTIONS = {
    'ABIN', 'ACTA', 'AFIX', 'ANIS', 'ANSC', 'ANSR', 'BASF', 'BIND', 'BLOC',
    'BOND', 'BUMP', 'CELL', 'CGLS', 'CHIV', 'CONF', 'CONN', 'DAMP', 'DANG',
    'DEFS', 'DELU', 'DFIX', 'DISP', 'EADP', 'END', 'EQIV', 'EXTI', 'EXYZ',
    'FEND', 'FLAT', 'FMAP', 'FRAG', 'FREE', 'FVAR', 'GRID', 'HFIX', 'HKLF',
    'HTAB', 'ISOR', 'LATT', 'LAUE', 'LIST', 'L.S.', 'MERG', 'MORE', 'MOVE',
    'MPLA', 'NCSY', 'NEUT', 'OMIT', 'PART', 'PLAN', 'PRIG', 'REM', 'RESI',
    'RIGU', 'RTAB', 'SADI', 'SAME', 'SFAC', 'SHEL', 'SIMU', 'SIZE', 'SPEC',
    'STIR', 'SUMP', 'SWAT', 'SYMM', 'TEMP', 'TITL', 'TWIN', 'TWST', 'UNIT',
    'WGHT', 'WIGL', 'WPDB', 'XNPD', 'ZERR', 'BEDE', 'LONE', 'SPIN', 'HOPE',
    'MOLE', 'ANIS', 'NOTE', 'ESEL', 'EGEN',
}

# lattice centering translations of LATT |n|
_CENTERING = {
    1: [(0, 0, 0)],
    2: [(0, 0, 0), (0.5, 0.5, 0.5)],
    3: [(0, 0, 0), (2 / 3, 1 / 3, 1 / 3), (1 / 3, 2 / 3, 2 / 3)],
    4: [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)],
    5: [(0, 0, 0), (0, 0.5, 0.5)],
    6: [(0, 0, 0), (0.5, 0, 0.5)],
    7: [(0, 0, 0), (0.5, 0.5, 0)],
}


def _logical_lines(filepath):
    """Instruction lines with comments dropped and '=' continuations
    joined, up to END."""
    with open(filepath, errors='replace') as fh:
        raw = fh.read().splitlines()
    lines, pending = [], ''
    for line in raw:
        line = line.split('!', 1)[0].rstrip()
        if pending:
            line = pending + ' ' + line.strip()
            pending = ''
        if not line.strip():
            continue
        if line.upper().startswith('REM'):
            continue
        if line.endswith('='):
            pending = line[:-1]
            continue
        lines.append(line)
        if line.split()[0].upper() == 'END':
            break
    if pending:
        lines.append(pending)
    return lines


def _symmetry_operation(text):
    """'1/2-X, 1/2+Y, 1/2-Z' as (W, t) in fractional coordinates."""
    rows = text.upper().replace(' ', '').split(',')
    if len(rows) != 3:
        raise ValueError(f'cannot read the symmetry operation {text!r}')
    W, t = np.zeros((3, 3)), np.zeros(3)
    for i, row in enumerate(rows):
        for sign, number, axis in re.findall(r'([+-]?)([0-9./]*)([XYZ]?)', row):
            if not number and not axis:
                continue
            if '/' in number:
                numerator, denominator = number.split('/')
                value = float(numerator) / float(denominator)
            else:
                value = float(number) if number else 1.0
            value = -value if sign == '-' else value
            if axis:
                W[i, 'XYZ'.index(axis)] += value
            else:
                t[i] += value
    return W, t


def _decode(value, free_variables):
    """A SHELXL parameter with its fixed / free-variable coding resolved."""
    magnitude = abs(value)
    if magnitude <= 5:
        return value
    m = int((magnitude + 1e-6) // 10)
    p = magnitude - 10 * m
    if m == 1:
        return math.copysign(p, value)
    fv = free_variables[m - 1] if m - 1 < len(free_variables) else 1.0
    return p * fv if value > 0 else p * (fv - 1)


def _element(symbol):
    symbol = re.sub(r'[^A-Za-z]', '', symbol)
    return symbol[:1].upper() + symbol[1:].lower()


def read_shelx(filepath):
    """Read a SHELX .res / .ins file into an ase.Atoms object of the full
    unit cell, with atom-site tags like a CIF read with store_tags=True
    (labels, fractional coordinates, U_iso and the aniso U table), so
    adp.atom_adps works on it as on a CIF."""
    from ase.spacegroup import crystal
    from ase.spacegroup.spacegroup import spacegroup_from_data

    cellpar = None
    latt = 1
    symm = []
    sfac = []
    free_variables = []
    sites = []   # (label, element, fractional, U list or Uiso)
    for line in _logical_lines(filepath):
        tokens = line.split()
        keyword = tokens[0].upper()
        if keyword == 'CELL':
            cellpar = [float(x) for x in tokens[2:8]]
        elif keyword == 'LATT':
            latt = int(tokens[1])
        elif keyword == 'SYMM':
            symm.append(_symmetry_operation(line.split(None, 1)[1]))
        elif keyword == 'SFAC':
            # short form: element symbols only; long form: one element
            # followed by its scattering factor numbers
            for token in tokens[1:]:
                try:
                    float(token)
                except ValueError:
                    sfac.append(_element(token))
        elif keyword == 'FVAR':
            # fv(1) is the overall scale factor, fv(2) the first free variable
            free_variables.extend(float(x) for x in tokens[1:])
        elif keyword in _INSTRUCTIONS or re.fullmatch(r'Q\d+', keyword):
            continue
        else:
            atom = _atom(tokens, sfac, free_variables)
            if atom is not None:
                sites.append(atom)

    if cellpar is None:
        raise ValueError(f'{filepath}: no CELL instruction')
    if not sites:
        raise ValueError(f'{filepath}: no atoms')

    from ase.cell import Cell
    cell = Cell.fromcellpar(cellpar)

    # riding hydrogens: -Uiso times U_eq of the previous non-hydrogen atom
    labels, uiso, aniso = [], [], {}
    parent_ueq = None
    seen = {}
    for label, element, _, U in sites:
        # labels repeat across residues (RESI); keep them unique
        seen[label] = seen.get(label, 0) + 1
        unique = label if seen[label] == 1 else f'{label}#{seen[label]}'
        labels.append(unique)
        if isinstance(U, list):
            tensor = np.array([[U[0], U[5], U[4]],
                               [U[5], U[1], U[3]],
                               [U[4], U[3], U[2]]])
            ueq = float(np.trace(cif_to_cartesian(tensor, cell)) / 3)
            aniso[unique] = U
        elif -5 < U < -0.5:
            ueq = -U * parent_ueq if parent_ueq is not None else np.nan
        else:
            ueq = U
        uiso.append(ueq)
        if element != 'H' and np.isfinite(ueq):
            parent_ueq = ueq

    rotations, translations = [np.eye(3)], [np.zeros(3)]
    for W, t in symm:
        rotations.append(W)
        translations.append(t)
    if latt > 0:
        rotations += [-W for W in rotations]
        translations += [-t for t in translations]
    centering = _CENTERING.get(abs(latt))
    if centering is None:
        raise ValueError(f'{filepath}: unknown LATT {latt}')
    all_rotations = np.array([W for _ in centering for W in rotations])
    all_translations = np.array([(t + c) % 1.0 for c in centering
                                 for t in translations])
    spacegroup = spacegroup_from_data(
        no=1, setting=1, centrosymmetric=False, subtrans=[[0, 0, 0]],
        sitesym=[], rotations=all_rotations, translations=all_translations)

    info = {
        '_atom_site_label': labels,
        '_atom_site_type_symbol': [element for _, element, _, _ in sites],
        '_atom_site_fract_x': [frac[0] for _, _, frac, _ in sites],
        '_atom_site_fract_y': [frac[1] for _, _, frac, _ in sites],
        '_atom_site_fract_z': [frac[2] for _, _, frac, _ in sites],
        '_atom_site_u_iso_or_equiv': uiso,
    }
    if aniso:
        info['_atom_site_aniso_label'] = list(aniso)
        # SHELX order: U11 U22 U33 U23 U13 U12
        for column, key in enumerate(('11', '22', '33', '23', '13', '12')):
            info[f'_atom_site_aniso_u_{key}'] = [U[column] for U in aniso.values()]

    atoms = crystal([element for _, element, _, _ in sites],
                    basis=[frac for _, _, frac, _ in sites],
                    spacegroup=spacegroup, cellpar=cellpar,
                    onduplicates='keep')
    atoms.info.update(info)
    return atoms


def _atom(tokens, sfac, free_variables):
    """(label, element, fractional coordinates, U) of an atom line, U being
    the six U_ij in SHELX order or a single (possibly riding) Uiso; None for
    a line that is not an atom."""
    if len(tokens) < 5:
        return None
    try:
        index = int(tokens[1])
        numbers = [float(x) for x in tokens[2:]]
    except ValueError:
        return None
    if not 1 <= index <= len(sfac):
        return None
    frac = [_decode(x, free_variables) for x in numbers[:3]]
    rest = numbers[3:]
    # rest: [sof] [Uiso | U11 U22 U33 U23 U13 U12]
    if len(rest) >= 7:
        U = [_decode(x, free_variables) for x in rest[1:7]]
    elif len(rest) >= 2:
        # a riding Uiso keeps its sign: -1.2 is 1.2 x U_eq of the parent
        U = rest[1] if -5 < rest[1] < -0.5 else _decode(rest[1], free_variables)
    else:
        U = 0.05
    return tokens[0], sfac[index - 1], frac, U
