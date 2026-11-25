"""Hückel pi molecular orbitals of conjugated chains and rings, solved and drawn.

The orbitals come from diagonalising the Hückel matrix (alpha on the diagonal, beta between
bonded atoms), so lobe sizes, phases, nodes, degeneracies and the electron filling all follow
from the topology. Each level is drawn as a row of p orbitals whose size is the coefficient and
whose shaded lobe shows the sign.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

from .canvas import Canvas
from .energy import MINUS, _result
from .orbitals import _spin_arrows, hund_fill

# name -> (atoms, cyclic, pi electrons)
SYSTEMS = {
    'ethene': (2, False, 2), 'allyl_cation': (3, False, 2), 'allyl_radical': (3, False, 3),
    'allyl_anion': (3, False, 4), 'butadiene': (4, False, 4), 'pentadienyl_cation': (5, False, 4),
    'pentadienyl_anion': (5, False, 6), 'hexatriene': (6, False, 6),
    'cyclopropenyl_cation': (3, True, 2), 'cyclobutadiene': (4, True, 4),
    'cyclopentadienyl_anion': (5, True, 6), 'cyclopentadienyl_cation': (5, True, 4),
    'benzene': (6, True, 6), 'tropylium': (7, True, 6), 'cyclooctatetraene': (8, True, 8),
}
TOL = 1e-6
# Streitwieser heteroatom parameters: alpha_X = alpha + h beta, beta_CX = k beta, and the pi electrons X gives.
HETEROATOMS = {
    'N_pyridine': {'symbol': 'N', 'h': 0.5, 'k': 1.0, 'electrons': 1},   # pyridine-type N (one pi electron)
    'N_pyrrole': {'symbol': 'N', 'h': 1.5, 'k': 0.8, 'electrons': 2},    # pyrrole-type N (lone pair in the pi system)
    'O_carbonyl': {'symbol': 'O', 'h': 1.0, 'k': 1.0, 'electrons': 1},   # carbonyl O (one pi electron)
    'O_furan': {'symbol': 'O', 'h': 2.0, 'k': 0.8, 'electrons': 2},      # furan-type O (lone pair in the pi system)
}
# name -> (atoms, cyclic, {index: heteroatom type})
HETERO_SYSTEMS = {
    'formaldehyde': (2, False, {1: 'O_carbonyl'}), 'acrolein': (4, False, {3: 'O_carbonyl'}),
    'pyridine': (6, True, {0: 'N_pyridine'}), 'pyrrole': (5, True, {0: 'N_pyrrole'}),
    'furan': (5, True, {0: 'O_furan'}), 'pyrimidine': (6, True, {0: 'N_pyridine', 2: 'N_pyridine'}),
}


def solve(atoms, cyclic=False, bonds=None, h_values=None, k_values=None):
    """Energies ``x`` (E = alpha + x beta, bonding x > 0) and coefficient vectors, highest x first.

    ``h_values``: ``{atom: h}``; ``k_values``: ``{(a, b): k}`` for heteroatoms (Streitwieser)."""
    import numpy as np
    n = int(atoms)
    if not 1 <= n <= 12:
        raise ValueError('A Hückel system here has 1 to 12 atoms')
    h = np.zeros((n, n))
    pairs = bonds if bonds is not None else [(i, i + 1) for i in range(n - 1)] + ([(n - 1, 0)] if cyclic and n > 2 else [])
    for a, b in pairs:
        a, b = int(a), int(b)
        if not (0 <= a < n and 0 <= b < n) or a == b:
            raise ValueError(f'bond {a}-{b} does not join two of the {n} atoms')
        h[a, b] = h[b, a] = (k_values or {}).get((a, b), (k_values or {}).get((b, a), 1.0))
    for atom, value in (h_values or {}).items():
        h[int(atom), int(atom)] = value
    values, vectors = np.linalg.eigh(h)
    order = list(range(n))[::-1]
    out = []
    for k in order:
        vector = vectors[:, k]
        first = next((v for v in vector if abs(v) > TOL), 1.0)
        out.append((float(values[k]), [float(v) * (1 if first > 0 else -1) for v in vector]))
    return out, pairs


def huckel_mo(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Hückel pi-MO diagram.

    ``{"system": "butadiene"}`` (see ``SYSTEMS``), or ``{"atoms": 5, "cyclic": true,
    "electrons": 6}``; ``"bonds": [[0, 1], ...]`` for another topology (atom indices from 0,
    drawn in a row). The report lists each orbital's energy (alpha + x beta), occupation and
    nodes, the HOMO and LUMO, the total pi energy and, for rings, the delocalisation energy
    against isolated ethene units.
    """
    hetero = {}
    if spec.get('system') is not None:
        key = str(spec['system']).lower().replace(' ', '_').replace('-', '_')
        if key in HETERO_SYSTEMS:
            atoms, cyclic, kinds = HETERO_SYSTEMS[key]
            hetero = {i: dict(HETEROATOMS[kind], type=kind) for i, kind in kinds.items()}
            electrons = sum(hetero[i]['electrons'] if i in hetero else 1 for i in range(atoms))
        elif key in SYSTEMS:
            atoms, cyclic, electrons = SYSTEMS[key]
        else:
            raise ValueError(f'Unknown system {spec["system"]!r}; available: {", ".join([*SYSTEMS, *HETERO_SYSTEMS])}')
        name = key
    else:
        atoms, cyclic = int(spec.get('atoms', 0)), bool(spec.get('cyclic', False))
        if spec.get('electrons') is None:
            raise ValueError('Give "electrons" (pi electrons) with "atoms"')
        electrons, name = int(spec['electrons']), str(spec.get('name', 'custom'))
    for item in spec.get('heteroatoms') or []:
        index = int(item.get('index', -1))
        if not 0 <= index < atoms:
            raise ValueError(f'heteroatom index {index} is not one of the {atoms} atoms')
        if item.get('type') is not None:
            if item['type'] not in HETEROATOMS:
                raise ValueError(f'heteroatom type must be one of {", ".join(HETEROATOMS)} (or give "h" and "k")')
            hetero[index] = dict(HETEROATOMS[item['type']], type=item['type'])
        else:
            hetero[index] = {'symbol': str(item.get('symbol', 'X')), 'h': float(item.get('h', 0.0)),
                             'k': float(item.get('k', 1.0)), 'type': 'custom'}
    if not 0 <= electrons <= 2 * atoms:
        raise ValueError(f'{electrons} pi electrons do not fit {atoms} pi orbitals')
    bond_list = spec.get('bonds')
    pairs_default = bond_list if bond_list is not None else [(i, i + 1) for i in range(atoms - 1)] + (
        [(atoms - 1, 0)] if cyclic and atoms > 2 else [])
    k_values = {}
    for a, b in pairs_default:
        for end in (int(a), int(b)):
            if end in hetero:
                k_values[(int(a), int(b))] = hetero[end]['k']
    levels, pairs = solve(atoms, cyclic, bond_list, {i: v['h'] for i, v in hetero.items()}, k_values)
    groups = []
    for energy, vector in levels:
        if groups and abs(groups[-1]['x'] - energy) < 1e-6:
            groups[-1]['vectors'].append(vector)
        else:
            groups.append({'x': energy, 'vectors': [vector]})
    remaining = electrons
    for group in groups:                                 # fill from the most bonding level
        take = min(remaining, 2 * len(group['vectors']))
        group['occupation'] = hund_fill(take, len(group['vectors']))
        remaining -= take

    # ---- drawing: chains as a row of p orbitals (side view), rings from above
    c = Canvas()
    spacing = float(spec.get('level_spacing', 40.0 if cyclic else 34.0))
    lobe = float(spec.get('lobe', 11.0))
    ring_r = 9.0 + 2.2 * atoms
    block_w = (atoms - 1) * 2.2 * lobe + 2 * lobe if not cyclic else 2 * ring_r + 2 * lobe
    xs = sorted({round(g['x'], 6) for g in groups})
    y_of = {x: (x - xs[0]) / ((xs[-1] - xs[0]) or 1) * spacing * (len(xs) - 1) for x in xs}
    # Keep the order but never let two drawn levels overlap (energies are then not to scale).
    room = 2 * ring_r * 0.95 + 14.0 if cyclic else 2 * 1.6 * lobe * 0.62 + 20.0
    for lower, upper in zip(xs, xs[1:]):
        if y_of[upper] - y_of[lower] < room:
            shift = room - (y_of[upper] - y_of[lower])
            for x in xs:
                if x >= upper:
                    y_of[x] += shift
    ring_labels = _ring_labels(groups[0], atoms, lobe, ring_r, block_w, hetero) if cyclic and hetero else {}
    if cyclic and hetero and len(xs) >= 2:
        y1 = y_of[xs[-1]]
        y2 = y_of[xs[-2]]
        label_tops = [c.measure(xy, hetero[index]['symbol'], size=8, align='Center', valign='middle')[1]
                      for index, xy in ring_labels.items()]
        floor = y2 + ring_r + lobe + 4.1 - min(label_tops)
        if floor > y1:
            y_of[xs[-1]] = floor
    report_orbitals = []
    psi = 1
    widest = 0.0
    for group in sorted(groups, key=lambda g: -g['x']):
        y = y_of[round(group['x'], 6)]
        count = len(group['vectors'])
        for k, (vector, occupation) in enumerate(zip(group['vectors'], group['occupation'])):
            x0 = 40.0 + k * (block_w + 30.0)
            if cyclic:
                _ring_orbital(c, x0, y, atoms, lobe, ring_r, vector, occupation)
                if psi == 1:
                    for index, (lx, ly) in ring_labels.items():
                        box = c.label((lx, y + ly), hetero[index]['symbol'],
                                      align='Center', valign='middle', size=8)
                        widest = max(widest, box[2])
            else:
                c.segment((x0 - 14.0, y), (x0 - 2.0, y))
                if occupation:
                    _spin_arrows(c, x0 - 8.0, y, occupation)
                left = x0 + 6.0
                for i, coefficient in enumerate(vector):
                    cx = left + lobe + i * 2.2 * lobe
                    if psi == 1 and i in hetero:
                        c.label((cx, y + 1.6 * lobe + 10.0), hetero[i]['symbol'], align='Center', size=8)
                    size = abs(coefficient) * lobe * 1.6
                    if size < 0.8:
                        c.circle((cx, y), 0.8, fill='#000000', outline=None)   # a node at this atom
                        continue
                    _p_orbital(c, cx, y, size, coefficient > 0)
            widest = max(widest, x0 + block_w + 6.0)
            report_orbitals.append({'psi': psi, 'x': round(group['x'], 4) + 0.0,
                                    'energy': 'alpha ' + _energy_text(group['x'], 'b'),
                                    'electrons': len(occupation),
                                    'nodes_between_atoms': sum(1 for a, b in pairs if vector[a] * vector[b] < -TOL),
                                    'coefficients': [round(v, 3) for v in vector],
                                    'y': round(y, 2)})
            psi += 1
        label = '\u03c8{sub:%d}' % (psi - count) + ('' if count == 1 else ', \u03c8{sub:%d}' % (psi - 1))
        c.label((22.0, y), label, align='Right', valign='middle')
        group['y'] = y
    for group in groups:
        c.label((widest + 8.0, group['y']), '\u03b1 ' + _energy_text(group['x'], '\u03b2'), valign='middle')
    top, bottom = min(y_of.values()) - (ring_r + 12.0 if cyclic else 24.0), max(y_of.values()) + 16.0
    axis_x = c.extent[0] - 10.0
    c.arrow((axis_x, bottom), (axis_x, top))
    c.label((axis_x - 4.0, top - 6.0), 'E')
    title = spec.get('title', name.replace('_', ' ') if name != 'custom' else None)
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')

    filled = [o for o in report_orbitals if o['electrons']]
    homo = filled[-1]['psi'] if filled else None
    lumo = next((o['psi'] for o in report_orbitals if o['electrons'] == 0), None)
    somo = [o['psi'] for o in report_orbitals if o['electrons'] == 1]
    total = sum(o['x'] * o['electrons'] for o in report_orbitals)
    report: dict[str, Any] = {'diagram': 'huckel_mo', 'system': name, 'atoms': atoms, 'cyclic': cyclic,
                              'pi_electrons': electrons, 'orbitals': report_orbitals, 'homo': homo, 'lumo': lumo, 'somo': somo,
                              'pi_energy_beta': round(total, 4),
                              'unpaired_electrons': sum(1 for o in report_orbitals if o['electrons'] == 1),
                              'energies': 'E = alpha + x beta; beta < 0, so x > 0 is bonding'}
    if hetero:
        report['heteroatoms'] = [{'index': i, 'symbol': v['symbol'], 'type': v['type'], 'h': v['h'], 'k': v['k']}
                                 for i, v in sorted(hetero.items())]
    if cyclic:
        # Against electrons/2 isolated ethene double bonds (2 beta each).
        if not hetero:
            report['delocalisation_energy_beta'] = round(total - 2.0 * (electrons // 2), 4)
        open_shell = report['unpaired_electrons'] > 0
        report['huckel_rule'] = ('aromatic (4n+2)' if electrons % 4 == 2 and not open_shell
                                 else 'antiaromatic (4n)' if electrons % 4 == 0 else 'neither')
    report['warnings'] = []
    return _result(c, output, report)


def _ring_orbital(c, x0, y, atoms, lobe, ring_r, vector, occupation):
    c.segment((x0 - 14.0, y), (x0 - 2.0, y))
    if occupation:
        _spin_arrows(c, x0 - 8.0, y, occupation)
    cx0 = x0 + 6.0 + lobe + ring_r
    corners = [(cx0 + ring_r * math.cos(-math.pi / 2 + 2 * math.pi * i / atoms),
                y + ring_r * math.sin(-math.pi / 2 + 2 * math.pi * i / atoms)) for i in range(atoms)]
    c.path(corners, closed=True)
    for (px, py), coefficient in zip(corners, vector):
        radius = abs(coefficient) * lobe
        if radius >= 0.8:
            c.circle((px, py), radius, fill='#000000' if coefficient > 0 else '#FFFFFF')
    return corners


def _radial_label_position(origin, direction, label_box, obstacles, gap=4.2):
    intervals = []
    for box in obstacles:
        lower = (box[0] - gap - label_box[2], box[1] - gap - label_box[3])
        upper = (box[2] + gap - label_box[0], box[3] + gap - label_box[1])
        start, end = 0.0, math.inf
        for i in (0, 1):
            if abs(direction[i]) < 1e-12:
                if not lower[i] <= origin[i] <= upper[i]:
                    break
            else:
                crossings = sorted(((lower[i] - origin[i]) / direction[i],
                                    (upper[i] - origin[i]) / direction[i]))
                start, end = max(start, crossings[0]), min(end, crossings[1])
        else:
            if start <= end:
                intervals.append((start, end))
    distance = 0.0
    for start, end in sorted(intervals):
        if start > distance:
            break
        if end >= distance:
            distance = end + 0.05
    return tuple(origin[i] + distance * direction[i] for i in (0, 1))


def _ring_labels(group, atoms, lobe, ring_r, block_w, hetero):
    layout = Canvas()
    first_corners = None
    for k, (vector, occupation) in enumerate(zip(group['vectors'], group['occupation'])):
        corners = _ring_orbital(layout, 40.0 + k * (block_w + 30.0), 0.0,
                                atoms, lobe, ring_r, vector, occupation)
        if k == 0:
            first_corners = corners
    obstacles = []
    for node in layout.page:
        if node.get('BoundingBox'):
            v = [float(q) for q in node.get('BoundingBox').split()]
            pad = 2.0 if node.tag == 'arrow' else 0.0
            obstacles.append((min(v[0], v[2]) - pad, min(v[1], v[3]) - pad,
                              max(v[0], v[2]) + pad, max(v[1], v[3]) + pad))
        elif node.get('CurvePoints'):
            v = [float(q) for q in node.get('CurvePoints').split()]
            obstacles.append((min(v[0::2]), min(v[1::2]), max(v[0::2]), max(v[1::2])))
    count = len(group['vectors'])
    psi_label = '\u03c8{sub:1}' + (', \u03c8{sub:%d}' % count if count > 1 else '')
    obstacles.append(layout.measure((22.0, 0.0), psi_label, align='Right', valign='middle'))
    positions = {}
    for index, info in sorted(hetero.items()):
        angle = -math.pi / 2 + 2 * math.pi * index / atoms
        box = layout.measure((0.0, 0.0), info['symbol'], size=8, align='Center', valign='middle')
        xy = _radial_label_position(first_corners[index], (math.cos(angle), math.sin(angle)), box, obstacles)
        positions[index] = xy
        obstacles.append(layout.measure(xy, info['symbol'], size=8, align='Center', valign='middle'))
    return positions


def _energy_text(x, beta):
    """``+ 1.618b``-style text for E = alpha + x beta; empty for x = 0 (a non-bonding level)."""
    if abs(x) < 1e-6:
        return ''
    value = f'{abs(x):.3f}'.rstrip('0').rstrip('.')
    return f'{"+" if x > 0 else MINUS} {"" if value == "1" else value}{beta}'


def _p_orbital(c, cx, cy, size, positive):
    """A vertical p orbital: two lobes, the upper one filled when the coefficient is positive."""
    width = size * 0.62
    half = [(width * math.sin(math.pi * k / 24), size * (1 - math.cos(math.pi * k / 24)) / 2) for k in range(25)]
    for direction, shaded in ((-1.0, positive), (1.0, not positive)):
        outline = [(cx + dx, cy + direction * dy) for dx, dy in half]
        outline += [(cx - dx, cy + direction * dy) for dx, dy in reversed(half[1:-1])]
        c.path(outline, closed=True, filled=shaded)
