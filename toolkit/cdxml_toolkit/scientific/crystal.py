"""Cubic unit cells in oblique projection, with the cell contents counted.

The spec names a structure type or lists fractional coordinates; the builder expands each
site to every equivalent position on the cell boundary, draws the atoms back to front and
reports the atoms per cell (corner 1/8, edge 1/4, face 1/2, interior 1), the formula and
each species' coordination number.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction
from functools import reduce
from typing import Any, Mapping

from .canvas import BLACK, BLUE, ORANGE, WHITE, Canvas, resolve_color
from .common import number
from .energy import _result

_FCC = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
_TETRA = [(0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75)]
_TETRA_ALL = [(x, y, z) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)]
# preset -> (default element per role, sites per role)
PRESETS = {
    'sc': ({'A': 'Po'}, {'A': [(0, 0, 0)]}),
    'bcc': ({'A': 'Fe'}, {'A': [(0, 0, 0), (0.5, 0.5, 0.5)]}),
    'fcc': ({'A': 'Cu'}, {'A': _FCC}),
    'diamond': ({'A': 'C'}, {'A': _FCC + _TETRA}),
    'cscl': ({'A': 'Cs', 'B': 'Cl'}, {'A': [(0.5, 0.5, 0.5)], 'B': [(0, 0, 0)]}),
    'nacl': ({'A': 'Na', 'B': 'Cl'}, {'A': [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5)], 'B': _FCC}),
    'zinc_blende': ({'A': 'Zn', 'B': 'S'}, {'A': _TETRA, 'B': _FCC}),
    'fluorite': ({'A': 'Ca', 'B': 'F'}, {'A': _FCC, 'B': _TETRA_ALL}),
    'perovskite': ({'A': 'Ca', 'B': 'Ti', 'X': 'O'},
                   {'A': [(0, 0, 0)], 'B': [(0.5, 0.5, 0.5)], 'X': [(0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]}),
}
_ALIASES = {'simple_cubic': 'sc', 'primitive': 'sc', 'body_centred': 'bcc', 'body_centered': 'bcc',
            'face_centred': 'fcc', 'face_centered': 'fcc', 'ccp': 'fcc', 'rock_salt': 'nacl', 'rocksalt': 'nacl',
            'zns': 'zinc_blende', 'sphalerite': 'zinc_blende', 'caf2': 'fluorite', 'catio3': 'perovskite'}
_DEFAULT_FILL = (BLUE, WHITE, ORANGE, BLACK)
EPS = 1e-6


def _expand(site):
    """All positions inside or on the cell that are equivalent to ``site`` by lattice translation."""
    base = [number(v, 'fractional coordinate', -1e-9, 1 + 1e-9) % 1.0 for v in site]
    options = [(0.0, 1.0) if abs(v) < EPS or abs(v - 1) < EPS else (v,) for v in base]
    return [tuple(combo) for combo in itertools.product(*options)]


def _share(position):
    on_boundary = sum(1 for v in position if abs(v) < EPS or abs(v - 1) < EPS)
    return Fraction(1, 2 ** on_boundary)


def _nearest(position, element, basis):
    """Nearest neighbours of ``position`` in the infinite crystal: (distance, counts per element)."""
    best = math.inf
    found = []
    for other_element, sites in basis.items():
        for site in sites:
            for shift in itertools.product((-1, 0, 1), repeat=3):
                d = math.dist(position, [s + t for s, t in zip(site, shift)])
                if d < EPS:
                    continue
                found.append((d, other_element))
                best = min(best, d)
    counts: dict[str, int] = {}
    for d, other_element in found:
        if d <= best * 1.02:
            counts[other_element] = counts.get(other_element, 0) + 1
    return best, counts


def unit_cell(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Cubic unit cell.

    ``{"preset": "NaCl"}`` (sc, bcc, fcc, diamond, CsCl, NaCl, zinc_blende, fluorite,
    perovskite), optionally with ``"elements": {"A": "K", "B": "Br"}`` to rename the
    species, or ``"sites": [{"element": "Na", "frac": [0.5, 0, 0]}, ...]`` with one entry
    per site of the basis. ``"bonds": true`` joins nearest unlike neighbours.
    """
    warnings: list[str] = []
    if spec.get('sites') is not None:
        basis: dict[str, list] = {}
        for item in spec['sites']:
            element = str(item.get('element', '')).strip()
            frac = item.get('frac')
            if not element or frac is None or len(frac) != 3:
                raise ValueError('Each site needs an "element" and "frac": [x, y, z]')
            basis.setdefault(element, []).append(tuple(float(v) % 1.0 for v in frac))
        name = str(spec.get('name', 'custom'))
    else:
        key = str(spec.get('preset', '')).lower().replace(' ', '_').replace('-', '_')
        key = _ALIASES.get(key, key)
        if key not in PRESETS:
            raise ValueError(f'preset must be one of: {", ".join(PRESETS)} (or give "sites")')
        roles, sites = PRESETS[key]
        names = {**roles, **{str(k): str(v) for k, v in (spec.get('elements') or {}).items()}}
        unknown = set(names) - set(roles)
        if unknown:
            raise ValueError(f'{key} has the roles {", ".join(roles)}; unknown: {", ".join(sorted(unknown))}')
        basis = {}
        for role, positions in sites.items():
            basis.setdefault(names[role], []).extend(positions)
        name = key
    for element, positions in basis.items():
        unique = {tuple(round(v, 6) for v in pos) for pos in positions}
        if len(unique) != len(positions):
            raise ValueError(f'{element}: the same site is listed twice')

    atoms = []
    contents: dict[str, Fraction] = {}
    for element, positions in basis.items():
        for site in positions:
            for position in _expand(site):
                atoms.append({'element': element, 'pos': position})
                contents[element] = contents.get(element, Fraction(0)) + _share(position)
    for element, count in contents.items():
        if count.denominator != 1:
            raise ValueError(f'{element} adds up to {count} atoms per cell; a basis must fill the cell with whole atoms')
    counts = {element: int(count) for element, count in contents.items()}
    divisor = reduce(math.gcd, counts.values())
    formula = ''.join(f'{element}{counts[element] // divisor if counts[element] // divisor > 1 else ""}'
                      for element in basis)
    coordination = {}
    for element, positions in basis.items():
        distance, neighbours = _nearest(positions[0], element, basis)
        coordination[element] = {'coordination_number': sum(neighbours.values()), 'neighbours': neighbours,
                                 'distance_in_cell_edges': round(distance, 4)}

    # ---- appearance
    edge = float(spec.get('edge', 92.0))
    if not 40 <= edge <= 400:
        raise ValueError('"edge" must be between 40 and 400 points')
    style = spec.get('style') or {}
    radii = {str(k): float(v) for k, v in (style.get('radius') or {}).items()}
    fills = {str(k): v for k, v in (style.get('fill') or {}).items()}
    look = {}
    for i, element in enumerate(basis):
        fill = fills.get(element, _DEFAULT_FILL[i % len(_DEFAULT_FILL)])
        look[element] = {'radius': radii.get(element, 5.0 if len(basis[element]) > 4 else 6.0),
                         'fill': None if fill in (None, 'none') else resolve_color(fill)}
    # Chosen so that no two sites of the built-in structures project closer than 0.16 of the edge.
    depth_x, depth_y = 0.50, 0.20

    def project(position):
        x, y, z = position
        return edge * (x + depth_x * y), -edge * (z + depth_y * y)

    c = Canvas()
    corners = list(itertools.product((0, 1), repeat=3))
    hidden = (0, 1, 0)
    back, front = [], []
    for a, b in itertools.combinations(corners, 2):
        if sum(abs(p - q) for p, q in zip(a, b)) != 1:
            continue
        if hidden in (a, b):
            c.segment(project(a), project(b), dashed=True)
        elif a[1] == 0 and b[1] == 0:
            front.append((a, b))
        else:
            back.append((a, b))
    for a, b in back:
        c.segment(project(a), project(b))

    bonds_drawn = 0
    if spec.get('bonds'):
        if len(basis) < 2 and name != 'diamond':
            warnings.append('"bonds" joins nearest neighbours; in a one-element lattice only diamond has them drawn')
        reach = min(info['distance_in_cell_edges'] for info in coordination.values())
        for a, b in itertools.combinations(atoms, 2):
            if abs(math.dist(a['pos'], b['pos']) - reach) < 1e-3 and (a['element'] != b['element'] or name == 'diamond'):
                c.segment(project(a['pos']), project(b['pos']))
                bonds_drawn += 1

    def paint(group):
        for atom in sorted(group, key=lambda item: (-item['pos'][1], item['pos'][2], item['pos'][0])):
            style_of = look[atom['element']]
            c.circle(project(atom['pos']), style_of['radius'], fill=style_of['fill'] or WHITE, outline=BLACK)

    paint([atom for atom in atoms if atom['pos'][1] > EPS])
    for a, b in front:
        c.segment(project(a), project(b))
    paint([atom for atom in atoms if atom['pos'][1] <= EPS])

    if spec.get('legend', True):
        x = c.extent[2] + 18.0
        y = c.extent[1] + 12.0
        for element, style_of in look.items():
            c.circle((x, y), style_of['radius'], fill=style_of['fill'] or WHITE, outline=BLACK)
            c.label((x + style_of['radius'] + 5.0, y), element, valign='middle')
            y += 2 * max(style_of['radius'], 6.0) + 6.0
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report = {'diagram': 'unit_cell', 'structure': name, 'atoms_per_cell': counts, 'formula': formula,
              'coordination': coordination, 'atoms_drawn': len(atoms), 'bonds_drawn': bonds_drawn,
              'warnings': warnings}
    return _result(c, output, report)
