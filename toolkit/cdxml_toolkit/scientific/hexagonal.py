"""Hexagonal structures drawn as the textbook hexagonal prism (three primitive cells).

Sites are given in the primitive hexagonal cell (a1 = a(1, 0), a2 = a(-1/2, sqrt3/2), c);
the builder finds every copy inside the prism, weights each by the share the prism owns
(vertex of the hexagon 1/3, hexagon edge 1/2, top or bottom face 1/2, so a prism corner is
1/6) and reports atoms per prism, formula, coordination numbers and, for one element, the
packing fraction.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction
from functools import reduce
from typing import Any, Mapping

from .canvas import BLACK, BLUE, ORANGE, WHITE, Canvas, resolve_color
from .energy import _result

IDEAL_C = math.sqrt(8.0 / 3.0)        # c/a of ideal close packing
# preset -> (default element per role, sites per role in the primitive cell, c/a)
HEX_PRESETS = {
    'hcp': ({'A': 'Mg'}, {'A': [(0, 0, 0), (1 / 3, 2 / 3, 0.5)]}, IDEAL_C),
    'wurtzite': ({'A': 'Zn', 'B': 'S'}, {'B': [(0, 0, 0), (1 / 3, 2 / 3, 0.5)],
                                         'A': [(0, 0, 0.375), (1 / 3, 2 / 3, 0.875)]}, IDEAL_C),
    'simple_hexagonal': ({'A': 'X'}, {'A': [(0, 0, 0)]}, 1.0),
}
HEX_ALIASES = {'hexagonal_close_packed': 'hcp', 'hexagonal_closest_packed': 'hcp', 'zno': 'wurtzite',
               'hexagonal': 'simple_hexagonal'}
EPS = 1e-6
_DEFAULT_FILL = (BLUE, WHITE, ORANGE, BLACK)
_SQ3 = math.sqrt(3.0)


def _cartesian(frac, c_ratio):
    u, v, w = frac
    return (u - v / 2.0, v * _SQ3 / 2.0, w * c_ratio)


def _hex_share(x, y):
    """Share of a point the hexagon (circumradius 1, vertices at 0, 60, ... degrees) owns; 0 outside."""
    distances = []
    for k in range(6):                         # signed distance to each edge, inside positive
        angle = math.radians(30 + 60 * k)
        distances.append(_SQ3 / 2.0 - (x * math.cos(angle) + y * math.sin(angle)))
    if min(distances) < -EPS:
        return Fraction(0)
    on = sum(1 for d in distances if abs(d) < EPS)
    return Fraction(1) if on == 0 else Fraction(1, 2) if on == 1 else Fraction(1, 3)


def hexagonal_cell(key: str, spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    roles, sites, c_ratio = HEX_PRESETS[key]
    c_ratio = float(spec.get('c_over_a', c_ratio))
    if not 0.5 <= c_ratio <= 3.0:
        raise ValueError('"c_over_a" must be between 0.5 and 3')
    names = {**roles, **{str(k): str(v) for k, v in (spec.get('elements') or {}).items()}}
    unknown = set(names) - set(roles)
    if unknown:
        raise ValueError(f'{key} has the roles {", ".join(roles)}; unknown: {", ".join(sorted(unknown))}')
    basis: dict[str, list] = {}
    for role, positions in sites.items():
        basis.setdefault(names[role], []).extend(positions)

    # ---- every copy inside the prism, with the share the prism owns
    atoms = []
    contents: dict[str, Fraction] = {}
    for element, positions in basis.items():
        for site in positions:
            for i, j, k in itertools.product(range(-2, 3), range(-2, 3), (0, 1)):
                frac = (site[0] + i, site[1] + j, site[2] + k)
                if frac[2] > 1 + EPS:
                    continue
                x, y, z = _cartesian(frac, c_ratio)
                share = _hex_share(x, y)
                if not share:
                    continue
                if abs(z) < EPS or abs(z - c_ratio) < EPS:
                    share /= 2
                atoms.append({'element': element, 'pos': (x, y, z)})
                contents[element] = contents.get(element, Fraction(0)) + share
    counts = {element: int(count) if count.denominator == 1 else float(count) for element, count in contents.items()}
    whole = {e: int(v) for e, v in counts.items()}
    divisor = reduce(math.gcd, whole.values())
    formula = ''.join(f'{e}{whole[e] // divisor if whole[e] // divisor > 1 else ""}' for e in basis)

    # ---- coordination in the infinite crystal
    lattice = []
    for element, positions in basis.items():
        for site in positions:
            for i, j, k in itertools.product((-1, 0, 1), repeat=3):
                lattice.append((element, _cartesian((site[0] + i, site[1] + j, site[2] + k), c_ratio)))
    coordination = {}
    for element, positions in basis.items():
        here = _cartesian(positions[0], c_ratio)
        found = [(math.dist(here, p), e) for e, p in lattice if math.dist(here, p) > EPS]
        nearest = min(d for d, _e in found)
        neighbours: dict[str, int] = {}
        for d, e in found:
            if d <= nearest * 1.02:
                neighbours[e] = neighbours.get(e, 0) + 1
        coordination[element] = {'coordination_number': sum(neighbours.values()), 'neighbours': neighbours,
                                 'distance_in_a': round(nearest, 4)}
    report: dict[str, Any] = {'diagram': 'unit_cell', 'structure': key, 'cell': 'hexagonal prism (3 primitive cells)',
                              'atoms_per_cell': whole, 'formula': formula, 'c_over_a': round(c_ratio, 4),
                              'coordination': coordination, 'atoms_drawn': len(atoms)}
    if len(basis) == 1:
        element = next(iter(basis))
        radius = coordination[element]['distance_in_a'] / 2
        volume = 3 * _SQ3 / 2 * c_ratio
        report['packing_fraction'] = round(whole[element] * 4 / 3 * math.pi * radius ** 3 / volume, 4)

    # ---- drawing: oblique view from slightly above
    edge = float(spec.get('edge', 56.0))
    if not 30 <= edge <= 300:
        raise ValueError('"edge" must be between 30 and 300 points')
    # Seen from slightly above and to the side: with this tilt and shear no two sites of the
    # built-in structures project closer than 0.22 a, so no atom hides another.
    tilt, shear = 0.36, 0.36

    def project(p):
        return edge * (p[0] + shear * p[1]), -edge * (p[2] - tilt * p[1])

    c = Canvas()
    corners = [(math.cos(math.radians(60 * k)), math.sin(math.radians(60 * k))) for k in range(6)]
    back = {k for k, (_x, y) in enumerate(corners) if y > EPS}
    for z in (0.0, c_ratio):
        for k in range(6):
            a, b = corners[k], corners[(k + 1) % 6]
            hidden = z == 0.0 and (k in back or (k + 1) % 6 in back)
            c.segment(project((*a, z)), project((*b, z)), dashed=hidden)
    for k, (x, y) in enumerate(corners):
        c.segment(project((x, y, 0.0)), project((x, y, c_ratio)), dashed=k in back)
    style = spec.get('style') or {}
    radii = {str(k): float(v) for k, v in (style.get('radius') or {}).items()}
    fills = {str(k): v for k, v in (style.get('fill') or {}).items()}
    look = {}
    for i, element in enumerate(basis):
        fill = fills.get(element, _DEFAULT_FILL[i % len(_DEFAULT_FILL)])
        look[element] = {'radius': radii.get(element, 4.5 if len(basis) > 1 else 5.0),
                         'fill': WHITE if fill in (None, 'none') else resolve_color(fill)}
    if spec.get('bonds'):
        reach = min(info['distance_in_a'] for info in coordination.values())
        for a, b in itertools.combinations(atoms, 2):
            if a['element'] != b['element'] and abs(math.dist(a['pos'], b['pos']) - reach) < 1e-3:
                c.segment(project(a['pos']), project(b['pos']))
    for atom in sorted(atoms, key=lambda item: (-item['pos'][1], item['pos'][2])):
        style_of = look[atom['element']]
        c.circle(project(atom['pos']), style_of['radius'], fill=style_of['fill'], outline=BLACK)
    if spec.get('legend', True) and len(basis) > 1:
        x, y = c.extent[2] + 16.0, c.extent[1] + 10.0
        for element, style_of in look.items():
            c.circle((x, y), style_of['radius'], fill=style_of['fill'], outline=BLACK)
            c.label((x + style_of['radius'] + 5.0, y), element, valign='middle')
            y += 2 * style_of['radius'] + 8.0
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report['warnings'] = []
    return _result(c, output, report)
