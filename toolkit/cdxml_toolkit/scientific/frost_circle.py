"""Frost polygon construction for cyclic HMO levels; redox Frost diagrams live in redox.py."""
from __future__ import annotations

import math
from typing import Any, Mapping

from .canvas import Canvas
from .energy import _result
from .orbitals import _spin_arrows, hund_fill

_RULES = (
    '1. Draw a circle of radius 2|\u03b2| about \u03b1.',
    '2. Inscribe a regular n-gon with one vertex down.',
    '3. Project the vertices; fill from the lowest level (Hund).',
)


def _integer(value, name, lower, upper):
    if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
        raise ValueError(f'{name} must be an integer from {lower} to {upper}')
    return value


def frost_circle(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Draw a cyclic pi HMO Frost circle, distinct from the redox ``frost`` operation.

    ``centres``: 3..12, ``electrons``: 0..2n; ``show_construction`` (default true)
    adds projections, centre, radius and rule ``lines``. ``variant`` is ``vertex-down``
    (default), ``vertex-up`` (odd n only) or ``flat-bottom``; the last two are marked
    incorrect and their projections are illustrative, not HMO predictions.
    """
    n = _integer(spec.get('centres'), 'centres', 3, 12)
    electrons = _integer(spec.get('electrons'), 'electrons', 0, 2 * n)
    show = spec.get('show_construction', True)
    if not isinstance(show, bool):
        raise ValueError('show_construction must be a boolean')
    variant = spec.get('variant', 'vertex-down')
    if variant not in ('vertex-down', 'vertex-up', 'flat-bottom'):
        raise ValueError('variant must be vertex-down, vertex-up or flat-bottom')
    if variant == 'vertex-up' and n % 2 == 0:
        raise ValueError('For even centres, vertex-up is also vertex-down; use flat-bottom for an incorrect example')
    lines = spec.get('lines', list(_RULES))
    if not isinstance(lines, list) or not all(isinstance(line, str) for line in lines):
        raise ValueError('lines must be a list of strings')
    valid = variant == 'vertex-down'
    rotation = {'vertex-down': 0.0, 'vertex-up': math.pi, 'flat-bottom': math.pi / n}[variant]
    radius = round(max(72.0, 20.0 / (1.0 - math.cos(2 * math.pi / n))), 2)
    vertices = [(radius * math.sin(rotation + 2 * math.pi * k / n),
                 radius * math.cos(rotation + 2 * math.pi * k / n)) for k in range(n)]
    groups = []
    for k in sorted(range(n), key=lambda k: -vertices[k][1]):
        coefficient = 2 * vertices[k][1] / radius
        if abs(coefficient) < 1e-10:
            coefficient = 0.0
        if groups and abs(groups[-1]['x'] - coefficient) < 1e-9:
            groups[-1]['vertices'].append(k)
        else:
            groups.append({'x': coefficient, 'vertices': [k]})
    remaining = electrons
    for group in groups:
        group['degeneracy'] = len(group['vertices'])
        take = min(remaining, 2 * group['degeneracy'])
        group['occupation'] = hund_fill(take, group['degeneracy'])
        remaining -= take

    c = Canvas()
    axis = radius + 42.0
    left, right = axis + 22.0, axis + 70.0
    c.circle((0.0, 0.0), radius)
    c.path(vertices, closed=True)
    c.segment((-radius - 8.0, 0.0), (axis + 6.0, 0.0))
    c.label((-radius - 12.0, 0.0), '\u03b1', align='Right', valign='middle')
    c.arrow((axis, radius + 10.0), (axis, -radius - 18.0), head_size=700)
    c.label((axis, -radius - 24.0), 'E', align='Center')
    if show:
        for vertex in vertices:
            c.segment(vertex, (left - 3.0, vertex[1]), dashed=True)
        c.circle((0.0, 0.0), 1.5, fill='black', outline=None)
        c.arrow((0.0, 0.0), (-radius, 0.0), tail_style='Full', head_size=600)
        c.label((-radius / 2, -8.0), 'r = 2|\u03b2|', align='Center')
    for group in groups:
        y = radius * group['x'] / 2
        c.segment((left, y), (right, y))
        for slot, occupation in enumerate(group['occupation']):
            x = left + (slot + 0.5) * (right - left) / group['degeneracy']
            _spin_arrows(c, x, y, occupation)
        text = f'x = {group["x"]:.6f}; d = {group["degeneracy"]}'
        c.label((right + 10.0, y), text, valign='middle')
    caption_y = radius + 32.0
    c.label((-radius, caption_y), f'n = {n}; \u03c0 electrons = {electrons}; \u03b2 < 0')
    c.label((-radius, caption_y + 16.0), '\u03b5 = \u03b1 + x\u03b2; y = y(\u03b1) + rx/2')
    if show:
        for index, line in enumerate(lines):
            c.label((-radius, caption_y + 36.0 + 16.0 * index), line)
    warnings = []
    if not valid:
        warnings.append(f'Incorrect {variant} construction: projected heights are not cyclic HMO energies')
        c.segment((-radius, -radius - 28.0), (-radius + 12.0, -radius - 16.0), color='#FF0000')
        c.segment((-radius, -radius - 16.0), (-radius + 12.0, -radius - 28.0), color='#FF0000')
        c.label((-radius + 20.0, -radius - 22.0), f'Incorrect: {variant}', valign='middle', color='#FF0000')
    if spec.get('title'):
        _, _, _, bottom = c.extent
        c.label((-radius, bottom + 20.0), str(spec['title']))
    rule = ('4q+2' if electrons >= 2 and electrons % 4 == 2 else
            '4q' if electrons > 0 and electrons % 4 == 0 else 'other') if valid else None
    report = {
        'diagram': 'frost_circle', 'centres': n, 'electrons': electrons,
        'variant': variant, 'valid_construction': valid, 'radius_pt': radius,
        'energy_convention': 'epsilon = alpha + x beta; beta < 0; circle radius = 2|beta|',
        'coordinate_mapping': 'y = y(alpha) + radius_pt * x / 2 (CDXML y points down)',
        'levels': groups, 'unpaired_electrons': sum(occupation == 'u' for group in groups
                                                   for occupation in group['occupation']),
        'huckel_rule': rule,
        'model': 'Ideal planar regular cyclic pi system; electron count alone does not establish molecular aromaticity',
        'warnings': warnings,
    }
    return _result(c, output, report)
