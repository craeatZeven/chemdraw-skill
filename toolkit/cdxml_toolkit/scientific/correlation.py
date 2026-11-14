"""Orbital correlation diagrams (Woodward-Hoffmann) and Walsh diagrams.

A correlation diagram joins reactant and product orbitals of the same symmetry, the lowest of
each symmetry to the lowest, never crossing within a symmetry (the non-crossing rule). The
reaction is thermally allowed when the occupied reactant orbitals correlate with the lowest
product orbitals; otherwise the ground state correlates with an excited state.
"""
from __future__ import annotations

from typing import Any, Mapping

from .canvas import BLUE, ORANGE, Canvas
from .energy import _result
from .orbitals import _spin_arrows

# name -> (element of symmetry, reactant levels low->high, product levels low->high, electrons)
_ELECTROCYCLIC = {
    'butadiene_conrotatory': ('C{sub:2}', [('\u03c8{sub:1}', 'A'), ('\u03c8{sub:2}', 'S'), ('\u03c8{sub:3}', 'A'), ('\u03c8{sub:4}', 'S')],
                              [('\u03c3', 'S'), ('\u03c0', 'A'), ('\u03c0*', 'S'), ('\u03c3*', 'A')], 4),
    'butadiene_disrotatory': ('\u03c3 (mirror)', [('\u03c8{sub:1}', 'S'), ('\u03c8{sub:2}', 'A'), ('\u03c8{sub:3}', 'S'), ('\u03c8{sub:4}', 'A')],
                              [('\u03c3', 'S'), ('\u03c0', 'S'), ('\u03c0*', 'A'), ('\u03c3*', 'A')], 4),
    'hexatriene_conrotatory': ('C{sub:2}', [('\u03c8{sub:%d}' % k, 'AS'[(k - 1) % 2]) for k in range(1, 7)],
                               [('\u03c3', 'S'), ('\u03c0{sub:1}', 'A'), ('\u03c0{sub:2}', 'S'), ('\u03c0{sub:3}*', 'A'),
                                ('\u03c0{sub:4}*', 'S'), ('\u03c3*', 'A')], 6),
    'hexatriene_disrotatory': ('\u03c3 (mirror)', [('\u03c8{sub:%d}' % k, 'SA'[(k - 1) % 2]) for k in range(1, 7)],
                               [('\u03c3', 'S'), ('\u03c0{sub:1}', 'S'), ('\u03c0{sub:2}', 'A'), ('\u03c0{sub:3}*', 'S'),
                                ('\u03c0{sub:4}*', 'A'), ('\u03c3*', 'A')], 6),
    'cycloaddition_2+2': ('\u03c3{sub:1}, \u03c3{sub:2}', [('\u03c0 + \u03c0', 'SS'), ('\u03c0 \u2212 \u03c0', 'SA'),
                                                       ('\u03c0* + \u03c0*', 'AS'), ('\u03c0* \u2212 \u03c0*', 'AA')],
                          [('\u03c3{sub:1} + \u03c3{sub:2}', 'SS'), ('\u03c3{sub:1} \u2212 \u03c3{sub:2}', 'AS'),
                           ('\u03c3{sub:1}* + \u03c3{sub:2}*', 'SA'), ('\u03c3{sub:1}* \u2212 \u03c3{sub:2}*', 'AA')], 4),
    # H\u00fcckel order: the diene's psi1 (alpha + 1.62 beta) lies below the ethene pi (alpha + beta).
    'cycloaddition_4+2': ('\u03c3 (mirror)', [('\u03c8{sub:1}', 'S'), ('\u03c0 (ethene)', 'S'), ('\u03c8{sub:2}', 'A'),
                                               ('\u03c8{sub:3}', 'S'), ('\u03c0* (ethene)', 'A'), ('\u03c8{sub:4}', 'A')],
                          [('\u03c3{sub:1}', 'S'), ('\u03c3{sub:2}', 'A'), ('\u03c0', 'S'), ('\u03c0*', 'A'),
                           ('\u03c3{sub:1}*', 'S'), ('\u03c3{sub:2}*', 'A')], 6),
}


def correlation_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Orbital correlation diagram.

    ``{"preset": "butadiene_conrotatory"}`` (also ``_disrotatory``, ``hexatriene_*``,
    ``cycloaddition_2+2``, ``cycloaddition_4+2``) or ``left``/``right``:
    ``[{"label", "symmetry"}]`` from the lowest orbital up, ``electrons`` and an optional
    ``symmetry_element`` caption. Energies are a schematic order.
    """
    if spec.get('preset') is not None:
        key = str(spec['preset']).lower().replace(' ', '_')
        if key not in _ELECTROCYCLIC:
            raise ValueError(f'Unknown preset {spec["preset"]!r}; available: {", ".join(_ELECTROCYCLIC)}')
        element, left, right, electrons = _ELECTROCYCLIC[key]
        left = [{'label': a, 'symmetry': b} for a, b in left]
        right = [{'label': a, 'symmetry': b} for a, b in right]
        spec = {'symmetry_element': element, 'left_title': spec.get('left_title', 'reactant'),
                'right_title': spec.get('right_title', 'product'), **spec, 'left': left, 'right': right,
                'electrons': electrons}
    left, right = spec.get('left'), spec.get('right')
    if not isinstance(left, list) or not isinstance(right, list) or not left or not right:
        raise ValueError('A correlation diagram needs "left" and "right" orbital lists')
    for side in (left, right):
        for item in side:
            if not item.get('symmetry'):
                raise ValueError('Every orbital needs a "symmetry" label')
    sym_l = sorted({str(o['symmetry']) for o in left})
    sym_r = sorted({str(o['symmetry']) for o in right})
    counts_l = {s: sum(1 for o in left if str(o['symmetry']) == s) for s in sym_l}
    counts_r = {s: sum(1 for o in right if str(o['symmetry']) == s) for s in sym_r}
    if counts_l != counts_r:
        raise ValueError(f'Orbital symmetry is conserved: the reactant has {counts_l} and the product {counts_r}')
    electrons = int(spec.get('electrons', 0))
    if not 0 <= electrons <= 2 * len(left) or electrons % 2:
        raise ValueError('"electrons" is an even number that fits the reactant orbitals')
    # Non-crossing rule: within each symmetry, k-th lowest joins k-th lowest.
    links = []
    for s in sym_l:
        a = [i for i, o in enumerate(left) if str(o['symmetry']) == s]
        b = [i for i, o in enumerate(right) if str(o['symmetry']) == s]
        links += list(zip(a, b))
    occupied = list(range(electrons // 2))
    targets = sorted(j for i, j in links if i in occupied)
    allowed = targets == occupied

    c = Canvas()
    step, gap, bar = 26.0, 130.0, 26.0
    for k, o in enumerate(left):
        o['y'] = -k * step
    for k, o in enumerate(right):
        o['y'] = -k * step
    for side, x0 in ((left, 0.0), (right, gap)):
        for k, o in enumerate(side):
            c.segment((x0, o['y']), (x0 + bar, o['y']))
            if side is left and k < len(occupied):
                _spin_arrows(c, x0 + bar / 2, o['y'], 'ud')
            if side is left:
                c.label((x0 - 5.0, o['y']), f'{o["label"]}', align='Right', valign='middle')
                c.label((x0 - 5.0 - 40.0, o['y']), str(o['symmetry']), align='Right', valign='middle', size=8)
            else:
                c.label((x0 + bar + 5.0, o['y']), f'{o["label"]}', valign='middle')
                c.label((x0 + bar + 45.0, o['y']), str(o['symmetry']), valign='middle', size=8)
    for i, j in links:
        bad = i in occupied and j not in occupied
        c.segment((bar + 2.0, left[i]['y']), (gap - 2.0, right[j]['y']), dashed=True, color=ORANGE if bad else BLUE)
    bottom = 18.0
    if spec.get('left_title'):
        c.label((bar / 2, bottom), str(spec['left_title']), align='Center')
    if spec.get('right_title'):
        c.label((gap + bar / 2, bottom), str(spec['right_title']), align='Center')
    verdict = 'thermally allowed' if allowed else 'thermally forbidden (photochemically allowed)'
    caption = (f'{spec["symmetry_element"]}: ' if spec.get('symmetry_element') else '') + verdict
    c.label(((gap + bar) / 2, bottom + 16.0), caption, align='Center')
    report = {'diagram': 'correlation', 'allowed_thermally': allowed, 'verdict': verdict,
              'links': [{'from': left[i]['label'], 'to': right[j]['label'], 'symmetry': left[i]['symmetry']}
                        for i, j in links],
              'ground_state_correlates_with': [right[j]['label'] for j in targets], 'warnings': []}
    return _result(c, output, report)


# ---------------------------------------------------------------------------- Walsh diagram
# AH2, bent (90 degrees) to linear (180 degrees): qualitative energies at 90, 105, 120, 135, 150, 165, 180.
_AH2 = {
    'x': [90, 105, 120, 135, 150, 165, 180],
    'orbitals': [
        {'left': '1a{sub:1}', 'right': '2\u03c3{sub:g}', 'energy': [-1.06, -1.05, -1.04, -1.03, -1.02, -1.01, -1.00]},
        {'left': '1b{sub:2}', 'right': '1\u03c3{sub:u}', 'energy': [-0.48, -0.55, -0.61, -0.66, -0.70, -0.73, -0.75]},
        {'left': '2a{sub:1}', 'right': '1\u03c0{sub:u}', 'energy': [-0.62, -0.58, -0.53, -0.48, -0.42, -0.36, -0.30]},
        {'left': '1b{sub:1}', 'right': '1\u03c0{sub:u}', 'energy': [-0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30]},
        {'left': '3a{sub:1}', 'right': '3\u03c3{sub:g}', 'energy': [0.30, 0.25, 0.20, 0.16, 0.13, 0.11, 0.10]},
        {'left': '2b{sub:2}', 'right': '2\u03c3{sub:u}', 'energy': [0.28, 0.32, 0.36, 0.40, 0.43, 0.45, 0.46]},
    ],
}


def walsh_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Walsh diagram: orbital energies against a geometric coordinate.

    ``{"preset": "AH2", "electrons": 8}`` (valence electrons; BeH2 4, BH2 5, CH2 6, NH2 7,
    H2O 8) or ``{"x": [...], "orbitals": [{"left", "right", "energy": [...]}], "electrons"}``.
    The sum of the occupied orbital energies along the coordinate gives the preferred end.
    Preset energies are qualitative.
    """
    data = dict(spec)
    if spec.get('preset') is not None:
        if str(spec['preset']).upper() != 'AH2':
            raise ValueError('The Walsh preset is "AH2"')
        data = {**_AH2, 'x_label': 'H\u2013A\u2013H angle / \u00b0', **{k: v for k, v in spec.items() if k != 'preset'}}
    xs = [float(v) for v in data.get('x') or []]
    orbitals = data.get('orbitals') or []
    if len(xs) < 2 or not orbitals or any(len(o.get('energy', [])) != len(xs) for o in orbitals):
        raise ValueError('Give "x" and, for every orbital, one "energy" per x value')
    electrons = int(data.get('electrons', 0))
    if not 0 <= electrons <= 2 * len(orbitals):
        raise ValueError(f'{electrons} electrons do not fit {len(orbitals)} orbitals')
    # Fill in energy order at each x (orbitals may swap order along the coordinate).
    totals = []
    for k in range(len(xs)):
        ranked = sorted(range(len(orbitals)), key=lambda i: orbitals[i]['energy'][k])
        left = electrons
        total = 0.0
        for i in ranked:
            take = min(2, left)
            total += take * orbitals[i]['energy'][k]
            left -= take
        totals.append(total)
    best = min(range(len(xs)), key=lambda k: totals[k])
    from .axes import Plot
    c = Canvas()
    energies = [e for o in orbitals for e in o['energy']]
    pad = (max(energies) - min(energies)) * 0.08
    plot = Plot(c, (min(xs), max(xs)), (min(energies) - pad, max(energies) + pad), 200.0, 170.0)
    plot.frame(data.get('x_label', 'coordinate'), data.get('y_label', 'orbital energy'), ticks=False)
    from .axes import nice_ticks, tick_text
    ticks, step = nice_ticks(min(xs), max(xs), 4)
    for value in ticks:
        px, _ = plot.at(value, plot.y0)
        c.label((px, plot.h + 10.0), tick_text(value, step), align='Center', size=8)
    end_rank = sorted(range(len(orbitals)), key=lambda i: orbitals[i]['energy'][0])
    occupied_left = set(end_rank[:(electrons + 1) // 2])
    for i, o in enumerate(orbitals):
        pts = []
        for k in range(len(xs) - 1):
            for t in range(10):
                f = t / 10
                pts.append((xs[k] + (xs[k + 1] - xs[k]) * f, o['energy'][k] + (o['energy'][k + 1] - o['energy'][k]) * f))
        pts.append((xs[-1], o['energy'][-1]))
        plot.polyline(pts, color=BLUE if i in occupied_left else None, dashed=i not in occupied_left)
        lx, ly = plot.at(xs[0], o['energy'][0])
        rx, ry = plot.at(xs[-1], o['energy'][-1])
        c.place([((lx - 4.0, ly + d), 'Right', 'middle') for d in (0, -7, 7)], str(o.get('left', '')), size=8)
        c.place([((rx + 4.0, ry + d), 'Left', 'middle') for d in (0, -7, 7, -14, 14)], str(o.get('right', '')), size=8)
    report = {'diagram': 'walsh', 'electrons': electrons, 'total_energy': [round(t, 4) for t in totals],
              'preferred_x': xs[best], 'preferred': ('the low end' if best == 0 else 'the high end' if best == len(xs) - 1
                                                     else 'in between'),
              'note': 'preset energies are qualitative; the trend, not the angle, is the result', 'warnings': []}
    if spec.get('preset') is not None:
        report['geometry'] = 'linear' if xs[best] >= 170 else 'bent'
    return _result(c, output, report)
