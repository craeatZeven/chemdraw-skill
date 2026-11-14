"""Redox and thermodynamic stability diagrams computed from the data in the spec.

* ``pourbaix``: potential-pH diagram of one element from standard Gibbs energies of formation;
  every field is the region where its species has the lowest free energy, found exactly by
  clipping the plot with half-planes (all boundaries are straight lines).
* ``latimer`` and ``frost``: from a chain of standard reduction potentials; non-adjacent
  couples are derived (sum of n E over sum of n) and species that disproportionate are found.
* ``ellingham``: Delta G = Delta H - T Delta S lines per mole of O2; crossings give the
  temperature above which one element reduces the other's oxide.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

from .axes import Plot, clip_polygon
from .canvas import BLUE, ORANGE, SERIES_COLORS, Canvas, boxes_overlap as _overlap, resolve_color
from .common import number
from .energy import MINUS, _result

F = 96485.332          # C/mol
R = 8.314462618        # J/(mol K)
G_WATER = -237.13      # kJ/mol, H2O(l)


# ---------------------------------------------------------------------------- Pourbaix
def pourbaix(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Potential-pH (Pourbaix) diagram of one element.

    ``species``: ``{"name", "dGf" (kJ/mol, per formula unit), "atoms" (of the element, default
    1), "O", "H", "charge", "aqueous"}``. The element itself (dGf 0, no O, H or charge) is
    added when it is missing. ``activity`` (default 1e-6) applies to dissolved species;
    ``T`` (default 298.15 K); ``pH`` and ``E`` ranges; ``points``: ``[{"pH", "E", "label"}]``
    are reported with their stable species. Water's stability limits are drawn dashed.
    """
    element = str(spec.get('element', 'M'))
    temperature = number(spec.get('T', 298.15), 'T', 200.0, 600.0)
    activity = number(spec.get('activity', 1e-6), 'activity', 1e-12, 1.0)
    slope = math.log(10) * R * temperature / F          # V per pH unit, 0.0592 at 298 K
    ph0, ph1 = (number(v, 'pH', -4.0, 18.0) for v in spec.get('pH', [0, 14]))
    e0, e1 = (number(v, 'E', -5.0, 5.0) for v in spec.get('E', [-1.5, 2.0]))
    raw = spec.get('species')
    if not isinstance(raw, list) or len(raw) < 2:
        raise ValueError('A Pourbaix diagram needs at least two "species"')
    species = []
    for i, item in enumerate(raw):
        name = str(item.get('name') or f'species {i + 1}')
        atoms = number(item.get('atoms', 1), f'{name} atoms', 0.01, 100.0)
        o, h = number(item.get('O', 0), f'{name} O', 0, 200) / atoms, number(item.get('H', 0), f'{name} H', 0, 400) / atoms
        z = number(item.get('charge', 0), f'{name} charge', -20, 20) / atoms
        dgf = number(item.get('dGf'), f'{name} dGf', -1e5, 1e5) * 1000 / atoms   # J per atom of the element
        aqueous = bool(item.get('aqueous', item.get('charge', 0) != 0))
        # M + o H2O -> species + p H+ + n e-   (per atom of M)
        p = 2 * o - h
        n = z + p
        dg = dgf - o * G_WATER * 1000
        if aqueous:
            dg += R * temperature * math.log(activity) / atoms
        # Free energy per atom of M, in volts: g(pH, E) = dg / F - p (ln 10) RT/F pH - n E
        # Released H+ lowers the free energy as pH rises (-p slope pH); released electrons at potential E add -n E.
        species.append({'name': name, 'a': -p * slope, 'b': -n, 'c': dg / F, 'n': n, 'p': p,
                        'aqueous': aqueous})
    if not any(abs(s['n']) < 1e-9 and abs(s['p']) < 1e-9 and abs(s['c']) < 1e-9 for s in species):
        species.append({'name': element, 'a': 0.0, 'b': 0.0, 'c': 0.0, 'n': 0.0, 'p': 0.0, 'aqueous': False})
    # g_i(pH, E) = c_i + a_i pH + b_i E ; stable species minimises g.
    box = [(ph0, e0), (ph1, e0), (ph1, e1), (ph0, e1)]
    fields = []
    for i, s in enumerate(species):
        polygon = list(box)
        for j, t in enumerate(species):
            if i == j or not polygon:
                continue
            # g_i <= g_j  ->  (a_i - a_j) pH + (b_i - b_j) E <= c_j - c_i
            polygon = clip_polygon(polygon, s['a'] - t['a'], s['b'] - t['b'], t['c'] - s['c'])
        area = 0.0
        for k in range(len(polygon)):
            (x1, y1), (x2, y2) = polygon[k], polygon[(k + 1) % len(polygon)]
            area += x1 * y2 - x2 * y1
        if len(polygon) >= 3 and abs(area) > 1e-6:
            fields.append({'species': s, 'polygon': polygon})

    c = Canvas()
    plot = Plot(c, (ph0, ph1), (e0, e1), float(spec.get('width', 220.0)), float(spec.get('height', 180.0)))
    plot.frame(spec.get('x_label', 'pH'), spec.get('y_label', 'E / V vs SHE'))
    edges = {}
    for field in fields:
        poly = field['polygon']
        for k in range(len(poly)):
            a, b = poly[k], poly[(k + 1) % len(poly)]
            on_frame = any(abs(a[k] - edge) < 1e-9 and abs(b[k] - edge) < 1e-9
                           for k, edge in ((0, ph0), (0, ph1), (1, e0), (1, e1)))
            if on_frame or math.dist(a, b) < 1e-9:
                continue
            key = tuple(sorted((tuple(round(v, 6) for v in a), tuple(round(v, 6) for v in b))))
            edges[key] = (a, b)
    for a, b in edges.values():
        plot.segment(a, b)
    if spec.get('water', True):
        # Water's limits at 1 bar O2 and H2 (E\u00b0 of O2/H2O taken at 298 K).
        for intercept, name in ((1.229, 'O{sub:2}/H{sub:2}O'), (0.0, 'H{sub:2}O/H{sub:2}')):
            ends = ((ph0, intercept - slope * ph0), (ph1, intercept - slope * ph1))
            if plot.segment(*ends, dashed=True, color=BLUE):
                px, py = plot.at(ph1, intercept - slope * ph1)
                if 0 <= py <= plot.h:
                    c.place([((px - 3.0, py - 4.0), 'Right', 'baseline'), ((px - 3.0, py + 10.0), 'Right', 'baseline')],
                            name, color=BLUE, size=8)
    warnings = []
    for field in fields:
        poly = field['polygon']
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        ring = [(cx + (p[0] - cx) * f, cy + (p[1] - cy) * f) for f in (0.3, 0.55) for p in poly]
        candidates = [(plot.at(cx, cy), 'Center', 'middle')] + [(plot.at(*p), 'Center', 'middle') for p in ring]
        name = field['species']['name']
        inside = next((cand for cand in candidates if not any(
            _overlap(c.measure(cand[0], name, None, cand[1], cand[2]), o) for o in c.obstacles)), None)
        if inside is not None:
            c.place([inside], name)
            continue
        # A sliver of a field: name it outside the frame, with a leader to the field's centre.
        fx, fy = plot.at(cx, cy)
        outside = ([((plot.w + 10.0, fy + dy), 'Left', 'middle') for dy in (0, -12, 12, -24, 24)]
                   + [((fx + dx, -8.0), 'Center', 'baseline') for dx in (0, 20, -20)])
        box, clear = c.place(outside, name)
        tx = box[0] - 2.0 if box[0] >= plot.w else (box[0] + box[2]) / 2
        ty = (box[1] + box[3]) / 2 if box[0] >= plot.w else box[3] + 1.0
        c.segment((fx, fy), (tx, ty))
        if not clear:
            warnings.append(f'the label {name} overlaps other ink')
    located = []
    for i, item in enumerate(spec.get('points') or []):
        ph, e = number(item.get('pH'), f'point {i + 1} pH', ph0, ph1), number(item.get('E'), f'point {i + 1} E', e0, e1)
        best = min(species, key=lambda s: s['c'] + s['a'] * ph + s['b'] * e)
        c.circle(plot.at(ph, e), 1.4, fill='#000000', outline=None)
        if item.get('label'):
            px, py = plot.at(ph, e)
            c.place([((px + 4.0, py - 4.0), 'Left', 'baseline'), ((px - 4.0, py - 4.0), 'Right', 'baseline')],
                    str(item['label']))
        located.append({'label': item.get('label'), 'pH': ph, 'E': e, 'stable': best['name']})
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report = {'diagram': 'pourbaix', 'element': element, 'T': temperature, 'activity': activity,
              'fields': [f['species']['name'] for f in fields],
              'absent': [s['name'] for s in species if s['name'] not in {f['species']['name'] for f in fields}],
              'points': located, 'warnings': warnings,
              'note': 'fields follow from the given Gibbs energies; check them against a data table'}
    return _result(c, output, report)


# ---------------------------------------------------------------------------- Latimer and Frost
def _redox_chain(spec):
    raw = spec.get('species')
    if not isinstance(raw, list) or len(raw) < 2:
        raise ValueError('Give at least two "species" from the highest oxidation state to the lowest')
    chain = []
    for i, item in enumerate(raw):
        name = str(item.get('name') or f'species {i + 1}')
        ox = number(item.get('oxidation_state'), f'{name} oxidation state', -10, 10)
        chain.append({'name': name, 'ox': ox})
    for a, b in zip(chain, chain[1:]):
        if not b['ox'] < a['ox']:
            raise ValueError('List the species from the highest oxidation state to the lowest')
    potentials = spec.get('potentials')
    if not isinstance(potentials, list) or len(potentials) != len(chain) - 1:
        raise ValueError(f'"potentials" needs {len(chain) - 1} values, one per adjacent couple')
    steps = [number(v, 'potential', -10, 10) for v in potentials]
    # Volt-equivalents N E (relative to the lowest state, then shifted to oxidation state 0 when present).
    ve = [0.0] * len(chain)
    for k in range(len(chain) - 2, -1, -1):
        ve[k] = ve[k + 1] + (chain[k]['ox'] - chain[k + 1]['ox']) * steps[k]
    zero = next((v for v, s in zip(ve, chain) if abs(s['ox']) < 1e-9), None)
    shift = zero if zero is not None else ve[-1]
    for s, v in zip(chain, ve):
        s['ve'] = v - shift
    unstable = []
    for k in range(1, len(chain) - 1):
        if steps[k] > steps[k - 1] + 1e-9:     # reduction of X is easier than reduction to X
            unstable.append(chain[k]['name'])
    couples = []
    for i in range(len(chain)):
        for j in range(i + 1, len(chain)):
            electrons = chain[i]['ox'] - chain[j]['ox']
            couples.append({'couple': f'{chain[i]["name"]}/{chain[j]["name"]}', 'n': electrons,
                            'E': round((chain[i]['ve'] - chain[j]['ve']) / electrons, 4)})
    return chain, steps, unstable, couples, zero is not None


def latimer(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Latimer diagram: ``{"species": [{"name": "MnO4{sup:\u2212}", "oxidation_state": 7}, ...],
    "potentials": [0.56, ...], "medium": "acid"}``, highest oxidation state first. ``derived``
    lists couples (``[["MnO4{sup:\u2212}", "MnO2"]]``) to add as arcs with their computed E\u00b0."""
    chain, steps, unstable, couples, _ = _redox_chain(spec)
    c = Canvas()
    x = 0.0
    centres = []
    for k, s in enumerate(chain):
        box = c.label((x, 0.0), s['name'], valign='middle', color=ORANGE if s['name'] in unstable else None)
        centres.append(((box[0] + box[2]) / 2, box[0], box[2]))
        if k < len(steps):
            start = box[2] + 4.0
            text = f'{steps[k]:+.2f}'.replace('-', MINUS)
            length = max(34.0, len(text) * 5.5 + 10.0)
            c.arrow((start, 0.0), (start + length, 0.0), head_size=800)
            c.label((start + length / 2, -5.0), text, align='Center', size=8)
            x = start + length + 4.0
    by_name = {s['name']: k for k, s in enumerate(chain)}
    level = 0
    derived = []
    for pair in spec.get('derived') or []:
        if len(pair) != 2 or pair[0] not in by_name or pair[1] not in by_name:
            raise ValueError(f'derived couple {pair!r} names unknown species')
        i, j = sorted((by_name[pair[0]], by_name[pair[1]]))
        if j - i < 2:
            raise ValueError(f'{pair[0]}/{pair[1]} is an adjacent couple, already drawn')
        value = next(cp['E'] for cp in couples if cp['couple'] == f'{chain[i]["name"]}/{chain[j]["name"]}')
        level += 1
        y = -8.0 - 14.0 * level
        a, b = centres[i][0], centres[j][0]
        c.path([(a, -7.0), (a, y), (b, y)])
        c.arrow((b, y), (b, -7.0), head_size=800)
        c.label(((a + b) / 2, y - 4.0), f'{value:+.2f}'.replace('-', MINUS), align='Center', size=8)
        derived.append({'couple': f'{chain[i]["name"]}/{chain[j]["name"]}', 'E': value})
    if spec.get('medium'):
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), f'{spec["medium"]}, E\u00b0 / V', align='Center')
    report = {'diagram': 'latimer', 'species': [s['name'] for s in chain], 'disproportionate': unstable,
              'derived': derived, 'all_couples': couples, 'warnings': []}
    return _result(c, output, report)


def frost(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Frost diagram (N E\u00b0 against oxidation state N) from the same spec as :func:`latimer`.
    Several chains go in ``"series": [{"label", "species", "potentials"}]`` (e.g. acid and base)."""
    series = spec.get('series') or [spec]
    chains = []
    for i, item in enumerate(series):
        chain, steps, unstable, couples, anchored = _redox_chain(item)
        chains.append({'label': item.get('label'), 'chain': chain, 'unstable': unstable, 'anchored': anchored,
                       'color': resolve_color(item.get('color'), SERIES_COLORS[i % 3])})
    oxs = [s['ox'] for ch in chains for s in ch['chain']]
    ves = [s['ve'] for ch in chains for s in ch['chain']]
    xpad = 0.6
    ypad = (max(ves) - min(ves)) * 0.12 or 1.0
    c = Canvas()
    plot = Plot(c, (min(oxs) - xpad, max(oxs) + xpad), (min(ves) - ypad, max(ves) + ypad),
                float(spec.get('width', 220.0)), float(spec.get('height', 170.0)),
                x_ticks=list(range(math.ceil(min(oxs) - xpad), math.floor(max(oxs) + xpad) + 1)))
    plot.frame(spec.get('x_label', 'oxidation state N'), spec.get('y_label', 'N E\u00b0 / V'))
    warnings = []
    for ch in chains:
        pts = [(s['ox'], s['ve']) for s in ch['chain']]
        plot.polyline(pts, color=ch['color'], dashed=len(chains) > 1 and chains.index(ch) > 0)
        for s in ch['chain']:
            px, py = plot.at(s['ox'], s['ve'])
            c.circle((px, py), 1.8, fill=ch['color'], outline=None)
            c.block((px - 2, py - 2, px + 2, py + 2))
        for s in ch['chain']:
            px, py = plot.at(s['ox'], s['ve'])
            near = [((px, py - 6.0), 'Center', 'baseline'), ((px, py + 13.0), 'Center', 'baseline'),
                    ((px + 5.0, py), 'Left', 'middle'), ((px - 5.0, py), 'Right', 'middle')]
            far = [((px + dx, py + dy), 'Left' if dx > 0 else 'Right', 'middle')
                   for dy in (-10.0, 10.0, -18.0, 18.0, -26.0, 26.0) for dx in (5.0, -5.0)]
            _b, clear = c.place(near + far,
                                s['name'], color=ch['color'] if len(chains) > 1 else None, size=8)
            if not clear:
                warnings.append(f'the label {s["name"]} overlaps other ink')
    if len(chains) > 1:
        x, y = plot.w + 12.0, 10.0
        for ch in chains:
            if ch['label']:
                c.segment((x, y), (x + 16.0, y), color=ch['color'], dashed=chains.index(ch) > 0)
                c.label((x + 20.0, y), str(ch['label']), valign='middle', color=ch['color'])
                y += 14.0
    report = {'diagram': 'frost', 'series': [{'label': ch['label'], 'points': [{'species': s['name'], 'N': s['ox'],
                                                                                'NE': round(s['ve'], 4)} for s in ch['chain']],
                                              'disproportionate': ch['unstable'],
                                              'reference': 'oxidation state 0' if ch['anchored'] else 'lowest state listed'}
                                             for ch in chains], 'warnings': warnings}
    return _result(c, output, report)


# ---------------------------------------------------------------------------- Ellingham
def ellingham(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Ellingham diagram: ``lines``: ``[{"label": "2Mg + O{sub:2} \u2192 2MgO", "dH": -1203.4, "dS": -216.8}]``
    (kJ/mol and J/(mol K), both per mole of O2 and taken as constant), optional ``segments``
    ``[{"from_T": 1380, "dH": ..., "dS": ...}]`` after a phase change; ``T`` range in K."""
    t0, t1 = (number(v, 'T', 0.0, 5000.0) for v in spec.get('T', [273, 2273]))
    if t1 <= t0:
        raise ValueError('The temperature range must increase')
    raw = spec.get('lines')
    if not isinstance(raw, list) or not raw:
        raise ValueError('An Ellingham diagram needs "lines"')
    lines = []
    for i, item in enumerate(raw):
        label = str(item.get('label', f'line {i + 1}'))
        pieces = [(t0, number(item.get('dH'), f'{label} dH', -1e5, 1e5), number(item.get('dS'), f'{label} dS', -1e5, 1e5))]
        for seg in item.get('segments') or []:
            start = number(seg.get('from_T'), f'{label} segment T', t0, t1)
            if start <= pieces[-1][0]:
                raise ValueError(f'{label}: segments must follow in increasing temperature')
            pieces.append((start, number(seg.get('dH'), 'dH', -1e5, 1e5), number(seg.get('dS'), 'dS', -1e5, 1e5)))

        def dg(t, pieces=pieces):
            piece = [p for p in pieces if p[0] <= t + 1e-9][-1]
            return piece[1] - t * piece[2] / 1000.0
        # Continuity check at each phase change (a kink, never a jump).
        for start, *_ in pieces[1:]:
            before = [p for p in pieces if p[0] < start][-1]
            jump = abs((before[1] - start * before[2] / 1000.0) - dg(start))
            if jump > 1.0:
                raise ValueError(f'{label}: Delta G jumps by {jump:.1f} kJ at {start:g} K; a phase change is a kink, '
                                 'so the segment\'s dH and dS must meet the previous line')
        lines.append({'label': label, 'fn': dg, 'pieces': pieces,
                      'color': resolve_color(item.get('color'), SERIES_COLORS[i % 3])})
    samples = [t0 + (t1 - t0) * k / 200 for k in range(201)]
    values = [ln['fn'](t) for ln in lines for t in samples]
    g0, g1 = min(values + [0.0]), max(values + [0.0])
    pad = (g1 - g0) * 0.06
    c = Canvas()
    plot = Plot(c, (t0, t1), (g0 - pad, g1 + pad), float(spec.get('width', 240.0)), float(spec.get('height', 180.0)))
    plot.frame(spec.get('x_label', 'T / K'), spec.get('y_label', '\u0394G\u00b0 / (kJ mol{sup:\u22121} O{sub:2})'))
    plot.segment((t0, 0.0), (t1, 0.0), dashed=True)
    warnings = []
    for ln in lines:
        plot.polyline([(t, ln['fn'](t)) for t in samples], color=ln['color'])
        for start, *_ in ln['pieces'][1:]:
            c.circle(plot.at(start, ln['fn'](start)), 1.4, fill=ln['color'], outline=None)
    for ln in lines:
        px, py = plot.at(t1, ln['fn'](t1))
        _b, clear = c.place([((px + 4.0, py + dy), 'Left', 'middle') for dy in (0, -8, 8, -16, 16, -24, 24)],
                            ln['label'], color=ln['color'] if len(lines) > 1 else None, size=8)
        if not clear:
            warnings.append(f'the label {ln["label"]} overlaps another')
    crossings = []
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            diff = [a['fn'](t) - b['fn'](t) for t in samples]
            for k in range(len(samples) - 1):
                if diff[k] == 0 or diff[k] * diff[k + 1] < 0:
                    lo, hi = samples[k], samples[k + 1]
                    for _ in range(60):
                        mid = (lo + hi) / 2
                        if (a['fn'](mid) - b['fn'](mid)) * diff[k] > 0:
                            lo = mid
                        else:
                            hi = mid
                    t = (lo + hi) / 2
                    above_low = a['label'] if a['fn'](t + 1) < b['fn'](t + 1) else b['label']
                    crossings.append({'lines': [a['label'], b['label']], 'T': round(t, 1),
                                      'lower_above': above_low,
                                      'meaning': f'above {t:.0f} K the element of "{above_low}" can reduce the other oxide'})
                    c.circle(plot.at(t, a['fn'](t)), 1.6, fill='#000000', outline=None)
    report = {'diagram': 'ellingham', 'crossings': crossings,
              'lines': [{'label': ln['label'], 'dG_at_T0': round(ln['fn'](t0), 1), 'dG_at_T1': round(ln['fn'](t1), 1)}
                        for ln in lines], 'warnings': warnings,
              'note': 'dH and dS are taken as constant between phase changes'}
    return _result(c, output, report)
