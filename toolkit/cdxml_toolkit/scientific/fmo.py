"""Frontier-orbital interaction (FMO) diagrams.

A filled donor level sits on the left, an empty acceptor level on the right, and the
mixing block between them: the donor and acceptor mix into two new levels that straddle
the donor (the filled combination stabilised below it, the empty combination
destabilised above it), joined by dashed correlation lines. Ground-state delocalisation
and vertical (Franck-Condon) excitation arrows can be overlaid. Coefficient lobes scale
with ``|coefficient|``; each interacting donor/acceptor pair is labelled with ``|ci cj|``.
The double-headed dE marker spans the donor HOMO and the lowest empty acceptor orbital.
"""
from __future__ import annotations

from itertools import product
from typing import Any, Mapping

from .canvas import BLUE, Canvas, caption_box, plain
from .common import number
from .energy import _result, _spread
from .orbital_lobes import lobe
from .orbitals import _spin_arrows, hund_fill

BAR = 20.0          # width of one orbital level
LOBES = {'donor': 180.0, 'acceptor': 0.0}   # coefficient lobe points away from the mixing block


def _level(item, column, index):
    """Normalise one ``{"label", "energy", "electrons", ...}`` entry."""
    label = str(item.get('label', ''))
    energy = number(item.get('energy'), f'{column} level {index + 1} energy', -1e9, 1e9)
    electrons = int(item.get('electrons', 0))
    if electrons not in (0, 1, 2):
        raise ValueError(f'{column} level {index + 1}: an orbital holds 0, 1 or 2 electrons, not {electrons}')
    coefficient = item.get('coefficient')
    if coefficient is not None:
        coefficient = number(coefficient, f'{column} level {index + 1} coefficient', 0.0, 1.0)
    phase = int(item.get('phase', 1))
    if phase not in (-1, 1):
        raise ValueError(f'{column} level {index + 1}: phase must be 1 or -1')
    return {
        'id': str(item.get('id') or plain(label) or f'{column}{index + 1}'),
        'label': label,
        'energy': energy,
        'electrons': electrons,
        'from': [str(v) for v in item.get('from', [])],
        'coefficient': coefficient,
        'phase': phase,
        'column': column,
    }


def _columns(spec):
    """Read the three columns; ``mixing`` entries add the dashed ``from`` links."""
    columns = {}
    for key in ('donor', 'acceptor', 'mixing'):
        raw = spec.get(key)
        if not isinstance(raw, list) or not raw:
            raise ValueError(f'"{key}" needs a list of levels')
        columns[key] = [_level(item, key, i) for i, item in enumerate(raw)]
    donor, acceptor, mixing = columns['donor'], columns['acceptor'], columns['mixing']
    if not any(lv['electrons'] for lv in donor):
        raise ValueError('the donor needs at least one occupied level')
    if not any(lv['electrons'] == 0 for lv in acceptor):
        raise ValueError('the acceptor needs at least one empty orbital for its LUMO')
    ids = {lv['id'] for lv in donor + acceptor}
    for lv in mixing:
        if not lv['from']:
            raise ValueError(f'mixing level {lv["id"]}: give "from" donor/acceptor level ids')
        for ref in lv['from']:
            if ref not in ids:
                raise ValueError(f'mixing level {lv["id"]}: "from" names an unknown level {ref!r}')
    return columns


def fmo_interaction(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Frontier-orbital interaction diagram.

    ``{"donor": [{"label", "energy", "electrons", "id", "coefficient", "phase"}],
    "acceptor": [...], "mixing": [{"label", "energy", "electrons", "from": [ids]}]}``.
    The mixing levels are the generated combinations; their ``from`` lists name the donor
    and acceptor levels they combine, and dashed lines join each source to each mixing
    level. Electrons are conserved across the three columns, and the generated levels must
    straddle the donor (filled below, empty above). Optional ``"arrows"`` draws a
    ``"delocalisation"`` arrow (donor HOMO to the stabilised level) and an ``"excitation"``
    vertical Franck-Condon arrow; ``"shade_stabilisation": true`` shades the stabilisation
    between the donor and the filled mixing level. ``"coefficients": true`` (or a per-level
    ``coefficient``) draws phase-coloured lobes sized by ``|coefficient|`` and labels each
    interacting donor/acceptor pair with its magnitude product. The dE marker and ``gap``
    report the acceptor LUMO minus donor HOMO energy; the acceptor must have an empty orbital.
    """
    columns = _columns(spec)
    donor, acceptor, mixing = columns['donor'], columns['acceptor'], columns['mixing']
    every = donor + acceptor + mixing

    # Energy axis: energies to y (lowest energy at the bottom of the page).
    height = float(spec.get('height', max(150.0, 30.0 * (len({round(lv['energy'], 6) for lv in every}) - 1))))
    warnings: list[str] = []
    if _spread(every, 'energy', height, 32.0) and not spec.get('schematic'):
        warnings.append('close levels were moved apart: the vertical spacing is not to scale')

    # Physics checks.
    donor_homo = max((lv for lv in donor if lv['electrons']), key=lambda lv: lv['energy'])
    acceptor_lumo = min((lv for lv in acceptor if lv['electrons'] == 0), key=lambda lv: lv['energy'])
    filled = [lv for lv in mixing if lv['electrons']]
    empty = [lv for lv in mixing if not lv['electrons']]
    straddles = bool(filled and empty and max(lv['energy'] for lv in filled) < donor_homo['energy']
                     and min(lv['energy'] for lv in empty) > donor_homo['energy'])
    if not straddles:
        raise ValueError('the generated levels must straddle the donor: the filled level below it and '
                         'the empty level above it (stabilised down, destabilised up)')
    total_before = sum(lv['electrons'] for lv in donor + acceptor)
    total_after = sum(lv['electrons'] for lv in mixing)
    if total_before != total_after and not spec.get('allow_mismatch'):
        raise ValueError(f'the donor and acceptor hold {total_before} electrons but the generated levels hold '
                         f'{total_after}; electrons are conserved (set "allow_mismatch": true for a charged fragment)')
    highest_filled = max(filled, key=lambda lv: lv['energy'])
    lowest_empty = min(empty, key=lambda lv: lv['energy'])
    gap = acceptor_lumo['energy'] - donor_homo['energy']
    stabilisation = donor_homo['energy'] - highest_filled['energy']

    # Geometry: donor left, mixing centre, acceptor right.
    column_gap = float(spec.get('column_gap', 118.0))
    centers = {'donor': 0.0, 'mixing': column_gap, 'acceptor': 2 * column_gap}
    for lv in every:
        if lv['column'] == 'mixing':
            lv['x1'] = centers[lv['column']] - BAR / 2
            lv['x2'] = centers[lv['column']] + BAR / 2
        else:
            lv['x1'] = centers[lv['column']]
            lv['x2'] = lv['x1'] + BAR

    c = Canvas()
    links: list[dict[str, str]] = []
    source_index = {lv['id']: lv for lv in donor + acceptor}

    # Stabilisation shading first, behind the dashed lines.
    if spec.get('shade_stabilisation', True):
        y_top = min(donor_homo['y'], highest_filled['y'])
        y_bottom = max(donor_homo['y'], highest_filled['y'])
        c.filled_rect(centers['mixing'] - 14.0, y_top + 2.0, centers['mixing'] + 14.0,
                      y_bottom - 2.0, fill_color='#808080')

    # Dashed correlation lines: each source named by a mixing level joins that level.
    for lv in mixing:
        for ref in lv['from']:
            source = source_index[ref]
            if source['column'] == 'donor':
                tail, head = (source['x2'], source['y']), (lv['x1'], lv['y'])
            else:
                tail, head = (source['x1'], source['y']), (lv['x2'], lv['y'])
            node = c.segment(tail, head, dashed=True)
            node.set('Tail3D', f'{tail[0]} {tail[1]} 0')
            node.set('Head3D', f'{head[0]} {head[1]} 0')
            links.append({'from': source['id'], 'to': lv['id'], 'node_id': node.get('id')})

    # Level lines, spin arrows and labels.
    for lv in every:
        occupation = hund_fill(lv['electrons'], 1)[0]
        lv['occupation'] = occupation
        line_node = c.segment((lv['x1'], lv['y']), (lv['x2'], lv['y']))
        lv['line_node_id'] = line_node.get('id')
        cx = (lv['x1'] + lv['x2']) / 2
        if occupation:
            _spin_arrows(c, cx, lv['y'], occupation)
        if lv['column'] in LOBES and (lv['coefficient'] is not None or spec.get('coefficients')):
            coeff = lv['coefficient'] if lv['coefficient'] is not None else 1.0
            lobe(c, (cx, lv['y'] - 14.0), angle_deg=LOBES[lv['column']], size=26.0, width=12.0,
                 coefficient=coeff, phase=lv['phase'], fill='filled')
        if not lv['label']:
            continue
        if lv['column'] == 'donor':
            c.label((lv['x1'] - 5.0, lv['y']), lv['label'], align='Right', valign='middle')
        elif lv['column'] == 'acceptor':
            c.label((lv['x2'] + 5.0, lv['y']), lv['label'], align='Left', valign='middle')
        else:
            box = caption_box(0.0, lv['y'], lv['label'])
            if lv is lowest_empty and any(a.get('kind') == 'excitation' for a in spec.get('arrows') or []):
                baseline = lv['y'] - 6.5 - (box[3] - lv['y'])
            else:
                baseline = lv['y'] + 6.5 + (lv['y'] - box[1])
            c.label((cx, baseline), lv['label'], align='Center')

    gap_x = c.extent[2] + 18.0
    gap_arrow = c.arrow((gap_x, donor_homo['y']), (gap_x, acceptor_lumo['y']),
                        tail_style='Full', head_size=800) if gap else None
    c.label((gap_x + 5.0, (donor_homo['y'] + acceptor_lumo['y']) / 2),
            f'dE = {gap:g}', size=8, valign='middle')

    coefficient_pairs = []
    seen_pairs = set()
    pair_y = max(lv['y'] for lv in every) + 52.0
    for lv in mixing:
        sources = [source_index[ref] for ref in lv['from']]
        for left, right in product((s for s in sources if s['column'] == 'donor'),
                                   (s for s in sources if s['column'] == 'acceptor')):
            pair = (left['id'], right['id'])
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            if not spec.get('coefficients') and left['coefficient'] is None and right['coefficient'] is None:
                continue
            ci = left['coefficient'] if left['coefficient'] is not None else 1.0
            cj = right['coefficient'] if right['coefficient'] is not None else 1.0
            magnitude = abs(ci * cj)
            c.label((centers['mixing'], pair_y + 14.0 * len(coefficient_pairs)),
                    f'|c({left["id"]})c({right["id"]})| = {magnitude:.2f}', size=8, align='Center')
            coefficient_pairs.append({'donor': left['id'], 'acceptor': right['id'],
                                      'magnitude_product': magnitude,
                                      'label_node_id': c.page[-1].get('id')})

    # Optional arrows.
    drawn_arrows: list[dict[str, Any]] = []
    for i, arrow in enumerate(spec.get('arrows') or []):
        kind = str(arrow.get('kind', 'delocalisation')).lower()
        if kind not in ('delocalisation', 'excitation'):
            raise ValueError(f'arrow {i + 1}: kind must be "delocalisation" or "excitation"')
        label = arrow.get('label')
        color = arrow.get('color')
        if kind == 'delocalisation':
            tail = (BAR + 2.0, donor_homo['y'])
            head = (centers['mixing'] - BAR / 2 - 2.0, highest_filled['y'])
            c.arrow(tail, head, color=color or BLUE, head_style='Full')
            mid = ((tail[0] + head[0]) / 2, max(tail[1], head[1]) + 16.0)
            drawn_arrows.append({'kind': kind, 'from': donor_homo['id'], 'to': highest_filled['id']})
        else:
            tail = (centers['mixing'], highest_filled['y'])
            head = (centers['mixing'], lowest_empty['y'])
            c.arrow(tail, head, color=color, dashed=True, head_style='Full')
            mid = (centers['mixing'] + 16.0, (tail[1] + head[1]) / 2)
            drawn_arrows.append({'kind': kind, 'from': highest_filled['id'], 'to': lowest_empty['id']})
        if label:
            c.label(mid, str(label), size=8, align='Left', valign='middle')

    bottom = max(lv['y'] for lv in every) + 32.0
    if spec.get('donor_title'):
        c.label((centers['donor'] + BAR / 2, bottom), str(spec['donor_title']), align='Center')
    if spec.get('acceptor_title'):
        c.label((centers['acceptor'] + BAR / 2, bottom), str(spec['acceptor_title']), align='Center')
    if spec.get('axis', True):
        left = c.extent[0] - 12.0
        top = min(lv['y'] for lv in every) - 16.0
        c.arrow((left, max(lv['y'] for lv in every) + 10.0), (left, top))
        c.label((left - 4.0, top - 6.0), spec.get('y_label', 'Energy'))

    report: dict[str, Any] = {
        'diagram': 'fmo_interaction',
        'donor_homo': {'id': donor_homo['id'], 'label': donor_homo['label'], 'energy': donor_homo['energy']},
        'acceptor_lumo': {'id': acceptor_lumo['id'], 'label': acceptor_lumo['label'], 'energy': acceptor_lumo['energy']},
        'gap': round(gap, 6),
        'gap_arrow_node_id': gap_arrow.get('id') if gap_arrow is not None else None,
        'coefficient_pairs': coefficient_pairs,
        'stabilisation': round(stabilisation, 6),
        'straddles_donor': straddles,
        'electrons': total_after,
        'occupation': [{'id': lv['id'], 'label': lv['label'], 'electrons': lv['electrons'],
                        'occupation': lv['occupation'], 'energy': lv['energy'],
                        'line_node_id': lv['line_node_id']}
                       for lv in every],
        'generated_levels': [{'id': lv['id'], 'label': lv['label'], 'energy': lv['energy'],
                              'electrons': lv['electrons'], 'occupation': lv['occupation'],
                              'line_node_id': lv['line_node_id']}
                             for lv in sorted(mixing, key=lambda item: item['energy'])],
        'links': links,
        'arrows': drawn_arrows,
        'warnings': warnings,
    }
    if spec.get('schematic'):
        report['energies'] = 'schematic order, not to scale'
    return _result(c, output, report)
