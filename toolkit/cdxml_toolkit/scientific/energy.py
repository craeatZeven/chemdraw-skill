"""Energy diagrams: reaction energy profiles, enthalpy cycles and Jablonski diagrams.

Each builder takes a JSON-like spec, writes one editable CDXML and returns a report. The
report carries what the drawing implies (barriers, overall change, step sums) and every
inconsistency found, so a figure cannot silently contradict its own numbers.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

from .canvas import BLACK, BLUE, ORANGE, SERIES_COLORS, Canvas, plain, resolve_color, text_width
from .common import number

DAGGER = '‡'
MINUS = '−'


def _fmt(value, pattern='{:g}'):
    """A number for a caption, with a true minus sign."""
    return pattern.format(value).replace('-', MINUS)


def _energy(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} needs a numeric energy')
    return number(value, name, -1e9, 1e9)


def _result(canvas, output, report, margin=12.0):
    saved, xml = canvas.write(output, margin)
    report['warnings'] = list(dict.fromkeys(report.get('warnings', [])))
    return {'output_path': saved, 'cdxml': xml, **report}


def _save_text(text, output):
    """Publish CDXML text the way :meth:`Drawing.save` does (never over an existing file)."""
    from ..mcp_runtime import artifact_safety
    target = artifact_safety.resolve_destination(source=None, output_path=str(output), tag='scientific', suffix='.cdxml')
    with artifact_safety.staging_file(target) as stage:
        Path(stage).write_text(text, encoding='utf-8')
        artifact_safety.publish_file(stage, target)
    return str(target)


# ---------------------------------------------------------------------------- energy profile
def _kinds(points, name):
    """Classify each stationary point; reject a transition state that is not a maximum."""
    warnings = []
    for i, pt in enumerate(points):
        left = points[i - 1]['energy'] if i else None
        right = points[i + 1]['energy'] if i + 1 < len(points) else None
        peak = left is not None and right is not None and pt['energy'] > left and pt['energy'] > right
        kind = pt.get('kind')
        if kind is None:
            pt['kind'] = 'ts' if peak else 'minimum'
        elif kind == 'ts':
            if not peak:
                raise ValueError(
                    f'{name}: {pt["name"]} is marked as a transition state but is not higher than both '
                    'neighbours; a transition state is an energy maximum along the reaction coordinate')
        elif kind in ('minimum', 'intermediate', 'reactant', 'product'):
            pt['kind'] = 'minimum'
            if peak:
                warnings.append(f'{name}: {pt["name"]} is drawn as a minimum but lies above both neighbours '
                                '(set "kind": "ts" if it is a transition state)')
        else:
            raise ValueError(f'{name}: unknown kind {kind!r} for {pt["name"]}')
    return warnings


def energy_profile(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Reaction energy profile (potential energy against reaction coordinate).

    ``style`` is ``"smooth"`` (one continuous curve through the stationary points) or
    ``"levels"`` (a bar per species joined by dashed lines). ``pathways`` is a list of
    ``{"label", "color", "dashed", "points": [{"label", "energy", "id", "kind", "x"}]}``;
    a point with no ``kind`` is a transition state when it lies above both neighbours.
    ``arrows`` mark energy differences: ``{"from", "to", "label", "pathway", "heads"}``.
    """
    style = str(spec.get('style', 'smooth')).lower()
    if style not in ('smooth', 'levels'):
        raise ValueError('style must be "smooth" or "levels"')
    raw_paths = spec.get('pathways')
    if raw_paths is None and spec.get('points') is not None:
        raw_paths = [{'points': spec['points'], 'label': spec.get('label')}]
    if not isinstance(raw_paths, list) or not raw_paths:
        raise ValueError('An energy profile needs "pathways" (or a single "points" list)')
    if len(raw_paths) > 6:
        raise ValueError('At most six pathways fit one energy profile')

    warnings: list[str] = []
    paths = []
    for p_idx, raw in enumerate(raw_paths):
        name = str(raw.get('label') or f'pathway {p_idx + 1}')
        raw_points = raw.get('points')
        if not isinstance(raw_points, list) or len(raw_points) < 2:
            raise ValueError(f'{name}: at least two points are needed')
        points = []
        for i, item in enumerate(raw_points):
            label = str(item.get('label', ''))
            pid = str(item.get('id') or plain(label) or f'p{i + 1}')
            points.append({'id': pid, 'name': pid, 'label': label, 'kind': item.get('kind'),
                           'energy': _energy(item.get('energy'), f'{name} point {pid}'),
                           'x': item.get('x'), 'side': item.get('label_side'), 'show_value': item.get('show_value'),
                           'structure': item.get('structure')})
        ids = [pt['id'] for pt in points]
        if len(set(ids)) != len(ids):
            raise ValueError(f'{name}: point ids must be unique (give "id" where labels repeat)')
        warnings += _kinds(points, name)
        paths.append({'name': name, 'label': raw.get('label'), 'points': points, 'dashed': bool(raw.get('dashed')),
                      'color': resolve_color(raw.get('color'), SERIES_COLORS[p_idx % len(SERIES_COLORS)])})

    # ---- structures drawn at stationary points (separate CDXML files, merged in at the end)
    structures = []
    for p in paths:
        for pt in p['points']:
            if pt['structure']:
                from ..cdxml_compose import content_box, parse_cdxml
                root = parse_cdxml(str(pt['structure']))
                box = content_box([child for page in root.findall('page') for child in page])
                pt['structure_size'] = (box[2] - box[0] + 6.0, box[3] - box[1] + 6.0)
                structures.append(pt)

    # ---- geometry
    count = max(len(p['points']) for p in paths)
    step = float(spec.get('step', 62.0 if style == 'smooth' else 70.0))
    if structures and spec.get('step') is None:
        step = max(step, max(pt['structure_size'][0] for pt in structures) + 10.0)
    lead = 16.0
    half = min(16.0, step * 0.26)
    span = float(spec.get('width', max(step * (count - 1), 170.0)))
    plot_h = float(spec.get('height', 170.0))
    if span < 60 or plot_h < 60:
        raise ValueError('width and height must be at least 60')
    energies = [pt['energy'] for p in paths for pt in p['points']]
    e_min, e_max = min(energies), max(energies)
    if e_max - e_min < 1e-12:
        raise ValueError('All energies are equal; there is no profile to draw')
    inner_top, inner_bottom = 30.0, plot_h - 30.0

    def y_of(energy):
        return inner_bottom - (energy - e_min) / (e_max - e_min) * (inner_bottom - inner_top)

    # Room beside the first and last species so that a centred caption clears the axis.
    reach = lead if style == 'smooth' else half
    pad_left = max([reach + 6.0] + [text_width(p['points'][0]['label']) / 2 + 6.0 for p in paths])
    pad_right = max([reach + 6.0] + [text_width(p['points'][-1]['label']) / 2 + 6.0 for p in paths])
    axis_x = 20.0
    plot_w = pad_left + span + pad_right
    for p in paths:
        n = len(p['points'])
        last = -math.inf
        for i, pt in enumerate(p['points']):
            if pt['x'] is None:
                frac = i / (n - 1)
            else:
                frac = number(pt['x'], f'{p["name"]} x', 0.0, 1.0)
            if frac <= last and i:
                raise ValueError(f'{p["name"]}: x positions must increase along the reaction coordinate')
            last = frac
            pt['px'] = axis_x + pad_left + frac * span
            pt['py'] = y_of(pt['energy'])

    c = Canvas()
    axis_y = plot_h
    show_axis = spec.get('axis', True)
    unit = spec.get('unit')
    if show_axis:
        c.arrow((axis_x, axis_y), (axis_x, 2.0))
        c.arrow((axis_x, axis_y), (axis_x + plot_w + 14.0, axis_y))
        c.block((axis_x, 2.0, axis_x, axis_y), 2.0)
        c.block((axis_x, axis_y, axis_x + plot_w + 14.0, axis_y), 2.0)
        y_label = spec.get('y_label', (f'Energy / ({unit})' if '/' in str(unit) else f'Energy / {unit}')
                           if unit else 'Energy')
        x_label = spec.get('x_label', 'Reaction coordinate')
        if y_label:
            box = c.label((axis_x - 6.0, -4.0), y_label, align='Left')
            c.block(box)
        if x_label:
            box = c.label((axis_x + plot_w + 14.0, axis_y + 13.0), x_label, align='Right')
            c.block(box)
        if spec.get('y_ticks'):
            span = e_max - e_min
            raw_step = span / 4
            mag = 10 ** math.floor(math.log10(raw_step))
            tick = min((m * mag for m in (1, 2, 2.5, 5, 10)), key=lambda v: abs(v - raw_step))
            value = math.ceil(e_min / tick) * tick
            while value <= e_max + 1e-9:
                ty = y_of(value)
                c.segment((axis_x - 3.0, ty), (axis_x, ty))
                c.block(c.label((axis_x - 5.0, ty), _fmt(value), align='Right', valign='middle'))
                value += tick

    # ---- curves first, so captions can avoid every pathway
    for p in paths:
        pts = p['points']
        if style == 'levels':
            for i, pt in enumerate(pts):
                seg = [(pt['px'] - half, pt['py']), (pt['px'] + half, pt['py'])]
                c.segment(*seg, color=p['color'])
                c.block_path(seg, 2.0)
                if i + 1 < len(pts):
                    nxt = pts[i + 1]
                    link = [(pt['px'] + half, pt['py']), (nxt['px'] - half, nxt['py'])]
                    c.segment(*link, color=p['color'], dashed=True)
                    c.block_path(_subdivide(link, 8), 1.0)
        else:
            samples = [(pts[0]['px'] - lead, pts[0]['py'])]
            for a, b in zip(pts, pts[1:]):
                for k in range(25):
                    t = k / 24
                    samples.append((a['px'] + (b['px'] - a['px']) * t,
                                    a['py'] + (b['py'] - a['py']) * (1 - math.cos(math.pi * t)) / 2))
            samples.append((pts[-1]['px'] + lead, pts[-1]['py']))
            c.path(samples, color=p['color'], dashed=p['dashed'])
            c.block_path(samples, 1.5)

    # ---- energy-difference arrows
    by_path = [{pt['id']: pt for pt in p['points']} for p in paths]
    arrow_report = []
    for a_idx, arrow in enumerate(spec.get('arrows') or []):
        p_idx = int(arrow.get('pathway', 0))
        if not 0 <= p_idx < len(paths):
            raise ValueError(f'arrow {a_idx + 1}: pathway {p_idx} does not exist')
        lookup = by_path[p_idx]
        try:
            src, dst = lookup[str(arrow['from'])], lookup[str(arrow['to'])]
        except KeyError as exc:
            raise ValueError(f'arrow {a_idx + 1}: no point {exc} in {paths[p_idx]["name"]}') from None
        delta = dst['energy'] - src['energy']
        if abs(delta) < 1e-12:
            raise ValueError(f'arrow {a_idx + 1}: {src["id"]} and {dst["id"]} have the same energy')
        pts = paths[p_idx]['points']
        color = resolve_color(arrow.get('color'), paths[p_idx]['color'])
        if dst is pts[-1]:
            ax = dst['px'] + reach + 8.0
        elif dst is pts[0]:
            ax = dst['px'] - reach - 8.0
        else:
            ax = dst['px']
            # Step aside when another species sits on the arrow's line.
            lo, hi = sorted((src['py'], dst['py']))
            others = [q for path in paths for q in path['points'] if q is not src and q is not dst]
            for shift in (0.0, 18.0, -18.0, 32.0, -32.0):
                if not any(abs(q['px'] - (ax + shift)) < 12.0 and lo - 16.0 < q['py'] < hi + 16.0 for q in others):
                    ax += shift
                    break
        ax += float(arrow.get('dx', 0.0))
        # Dashed guides carry both energy levels to the arrow.
        for pt in (src, dst):
            side = 1.0 if ax >= pt['px'] else -1.0
            edge = pt['px'] + side * (half if style == 'levels' else (lead if pt in (pts[0], pts[-1]) else 0.0))
            if (ax - edge) * side > 3.0:
                guide = [(edge, pt['py']), (ax + side * 4.0, pt['py'])]
                c.segment(*guide, color=color, dashed=True)
                c.block_path(guide, 1.5)
        heads = str(arrow.get('heads', 'both')).lower()
        if heads not in ('both', 'to', 'none'):
            raise ValueError(f'arrow {a_idx + 1}: heads must be "both", "to" or "none"')
        gap = 1.2
        y1 = src['py'] + (gap if dst['py'] > src['py'] else -gap)
        y2 = dst['py'] - (gap if dst['py'] > src['py'] else -gap)
        c.arrow((ax, y1), (ax, y2), color=color, head_style=None if heads == 'none' else 'Full',
                tail_style='Full' if heads == 'both' else None, head_size=800)
        c.block((ax, y1, ax, y2), 3.0)
        label = arrow.get('label')
        if label is None and arrow.get('show_value'):
            label = _fmt(delta, str(spec.get('value_format', '{:+g}')))
        if label:
            candidates = []
            upper, lower = sorted((y1, y2))
            for frac in (0.5, 0.7, 0.85, 0.3, 0.15):  # a hump is widest near its base
                at = upper + (lower - upper) * frac
                candidates += [((ax + 5.0, at), 'Left', 'middle'), ((ax - 5.0, at), 'Right', 'middle')]
            _, clear = c.place(candidates, label, color=color)
            if not clear:
                warnings.append(f'the label of arrow {a_idx + 1} ({plain(label)}) overlaps other ink')
        arrow_report.append({'from': src['id'], 'to': dst['id'], 'pathway': paths[p_idx]['name'], 'difference': delta})

    # ---- species captions
    show_values = bool(spec.get('show_values', False))
    pattern = str(spec.get('value_format', '{:g}'))
    ts_marker = spec.get('ts_marker', True)
    drawn = set()
    for p in paths:
        for pt in p['points']:
            text = pt['label']
            if text and pt['kind'] == 'ts' and ts_marker and DAGGER not in text:
                text += '{sup:%s}' % DAGGER
            show = show_values if pt['show_value'] is None else bool(pt['show_value'])
            if show:
                value = _fmt(pt['energy'], pattern)
                text = f'{text}\n({value})' if text else value
            if not text:
                continue
            key = (round(pt['px'], 1), round(pt['py'], 1), text)
            if key in drawn:
                continue
            drawn.add(key)
            lines = text.count('\n') + 1
            above = ((pt['px'], pt['py'] - 6.0 - (lines - 1) * 12.0), 'Center', 'baseline')
            below = ((pt['px'], pt['py'] + 13.0), 'Center', 'baseline')
            shift = half + 4.0 if style == 'levels' else 8.0
            right = ((pt['px'] + shift, pt['py']), 'Left', 'middle')
            left = ((pt['px'] - shift, pt['py']), 'Right', 'middle')
            order = [above, below, right, left] if pt['kind'] == 'ts' else [below, above, right, left]
            named = {'above': above, 'below': below, 'right': right, 'left': left}
            if pt['side'] is not None:
                if pt['side'] not in named:
                    raise ValueError(f'label_side must be one of {", ".join(named)}')
                order = [named[pt['side']]] + [cand for cand in order if cand is not named[pt['side']]]
            # Then the same sides further out and shifted sideways, for crowded crossings.
            extra = []
            for dist in (9.0, 18.0):
                for (x, y), align, valign in order:
                    vertical = align == 'Center'
                    sign = -1.0 if y < pt['py'] else 1.0
                    if vertical:
                        for dx in (0.0, 14.0, -14.0):
                            extra.append(((x + dx, y + sign * dist), align, valign))
                    else:
                        side = 1.0 if x > pt['px'] else -1.0
                        for dy in (-dist, dist):
                            extra.append(((x + side * 4.0, y + dy), align, valign))
            box, clear = c.place(order + extra, text, color=p['color'] if len(paths) > 1 else None)
            if not clear:
                warnings.append(f'the label of {pt["id"]} in {p["name"]} overlaps other ink; set "label_side" '
                                'or move the point with "x"')
            pt['label_box'] = box
    for pt in structures:
        # The structure continues away from the curve, beyond the species caption.
        width, height = pt['structure_size']
        box = pt.get('label_box') or (pt['px'], pt['py'] - 4.0, pt['px'], pt['py'] + 4.0)
        above = (box[1] + box[3]) / 2 < pt['py']
        cx = (box[0] + box[2]) / 2
        if above:
            target = (cx - width / 2, box[1] - 3.0 - height, cx + width / 2, box[1] - 3.0)
            pt['structure_at'], pt['structure_anchor'] = (cx, box[1] - 6.0), 'bottom'
        else:
            target = (cx - width / 2, box[3] + 3.0, cx + width / 2, box[3] + 3.0 + height)
            pt['structure_at'], pt['structure_anchor'] = (cx, box[3] + 6.0), 'top'
        if any(a[0] < target[2] and target[0] < a[2] and a[1] < target[3] and target[1] < a[3] for a in c.obstacles):
            warnings.append(f'the structure of {pt["id"]} overlaps other ink; widen "step" or change "label_side"')
        c.block(target)
        c.see_box(target)

    # ---- legend for several pathways
    legend = [p for p in paths if p['label']]
    if len(paths) > 1 and legend and spec.get('legend', True):
        lx, ly = axis_x + plot_w + 26.0, 20.0
        for p in legend:
            c.segment((lx, ly), (lx + 18.0, ly), color=p['color'], dashed=p['dashed'])
            c.label((lx + 23.0, ly), str(p['label']), valign='middle', color=p['color'])
            ly += 15.0

    # ---- what the drawing says
    report_paths = []
    for p in paths:
        pts = p['points']
        barriers = []
        for i, pt in enumerate(pts):
            if pt['kind'] != 'ts':
                continue
            before = [q for q in pts[:i] if q['kind'] == 'minimum']
            after = [q for q in pts[i + 1:] if q['kind'] == 'minimum']
            entry = {'ts': pt['id']}
            if before:
                entry.update({'from': before[-1]['id'], 'forward_barrier': pt['energy'] - before[-1]['energy']})
            if after:
                entry.update({'to': after[0]['id'], 'reverse_barrier': pt['energy'] - after[0]['energy']})
            barriers.append(entry)
        forward = [b for b in barriers if 'forward_barrier' in b]
        report_paths.append({
            'label': p['name'],
            'points': [{'id': pt['id'], 'energy': pt['energy'], 'kind': pt['kind']} for pt in pts],
            'barriers': barriers,
            'overall_change': pts[-1]['energy'] - pts[0]['energy'],
            'thermodynamics': 'exothermic' if pts[-1]['energy'] < pts[0]['energy'] else
                              'endothermic' if pts[-1]['energy'] > pts[0]['energy'] else 'thermoneutral',
            'largest_forward_barrier': max(forward, key=lambda b: b['forward_barrier'])['ts'] if forward else None,
        })
    report = {'diagram': 'energy_profile', 'style': style, 'unit': unit, 'pathways': report_paths,
              'arrows': arrow_report, 'warnings': warnings}
    if not structures:
        return _result(c, output, report)
    from ..cdxml_compose import content_box, merge_cdxml, parse_cdxml
    margin = 12.0
    x1, y1, _x2, _y2 = c.extent
    c.finish(margin)
    shift = (margin - x1, margin - y1)
    base = c.to_string()
    origin = content_box(list(parse_cdxml(base).find('page')))
    parts = [{'cdxml': base, 'at': [origin[0], origin[1]], 'anchor': 'top-left'}]
    for pt in structures:
        parts.append({'file': str(pt['structure']), 'anchor': pt['structure_anchor'],
                      'at': [pt['structure_at'][0] + shift[0], pt['structure_at'][1] + shift[1]]})
    merged = merge_cdxml(parts, margin=margin)
    report['structures'] = [pt['id'] for pt in structures]
    report['warnings'] = list(dict.fromkeys(warnings + [w for w in merged['warnings'] if 'overlap' not in w]))
    saved = _save_text(merged['cdxml'], output) if output is not None else None
    return {'output_path': saved, 'cdxml': merged['cdxml'], **report}


def _subdivide(segment, parts):
    (x1, y1), (x2, y2) = segment
    return [(x1 + (x2 - x1) * i / parts, y1 + (y2 - y1) * i / parts) for i in range(parts + 1)]


# ---------------------------------------------------------------------------- enthalpy cycle
def _spread(levels, key, height, min_gap):
    """Map energies to y so that the order is kept and neighbours stay ``min_gap`` apart."""
    ordered = sorted(levels, key=lambda lv: lv[key])
    low, high = ordered[0][key], ordered[-1][key]
    span = (high - low) or 1.0
    adjusted = False
    previous = None
    for lv in ordered:
        y = height - (lv[key] - low) / span * height
        if previous is not None and abs(lv[key] - previous[key]) > 1e-12 and previous['y'] - y < min_gap:
            y = previous['y'] - min_gap
            adjusted = True
        elif previous is not None and abs(lv[key] - previous[key]) <= 1e-12:
            y = previous['y']
        lv['y'] = y
        previous = lv
    return adjusted


def energy_cycle(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Enthalpy-level diagram with stepwise arrows: Born-Haber and Hess cycles.

    ``levels`` is a list of ``{"id", "label", "energy"}``; ``steps`` is a list of
    ``{"from", "to", "label", "value"}`` drawn left to right in the given order. Level
    energies may be omitted when the step values determine them; values that disagree
    with the level energies, or that do not close around a cycle, are an error.
    """
    raw_levels = spec.get('levels')
    raw_steps = spec.get('steps')
    if not isinstance(raw_levels, list) or len(raw_levels) < 2:
        raise ValueError('An energy cycle needs at least two "levels"')
    if not isinstance(raw_steps, list) or not raw_steps:
        raise ValueError('An energy cycle needs "steps"')
    levels = {}
    for i, item in enumerate(raw_levels):
        lid = str(item.get('id') or f'L{i + 1}')
        if lid in levels:
            raise ValueError(f'Duplicate level id {lid!r}')
        energy = item.get('energy')
        levels[lid] = {'id': lid, 'label': str(item.get('label', lid)),
                       'energy': None if energy is None else _energy(energy, f'level {lid}')}
    steps = []
    for i, item in enumerate(raw_steps):
        src, dst = str(item.get('from')), str(item.get('to'))
        if src not in levels or dst not in levels:
            raise ValueError(f'step {i + 1}: unknown level in {src!r} -> {dst!r}')
        if src == dst:
            raise ValueError(f'step {i + 1}: a step joins two different levels')
        value = item.get('value')
        steps.append({'from': src, 'to': dst, 'label': item.get('label'),
                      'value': None if value is None else _energy(value, f'step {i + 1}'),
                      'color': item.get('color'), 'dashed': bool(item.get('dashed'))})

    tolerance = float(spec.get('tolerance', 0.5))
    # Fill in level energies from step values, then verify every step (Hess's law).
    if all(lv['energy'] is None for lv in levels.values()):
        next(iter(levels.values()))['energy'] = 0.0
    changed = True
    while changed:
        changed = False
        for st in steps:
            a, b = levels[st['from']], levels[st['to']]
            if st['value'] is None:
                continue
            if a['energy'] is not None and b['energy'] is None:
                b['energy'] = a['energy'] + st['value']; changed = True
            elif b['energy'] is not None and a['energy'] is None:
                a['energy'] = b['energy'] - st['value']; changed = True
    missing = [lv['id'] for lv in levels.values() if lv['energy'] is None]
    if missing:
        raise ValueError(f'The energy of {", ".join(missing)} is neither given nor fixed by the step values')
    step_report = []
    for i, st in enumerate(steps):
        delta = levels[st['to']]['energy'] - levels[st['from']]['energy']
        if abs(delta) < 1e-12:
            raise ValueError(f'step {i + 1}: {st["from"]} and {st["to"]} have the same energy')
        if st['value'] is not None and abs(delta - st['value']) > tolerance:
            raise ValueError(
                f'step {i + 1} ({st["from"]} -> {st["to"]}) is given as {st["value"]:g} but the level energies '
                f'differ by {delta:g}; the cycle does not close (Hess\'s law)')
        st['delta'] = delta
        step_report.append({'from': st['from'], 'to': st['to'], 'difference': delta, 'given': st['value']})

    ordered = list(levels.values())
    height = float(spec.get('height', max(150.0, 30.0 * (len({lv['energy'] for lv in ordered}) - 1))))
    warnings = []
    if _spread(ordered, 'energy', height, 26.0):
        warnings.append('close levels were moved apart: the vertical spacing is not to scale')

    unit = spec.get('unit')
    show_values = bool(spec.get('show_values', True))
    pattern = str(spec.get('value_format', '{:+g}'))
    captions = []
    for st in steps:
        text = str(st['label']) if st['label'] else ''
        if show_values and (st['value'] is not None or not text):
            value = _fmt(st['delta'], pattern)
            text = f'{text} = {value}' if text else value
        captions.append(text)
    column = float(spec.get('column', max(46.0, max(text_width(t) for t in captions) + 16.0)))
    half = 15.0
    for i, st in enumerate(steps):
        st['x'] = 30.0 + i * column
    for lv in ordered:
        touching = [st['x'] for st in steps if lv['id'] in (st['from'], st['to'])]
        if not touching:
            raise ValueError(f'level {lv["id"]} is not used by any step')
        lv['x1'], lv['x2'] = min(touching) - half, max(touching) + half

    c = Canvas()
    for lv in ordered:
        c.segment((lv['x1'], lv['y']), (lv['x2'], lv['y']))
        c.block((lv['x1'], lv['y'], lv['x2'], lv['y']), 1.5)
    for st in steps:
        a, b = levels[st['from']], levels[st['to']]
        sign = 1.0 if b['y'] > a['y'] else -1.0
        color = resolve_color(st['color']) if st['color'] else None
        c.arrow((st['x'], a['y'] + sign * 1.0), (st['x'], b['y'] - sign * 1.2), color=color, dashed=st['dashed'],
                head_size=800)
        c.block((st['x'], a['y'], st['x'], b['y']), 3.0)
    for st, text in zip(steps, captions):
        if not text:
            continue
        mid = (levels[st['from']]['y'] + levels[st['to']]['y']) / 2
        _, clear = c.place([((st['x'] + 5.0, mid), 'Left', 'middle'), ((st['x'] - 5.0, mid), 'Right', 'middle'),
                            ((st['x'] + 5.0, mid - 12.0), 'Left', 'middle'), ((st['x'] + 5.0, mid + 12.0), 'Left', 'middle')],
                           text, color=resolve_color(st['color']) if st['color'] else None)
        if not clear:
            warnings.append(f'the label of step {st["from"]} -> {st["to"]} overlaps other ink; widen "column"')
    for lv in sorted(ordered, key=lambda item: -item['y']):
        candidates = [((lv['x2'] + 5.0, lv['y']), 'Left', 'middle'), ((lv['x1'] - 5.0, lv['y']), 'Right', 'middle'),
                      ((lv['x1'] + 1.0, lv['y'] - 4.0), 'Left', 'baseline'), ((lv['x1'] + 1.0, lv['y'] + 11.0), 'Left', 'baseline')]
        _, clear = c.place(candidates, lv['label'])
        if not clear:
            warnings.append(f'the label of level {lv["id"]} overlaps other ink; widen "column"')
    if spec.get('axis', True):
        ax = min(lv['x1'] for lv in ordered) - 14.0
        left = min([ax] + [box[0] for box in c.text_boxes]) - 10.0
        bottom, top = max(lv['y'] for lv in ordered) + 8.0, min(lv['y'] for lv in ordered) - 18.0
        c.arrow((left, bottom), (left, top))
        c.label((left - 4.0, top - 6.0), spec.get('y_label', f'Enthalpy / {unit}' if unit else 'Enthalpy'))
    report = {'diagram': 'energy_cycle', 'unit': unit, 'steps': step_report,
              'levels': [{'id': lv['id'], 'energy': lv['energy']} for lv in levels.values()], 'warnings': warnings}
    return _result(c, output, report)


# ---------------------------------------------------------------------------- Jablonski diagram
_TRANSITIONS = {
    'absorption': {'radiative': True, 'direction': 'up', 'spin': 'same', 'color': BLUE, 'label': 'A'},
    'fluorescence': {'radiative': True, 'direction': 'down', 'spin': 'same', 'color': ORANGE, 'label': 'F'},
    'phosphorescence': {'radiative': True, 'direction': 'down', 'spin': 'different', 'color': ORANGE, 'label': 'P'},
    'ic': {'radiative': False, 'direction': 'down', 'spin': 'same', 'color': BLACK, 'label': 'IC'},
    'isc': {'radiative': False, 'direction': 'any', 'spin': 'different', 'color': BLACK, 'label': 'ISC'},
    'vr': {'radiative': False, 'direction': 'down', 'spin': 'same', 'color': BLACK, 'label': 'VR'},
}
_JABLONSKI_PRESET = {
    'states': [
        {'id': 'S0', 'label': 'S{sub:0}', 'energy': 0.0, 'multiplicity': 'singlet'},
        {'id': 'S1', 'label': 'S{sub:1}', 'energy': 3.0, 'multiplicity': 'singlet'},
        {'id': 'S2', 'label': 'S{sub:2}', 'energy': 4.6, 'multiplicity': 'singlet'},
        {'id': 'T1', 'label': 'T{sub:1}', 'energy': 2.1, 'multiplicity': 'triplet'},
    ],
    'transitions': [
        {'from': 'S0', 'to': 'S1', 'kind': 'absorption'},
        {'from': 'S0', 'to': 'S2', 'kind': 'absorption'},
        {'from': 'S2', 'to': 'S1', 'kind': 'ic'},
        {'from': 'S1', 'to': 'S0', 'kind': 'fluorescence'},
        {'from': 'S1', 'to': 'T1', 'kind': 'isc'},
        {'from': 'T1', 'to': 'S0', 'kind': 'phosphorescence'},
    ],
}


def jablonski(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Jablonski diagram: electronic states with vibrational levels and their transitions.

    ``states``: ``{"id", "label", "energy", "multiplicity": "singlet"|"triplet", "vib"}``.
    ``transitions``: ``{"from", "to", "kind", "label"}`` with ``kind`` one of absorption,
    fluorescence, phosphorescence, ic, isc, vr. Radiative transitions are straight arrows,
    non-radiative ones wavy. ``{"preset": "standard"}`` draws S0, S1, S2 and T1.
    """
    if spec.get('preset') is not None:
        if str(spec['preset']).lower() != 'standard':
            raise ValueError('The only Jablonski preset is "standard"')
        spec = {**_JABLONSKI_PRESET, **{k: v for k, v in spec.items() if k != 'preset'}}
    raw_states = spec.get('states')
    if not isinstance(raw_states, list) or len(raw_states) < 2:
        raise ValueError('A Jablonski diagram needs at least two "states"')
    default_vib = int(spec.get('vib', 3))
    states = {}
    for item in raw_states:
        sid = str(item.get('id'))
        if sid in states or sid == 'None':
            raise ValueError(f'Every state needs a unique "id" ({sid!r})')
        mult = str(item.get('multiplicity', 'triplet' if sid.upper().startswith('T') else 'singlet')).lower()
        if mult not in ('singlet', 'triplet'):
            raise ValueError(f'state {sid}: multiplicity must be "singlet" or "triplet"')
        vib = int(item.get('vib', default_vib))
        if not 0 <= vib <= 6:
            raise ValueError(f'state {sid}: "vib" must be 0 to 6')
        states[sid] = {'id': sid, 'label': str(item.get('label', sid)), 'mult': mult, 'vib': vib,
                       'energy': _energy(item.get('energy'), f'state {sid}')}
    ground = min(states.values(), key=lambda s: s['energy'])
    warnings = []
    transitions = []
    for i, item in enumerate(spec.get('transitions') or []):
        kind = str(item.get('kind', '')).lower()
        if kind not in _TRANSITIONS:
            raise ValueError(f'transition {i + 1}: kind must be one of {", ".join(_TRANSITIONS)}')
        rule = _TRANSITIONS[kind]
        src_id, dst_id = str(item.get('from')), str(item.get('to', item.get('from')))
        if src_id not in states or dst_id not in states:
            raise ValueError(f'transition {i + 1}: unknown state in {src_id!r} -> {dst_id!r}')
        src, dst = states[src_id], states[dst_id]
        name = f'{kind} {src_id} -> {dst_id}'
        if kind == 'vr':
            if src is not dst:
                raise ValueError(f'{name}: vibrational relaxation stays within one electronic state')
            if src['vib'] == 0:
                raise ValueError(f'{name}: the state has no vibrational levels to relax through')
        else:
            if src is dst:
                raise ValueError(f'{name}: the two states must differ')
            if rule['direction'] == 'up' and dst['energy'] <= src['energy']:
                raise ValueError(f'{name}: absorption must end in a higher state')
            if rule['direction'] == 'down' and dst['energy'] >= src['energy']:
                raise ValueError(f'{name}: {kind} must end in a lower state')
            same = src['mult'] == dst['mult']
            if rule['spin'] == 'same' and not same:
                raise ValueError(f'{name}: {kind} keeps the spin multiplicity; between a singlet and a triplet '
                                 'use "isc" (non-radiative) or "phosphorescence" (radiative)')
            if rule['spin'] == 'different' and same:
                raise ValueError(f'{name}: {kind} changes the spin multiplicity; both states are {src["mult"]}s')
            if kind == 'isc' and dst['energy'] > src['energy']:
                warnings.append(f'{name} goes uphill (reverse intersystem crossing, thermally activated)')
        transitions.append({'kind': kind, 'src': src, 'dst': dst, 'label': item.get('label', rule['label']),
                            'color': resolve_color(item.get('color'), rule['color'])})

    height = float(spec.get('height', 190.0))
    ordered = list(states.values())
    if _spread(ordered, 'energy', height, 36.0):
        warnings.append('close states were moved apart: the vertical spacing is not to scale')
    vib_gap = 5.0
    singlet_n = sum(1 for t in transitions if t['kind'] != 'isc' and t['src']['mult'] == 'singlet'
                    and t['dst']['mult'] == 'singlet')
    triplet_n = sum(1 for t in transitions if t['kind'] != 'isc' and 'triplet' in (t['src']['mult'], t['dst']['mult']))
    slot = 20.0
    s_width = max(70.0, 18.0 + slot * max(singlet_n, 1))
    has_triplet = any(s['mult'] == 'triplet' for s in ordered)
    t_left = s_width + 34.0
    t_width = max(46.0, 18.0 + slot * max(triplet_n, 1))
    c = Canvas()
    for s in ordered:
        if s['mult'] == 'singlet':
            s['x1'], s['x2'] = 0.0, s_width
        else:
            s['x1'], s['x2'] = t_left, t_left + t_width
        if s is ground and has_triplet:
            s['x2'] = t_left + t_width
        c.segment((s['x1'], s['y']), (s['x2'], s['y']))
        c.block((s['x1'], s['y'], s['x2'], s['y']), 1.5)
        for v in range(1, s['vib'] + 1):
            c.segment((s['x1'], s['y'] - v * vib_gap), (s['x2'], s['y'] - v * vib_gap))
        if s['mult'] == 'singlet':
            c.block(c.label((s['x1'] - 6.0, s['y']), s['label'], align='Right', valign='middle'))
        else:
            c.block(c.label((s['x2'] + 6.0, s['y']), s['label'], align='Left', valign='middle'))

    def top(state):
        return state['y'] - state['vib'] * vib_gap

    next_slot = {'singlet': 14.0, 'triplet': t_left + 14.0}
    drawn = []
    for t in transitions:
        src, dst, kind = t['src'], t['dst'], t['kind']
        if kind == 'isc':
            left, right = (src, dst) if src['x1'] <= dst['x1'] else (dst, src)
            a = (src['x2'] if src is left else src['x1'], src['y'])
            b = (dst['x1'] if dst is right else dst['x2'], top(dst) if dst['energy'] < src['energy'] else dst['y'])
            c.wavy(a, b, color=t['color'])
            anchor = ((a[0] + b[0]) / 2, min(a[1], b[1]) - 5.0)
            candidates = [(anchor, 'Center', 'baseline'), ((anchor[0], max(a[1], b[1]) + 12.0), 'Center', 'baseline')]
        else:
            column = 'triplet' if 'triplet' in (src['mult'], dst['mult']) else 'singlet'
            x = next_slot[column]
            next_slot[column] += slot
            if kind == 'vr':
                a, b = (x, top(src) ), (x, src['y'] - 1.0)
            elif kind == 'absorption':
                a, b = (x, src['y'] - 1.0), (x, top(dst) + 1.0)
            elif kind == 'ic':
                a, b = (x, src['y'] + 1.0), (x, top(dst) - 1.0)
            else:  # emission ends on an excited vibrational level of the lower state
                a, b = (x, src['y'] + 1.0), (x, dst['y'] - min(dst['vib'], 2) * vib_gap - 1.0)
            if _TRANSITIONS[kind]['radiative']:
                c.arrow(a, b, color=t['color'], head_size=800)
            else:
                c.wavy(a, b, color=t['color'])
            c.block((x, a[1], x, b[1]), 3.0)
            mid = (a[1] + b[1]) / 2
            candidates = [((x + 5.0, mid), 'Left', 'middle'), ((x - 5.0, mid), 'Right', 'middle'),
                          ((x + 5.0, mid - 12.0), 'Left', 'middle'), ((x + 5.0, mid + 12.0), 'Left', 'middle')]
        if t['label']:
            _, clear = c.place(candidates, str(t['label']), color=t['color'])
            if not clear:
                warnings.append(f'the label of {kind} {src["id"]} -> {dst["id"]} overlaps other ink')
        drawn.append({'kind': kind, 'from': src['id'], 'to': dst['id'],
                      'radiative': _TRANSITIONS[kind]['radiative'], 'energy_change': dst['energy'] - src['energy']})
    if spec.get('axis', True):
        left = min(box[0] for box in c.text_boxes) - 12.0
        c.arrow((left, ground['y'] + 6.0), (left, min(top(s) for s in ordered) - 14.0))
        c.label((left - 4.0, min(top(s) for s in ordered) - 20.0), spec.get('y_label', 'Energy'))
    report = {'diagram': 'jablonski', 'states': [{'id': s['id'], 'energy': s['energy'], 'multiplicity': s['mult']}
                                                  for s in states.values()],
              'transitions': drawn, 'warnings': warnings}
    return _result(c, output, report)
