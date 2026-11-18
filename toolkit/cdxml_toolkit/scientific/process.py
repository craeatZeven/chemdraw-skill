"""Process diagrams: flow charts of a preparation or separation, and electrochemical cells.

A flow chart is laid out from its connections, so adding a step never needs coordinates.
An electrochemical cell is drawn from the two half-cells and the cell type; the polarity,
the direction of electron flow and the ion migration are derived, not stated, so the
drawing cannot contradict itself.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

from .canvas import Canvas, plain, resolve_color, text_height, text_width
from .energy import MINUS, _result

SHAPES = ('box', 'round', 'ellipse', 'plain')


# ---------------------------------------------------------------------------- flow chart
def flowchart(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Left-to-right flow chart.

    ``nodes``: ``{"id", "text", "shape": box|round|ellipse|plain, "row", "col"}`` (``row`` and
    ``col`` are optional). ``edges``: ``{"from", "to", "label", "below", "dashed"}``; an edge
    that points back to an earlier column is drawn as a return loop. ``inputs`` and
    ``outputs`` attach a short arrow with a caption to the top or bottom of a node:
    ``{"to"|"from": id, "text", "side": "top"|"bottom"}``.
    """
    raw_nodes = spec.get('nodes')
    if not isinstance(raw_nodes, list) or not raw_nodes:
        raise ValueError('A flow chart needs "nodes"')
    nodes = {}
    for i, item in enumerate(raw_nodes):
        nid = str(item.get('id') or f'n{i + 1}')
        if nid in nodes:
            raise ValueError(f'Duplicate node id {nid!r}')
        shape = str(item.get('shape', 'box')).lower()
        if shape not in SHAPES:
            raise ValueError(f'node {nid}: shape must be one of {", ".join(SHAPES)}')
        text = str(item.get('text', nid))
        pad = 3.0 if shape == 'plain' else 7.0
        nodes[nid] = {'id': nid, 'text': text, 'shape': shape, 'row': item.get('row'), 'col': item.get('col'),
                      'w': max(text_width(text) + 2 * pad, 0.0 if shape == 'plain' else 30.0) * (1.25 if shape == 'ellipse' else 1.0),
                      'h': max(text_height(text) + (2.0 if shape == 'plain' else 8.0), 20.0),
                      'color': item.get('color'), 'order': i}
    edges = []
    for i, item in enumerate(spec.get('edges') or []):
        src, dst = str(item.get('from')), str(item.get('to'))
        if src not in nodes or dst not in nodes:
            raise ValueError(f'edge {i + 1}: unknown node in {src!r} -> {dst!r}')
        if src == dst:
            raise ValueError(f'edge {i + 1}: an edge joins two different nodes')
        route = item.get('route')
        if route not in (None, 'top', 'bottom'):
            raise ValueError(f'edge {i + 1}: route must be "top" or "bottom"')
        edges.append({'from': src, 'to': dst, 'label': item.get('label'), 'below': item.get('below'),
                      'dashed': bool(item.get('dashed')), 'color': item.get('color'), 'back': bool(item.get('back')),
                      'route': route})

    # ---- columns: an edge back to a node still being explored closes a loop; the rest
    # give each node the length of the longest path that reaches it.
    state = {}

    def explore(nid):
        state[nid] = 'open'
        for edge in edges:
            if edge['from'] != nid or edge['back']:
                continue
            if state.get(edge['to']) == 'open':
                edge['back'] = True
            elif edge['to'] not in state:
                explore(edge['to'])
        state[nid] = 'done'

    for nid in nodes:
        if nid not in state:
            explore(nid)
    depth = {nid: 0 for nid in nodes}
    for _ in nodes:
        for edge in edges:
            if not edge['back']:
                depth[edge['to']] = max(depth[edge['to']], depth[edge['from']] + 1)
    used = {}
    for nid, node in nodes.items():
        if node['col'] is None:
            node['col'] = depth[nid]
        node['col'] = int(node['col'])
    for nid, node in sorted(nodes.items(), key=lambda kv: kv[1]['order']):
        if node['row'] is None:
            row = 0
            parents = [nodes[e['from']] for e in edges if e['to'] == nid and not e['back']
                       and nodes[e['from']]['row'] is not None]
            if parents:
                row = int(parents[0]['row'])
            while (node['col'], row) in used:
                row += 1
            node['row'] = row
        node['row'] = int(node['row'])
        if (node['col'], node['row']) in used:
            raise ValueError(f'nodes {used[(node["col"], node["row"])]} and {nid} share row {node["row"]}, '
                             f'column {node["col"]}')
        used[(node['col'], node['row'])] = nid

    cols = sorted({n['col'] for n in nodes.values()})
    rows = sorted({n['row'] for n in nodes.values()})
    col_w = {col: max(n['w'] for n in nodes.values() if n['col'] == col) for col in cols}
    row_h = {row: max(n['h'] for n in nodes.values() if n['row'] == row) for row in rows}
    gap_after = {col: float(spec.get('gap', 34.0)) for col in cols}
    for edge in edges:
        if edge['back']:
            continue
        a, b = nodes[edge['from']], nodes[edge['to']]
        if b['col'] > a['col']:
            need = max(text_width(str(edge[k])) for k in ('label', 'below') if edge[k]) + 16.0 if (
                edge['label'] or edge['below']) else 0.0
            span = [col for col in cols if a['col'] <= col < b['col']]
            share = need / len(span)
            for col in span:
                gap_after[col] = max(gap_after[col], share)
    side_notes = []
    for kind in ('inputs', 'outputs'):
        for i, item in enumerate(spec.get(kind) or []):
            target = str(item.get('to' if kind == 'inputs' else 'from'))
            if target not in nodes:
                raise ValueError(f'{kind} {i + 1}: unknown node {target!r}')
            side = str(item.get('side', 'top' if kind == 'inputs' else 'bottom')).lower()
            if side not in ('top', 'bottom'):
                raise ValueError(f'{kind} {i + 1}: side must be "top" or "bottom"')
            side_notes.append({'node': target, 'text': str(item.get('text', '')), 'side': side,
                               'incoming': kind == 'inputs', 'dashed': bool(item.get('dashed'))})
    above = {row: max([0.0] + [text_height(n['text']) + 22.0 for n in side_notes
                               if n['side'] == 'top' and nodes[n['node']]['row'] == row]) for row in rows}
    below = {row: max([0.0] + [text_height(n['text']) + 22.0 for n in side_notes
                               if n['side'] == 'bottom' and nodes[n['node']]['row'] == row]) for row in rows}
    loops = [e for e in edges if e['back']]
    for edge in loops:       # a loop keeps clear of the side arrows of the nodes it joins
        if edge['route'] is None:
            ends = (edge['from'], edge['to'])
            top_busy = any(n['node'] in ends and n['side'] == 'top' for n in side_notes)
            bottom_busy = any(n['node'] in ends and n['side'] == 'bottom' for n in side_notes)
            edge['route'] = 'bottom' if top_busy and not bottom_busy else 'top'
    top_loops = [e for e in loops if e['route'] == 'top']
    bottom_loops = [e for e in loops if e['route'] == 'bottom']
    x_of, y_of = {}, {}
    x = 0.0
    for col in cols:
        x_of[col] = x + col_w[col] / 2
        x += col_w[col] + gap_after[col]
    y = 0.0
    for row in rows:
        y += above[row]
        y_of[row] = y + row_h[row] / 2
        y += row_h[row] + below[row] + float(spec.get('row_gap', 26.0))
    chart_top = -above[rows[0]] - 8.0
    chart_bottom = y - float(spec.get('row_gap', 26.0)) + 8.0
    for node in nodes.values():
        node['x'], node['y'] = x_of[node['col']], y_of[node['row']]

    c = Canvas()
    warnings: list[str] = []
    for node in nodes.values():
        x1, y1 = node['x'] - node['w'] / 2, node['y'] - node['h'] / 2
        x2, y2 = node['x'] + node['w'] / 2, node['y'] + node['h'] / 2
        color = resolve_color(node['color']) if node['color'] else None
        if node['shape'] in ('box', 'round'):
            c.rect(x1, y1, x2, y2, color=color, rounded=node['shape'] == 'round')
        elif node['shape'] == 'ellipse':
            steps = 48
            c.path([(node['x'] + node['w'] / 2 * math.cos(2 * math.pi * k / steps),
                     node['y'] + node['h'] / 2 * math.sin(2 * math.pi * k / steps)) for k in range(steps)],
                   color=color, closed=True)
        c.label((node['x'], node['y']), node['text'], align='Center', valign='middle', color=color)
        c.block((x1, y1, x2, y2), 1.0)

    for edge in edges:
        if edge['back']:
            continue
        a, b = nodes[edge['from']], nodes[edge['to']]
        color = resolve_color(edge['color']) if edge['color'] else None
        if a['col'] == b['col']:
            down = b['y'] > a['y']
            start = (a['x'], a['y'] + (a['h'] / 2 if down else -a['h'] / 2))
            end = (b['x'], b['y'] - (b['h'] / 2 if down else -b['h'] / 2))
            c.arrow(start, end, color=color, dashed=edge['dashed'])
            mid = (start[1] + end[1]) / 2
            if edge['label']:
                c.label((a['x'] + 5.0, mid), str(edge['label']), valign='middle', color=color)
            if edge['below']:
                c.label((a['x'] - 5.0, mid), str(edge['below']), align='Right', valign='middle', color=color)
            continue
        forward = b['col'] > a['col']
        start = (a['x'] + (a['w'] / 2 if forward else -a['w'] / 2), a['y'])
        end = (b['x'] - (b['w'] / 2 if forward else -b['w'] / 2), b['y'])
        if abs(a['y'] - b['y']) < 1e-6:
            blockers = [n['id'] for n in nodes.values() if n is not a and n is not b and n['row'] == a['row']
                        and min(a['col'], b['col']) < n['col'] < max(a['col'], b['col'])]
            if blockers:
                warnings.append(f'the arrow {a["id"]} -> {b["id"]} passes through {", ".join(blockers)}; '
                                'put one of them on another row')
            c.arrow(start, end, color=color, dashed=edge['dashed'])
            label_x, label_y = (start[0] + end[0]) / 2, a['y']
        else:
            # Leave horizontally, turn in the gap before the target column, enter horizontally.
            turn = end[0] - (gap_after[cols[cols.index(b['col']) - 1]] / 2 if forward else -17.0)
            c.path([start, (turn, start[1]), (turn, end[1])], color=color, dashed=edge['dashed'])
            c.arrow((turn, end[1]), end, color=color, dashed=edge['dashed'])
            label_x, label_y = (turn + end[0]) / 2, end[1]
            if abs(end[0] - turn) < text_width(str(edge['label'] or edge['below'] or '')) + 6.0:
                label_x, label_y = (start[0] + turn) / 2, start[1]
        if edge['label']:
            c.label((label_x, label_y - 4.0 - (str(edge['label']).count('\n')) * 12.0), str(edge['label']),
                    align='Center', color=color)
        if edge['below']:
            c.label((label_x, label_y + 12.0), str(edge['below']), align='Center', color=color)

    for group, sign, base in ((top_loops, -1.0, chart_top), (bottom_loops, 1.0, chart_bottom)):
        for k, edge in enumerate(group):       # return loops run outside the chart
            a, b = nodes[edge['from']], nodes[edge['to']]
            color = resolve_color(edge['color']) if edge['color'] else None
            level = base + sign * 16.0 * k
            ax, bx = a['x'] + min(8.0, a['w'] / 4), b['x'] - min(8.0, b['w'] / 4)
            c.path([(ax, a['y'] + sign * a['h'] / 2), (ax, level), (bx, level)], color=color, dashed=edge['dashed'])
            c.arrow((bx, level), (bx, b['y'] + sign * b['h'] / 2), color=color, dashed=edge['dashed'])
            if edge['label']:
                c.label(((ax + bx) / 2, level - 4.0 if sign < 0 else level + 12.0), str(edge['label']),
                        align='Center', color=color)

    for note in side_notes:
        node = nodes[note['node']]
        offsets = [n for n in side_notes if n['node'] == note['node'] and n['side'] == note['side']]
        shift = (offsets.index(note) - (len(offsets) - 1) / 2) * min(26.0, node['w'] / max(len(offsets), 1))
        px = node['x'] + shift
        align = 'Center'
        # A return loop on the same side takes one corner; the caption moves to the other.
        for edge in loops:
            if edge['route'] != note['side'] or note['node'] not in (edge['from'], edge['to']):
                continue
            away = 1.0 if note['node'] == edge['to'] else -1.0
            px = node['x'] + away * min(8.0, node['w'] / 4)
            align = 'Left' if away > 0 else 'Right'
        if note['side'] == 'top':
            edge_y, far_y = node['y'] - node['h'] / 2, node['y'] - node['h'] / 2 - 18.0
            text_at, valign = (px, far_y - 3.0 - (note['text'].count('\n')) * 12.0), 'baseline'
        else:
            edge_y, far_y = node['y'] + node['h'] / 2, node['y'] + node['h'] / 2 + 18.0
            text_at, valign = (px, far_y + 2.0), 'top'
        if note['incoming']:
            c.arrow((px, far_y), (px, edge_y), dashed=note['dashed'])
        else:
            c.arrow((px, edge_y), (px, far_y), dashed=note['dashed'])
        if note['text']:
            anchor = text_at if align == 'Center' else (text_at[0] + (-3.0 if align == 'Left' else 3.0), text_at[1])
            c.label(anchor, note['text'], align=align, valign=valign)
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    for a_text, b_text in c.overlapping_labels():
        warnings.append(f'the captions "{a_text}" and "{b_text}" overlap; increase "gap" or "row_gap"')
    if c.extent[2] - c.extent[0] > 504.0:
        warnings.append(f'the chart is {c.extent[2] - c.extent[0]:.0f} pt wide, more than a double-column page '
                        '(504 pt); move later steps to a second row with "row" and "col"')
    report = {'diagram': 'flowchart', 'nodes': len(nodes), 'edges': len(edges),
              'return_loops': [f'{e["from"]} -> {e["to"]}' for e in loops],
              'layout': {nid: {'row': n['row'], 'col': n['col']} for nid, n in nodes.items()},
              'warnings': warnings}
    return _result(c, output, report)


# ---------------------------------------------------------------------------- electrochemical cell
_DEFAULT_LABELS = {
    'galvanic': {'anode': 'anode (' + MINUS + ')', 'cathode': 'cathode (+)'},
    'electrolytic': {'anode': 'anode (+)', 'cathode': 'cathode (' + MINUS + ')'},
}


def electrochemical_cell(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Galvanic or electrolytic cell.

    ``{"type": "galvanic"|"electrolytic", "anode": {"electrode", "solution", "reaction"},
    "cathode": {...}, "anode_side": "left"|"right", "bridge": "salt_bridge"|"membrane"|"none",
    "meter": "V"}``. Oxidation happens at the anode in both cell types; the anode is the
    negative pole of a galvanic cell and the positive pole of an electrolytic cell, electrons
    leave the anode through the wire, anions migrate to the anode and cations to the cathode.
    """
    kind = str(spec.get('type', 'galvanic')).lower()
    if kind not in ('galvanic', 'electrolytic'):
        raise ValueError('type must be "galvanic" or "electrolytic"')
    anode, cathode = spec.get('anode') or {}, spec.get('cathode') or {}
    if not anode.get('electrode') or not cathode.get('electrode'):
        raise ValueError('Give the "electrode" material of both "anode" and "cathode"')
    side = str(spec.get('anode_side', 'left')).lower()
    if side not in ('left', 'right'):
        raise ValueError('anode_side must be "left" or "right"')
    bridge = str(spec.get('bridge', 'salt_bridge' if kind == 'galvanic' else 'none')).lower()
    if bridge not in ('salt_bridge', 'membrane', 'none'):
        raise ValueError('bridge must be "salt_bridge", "membrane" or "none"')
    labels = {**_DEFAULT_LABELS[kind], **(spec.get('labels') or {})}
    polarity = {'anode': MINUS if kind == 'galvanic' else '+', 'cathode': '+' if kind == 'galvanic' else MINUS}
    halves = {'left': ('anode', anode) if side == 'left' else ('cathode', cathode),
              'right': ('cathode', cathode) if side == 'left' else ('anode', anode)}

    c = Canvas()
    beaker_w, beaker_h, top = 74.0, 74.0, 0.0
    two = bridge == 'salt_bridge'
    if two:
        boxes = {'left': (0.0, beaker_w), 'right': (beaker_w + 44.0, 2 * beaker_w + 44.0)}
    else:
        wide = 150.0
        boxes = {'left': (0.0, wide / 2), 'right': (wide / 2, wide)}
    surface = top + 16.0
    if two:
        for x1, x2 in boxes.values():
            c.path([(x1, top), (x1, top + beaker_h), (x2, top + beaker_h), (x2, top)])
            c.segment((x1, surface), (x2, surface))
    else:
        x1, x2 = boxes['left'][0], boxes['right'][1]
        c.path([(x1, top), (x1, top + beaker_h), (x2, top + beaker_h), (x2, top)])
        c.segment((x1, surface), (x2, surface))
        if bridge == 'membrane':
            mx = (x1 + x2) / 2
            c.segment((mx, surface - 6.0), (mx, top + beaker_h), dashed=True)
            c.label((mx, top + beaker_h + 13.0), spec.get('bridge_label', 'membrane'), align='Center')
    electrode_x = {}
    for name, (x1, x2) in boxes.items():
        # Electrodes sit toward the outer wall, leaving the middle for the bridge or membrane.
        centre = (x1 + x2) / 2 + (-1.0 if name == 'left' else 1.0) * (12.0 if two else 14.0)
        electrode_x[name] = centre
        c.rect(centre - 4.0, top - 16.0, centre + 4.0, top + beaker_h - 14.0)
    ions = spec.get('ions') or {}
    wire_y = top - (72.0 if two and ions else 44.0)
    lx, rx = electrode_x['left'], electrode_x['right']
    mid = (lx + rx) / 2
    meter = spec.get('meter', 'V' if kind == 'galvanic' else 'battery')
    if meter == 'battery' or kind == 'electrolytic' and meter in (None, 'battery'):
        # Long plate = positive pole, on the anode side of an electrolytic cell.
        positive_left = (halves['left'][0] == 'anode') == (kind == 'electrolytic')
        long_x, short_x = (mid - 3.0, mid + 3.0) if positive_left else (mid + 3.0, mid - 3.0)
        c.path([(lx, top - 16.0), (lx, wire_y), (min(long_x, short_x), wire_y)])
        c.path([(max(long_x, short_x), wire_y), (rx, wire_y), (rx, top - 16.0)])
        c.segment((long_x, wire_y - 8.0), (long_x, wire_y + 8.0))
        c.segment((short_x, wire_y - 4.0), (short_x, wire_y + 4.0))
        c.label((long_x + (-6.0 if positive_left else 6.0), wire_y - 9.0), '+', align='Center')
        c.label((short_x + (6.0 if positive_left else -6.0), wire_y - 9.0), MINUS, align='Center')
    else:
        radius = 9.0
        c.path([(lx, top - 16.0), (lx, wire_y), (mid - radius, wire_y)])
        c.path([(mid + radius, wire_y), (rx, wire_y), (rx, top - 16.0)])
        if meter:
            c.circle((mid, wire_y), radius)
            c.label((mid, wire_y), str(meter), align='Center', valign='middle')
        else:
            c.segment((mid - radius, wire_y), (mid + radius, wire_y))
    if two:
        bx1, bx2 = boxes['left'][1] - 16.0, boxes['right'][0] + 16.0
        y_in, y_top = top + 44.0, top - 14.0
        c.path([(bx1, y_in), (bx1, y_top), (bx2, y_top), (bx2, y_in)])
        c.path([(bx1 + 8.0, y_in), (bx1 + 8.0, y_top + 8.0), (bx2 - 8.0, y_top + 8.0), (bx2 - 8.0, y_in)])
        c.segment((bx1, y_in), (bx1 + 8.0, y_in))
        c.segment((bx2 - 8.0, y_in), (bx2, y_in))
        c.label(((bx1 + bx2) / 2, y_top - 5.0), spec.get('bridge_label', 'salt bridge'), align='Center')

    electrons_right = halves['left'][0] == 'anode'          # electrons leave the anode
    if spec.get('show_electron_flow', True):
        y = wire_y - 7.0
        a, b = (lx + 8.0, y), (mid - 16.0, y)
        c.arrow(a if electrons_right else b, b if electrons_right else a, head_size=800)
        c.label(((a[0] + b[0]) / 2, y - 5.0), 'e{sup:%s}' % MINUS, align='Center')
    for name, (role, half) in halves.items():
        x1, x2 = boxes[name]
        ex = electrode_x[name]
        outward = -1.0 if name == 'left' else 1.0
        c.label((ex + outward * 8.0, top - 24.0), str(half['electrode']), align='Right' if name == 'left' else 'Left')
        lines = [labels[role]] if spec.get('show_polarity', True) else []
        if half.get('solution'):
            lines.insert(0, str(half['solution']))
        if half.get('reaction'):
            lines.append(str(half['reaction']))
        centre = (x1 + x2) / 2 if two else (x1 + x2) / 2 + outward * 4.0
        y = top + beaker_h + (13.0 if two or bridge == 'none' else 27.0)
        for line in lines:
            c.label((centre if two else (boxes['left'][0] if name == 'left' else boxes['right'][1]), y), line,
                    align='Center' if two else ('Left' if name == 'left' else 'Right'))
            y += 13.0
    anode_side_name = 'left' if halves['left'][0] == 'anode' else 'right'
    if ions:
        # Anions migrate to the anode, cations to the cathode: above the bridge, or in the solution.
        base = (boxes['left'][0] + boxes['right'][1]) / 2
        y = top - 34.0 if two else top + 40.0
        for key, toward in (('anions', anode_side_name), ('cations', 'right' if anode_side_name == 'left' else 'left')):
            if not ions.get(key):
                continue
            direction = -1.0 if toward == 'left' else 1.0
            c.arrow((base - direction * 11.0, y), (base + direction * 11.0, y), head_size=700)
            c.label((base + direction * 14.0, y), str(ions[key]), align='Left' if direction > 0 else 'Right',
                    valign='middle')
            y += -12.0 if two else 14.0
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    warnings = [f'the captions "{a}" and "{b}" overlap' for a, b in c.overlapping_labels()]
    report = {'diagram': 'electrochemical_cell', 'type': kind,
              'anode': {'side': anode_side_name, 'electrode': plain(str(anode['electrode'])), 'process': 'oxidation',
                        'polarity': '-' if kind == 'galvanic' else '+'},
              'cathode': {'side': 'right' if anode_side_name == 'left' else 'left',
                          'electrode': plain(str(cathode['electrode'])), 'process': 'reduction',
                          'polarity': '+' if kind == 'galvanic' else '-'},
              'electron_flow': 'anode to cathode through the external circuit',
              'ion_migration': {'anions': 'toward the anode', 'cations': 'toward the cathode'},
              'warnings': warnings}
    return _result(c, output, report)
