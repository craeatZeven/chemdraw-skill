"""Periodic table with highlighted elements and trend arrows."""
from __future__ import annotations

from typing import Any, Mapping

from .canvas import Canvas, resolve_color
from .elements import NOBLE_GASES, SYMBOLS, atomic_number, group_period
from .energy import _result

_BLOCK_OF_GROUP = {**{g: 's' for g in (1, 2)}, **{g: 'd' for g in range(3, 13)}, **{g: 'p' for g in range(13, 19)}}


def table_position(z):
    """``(row, column, block)`` of element ``z`` in the 18-column table; the f block sits in two rows below."""
    if 57 <= z <= 71:
        return 8, z - 57 + 3, 'f'
    if 89 <= z <= 103:
        return 9, z - 89 + 3, 'f'
    group, period = group_period(z)
    block = 's' if z == 2 else _BLOCK_OF_GROUP[group]
    return period, group, block


def periodic_table(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Periodic table.

    ``periods``: ``[first, last]`` (default 1 to 7); ``f_block``: draw the lanthanoids and
    actinoids (default true when period 6 is shown); ``show_number``: atomic numbers;
    ``highlight``: a list of ``{"elements": [...]}``, ``{"group": 17}``, ``{"period": 3}``
    or ``{"block": "d"}`` with ``"color"`` and an optional legend ``"label"``;
    ``trends``: ``{"direction": "right"|"left"|"up"|"down", "label"}`` arrows beside the table.
    """
    periods = spec.get('periods', [1, 7])
    if len(periods) != 2 or not 1 <= int(periods[0]) <= int(periods[1]) <= 7:
        raise ValueError('"periods" is [first, last] within 1 to 7')
    first, last = int(periods[0]), int(periods[1])
    f_block = bool(spec.get('f_block', last >= 6))
    cell = float(spec.get('cell', 24.0))
    if not 16 <= cell <= 60:
        raise ValueError('"cell" must be between 16 and 60 points')
    show_number = bool(spec.get('show_number', False))
    marks = {}
    legend = []
    for i, item in enumerate(spec.get('highlight') or []):
        color = resolve_color(item.get('color'), '#FF6600')
        chosen = set()
        if item.get('elements'):
            chosen |= {atomic_number(e) for e in item['elements']}
        for z in range(1, len(SYMBOLS) + 1):
            row, column, block = table_position(z)
            group, period = group_period(z)
            if item.get('group') is not None and row <= 7 and column == int(item['group']) and block != 'f':
                chosen.add(z)
            if item.get('period') is not None and period == int(item['period']):
                chosen.add(z)
            if item.get('block') is not None and block == str(item['block']).lower():
                chosen.add(z)
        if not chosen:
            raise ValueError(f'highlight {i + 1} selects no element')
        for z in chosen:
            marks[z] = color
        if item.get('label'):
            legend.append((color, str(item['label'])))

    c = Canvas()
    drawn = []
    rows_used = set()
    highlighted = []
    for z in range(1, len(SYMBOLS) + 1):
        row, column, block = table_position(z)
        _group, period = group_period(z)
        if not first <= period <= last or (block == 'f' and not f_block):
            continue
        y_row = (row - first) if row <= 7 else (last - first + 1.5 + (row - 8))
        x1, y1 = (column - 1) * cell, y_row * cell
        rows_used.add(row)
        drawn.append(z)
        color = marks.get(z)
        if color is None:
            c.rect(x1, y1, x1 + cell, y1 + cell)
        else:
            highlighted.append((z, x1, y1, color))
        c.label((x1 + cell / 2, y1 + cell / 2 + (2.0 if show_number else 0.0)), SYMBOLS[z - 1], align='Center',
                valign='middle', color=color, bold=color is not None)
        if show_number:
            c.label((x1 + 2.0, y1 + 7.0), str(z), size=6)
    for z, x1, y1, color in highlighted:        # coloured frames go on top of their neighbours' edges
        c.rect(x1 + 0.6, y1 + 0.6, x1 + cell - 0.6, y1 + cell - 0.6, color=color)
        c.rect(x1, y1, x1 + cell, y1 + cell, color=color)
    for period, marker in ((6, '*'), (7, '**')):
        if first <= period <= last:
            x1, y1 = 2 * cell, (period - first) * cell
            c.rect(x1, y1, x1 + cell, y1 + cell)
            c.label((x1 + cell / 2, y1 + cell / 2), marker, align='Center', valign='middle')
            if f_block:
                c.label((2 * cell - 5.0, (last - first + 1.5 + (period - 6)) * cell + cell / 2), marker, align='Right',
                        valign='middle')
    if spec.get('labels', True):
        for group in range(1, 19):
            top_row = min((table_position(z)[0] for z in drawn if table_position(z)[1] == group
                           and table_position(z)[2] != 'f'), default=None)
            if top_row is not None:
                c.label(((group - 0.5) * cell, (top_row - first) * cell - 4.0), str(group), align='Center', size=8)
        for period in range(first, last + 1):
            c.label((-6.0, (period - first + 0.5) * cell), str(period), align='Right', valign='middle', size=8)
    x1, y1, x2, y2 = c.extent
    for i, trend in enumerate(spec.get('trends') or []):
        direction = str(trend.get('direction', 'right')).lower()
        label = str(trend.get('label', ''))
        if direction in ('right', 'left'):
            y = y1 - 14.0 - 22.0 * sum(1 for t in (spec.get('trends') or [])[:i]
                                       if str(t.get('direction', 'right')).lower() in ('right', 'left'))
            a, b = (2 * cell, y), (16 * cell, y)
            c.arrow(a if direction == 'right' else b, b if direction == 'right' else a)
            c.label((9 * cell, y - 4.0), label, align='Center')
        elif direction in ('up', 'down'):
            x = x1 - 16.0 - 22.0 * sum(1 for t in (spec.get('trends') or [])[:i]
                                       if str(t.get('direction', 'right')).lower() in ('up', 'down'))
            top, bottom = 0.5 * cell, (last - first + 0.5) * cell
            c.arrow((x, top) if direction == 'down' else (x, bottom), (x, bottom) if direction == 'down' else (x, top))
            c.label((x - 5.0, (top + bottom) / 2), label, align='Right', valign='middle')
        else:
            raise ValueError('trend direction must be right, left, up or down')
    if legend:
        x = 0.0
        y = c.extent[3] + 14.0
        for color, label in legend:
            c.rect(x, y - 8.0, x + 10.0, y + 2.0, color=color)
            box = c.label((x + 14.0, y), label, color=color)
            x = box[2] + 16.0
    title = spec.get('title')
    if title:
        ex1, _ey1, ex2, ey2 = c.extent
        c.label(((ex1 + ex2) / 2, ey2 + 16.0), str(title), align='Center')
    report = {'diagram': 'periodic_table', 'elements_drawn': len(drawn),
              'highlighted': [SYMBOLS[z - 1] for z in sorted(marks) if z in drawn], 'warnings': []}
    return _result(c, output, report)
