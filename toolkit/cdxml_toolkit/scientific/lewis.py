"""Lewis (electron-dot) structures drawn on a square grid.

The spec states the atoms, their grid positions and the bonds; the non-bonding electrons
follow from the valence count (``valence - bonds - formal charge``), so a structure cannot
gain or lose electrons by a drawing slip. The report lists every atom's electron count and
flags atoms that miss or exceed an octet.
"""
from __future__ import annotations

from typing import Any, Mapping

from .canvas import BLACK, Canvas, text_width
from .elements import SYMBOLS, atomic_number, group_period, valence_electrons
from .energy import MINUS, _result

SIDES = {'up': (0, -1), 'right': (1, 0), 'down': (0, 1), 'left': (-1, 0)}
_OPPOSITE = {'up': 'down', 'down': 'up', 'left': 'right', 'right': 'left'}
DOT = 0.9          # radius of an electron dot
PAIR = 2.3         # half the distance between the two electrons of a pair


def _mark(c, x, y, cross):
    if cross:
        c.segment((x - 1.3, y - 1.3), (x + 1.3, y + 1.3))
        c.segment((x - 1.3, y + 1.3), (x + 1.3, y - 1.3))
    else:
        c.circle((x, y), DOT, fill=BLACK, outline=None)


def _species(raw, index, style):
    atoms = {}
    for i, item in enumerate(raw.get('atoms') or []):
        symbol = SYMBOLS[atomic_number(item.get('element')) - 1]
        aid = str(item.get('id') or f'{symbol}{i + 1}')
        if aid in atoms:
            raise ValueError(f'species {index}: duplicate atom id {aid!r}')
        at = item.get('at', [i, 0])
        if len(at) != 2 or any(int(v) != v for v in at):
            raise ValueError(f'species {index}: atom {aid} needs integer grid coordinates "at": [column, row]')
        valence = item.get('valence')
        atoms[aid] = {'id': aid, 'symbol': symbol, 'col': int(at[0]), 'row': int(at[1]), 'bonds': {}, 'order': i,
                      'valence': int(valence) if valence is not None else valence_electrons(symbol),
                      'fc': item.get('formal_charge'), 'lone_pairs': item.get('lone_pairs'),
                      'radical': item.get('radical')}
    if not atoms:
        raise ValueError(f'species {index}: no atoms')
    if len({(a['col'], a['row']) for a in atoms.values()}) != len(atoms):
        raise ValueError(f'species {index}: two atoms share a grid position')
    bonds = []
    for item in raw.get('bonds') or []:
        if isinstance(item, Mapping):
            a_id, b_id, order = str(item.get('a')), str(item.get('b')), int(item.get('order', 1))
        else:
            a_id, b_id = str(item[0]), str(item[1])
            order = int(item[2]) if len(item) > 2 else 1
        if a_id not in atoms or b_id not in atoms:
            raise ValueError(f'species {index}: bond {a_id}-{b_id} names an unknown atom')
        if order not in (1, 2, 3):
            raise ValueError(f'species {index}: bond {a_id}-{b_id} must have order 1, 2 or 3')
        a, b = atoms[a_id], atoms[b_id]
        step = (b['col'] - a['col'], b['row'] - a['row'])
        side = next((name for name, vec in SIDES.items() if vec == step), None)
        if side is None:
            raise ValueError(f'species {index}: bonded atoms {a_id} and {b_id} must be neighbours on the grid '
                             '(one step left, right, up or down)')
        if side in a['bonds'] or _OPPOSITE[side] in b['bonds']:
            raise ValueError(f'species {index}: {a_id} and {b_id} are bonded twice; use "order" for multiple bonds')
        a['bonds'][side] = (b, order)
        b['bonds'][_OPPOSITE[side]] = (a, order)
        bonds.append((a, b, order, side))

    warnings, rows = [], []
    for a in atoms.values():
        bonded = sum(order for _other, order in a['bonds'].values())
        if a['lone_pairs'] is not None or a['radical'] is not None:
            nonbonding = 2 * int(a['lone_pairs'] or 0) + int(a['radical'] or 0)
            implied = a['valence'] - bonded - nonbonding
            if a['fc'] is not None and int(a['fc']) != implied:
                raise ValueError(
                    f'species {index}: {a["id"]} with {bonded} bond(s) and {nonbonding} non-bonding electron(s) '
                    f'has formal charge {implied:+d}, not {int(a["fc"]):+d}')
            a['fc'] = implied
        else:
            a['fc'] = int(a['fc'] or 0)
            nonbonding = a['valence'] - bonded - a['fc']
            if nonbonding < 0:
                raise ValueError(f'species {index}: {a["id"]} has {a["valence"]} valence electrons, which is not '
                                 f'enough for {bonded} bond(s) at formal charge {a["fc"]:+d}')
        a['nonbonding'] = nonbonding
        a['shell'] = nonbonding + 2 * bonded
        _group, period = group_period(a['symbol'])
        target = 2 if a['symbol'] in ('H', 'He') else 8
        status = 'complete'
        if a['shell'] == 0 and not a['bonds']:
            status = 'cation'                       # a bare cation shows its empty valence shell
        elif a['shell'] < target:
            status = 'incomplete'
            warnings.append(f'{a["id"]} has {a["shell"]} electrons around it (fewer than {target})')
        elif a['shell'] > target:
            status = 'expanded'
            if period < 3:
                warnings.append(f'{a["id"]} has {a["shell"]} electrons around it; a period-{period} atom cannot '
                                'exceed an octet')
        if nonbonding % 2:
            warnings.append(f'{a["id"]} carries an unpaired electron (radical)')
        rows.append({'atom': a['id'], 'element': a['symbol'], 'bonding_electrons': 2 * bonded,
                     'nonbonding_electrons': nonbonding, 'formal_charge': a['fc'], 'shell': a['shell'],
                     'octet': status})
    charge = sum(a['fc'] for a in atoms.values())
    if raw.get('charge') is not None and int(raw['charge']) != charge:
        raise ValueError(f'species {index}: the formal charges add up to {charge:+d}, not the stated charge '
                         f'{int(raw["charge"]):+d}')
    total = sum(a['valence'] for a in atoms.values()) - charge
    bare = len(atoms) == 1 and next(iter(atoms.values()))['nonbonding'] == 0
    bracket = raw.get('bracket')
    if bracket is None:
        bracket = charge != 0 and not bare
    return {'atoms': atoms, 'bonds': bonds, 'charge': charge, 'bracket': bool(bracket), 'bare': bare,
            'report': {'charge': charge, 'valence_electrons': total, 'atoms': rows}, 'warnings': warnings}


def _groups(atom):
    """Non-bonding electrons per free side: ``{side: 1 or 2}``."""
    free = [s for s in ('up', 'down', 'right', 'left') if s not in atom['bonds']]
    # Opposite a bond first, so that lone pairs sit away from the bonded neighbour.
    free.sort(key=lambda s: (_OPPOSITE[s] not in atom['bonds'], ('up', 'down', 'right', 'left').index(s)))
    n = atom['nonbonding']
    if not atom['bonds'] and n < 8:
        # A free atom: one electron per side first, then pair (the element's dot symbol).
        order = ['right', 'up', 'left', 'down']
        counts = {s: 0 for s in order}
        for k in range(n):
            counts[order[k % 4]] += 1
        return {s: v for s, v in counts.items() if v}
    pairs, single = divmod(n, 2)
    if len(atom['bonds']) == 1 and pairs + single == 2:
        # Two groups on a terminal atom sit either side of the bond axis (:O=C=O: is symmetric).
        axis = next(iter(atom['bonds']))
        free = [s for s in free if s not in (axis, _OPPOSITE[axis])] + [_OPPOSITE[axis]]
    if pairs + single > len(free):
        raise ValueError(f'{atom["id"]} needs {pairs + single} places for its non-bonding electrons but only '
                         f'{len(free)} side(s) are free; give it fewer neighbours on the grid')
    out = {side: 2 for side in free[:pairs]}
    if single:
        out[free[pairs]] = 1
    return out


def lewis_structure(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Lewis structure of one or more species.

    ``{"atoms": [{"id", "element", "at": [col, row], "formal_charge"}], "bonds": [[a, b, order]]}``
    for one species, or ``"species": [...]`` for several drawn left to right (an ionic
    compound). ``style`` is ``"dots"`` (shared pairs as dots), ``"lines"`` (bonds as lines)
    or ``"dot_cross"`` (alternate atoms drawn with dots and crosses).
    """
    style = str(spec.get('style', 'dots')).lower()
    if style not in ('dots', 'lines', 'dot_cross'):
        raise ValueError('style must be "dots", "lines" or "dot_cross"')
    raw_species = spec.get('species') or [spec]
    spacing = float(spec.get('spacing', 22.0 if style != 'lines' else 27.0))
    if spacing < 20:
        raise ValueError('spacing must be at least 20')
    size = float(spec.get('font_size', 10.0))
    c = Canvas(size)
    warnings: list[str] = []
    reports = []
    cursor = 0.0
    for s_idx, raw in enumerate(raw_species, start=1):
        sp = _species(raw, s_idx, style)
        warnings += sp['warnings']
        atoms = sp['atoms']
        min_col = min(a['col'] for a in atoms.values())
        min_row = min(a['row'] for a in atoms.values())
        left_pad = 12.0 if sp['bracket'] else 9.0
        if sp['bare']:
            left_pad = 4.0
        for a in atoms.values():
            a['x'] = cursor + left_pad + (a['col'] - min_col) * spacing
            a['y'] = (a['row'] - min_row) * spacing
            a['half'] = text_width(a['symbol'], size) / 2
        box = [1e9, 1e9, -1e9, -1e9]

        def grow(x, y):
            box[0] = min(box[0], x); box[1] = min(box[1], y); box[2] = max(box[2], x); box[3] = max(box[3], y)

        for a in atoms.values():
            text = a['symbol']
            if sp['bare'] and sp['charge']:
                mag = '' if abs(sp['charge']) == 1 else str(abs(sp['charge']))
                text += '{sup:%s%s}' % (mag, '+' if sp['charge'] > 0 else MINUS)
                c.label((a['x'] - a['half'], a['y']), text, valign='middle')
                grow(a['x'] - a['half'], a['y'] - 5.0); grow(a['x'] + a['half'] + 9.0, a['y'] + 5.0)
            else:
                c.label((a['x'], a['y']), text, align='Center', valign='middle')
                grow(a['x'] - a['half'], a['y'] - 5.0); grow(a['x'] + a['half'], a['y'] + 5.0)
            cross = style == 'dot_cross' and a['order'] % 2 == 1
            for side, count in _groups(a).items():
                dx, dy = SIDES[side]
                reach = (a['half'] + 3.6) if dx else 8.0
                px, py = a['x'] + dx * reach, a['y'] + dy * reach
                offsets = [0.0] if count == 1 else [-PAIR, PAIR]
                for off in offsets:
                    ex, ey = (px, py + off) if dx else (px + off, py)
                    _mark(c, ex, ey, cross)
                    grow(ex - DOT, ey - DOT); grow(ex + DOT, ey + DOT)
        for a, b, order, side in sp['bonds']:
            dx, dy = SIDES[side]
            if dx:
                start, end = a['x'] + dx * a['half'], b['x'] - dx * b['half']
                mx, my = (start + end) / 2, a['y']
            else:
                start, end = a['y'] + dy * 5.0, b['y'] - dy * 5.0
                mx, my = a['x'], (start + end) / 2
            if style == 'lines':
                gap = 2.2
                for k in range(order):
                    off = (k - (order - 1) / 2) * 2.6
                    if dx:
                        c.segment((start + dx * gap, my + off), (end - dx * gap, my + off))
                    else:
                        c.segment((mx + off, start + dy * gap), (mx + off, end - dy * gap))
                continue
            for k in range(order):
                along = (k - (order - 1) / 2) * 2 * PAIR
                for n, off in enumerate((-PAIR, PAIR)):
                    ex, ey = (mx + along, my + off) if dx else (mx + off, my + along)
                    cross = style == 'dot_cross' and (a if n == 0 else b)['order'] % 2 == 1
                    _mark(c, ex, ey, cross)
        if sp['bracket']:
            x1, y1, x2, y2 = box[0] - 5.0, box[1] - 5.0, box[2] + 5.0, box[3] + 5.0
            for x, lip in ((x1, 3.5), (x2, -3.5)):
                c.segment((x, y1), (x, y2))
                c.segment((x, y1), (x + lip, y1))
                c.segment((x, y2), (x + lip, y2))
            mag = '' if abs(sp['charge']) == 1 else str(abs(sp['charge']))
            sign = '+' if sp['charge'] > 0 else MINUS
            if sp['charge']:
                tag = c.label((x2 + 2.0, y1 + 3.0), f'{mag}{sign}', size=size * 0.8)
                box[2] = tag[2]
            else:
                box[2] = x2
        cursor = box[2] + float(spec.get('gap', 12.0))
        reports.append(sp['report'])
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    return _result(c, output, {'diagram': 'lewis_structure', 'style': style, 'species': reports,
                               'warnings': warnings}, margin=10.0)
