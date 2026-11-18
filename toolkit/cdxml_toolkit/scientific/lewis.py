"""Lewis (electron-dot) structures drawn on a square grid (90 degree bonds) or a hexagonal
grid (120 degree bonds).

The spec states the atoms, their grid positions and the bonds; the non-bonding electrons
follow from the valence count (``valence - bonds - formal charge``), so a structure cannot
gain or lose electrons by a drawing slip. The report lists every atom's electron count and
flags atoms that miss or exceed an octet.
"""
from __future__ import annotations

import itertools
import math
from typing import Any, Mapping

from .canvas import BLACK, Canvas, text_width
from .elements import SYMBOLS, atomic_number, group_period, valence_electrons
from .energy import MINUS, _result

# Grid directions as angles in degrees, measured clockwise from "right" (page y points down).
GRIDS = {'square': (0, 90, 180, 270), 'hex': (0, 60, 120, 180, 240, 300)}
_H3 = math.sqrt(3) / 2


def grid_xy(col, row, grid):
    """Page position of a grid point in bond lengths: square, or axial hexagonal coordinates."""
    return (col, row) if grid == 'square' else (col + 0.5 * row, row * _H3)


def _unit(angle):
    return math.cos(math.radians(angle)), math.sin(math.radians(angle))


def _gap(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)
DOT = 0.9          # radius of an electron dot
PAIR = 2.3         # half the distance between the two electrons of a pair


def _mark(c, x, y, cross):
    if cross:
        c.segment((x - 1.3, y - 1.3), (x + 1.3, y + 1.3))
        c.segment((x - 1.3, y + 1.3), (x + 1.3, y - 1.3))
    else:
        c.circle((x, y), DOT, fill=BLACK, outline=None)


def _species(raw, index, style, grid='square'):
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
        (ax, ay), (bx, by) = grid_xy(a['col'], a['row'], grid), grid_xy(b['col'], b['row'], grid)
        if abs(math.hypot(bx - ax, by - ay) - 1.0) > 1e-6:
            raise ValueError(f'species {index}: bonded atoms {a_id} and {b_id} must be neighbours on the '
                             f'{grid} grid (one step apart)')
        side = round(math.degrees(math.atan2(by - ay, bx - ax))) % 360
        if side in a['bonds'] or (side + 180) % 360 in b['bonds']:
            raise ValueError(f'species {index}: {a_id} and {b_id} are bonded twice; use "order" for multiple bonds')
        a['bonds'][side] = (b, order)
        b['bonds'][(side + 180) % 360] = (a, order)
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


def _groups(atom, grid='square'):
    """Non-bonding electrons per free direction: ``{angle: 1 or 2}``.

    The groups take the free grid directions that keep everything around the atom furthest
    apart (largest smallest angle, then the most balanced, then mirror-symmetric about a lone
    bond), so H2O has its pairs above and below and :O=C=O: is symmetric."""
    n = atom['nonbonding']
    if not atom['bonds'] and n < 8:
        # A free atom: one electron per side first, then pair (the element's dot symbol).
        order = [0, 270, 180, 90]
        counts = {s: 0 for s in order}
        for k in range(n):
            counts[order[k % 4]] += 1
        return {s: v for s, v in counts.items() if v}
    pairs, single = divmod(n, 2)
    groups = pairs + single
    free = [d for d in GRIDS[grid] if d not in atom['bonds']] if atom['bonds'] else list(GRIDS['square'])
    if groups > len(free):
        raise ValueError(f'{atom["id"]} needs {groups} places for its non-bonding electrons but only '
                         f'{len(free)} direction(s) are free; give it fewer neighbours on the grid')
    if not groups:
        return {}
    bonds = list(atom['bonds'])

    def score(chosen):
        every = bonds + list(chosen)
        smallest = min(_gap(x, y) for x, y in itertools.combinations(every, 2)) if len(every) > 1 else 360
        vx = sum(_unit(d)[0] for d in every)
        vy = sum(_unit(d)[1] for d in every)
        symmetric = len(bonds) == 1 and {(2 * bonds[0] - d) % 360 for d in chosen} == set(chosen)
        return (smallest, -round(math.hypot(vx, vy), 6), symmetric)

    best = max(itertools.combinations(free, groups), key=score)
    # Pairs go where they are furthest from the bonds; a lone electron takes the remaining place.
    ranked = sorted(best, key=lambda d: -min((_gap(d, b) for b in bonds), default=0))
    out = {d: 2 for d in ranked[:pairs]}
    if single:
        out[ranked[pairs]] = 1
    return out


def _edge(half_w, half_h, ux, uy):
    """Distance from a label's centre to its (elliptical) edge in direction (ux, uy)."""
    return 1.0 / math.sqrt((ux / half_w) ** 2 + (uy / half_h) ** 2)


def lewis_structure(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Lewis structure of one or more species.

    ``{"atoms": [{"id", "element", "at": [col, row], "formal_charge"}], "bonds": [[a, b, order]]}``
    for one species, or ``"species": [...]`` for several drawn left to right (an ionic
    compound). ``style`` is ``"dots"`` (shared pairs as dots), ``"lines"`` (bonds as lines)
    or ``"dot_cross"`` (alternate atoms drawn with dots and crosses). ``grid`` is ``"square"``
    (``at`` is [column, row], bonds at 90 degrees) or ``"hex"`` (``at`` is [column, row] in axial
    hexagonal coordinates, bonds at 60 and 120 degrees: trigonal planar centres, zigzag chains).
    """
    grid = str(spec.get('grid', 'square')).lower()
    if grid not in GRIDS:
        raise ValueError('grid must be "square" or "hex"')
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
        sp = _species(raw, s_idx, style, grid)
        warnings += sp['warnings']
        atoms = sp['atoms']
        places = {aid: grid_xy(a['col'], a['row'], grid) for aid, a in atoms.items()}
        min_x = min(x for x, _y in places.values())
        min_y = min(y for _x, y in places.values())
        left_pad = 12.0 if sp['bracket'] else 9.0
        if sp['bare']:
            left_pad = 4.0
        for aid, a in atoms.items():
            a['x'] = cursor + left_pad + (places[aid][0] - min_x) * spacing
            a['y'] = (places[aid][1] - min_y) * spacing
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
            for side, count in _groups(a, grid).items():
                ux, uy = _unit(side)
                # On the hexagonal grid neighbouring groups can be 60 degrees apart: sit them further out.
                reach = _edge(a['half'] + 3.6, 8.0, ux, uy) + (1.6 if grid == 'hex' else 0.0)
                px, py = a['x'] + ux * reach, a['y'] + uy * reach
                offsets = [0.0] if count == 1 else [-PAIR, PAIR]
                for off in offsets:
                    ex, ey = px - uy * off, py + ux * off
                    _mark(c, ex, ey, cross)
                    grow(ex - DOT, ey - DOT); grow(ex + DOT, ey + DOT)
        for a, b, order, side in sp['bonds']:
            ux, uy = _unit(side)
            ra, rb = _edge(a['half'], 5.0, ux, uy), _edge(b['half'], 5.0, ux, uy)
            start = (a['x'] + ux * ra, a['y'] + uy * ra)
            end = (b['x'] - ux * rb, b['y'] - uy * rb)
            mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
            if style == 'lines':
                gap = 2.2
                for k in range(order):
                    off = (k - (order - 1) / 2) * 2.6
                    c.segment((start[0] + ux * gap - uy * off, start[1] + uy * gap + ux * off),
                              (end[0] - ux * gap - uy * off, end[1] - uy * gap + ux * off))
                continue
            for k in range(order):
                along = (k - (order - 1) / 2) * 2 * PAIR
                for n, off in enumerate((-PAIR, PAIR)):
                    ex, ey = mx + ux * along - uy * off, my + uy * along + ux * off
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
