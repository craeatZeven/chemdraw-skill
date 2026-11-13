"""Orbital diagrams: electron-configuration boxes, molecular-orbital diagrams and
crystal-field splitting.

Electrons are placed by rule (Aufbau order, Pauli limit, Hund's rule) rather than drawn by
hand, and a spec that states occupations is checked against the same rules. The report
returns the configuration, the number of unpaired electrons and, where it is defined, the
bond order or the crystal-field stabilisation energy.
"""
from __future__ import annotations

import math
import re
from typing import Any, Mapping

from .canvas import Canvas, plain
from .common import number
from .elements import (CAPACITY, L_LETTERS, MADELUNG, SYMBOLS, atomic_number, d_electron_count,
                       electron_configuration, noble_gas_core, parse_ion, subshell_name)
from .energy import MINUS, _result, _spread

BAR = 20.0          # width of one orbital level
BAR_GAP = 6.0       # gap between degenerate levels
BOX = 15.0          # side of one orbital box


def _charge_markup(charge):
    if not charge:
        return ''
    size = '' if abs(charge) == 1 else str(abs(charge))
    return '{sup:%s%s}' % (size, '+' if charge > 0 else MINUS)


def hund_fill(electrons, orbitals):
    """Occupation of ``orbitals`` degenerate orbitals: singly with parallel spins first, then paired.

    Returns a list of ``'' | 'u' | 'ud'`` per orbital.
    """
    if not 0 <= electrons <= 2 * orbitals:
        raise ValueError(f'{electrons} electrons do not fit {orbitals} orbital(s)')
    singles = min(electrons, orbitals)
    pairs = electrons - singles
    return ['ud' if i < pairs else 'u' if i < singles else '' for i in range(orbitals)]


def _check_boxes(boxes, name, allow):
    """Pauli and Hund checks of an explicit occupation; returns the problems found."""
    problems = []
    for box in boxes:
        if box not in ('', 'u', 'd', 'ud', 'du', 'uu', 'dd'):
            raise ValueError(f'{name}: an orbital holds "", "u", "d" or "ud", not {box!r}')
        if box in ('uu', 'dd'):
            problems.append(f'{name}: two electrons with the same spin in one orbital violate the Pauli principle')
    singles = [b for b in boxes if len(b) == 1]
    if any(len(b) == 2 for b in boxes) and any(b == '' for b in boxes):
        problems.append(f'{name}: an orbital is doubly occupied while another is empty (Hund\'s rule)')
    if len(set(singles)) > 1:
        problems.append(f'{name}: the unpaired electrons do not have parallel spins (Hund\'s rule)')
    if problems and not allow:
        raise ValueError('; '.join(problems) + '. Set "allow_violations": true to draw it anyway')
    return problems


def _spin_arrows(c, x, y, occupation, top=-9.0, bottom=5.0):
    """Draw the electrons of one orbital centred at ``x`` on the level ``y``."""
    slots = [0.0] if len(occupation) == 1 else [-3.0, 3.0]
    for dx, ch in zip(slots, occupation):
        upper, lower = (x + dx, y + top), (x + dx, y + bottom)
        if ch == 'u':
            c.arrow(lower, upper, head_style='HalfLeft', head_size=600)
        else:
            c.arrow(upper, lower, head_style='HalfLeft', head_size=600)


def _unpaired(occupations):
    return sum(1 for box in occupations if len(box) == 1)


# ---------------------------------------------------------------------------- configuration boxes
def _config_text(config, core=None):
    body = ' '.join(f'{name}{count}' for name, count in config.items())
    return f'[{core}] {body}'.strip() if core else body


def electron_configuration_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Orbital-box diagram of an electron configuration.

    ``{"element": "Fe", "charge": 2, "core": "auto"|"full"|"noble"|"valence"}`` fills the
    boxes from the ground-state configuration. ``"subshells": [{"label": "2p", "electrons": 4}]``
    or ``{"label": "2p", "boxes": ["ud", "u", "u"]}`` draws a stated occupation, which must
    obey the Pauli principle and Hund's rule unless ``"allow_violations"`` is true.
    ``"layout"`` is ``"row"`` or ``"energy"`` (subshells stacked by energy).
    """
    allow = bool(spec.get('allow_violations', False))
    warnings: list[str] = []
    core_symbol = None
    report: dict[str, Any] = {'diagram': 'electron_configuration'}
    groups = []
    if spec.get('subshells') is not None:
        for i, item in enumerate(spec['subshells']):
            label = str(item.get('label', ''))
            match = re.fullmatch(r'(\d)([spdf])', label)
            if 'boxes' in item:
                boxes = [str(b) for b in item['boxes']]
                if match and len(boxes) != CAPACITY[match.group(2)] // 2:
                    raise ValueError(f'{label}: a {match.group(2)} subshell has {CAPACITY[match.group(2)] // 2} orbitals')
                warnings += _check_boxes(boxes, label or f'subshell {i + 1}', allow)
                boxes = ['ud' if b == 'du' else b for b in boxes]
            else:
                orbitals = int(item.get('orbitals', CAPACITY[match.group(2)] // 2 if match else 1))
                boxes = hund_fill(int(item.get('electrons', 0)), orbitals)
            groups.append({'label': label, 'boxes': boxes})
        if not groups:
            raise ValueError('"subshells" is empty')
        report['configuration'] = ' '.join(f'{g["label"]}{sum(len(b) for b in g["boxes"])}' for g in groups)
    else:
        if spec.get('ion') is not None:
            symbol, charge = parse_ion(spec['ion'])
        elif spec.get('element') is not None:
            symbol, charge = SYMBOLS[atomic_number(spec['element']) - 1], int(spec.get('charge', 0))
        else:
            raise ValueError('Give "element" (with an optional "charge"), "ion" or explicit "subshells"')
        config = electron_configuration(symbol, charge)
        core_name, core_shells = noble_gas_core(symbol, charge)
        mode = str(spec.get('core', 'auto')).lower()
        if mode not in ('auto', 'full', 'noble', 'valence'):
            raise ValueError('"core" must be auto, full, noble or valence')
        electrons = sum(config.values())
        if mode == 'auto':
            mode = 'full' if electrons <= 18 else 'noble'
        shown = dict(config)
        if mode in ('noble', 'valence') and core_name and all(config.get(k) == v for k, v in core_shells.items()):
            shown = {k: v for k, v in config.items() if k not in core_shells}
            core_symbol = core_name if mode == 'noble' else None
            if mode == 'valence':
                # A filled inner d or f subshell is not part of the valence shell of a p-block atom.
                outer = max(int(name[0]) for name in shown)
                if any(name[1] == 'p' and int(name[0]) == outer for name in shown):
                    shown = {k: v for k, v in shown.items() if int(k[0]) == outer}
        for name, count in shown.items():
            groups.append({'label': name, 'boxes': hund_fill(count, CAPACITY[name[1]] // 2)})
        if not groups:
            raise ValueError('There are no electrons left to draw outside the noble-gas core')
        report.update({'species': symbol + (f'{abs(charge) if abs(charge) > 1 else ""}{"+" if charge > 0 else "-"}'
                                            if charge else ''),
                       'electrons': electrons, 'configuration': _config_text(config),
                       'abbreviated': _config_text({k: v for k, v in config.items() if k not in core_shells}
                                                   if core_name else config, core_name),
                       'shown': mode})
    unpaired = sum(_unpaired(g['boxes']) for g in groups)
    report['unpaired_electrons'] = unpaired

    layout = str(spec.get('layout', 'row')).lower()
    if layout not in ('row', 'energy'):
        raise ValueError('"layout" must be "row" or "energy"')
    c = Canvas()

    def draw_group(x, y, group):
        for k, occupation in enumerate(group['boxes']):
            c.rect(x + k * BOX, y - BOX / 2, x + (k + 1) * BOX, y + BOX / 2)
            if occupation:
                _spin_arrows(c, x + (k + 0.5) * BOX, y, occupation, top=-5.5, bottom=5.5)
        return x + len(group['boxes']) * BOX

    if layout == 'row':
        x = 0.0
        if core_symbol:
            box = c.label((0.0, 0.0), f'[{core_symbol}]', valign='middle')
            x = box[2] + 8.0
        for group in groups:
            end = draw_group(x, 0.0, group)
            if group['label']:
                c.label(((x + end) / 2, BOX / 2 + 12.0), group['label'], align='Center')
            x = end + 10.0
    else:
        order = {subshell_name(n, l): i for i, (n, l) in enumerate(MADELUNG)}
        ranked = sorted(groups, key=lambda g: order.get(g['label'], len(order)))
        for level, group in enumerate(ranked):
            l_index = L_LETTERS.index(group['label'][1]) if re.fullmatch(r'\d[spdf]', group['label']) else 0
            x, y = 34.0 + l_index * 52.0, -level * 24.0
            draw_group(x, y, group)
            if group['label']:
                c.label((x - 5.0, y), group['label'], align='Right', valign='middle')
        top = -(len(ranked) - 1) * 24.0 - 20.0
        c.arrow((0.0, 12.0), (0.0, top))
        c.label((-4.0, top - 6.0), spec.get('y_label', 'Energy'))
        if core_symbol:
            c.label((34.0, 24.0), f'[{core_symbol}] core', align='Left')
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report['warnings'] = warnings
    return _result(c, output, report)


# ---------------------------------------------------------------------------- MO diagram
def _diatomic_preset(name, charge):
    match = re.fullmatch(r'([A-Z][a-z]?)2', str(name).strip())
    if not match:
        raise ValueError('An MO preset is a homonuclear diatomic of period 1 or 2, such as "O2" or "N2"')
    symbol = match.group(1)
    z = atomic_number(symbol)
    if z > 10:
        raise ValueError('MO presets cover H2 to Ne2; give "left", "right" and "center" for other molecules')
    charge = int(charge)

    def side(valence, title):
        if z <= 2:
            return {'title': title, 'levels': [{'id': '1s', 'label': '1s', 'energy': 0.0, 'electrons': valence}]}
        levels = [{'id': '2s', 'label': '2s', 'energy': -4.0, 'electrons': min(2, valence)}]
        if z >= 5:
            levels.append({'id': '2p', 'label': '2p', 'energy': 0.0, 'degeneracy': 3, 'electrons': max(0, valence - 2)})
        return {'title': title, 'levels': levels}

    valence = z if z <= 2 else z - 2
    if z <= 2:
        center = [{'label': 'σ{sub:1s}', 'energy': -1.2, 'character': 'bonding', 'from': ['1s']},
                  {'label': 'σ{sub:1s}{sup:*}', 'energy': 1.2, 'character': 'antibonding', 'from': ['1s']}]
    else:
        center = [{'label': 'σ{sub:2s}', 'energy': -5.0, 'character': 'bonding', 'from': ['2s']},
                  {'label': 'σ{sub:2s}{sup:*}', 'energy': -3.0, 'character': 'antibonding', 'from': ['2s']}]
        if z >= 5:
            sigma = {'label': 'σ{sub:2p}', 'character': 'bonding', 'from': ['2p']}
            pi = {'label': 'π{sub:2p}', 'degeneracy': 2, 'character': 'bonding', 'from': ['2p']}
            # s-p mixing puts sigma(2p) above pi(2p) up to N2; from O2 on it lies below.
            if z <= 7:
                center += [{**pi, 'energy': -1.4}, {**sigma, 'energy': -0.4}]
            else:
                center += [{**sigma, 'energy': -1.9}, {**pi, 'energy': -0.9}]
            center += [{'label': 'π{sub:2p}{sup:*}', 'energy': 1.0, 'degeneracy': 2, 'character': 'antibonding',
                        'from': ['2p']},
                       {'label': 'σ{sub:2p}{sup:*}', 'energy': 2.1, 'character': 'antibonding', 'from': ['2p']}]
    right_valence = valence - charge
    capacity = 2 if z <= 2 else (2 if z < 5 else 8)
    if not 0 <= right_valence <= capacity:
        raise ValueError(f'{symbol}2 cannot carry a charge of {charge:+d} in this valence diagram')
    return {'left': side(valence, symbol), 'right': side(right_valence, symbol + _charge_markup(charge)),
            'center': {'title': symbol + '{sub:2}' + _charge_markup(charge), 'levels': center}, 'fill': True,
            'schematic': True}


def mo_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Molecular-orbital energy-level diagram: atomic orbitals left and right, MOs between.

    ``{"preset": "O2", "charge": 0}`` draws a period-1/2 homonuclear diatomic with the
    correct sigma/pi order. Otherwise give ``left``, ``right`` and ``center``, each
    ``{"title", "levels": [{"id", "label", "energy", "degeneracy", "electrons"}]}``;
    centre levels add ``"from": [atomic level ids]`` and ``"character"``
    (bonding, antibonding or nonbonding). With ``"fill": true`` the MOs are filled from
    the atomic electron count; stated MO electrons must add up to it.
    """
    if spec.get('preset') is not None:
        preset = _diatomic_preset(spec['preset'], spec.get('charge', 0))
        spec = {**preset, **{k: v for k, v in spec.items() if k not in ('preset', 'charge')}}
    warnings: list[str] = []
    columns = {}
    for key in ('left', 'center', 'right'):
        raw = spec.get(key)
        if raw is None and key != 'center':
            columns[key] = {'title': '', 'levels': []}
            continue
        if not isinstance(raw, Mapping) or not raw.get('levels'):
            raise ValueError(f'"{key}" needs a "levels" list')
        levels = []
        for i, item in enumerate(raw['levels']):
            degeneracy = int(item.get('degeneracy', 1))
            if not 1 <= degeneracy <= 5:
                raise ValueError(f'{key} level {i + 1}: degeneracy must be 1 to 5')
            electrons = item.get('electrons')
            level = {'id': str(item.get('id') or plain(item.get('label', '')) or f'{key}{i + 1}'),
                     'label': str(item.get('label', '')), 'degeneracy': degeneracy,
                     'energy': number(item.get('energy'), f'{key} level {i + 1} energy', -1e6, 1e6),
                     'electrons': None if electrons is None else int(electrons),
                     'from': [str(v) for v in item.get('from', [])], 'character': item.get('character'),
                     'column': key}
            if level['electrons'] is not None and not 0 <= level['electrons'] <= 2 * degeneracy:
                raise ValueError(f'{key} level {level["id"]}: {level["electrons"]} electrons exceed the Pauli limit '
                                 f'of {2 * degeneracy}')
            if level['character'] not in (None, 'bonding', 'antibonding', 'nonbonding'):
                raise ValueError(f'{key} level {level["id"]}: character must be bonding, antibonding or nonbonding')
            levels.append(level)
        columns[key] = {'title': str(raw.get('title', '')), 'levels': levels}
    atomic = columns['left']['levels'] + columns['right']['levels']
    molecular = columns['center']['levels']
    total = sum(lv['electrons'] or 0 for lv in atomic)
    if spec.get('fill') or all(lv['electrons'] is None for lv in molecular):
        if not atomic:
            raise ValueError('Give the electrons of each MO, or atomic levels to fill them from')
        remaining = total
        for lv in sorted(molecular, key=lambda item: item['energy']):
            lv['electrons'] = min(remaining, 2 * lv['degeneracy'])
            remaining -= lv['electrons']
        if remaining:
            raise ValueError(f'{remaining} electron(s) do not fit the molecular orbitals given')
    else:
        for lv in molecular:
            lv['electrons'] = lv['electrons'] or 0
        placed = sum(lv['electrons'] for lv in molecular)
        if atomic and placed != total and not spec.get('allow_mismatch'):
            raise ValueError(f'The atomic orbitals hold {total} electrons but the molecular orbitals hold {placed}; '
                             'electrons are conserved (set "allow_mismatch": true for a charged fragment)')
        ordered = sorted(molecular, key=lambda item: item['energy'])
        for lower, upper in zip(ordered, ordered[1:]):
            if upper['electrons'] and lower['electrons'] < 2 * lower['degeneracy'] and upper['energy'] > lower['energy']:
                warnings.append(f'{upper["id"]} is occupied while the lower {lower["id"]} is not full '
                                '(not the ground state)')
    for lv in atomic:
        lv['electrons'] = lv['electrons'] or 0

    # ---- geometry
    every = atomic + molecular
    height = float(spec.get('height', max(150.0, 30.0 * (len({round(lv["energy"], 6) for lv in every}) - 1))))
    if _spread(every, 'energy', height, 28.0) and not spec.get('schematic'):
        warnings.append('close levels were moved apart: the vertical spacing is not to scale')
    gap = float(spec.get('column_gap', 118.0))
    centers = {'left': 0.0, 'center': gap, 'right': 2 * gap}
    c = Canvas()
    for lv in every:
        width = lv['degeneracy'] * BAR + (lv['degeneracy'] - 1) * BAR_GAP
        lv['x1'] = centers[lv['column']] - width / 2
        lv['x2'] = lv['x1'] + width
    index = {}
    for lv in atomic:
        index.setdefault(lv['id'], []).append(lv)
    for lv in molecular:
        for ref in lv['from']:
            if ref not in index:
                raise ValueError(f'MO {lv["id"]}: "from" names an unknown atomic level {ref!r}')
            for source in index[ref]:
                if source['column'] == 'left':
                    c.segment((source['x2'] + 2.0, source['y']), (lv['x1'] - 2.0, lv['y']), dashed=True)
                else:
                    c.segment((source['x1'] - 2.0, source['y']), (lv['x2'] + 2.0, lv['y']), dashed=True)
    for lv in every:
        occupation = hund_fill(lv['electrons'], lv['degeneracy'])
        lv['occupation'] = occupation
        for k in range(lv['degeneracy']):
            x = lv['x1'] + k * (BAR + BAR_GAP)
            c.segment((x, lv['y']), (x + BAR, lv['y']))
            if occupation[k]:
                _spin_arrows(c, x + BAR / 2, lv['y'], occupation[k])
        if not lv['label']:
            continue
        if lv['column'] == 'left':
            c.label((lv['x1'] - 5.0, lv['y']), lv['label'], align='Right', valign='middle')
        elif lv['column'] == 'right':
            c.label((lv['x2'] + 5.0, lv['y']), lv['label'], align='Left', valign='middle')
        else:
            c.label(((lv['x1'] + lv['x2']) / 2, lv['y'] + 12.0), lv['label'], align='Center')
    bottom = max(lv['y'] for lv in every) + 32.0
    for key in ('left', 'center', 'right'):
        if columns[key]['title']:
            c.label((centers[key], bottom), columns[key]['title'], align='Center')
    if spec.get('axis', True):
        left = c.extent[0] - 12.0
        top = min(lv['y'] for lv in every) - 16.0
        c.arrow((left, max(lv['y'] for lv in every) + 10.0), (left, top))
        c.label((left - 4.0, top - 6.0), spec.get('y_label', 'Energy'))

    report: dict[str, Any] = {'diagram': 'mo_diagram', 'electrons': sum(lv['electrons'] for lv in molecular),
                              'occupation': [{'id': lv['id'], 'electrons': lv['electrons'],
                                              'character': lv['character']}
                                             for lv in sorted(molecular, key=lambda item: item['energy'])]}
    unpaired = sum(_unpaired(lv['occupation']) for lv in molecular)
    report['unpaired_electrons'] = unpaired
    report['magnetism'] = 'paramagnetic' if unpaired else 'diamagnetic'
    if all(lv['character'] for lv in molecular):
        bonding = sum(lv['electrons'] for lv in molecular if lv['character'] == 'bonding')
        anti = sum(lv['electrons'] for lv in molecular if lv['character'] == 'antibonding')
        report['bond_order'] = (bonding - anti) / 2
        if report['bond_order'] <= 0:
            warnings.append('the bond order is zero or negative: the species is not bound')
    if spec.get('schematic'):
        report['energies'] = 'schematic order, not to scale'
    report['warnings'] = warnings
    return _result(c, output, report)


# ---------------------------------------------------------------------------- crystal field
# geometry -> sets from lowest to highest energy: (name, orbital labels, energy in units of the splitting)
_FIELDS = {
    'octahedral': {'symbol': 'Δ{sub:o}', 'default_spin': 'high', 'sets': [
        ('t{sub:2g}', ['d{sub:xy}', 'd{sub:xz}', 'd{sub:yz}'], -0.4),
        ('e{sub:g}', ['d{sub:z²}', 'd{sub:x²−y²}'], 0.6)]},
    'tetrahedral': {'symbol': 'Δ{sub:t}', 'default_spin': 'high', 'sets': [
        ('e', ['d{sub:z²}', 'd{sub:x²−y²}'], -0.6),
        ('t{sub:2}', ['d{sub:xy}', 'd{sub:xz}', 'd{sub:yz}'], 0.4)]},
    'square_planar': {'symbol': None, 'default_spin': 'low', 'sets': [
        ('', ['d{sub:xz}', 'd{sub:yz}'], -0.51),
        ('', ['d{sub:z²}'], -0.43),
        ('', ['d{sub:xy}'], 0.23),
        ('', ['d{sub:x²−y²}'], 1.23)]},
}


def crystal_field_fill(geometry, d_electrons, spin=None):
    """Occupation of each orbital set: ``[['ud', 'u', 'u'], ['u', 'u']]`` from the lowest set up."""
    field = _FIELDS[geometry]
    spin = spin or field['default_spin']
    sizes = [len(orbitals) for _name, orbitals, _e in field['sets']]
    if not 0 <= d_electrons <= 10:
        raise ValueError('A d subshell holds 0 to 10 electrons')
    counts = [0] * len(sizes)
    remaining = d_electrons
    if spin == 'high':
        for i, size in enumerate(sizes):       # one electron per orbital, lowest set first
            take = min(size, remaining); counts[i] += take; remaining -= take
        for i, size in enumerate(sizes):       # then pair up, lowest set first
            take = min(size, remaining); counts[i] += take; remaining -= take
    elif spin == 'low':
        for i, size in enumerate(sizes):       # fill each set completely before the next
            take = min(2 * size, remaining); counts[i] += take; remaining -= take
    else:
        raise ValueError('spin must be "high" or "low"')
    return [hund_fill(count, size) for count, size in zip(counts, sizes)]


def crystal_field(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Crystal-field splitting diagram of the d orbitals with the electrons filled in.

    ``{"geometry": "octahedral"|"tetrahedral"|"square_planar", "ion": "Fe3+" or
    "d_electrons": 5, "spin": "high"|"low"}``. The report gives the configuration, the
    unpaired electrons, the spin-only magnetic moment and the stabilisation energy.
    """
    geometry = str(spec.get('geometry', 'octahedral')).lower().replace(' ', '_').replace('-', '_')
    if geometry not in _FIELDS:
        raise ValueError(f'geometry must be one of: {", ".join(_FIELDS)}')
    field = _FIELDS[geometry]
    warnings: list[str] = []
    species = None
    if spec.get('d_electrons') is not None:
        d = int(spec['d_electrons'])
    elif spec.get('ion') is not None:
        symbol, charge = parse_ion(spec['ion'])
        d = d_electron_count(symbol, charge)
        species = f'{symbol}{_charge_markup(charge)}'
    else:
        raise ValueError('Give "d_electrons" or "ion" (for example "Fe3+")')
    spin = spec.get('spin')
    if spin is not None:
        spin = str(spin).lower()
    occupations = crystal_field_fill(geometry, d, spin)
    used_spin = spin or field['default_spin']
    other = crystal_field_fill(geometry, d, 'low' if used_spin == 'high' else 'high')
    spin_matters = other != occupations
    if geometry == 'tetrahedral' and used_spin == 'low' and spin_matters:
        warnings.append('tetrahedral complexes are almost always high spin (the splitting is small)')
    if geometry == 'square_planar' and used_spin == 'high' and spin_matters:
        warnings.append('square-planar complexes are almost always low spin (the splitting is large)')
    if spin is None and spin_matters and geometry == 'octahedral':
        warnings.append(f'd{d} octahedral complexes exist in high- and low-spin forms; high spin was drawn, '
                        'set "spin" to choose')

    scale = float(spec.get('splitting', 64.0))
    sets = []
    for (name, orbitals, energy), occupation in zip(field['sets'], occupations):
        sets.append({'name': name, 'orbitals': orbitals, 'energy': energy, 'occupation': occupation})
    _spread(sets, 'energy', (sets[-1]['energy'] - sets[0]['energy']) * scale, 28.0)
    show_orbitals = bool(spec.get('orbital_labels', geometry == 'square_planar'))
    if show_orbitals:   # room for the orbital names under each level
        for k, s in enumerate(sets):
            s['y'] += -6.0 * k
    c = Canvas()
    right = float(spec.get('column_gap', 96.0))
    lowest, highest = max(s['y'] for s in sets), min(s['y'] for s in sets)
    barycentre = lowest + (0.0 - sets[0]['energy']) / (sets[-1]['energy'] - sets[0]['energy']) * (highest - lowest)
    free_width = 5 * BAR + 4 * BAR_GAP
    if spec.get('free_ion', True):
        for k in range(5):
            x = -free_width / 2 + k * (BAR + BAR_GAP)
            c.segment((x, barycentre), (x + BAR, barycentre))
        c.label((0.0, barycentre + 13.0), spec.get('free_ion_label', 'free ion'), align='Center')
    edge = 0.0
    for s in sets:
        n = len(s['orbitals'])
        width = n * BAR + (n - 1) * BAR_GAP
        x1 = right + free_width / 2 - width / 2 + 30.0
        s['x1'], s['x2'] = x1, x1 + width
        edge = max(edge, s['x2'])
        if spec.get('free_ion', True):
            c.segment((free_width / 2 + 3.0, barycentre), (x1 - 3.0, s['y']), dashed=True)
        for k in range(n):
            x = x1 + k * (BAR + BAR_GAP)
            c.segment((x, s['y']), (x + BAR, s['y']))
            if s['occupation'][k]:
                _spin_arrows(c, x + BAR / 2, s['y'], s['occupation'][k])
            if show_orbitals:
                c.label((x + BAR / 2, s['y'] + 12.0), s['orbitals'][k], align='Center', size=8)
    for s in sets:
        if s['name']:
            c.label((edge + 6.0, s['y']), s['name'], align='Left', valign='middle')
    if field['symbol'] and spec.get('splitting_arrow', True):
        ax = edge + 30.0
        c.segment((edge + 22.0, sets[0]['y']), (ax + 4.0, sets[0]['y']), dashed=True)
        c.segment((edge + 22.0, sets[1]['y']), (ax + 4.0, sets[1]['y']), dashed=True)
        c.arrow((ax, sets[0]['y'] - 1.2), (ax, sets[1]['y'] + 1.2), tail_style='Full', head_size=800)
        c.label((ax + 5.0, (sets[0]['y'] + sets[1]['y']) / 2), spec.get('splitting_label', field['symbol']),
                valign='middle')
    title = spec.get('title')
    if title is None and species:
        title = f'{species}, d{{sup:{d}}}' + (f', {used_spin} spin' if spin_matters else '')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    if spec.get('axis', True):
        left = c.extent[0] - 12.0
        c.arrow((left, lowest + 10.0), (left, highest - 18.0))
        c.label((left - 4.0, highest - 24.0), spec.get('y_label', 'Energy'))

    counts = [sum(len(b) for b in s['occupation']) for s in sets]
    unpaired = sum(_unpaired(s['occupation']) for s in sets)
    report: dict[str, Any] = {
        'diagram': 'crystal_field', 'geometry': geometry, 'd_electrons': d, 'spin': used_spin,
        'spin_matters': spin_matters, 'unpaired_electrons': unpaired,
        'spin_only_moment_bm': round(math.sqrt(unpaired * (unpaired + 2)), 2),
        'configuration': ' '.join(f'{plain(s["name"]) or "(" + ",".join(plain(o) for o in s["orbitals"]) + ")"}^{n}'
                                  for s, n in zip(sets, counts)),
    }
    if geometry in ('octahedral', 'tetrahedral'):
        stabilisation = sum(s['energy'] * n for s, n in zip(sets, counts))
        pairs = sum(1 for s in sets for b in s['occupation'] if len(b) == 2)
        extra_pairs = pairs - max(0, d - 5)
        unit = 'Delta_o' if geometry == 'octahedral' else 'Delta_t'
        report['cfse'] = {'splitting_units': round(stabilisation, 2), 'unit': unit, 'extra_pairing_energies': extra_pairs}
    report['warnings'] = warnings
    return _result(c, output, report)
