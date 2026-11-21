"""CDXML geometry, occupation and chemistry contracts for frontier-orbital interaction diagrams."""
from __future__ import annotations

from copy import deepcopy
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.scientific.diagrams import DIAGRAM_OPERATIONS, build
from cdxml_toolkit.scientific.fmo import fmo_interaction

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples' / 'fmo'
ORBITAL_PALETTE = ['#000000', '#FFFFFF', '#FF6600', '#336699', '#808080', '#C0392B', '#0072B2']


def _root(result):
    return ET.fromstring(result['cdxml'])


def _solid_lines(root):
    return [node for node in root.iter('graphic')
            if node.get('GraphicType') == 'Line' and node.get('LineType') is None]


def _dashed_lines(root):
    return [node for node in root.iter('graphic')
            if node.get('GraphicType') == 'Line' and node.get('LineType') == 'Dashed']


def _half_arrows(root):
    return [node for node in root.iter('arrow') if node.get('ArrowheadHead') in ('HalfLeft', 'HalfRight')]


def _node(root, node_id):
    return next(node for node in root.iter() if node.get('id') == node_id)


def _coords(node, attr):
    return tuple(float(value) for value in node.get(attr).split())


def _line_y(root, node_id):
    return _coords(_node(root, node_id), 'BoundingBox')[1]


def _assert_occupations(root, levels, expected):
    spins = {level_id: [] for level_id in levels}
    for node in _half_arrows(root):
        tail, head = _coords(node, 'Tail3D'), _coords(node, 'Head3D')
        assert tail[0] == pytest.approx(head[0], abs=0.5)
        owners = []
        for level_id, line in levels.items():
            x1, y1, x2, y2 = _coords(line, 'BoundingBox')
            assert y1 == pytest.approx(y2, abs=0.5)
            if x1 <= tail[0] <= x2 and min(tail[1], head[1]) < y1 < max(tail[1], head[1]):
                owners.append(level_id)
        assert len(owners) == 1, f'electron {node.get("id")} must cross exactly one level'
        spins[owners[0]].append('u' if head[1] < tail[1] else 'd')
    for level_id, occupation in expected.items():
        assert sorted(spins[level_id]) == sorted(occupation), f'occupation at {level_id}: {spins[level_id]}'


def _assert_links(root, result, levels):
    expected = {(source['id'], mixed['id'])
                for source in _ONE_DONOR_ONE_ACCEPTOR['donor'] + _ONE_DONOR_ONE_ACCEPTOR['acceptor']
                for mixed in _ONE_DONOR_ONE_ACCEPTOR['mixing']}
    assert len(result['links']) == len(_dashed_lines(root)) == len(expected)
    assert {(link['from'], link['to']) for link in result['links']} == expected
    assert {link['node_id'] for link in result['links']} == {node.get('id') for node in _dashed_lines(root)}
    for link in result['links']:
        node = _node(root, link['node_id'])
        source, mixed = levels[link['from']], levels[link['to']]
        sx1, sy, sx2, _ = _coords(source, 'BoundingBox')
        mx1, my, mx2, _ = _coords(mixed, 'BoundingBox')
        tail = (sx2 if sx2 < mx1 else sx1, sy, 0.0)
        head = (mx1 if sx2 < mx1 else mx2, my, 0.0)
        assert _coords(node, 'Tail3D') == pytest.approx(tail, abs=0.5), f'link tail {link["from"]}->{link["to"]}'
        assert _coords(node, 'Head3D') == pytest.approx(head, abs=0.5), f'link head {link["from"]}->{link["to"]}'
        assert _coords(node, 'BoundingBox') == pytest.approx((*tail[:2], *head[:2]), abs=0.5)


def _colors(root):
    table = {0: '#000000', 1: '#FFFFFF'}
    for index, color in enumerate(root.findall('./colortable/color'), start=2):
        rgb = [round(float(color.get(channel)) * 255) for channel in 'rgb']
        table[index] = '#' + ''.join(f'{value:02X}' for value in rgb)
    return table


_ONE_DONOR_ONE_ACCEPTOR = {
    'donor': [{'id': 'nCl', 'label': 'n(Cl)', 'energy': 0.0, 'electrons': 2}],
    'acceptor': [{'id': 'pi*', 'label': '\u03c0*(C=C)', 'energy': 6.0, 'electrons': 0}],
    'mixing': [
        {'id': 'psi-', 'label': '\u03c8{sub:\u2212}', 'energy': -2.0, 'electrons': 2, 'from': ['nCl', 'pi*']},
        {'id': 'psi+', 'label': '\u03c8{sub:+}', 'energy': 8.0, 'electrons': 0, 'from': ['nCl', 'pi*']},
    ],
}


def test_one_donor_one_acceptor_straddles_and_carries_two_electrons():
    result = fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR)
    assert result['diagram'] == 'fmo_interaction'
    assert result['straddles_donor'] is True
    generated = result['generated_levels']
    assert [level['electrons'] for level in generated] == [2, 0]
    assert [level['occupation'] for level in generated] == ['ud', '']
    assert generated[0]['energy'] < result['donor_homo']['energy'] < generated[1]['energy']
    donor = next(level for level in result['occupation'] if level['id'] == 'nCl')
    assert donor['electrons'] == 2 and donor['occupation'] == 'ud'
    assert result['gap'] == pytest.approx(6.0) and result['stabilisation'] == pytest.approx(2.0)

    root = _root(result)
    assert len(_solid_lines(root)) == 4
    levels = {level['id']: _node(root, level['line_node_id']) for level in result['occupation']}
    y = {level_id: _coords(line, 'BoundingBox')[1] for level_id, line in levels.items()}
    assert y['psi-'] > y['nCl'] > y['psi+']
    assert y['nCl'] > y['pi*']
    _assert_occupations(root, levels, {'nCl': 'ud', 'pi*': '', 'psi-': 'ud', 'psi+': ''})
    _assert_links(root, result, levels)
    assert check_cdxml(result['cdxml'])['ok']


def _translate(node, dx=0.0, dy=0.0):
    for attr, stride in (('BoundingBox', 2), ('Tail3D', 3), ('Head3D', 3)):
        if node.get(attr) is None:
            continue
        values = [float(value) for value in node.get(attr).split()]
        for index in range(0, len(values), stride):
            values[index] += dx
            values[index + 1] += dy
        node.set(attr, ' '.join(map(str, values)))


def _move_filled_electrons(result):
    root = _root(result)
    levels = {level['id']: level for level in result['occupation']}
    low_y = _line_y(root, levels['psi-']['line_node_id'])
    high_y = _line_y(root, levels['psi+']['line_node_id'])
    moved = 0
    for node in _half_arrows(root):
        tail_y = float(node.get('Tail3D').split()[1])
        head_y = float(node.get('Head3D').split()[1])
        if min(tail_y, head_y) < low_y < max(tail_y, head_y):
            _translate(node, dy=high_y - low_y)
            moved += 1
    assert moved == 2
    return {**result, 'cdxml': ET.tostring(root, encoding='unicode')}


def _shift_correlation_lines(result):
    root = _root(result)
    lines = _dashed_lines(root)
    assert len(lines) == 4
    for node in lines:
        _translate(node, dx=8.0, dy=8.0)
    return {**result, 'cdxml': ET.tostring(root, encoding='unicode')}


def test_acceptance_rejects_electrons_moved_to_empty_level(monkeypatch):
    mutated = _move_filled_electrons(fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR))
    monkeypatch.setitem(globals(), 'fmo_interaction', lambda spec: mutated)
    with pytest.raises(AssertionError, match='occupation at psi-') as failure:
        test_one_donor_one_acceptor_straddles_and_carries_two_electrons()
    print(f'Moved-electron CDXML: FAIL ({failure.value})')


def test_acceptance_rejects_correlation_lines_shifted_eight_points(monkeypatch):
    mutated = _shift_correlation_lines(fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR))
    monkeypatch.setitem(globals(), 'fmo_interaction', lambda spec: mutated)
    with pytest.raises(AssertionError, match='link tail') as failure:
        test_one_donor_one_acceptor_straddles_and_carries_two_electrons()
    print(f'Shifted-link CDXML: FAIL ({failure.value})')


def test_generated_levels_must_straddle_the_donor():
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['mixing'][0]['energy'] = 1
    with pytest.raises(ValueError, match='straddle'):
        fmo_interaction(spec)


def test_electrons_are_conserved_across_the_columns():
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['mixing'][0]['electrons'] = 1
    with pytest.raises(ValueError, match='conserved'):
        fmo_interaction(spec)


def test_mixing_from_names_an_unknown_level():
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['mixing'][0]['from'] = ['nCl', 'unknown']
    with pytest.raises(ValueError, match='unknown level'):
        fmo_interaction(spec)


@pytest.mark.parametrize('ci,cj,product', [(0.5, 0.4, 0.20), (0.6, 0.3, 0.18)])
def test_coefficient_lobes_reuse_the_s2_colors_and_scale_with_magnitude(ci, cj, product):
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['donor'][0]['coefficient'] = ci
    spec['acceptor'][0].update(coefficient=cj, phase=-1)
    result = fmo_interaction(spec)
    root = _root(result)
    ovals = [node for node in root.iter('graphic') if node.get('GraphicType') == 'Oval']
    assert len(ovals) == 2
    colors = _colors(root)
    assert {colors[int(node.get('color'))] for node in ovals} == {'#C0392B', '#0072B2'}
    lengths = []
    for node in ovals:
        center = _coords(node, 'Center3D')[:2]
        major = _coords(node, 'MajorAxisEnd3D')[:2]
        lengths.append(2 * ((major[0] - center[0]) ** 2 + (major[1] - center[1]) ** 2) ** 0.5)
    assert lengths == pytest.approx([26.0 * ci, 26.0 * cj], abs=0.05)
    pair, = result['coefficient_pairs']
    assert (pair['donor'], pair['acceptor']) == ('nCl', 'pi*')
    assert pair['magnitude_product'] == pytest.approx(product)
    label = _node(root, pair['label_node_id'])
    assert label.tag == 't'
    assert ''.join(label.itertext()) == f'|c(nCl)c(pi*)| = {product:.2f}'
    assert check_cdxml(result['cdxml'], palette=ORBITAL_PALETTE)['ok']


def _assert_gap(result):
    root = _root(result)
    levels = {level['id']: level for level in result['occupation']}
    marker = _node(root, result['gap_arrow_node_id'])
    tail, head = _coords(marker, 'Tail3D'), _coords(marker, 'Head3D')
    assert marker.get('ArrowheadHead') == marker.get('ArrowheadTail') == 'Full'
    assert marker.get('LineType') is None
    assert tail[0] == pytest.approx(head[0], abs=0.5)
    assert tail[1] == pytest.approx(_line_y(root, levels[result['donor_homo']['id']]['line_node_id']), abs=0.5)
    assert head[1] == pytest.approx(_line_y(root, levels[result['acceptor_lumo']['id']]['line_node_id']), abs=0.5)
    assert f'dE = {result["gap"]:g}' in [''.join(node.itertext()) for node in root.iter('t')]


def test_gap_is_marked_between_donor_homo_and_acceptor_lumo():
    _assert_gap(fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR))


def test_lumo_ignores_occupied_acceptor_level():
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['acceptor'].insert(0, {'id': 'core', 'label': 'core', 'energy': -10, 'electrons': 2})
    spec['mixing'].insert(0, {'id': 'core-mixed', 'label': 'core', 'energy': -11,
                               'electrons': 2, 'from': ['core']})
    result = fmo_interaction(spec)
    assert result['donor_homo']['energy'] == 0
    assert result['acceptor_lumo']['id'] == 'pi*'
    assert result['acceptor_lumo']['energy'] == 6
    assert result['gap'] == 6
    _assert_gap(result)
    print('Donor=0; acceptor=-10 (2e), +6 (0e): LUMO=+6, gap=+6')


def test_acceptor_without_empty_orbitals_is_rejected():
    spec = deepcopy(_ONE_DONOR_ONE_ACCEPTOR)
    spec['acceptor'][0]['electrons'] = 2
    with pytest.raises(ValueError, match='acceptor needs at least one empty orbital'):
        fmo_interaction(spec)


def test_delocalisation_is_solid_and_excitation_is_a_vertical_dash():
    spec = {**_ONE_DONOR_ONE_ACCEPTOR,
            'arrows': [{'kind': 'delocalisation'}, {'kind': 'excitation', 'label': 'h\u03bd'}]}
    result = fmo_interaction(spec)
    root = _root(result)
    assert result['arrows'] == [{'kind': 'delocalisation', 'from': 'nCl', 'to': 'psi-'},
                                {'kind': 'excitation', 'from': 'psi-', 'to': 'psi+'}]
    full = []
    for node in root.iter('arrow'):
        if node.get('ArrowheadHead') != 'Full':
            continue
        tail = _coords(node, 'Tail3D')[:2]
        head = _coords(node, 'Head3D')[:2]
        full.append({'dashed': node.get('LineType') == 'Dashed', 'vertical': tail[0] == head[0]})
    assert any(a['dashed'] and a['vertical'] for a in full)
    assert any(not a['dashed'] and not a['vertical'] for a in full)
    assert 'h\u03bd' in ''.join(''.join(node.itertext()) for node in root.iter('t'))
    assert check_cdxml(result['cdxml'])['ok']


def test_labels_support_superscript_and_subscript_markup():
    spec = {'donor': [{'id': 'd', 'label': 'n(Cl)', 'energy': 0, 'electrons': 2}],
            'acceptor': [{'id': 'a', 'label': '\u03c0{sup:*}(C=C)', 'energy': 6, 'electrons': 0}],
            'mixing': [{'id': 'low', 'label': '\u03c8{sub:1}', 'energy': -2, 'electrons': 2, 'from': ['d', 'a']},
                       {'id': 'high', 'label': '\u03c8{sub:2}{sup:*}', 'energy': 8, 'electrons': 0, 'from': ['d', 'a']}]}
    result = fmo_interaction(spec)
    root = _root(result)
    faces = {int(node.get('face', '0')) & 96 for node in root.iter('s')}
    assert 32 in faces and 64 in faces
    assert check_cdxml(result['cdxml'])['ok']


def test_writes_a_new_file_and_never_overwrites(tmp_path):
    out = tmp_path / 'fmo.cdxml'
    result = fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR, out)
    assert Path(result['output_path']).exists()
    with pytest.raises(ValueError, match='overwrite'):
        fmo_interaction(_ONE_DONOR_ONE_ACCEPTOR, out)


def test_operation_is_registered_in_the_diagram_builders():
    assert 'fmo_interaction' in DIAGRAM_OPERATIONS
    result = build('fmo_interaction', _ONE_DONOR_ONE_ACCEPTOR)
    assert result['diagram'] == 'fmo_interaction'


@pytest.mark.parametrize('name,palette', [
    ('fmo-interaction', None),
    ('fmo-coefficients', ORBITAL_PALETTE),
])
def test_bundled_examples_build_and_meet_the_house_rules(name, palette):
    spec_path = EXAMPLES / f'{name}.json'
    if not spec_path.exists():
        pytest.skip('the fmo examples are not beside this toolkit checkout')
    result = build('fmo_interaction', json.loads(spec_path.read_text(encoding='utf-8')))
    _assert_gap(result)
    assert check_cdxml(result['cdxml'], palette=palette)['ok']
