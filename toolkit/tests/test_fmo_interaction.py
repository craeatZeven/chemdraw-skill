"""CDXML geometry, occupation and chemistry contracts for frontier-orbital interaction diagrams."""
from __future__ import annotations

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
    return [node for node in root.iter('arrow') if node.get('ArrowheadHead') == 'HalfLeft']


def _node(root, node_id):
    return next(node for node in root.iter() if node.get('id') == node_id)


def _line_y(root, node_id):
    return tuple(float(value) for value in _node(root, node_id).get('BoundingBox').split())[1]


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

    # The numbers the acceptance asks for are in the CDXML itself.
    root = _root(result)
    assert len(_solid_lines(root)) == 4
    assert len(_dashed_lines(root)) == len(result['links']) == 4
    assert {link['node_id'] for link in result['links']} == {node.get('id') for node in _dashed_lines(root)}
    y = {level['id']: _line_y(root, level['line_node_id']) for level in result['occupation']}
    assert y['psi-'] > y['nCl'] > y['psi+']          # filled below the donor, empty above it
    assert y['nCl'] > y['pi*']                        # the acceptor sits above the donor
    assert len(_half_arrows(root)) == 4               # two electrons on the donor, two on the filled level
    assert check_cdxml(result['cdxml'])['ok']


def test_generated_levels_must_straddle_the_donor():
    spec = {'donor': [{'id': 'd', 'label': 'd', 'energy': 0, 'electrons': 2}],
            'acceptor': [{'id': 'a', 'label': 'a', 'energy': 6, 'electrons': 0}],
            'mixing': [{'id': 'low', 'label': 'low', 'energy': 1, 'electrons': 2, 'from': ['d', 'a']},
                       {'id': 'high', 'label': 'high', 'energy': 8, 'electrons': 0, 'from': ['d', 'a']}]}
    with pytest.raises(ValueError, match='straddle'):
        fmo_interaction(spec)


def test_electrons_are_conserved_across_the_columns():
    spec = {'donor': [{'id': 'd', 'label': 'd', 'energy': 0, 'electrons': 2}],
            'acceptor': [{'id': 'a', 'label': 'a', 'energy': 6, 'electrons': 0}],
            'mixing': [{'id': 'low', 'label': 'low', 'energy': -2, 'electrons': 1, 'from': ['d', 'a']},
                       {'id': 'high', 'label': 'high', 'energy': 8, 'electrons': 0, 'from': ['d', 'a']}]}
    with pytest.raises(ValueError, match='conserved'):
        fmo_interaction(spec)


def test_mixing_from_names_an_unknown_level():
    spec = {'donor': [{'id': 'd', 'label': 'd', 'energy': 0, 'electrons': 2}],
            'acceptor': [{'id': 'a', 'label': 'a', 'energy': 6, 'electrons': 0}],
            'mixing': [{'id': 'low', 'label': 'low', 'energy': -2, 'electrons': 2, 'from': ['d', 'x']},
                       {'id': 'high', 'label': 'high', 'energy': 8, 'electrons': 0, 'from': ['d', 'a']}]}
    with pytest.raises(ValueError, match='unknown level'):
        fmo_interaction(spec)


def test_coefficient_lobes_reuse_the_s2_colors_and_scale_with_magnitude():
    spec = {'donor': [{'id': 'd', 'label': 'd', 'energy': 0, 'electrons': 2, 'coefficient': 1.0}],
            'acceptor': [{'id': 'a', 'label': 'a', 'energy': 6, 'electrons': 0, 'coefficient': 0.5, 'phase': -1}],
            'mixing': [{'id': 'low', 'label': 'low', 'energy': -2, 'electrons': 2, 'from': ['d', 'a']},
                       {'id': 'high', 'label': 'high', 'energy': 8, 'electrons': 0, 'from': ['d', 'a']}]}
    result = fmo_interaction(spec)
    root = _root(result)
    ovals = [node for node in root.iter('graphic') if node.get('GraphicType') == 'Oval']
    assert len(ovals) == 2
    colors = _colors(root)
    assert {colors[int(node.get('color'))] for node in ovals} == {'#C0392B', '#0072B2'}
    lengths = []
    for node in ovals:
        center = tuple(float(v) for v in node.get('Center3D').split())[:2]
        major = tuple(float(v) for v in node.get('MajorAxisEnd3D').split())[:2]
        lengths.append(2 * ((major[0] - center[0]) ** 2 + (major[1] - center[1]) ** 2) ** 0.5)
    assert lengths[0] == pytest.approx(26.0, abs=0.05)      # |c| = 1.0
    assert lengths[1] == pytest.approx(13.0, abs=0.05)      # |c| = 0.5
    assert check_cdxml(result['cdxml'], palette=ORBITAL_PALETTE)['ok']


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
        tail = tuple(float(v) for v in node.get('Tail3D').split())[:2]
        head = tuple(float(v) for v in node.get('Head3D').split())[:2]
        full.append({'dashed': node.get('LineType') == 'Dashed', 'vertical': tail[0] == head[0]})
    assert any(a['dashed'] and a['vertical'] for a in full)       # the Franck-Condon excitation
    assert any(not a['dashed'] and not a['vertical'] for a in full)   # the delocalisation arrow
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
    assert 32 in faces and 64 in faces     # sub and sup runs both present
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
    assert check_cdxml(result['cdxml'], palette=palette)['ok']
