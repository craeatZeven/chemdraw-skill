"""Diagram generators: the chemistry each report states, the errors a self-contradicting spec
raises, and the house rules every bundled example meets. None of this needs ChemDraw or RDKit."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.cdxml_compose import add_ring_circles, merge_cdxml
from cdxml_toolkit.scientific import elements as E
from cdxml_toolkit.scientific.crystal import PRESETS, unit_cell
from cdxml_toolkit.scientific.diagrams import DIAGRAM_OPERATIONS, build
from cdxml_toolkit.scientific.energy import energy_cycle, energy_profile, jablonski
from cdxml_toolkit.scientific.lewis import lewis_structure
from cdxml_toolkit.scientific.orbitals import crystal_field, electron_configuration_diagram, hund_fill, mo_diagram
from cdxml_toolkit.scientific.periodic import periodic_table
from cdxml_toolkit.scientific.process import electrochemical_cell, flowchart

EXAMPLES = Path(__file__).resolve().parents[1] / 'docs' / 'examples' / 'scientific'


def _operation(name):
    return 'energy_profile' if name.startswith('energy_profile') else name


@pytest.mark.parametrize('example', sorted(p.stem for p in EXAMPLES.glob('*.json')
                                           if _operation(p.stem) in DIAGRAM_OPERATIONS))
def test_every_example_meets_the_house_rules_without_layout_warnings(example):
    spec = json.loads((EXAMPLES / f'{example}.json').read_text(encoding='utf-8'))
    result = build(_operation(example), spec)
    check = check_cdxml(result['cdxml'], allow_sizes=[8, 6])
    assert check['ok'], check['errors']
    # Only the honest "not to scale" note is allowed; no caption may sit on other ink.
    assert all('not to scale' in w for w in result['warnings']), result['warnings']


def test_every_example_is_one_of_the_documented_operations():
    for path in EXAMPLES.glob('*.json'):
        assert _operation(path.stem) in DIAGRAM_OPERATIONS or path.stem in (
            'nmr', 'apparatus', 'tlc', 'mechanism'), path.name


def test_diagram_writes_a_new_file_and_never_overwrites(tmp_path):
    out = tmp_path / 'p.cdxml'
    result = build('periodic_table', {'periods': [1, 2]}, out)
    assert Path(result['output_path']).exists() and 'cdxml' not in result
    with pytest.raises(ValueError, match='overwrite'):
        build('periodic_table', {'periods': [1, 2]}, out)


# ---------------------------------------------------------------- elements
def test_configurations_conserve_electrons_and_know_the_exceptions():
    for z in range(1, 119):
        assert sum(E.electron_configuration(z).values()) == z
    assert E.electron_configuration('Cr') == {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '3d': 5, '4s': 1}
    assert E.electron_configuration('Cu')['3d'] == 10 and E.electron_configuration('Cu')['4s'] == 1
    assert E.electron_configuration('Fe', 3)['3d'] == 5 and '4s' not in E.electron_configuration('Fe', 3)
    assert E.electron_configuration('Cl', -1)['3p'] == 6


def test_table_positions_valence_and_d_counts():
    assert E.group_period('He') == (18, 1) and E.group_period('Al') == (13, 3) and E.group_period('Hg') == (12, 6)
    assert E.group_period('La') == (None, 6) and E.group_period('Lu') == (3, 6)
    assert [E.valence_electrons(s) for s in ('H', 'C', 'N', 'O', 'Cl', 'Xe', 'Mg')] == [1, 4, 5, 6, 7, 8, 2]
    assert E.d_electron_count('Fe', 3) == 5 and E.d_electron_count('Pt', 2) == 8 and E.d_electron_count('Cu', 2) == 9
    assert E.parse_ion('Fe(III)') == ('Fe', 3) and E.parse_ion('O2-') == ('O', -2) and E.parse_ion('Cl-') == ('Cl', -1)
    with pytest.raises(ValueError):
        E.valence_electrons('Fe')


# ---------------------------------------------------------------- energy diagrams
def test_energy_profile_reports_barriers_and_rejects_a_false_transition_state():
    result = energy_profile({'points': [{'label': 'A', 'energy': 0}, {'label': 'TS1', 'energy': 20},
                                        {'label': 'I', 'energy': 5}, {'label': 'TS2', 'energy': 25},
                                        {'label': 'B', 'energy': -10}]})
    path = result['pathways'][0]
    assert [b['forward_barrier'] for b in path['barriers']] == [20, 20]
    assert path['barriers'][1]['reverse_barrier'] == 35
    assert path['overall_change'] == -10 and path['thermodynamics'] == 'exothermic'
    assert '‡' in result['cdxml']                     # transition states are marked
    with pytest.raises(ValueError, match='maximum'):
        energy_profile({'points': [{'label': 'A', 'energy': 0}, {'label': 'X', 'energy': -5, 'kind': 'ts'},
                                   {'label': 'B', 'energy': 3}]})


def test_energy_profile_places_stationary_point_structures(tmp_path):
    structure = tmp_path / 'ts.cdxml'
    structure.write_text(build('lewis', {'atoms': [{'element': 'O'}]})['cdxml'], encoding='utf-8')
    result = energy_profile({'points': [{'label': 'A', 'energy': 0}, {'label': 'TS', 'energy': 20,
                                                                      'structure': str(structure)},
                                        {'label': 'B', 'energy': -5}]})
    assert result['structures'] == ['TS']
    assert check_cdxml(result['cdxml'])['ok']


def test_energy_cycle_closes_by_hess_law():
    result = energy_cycle({'levels': [{'id': 'a'}, {'id': 'b'}, {'id': 'c'}],
                           'steps': [{'from': 'a', 'to': 'b', 'value': 10}, {'from': 'b', 'to': 'c', 'value': 5},
                                     {'from': 'a', 'to': 'c', 'value': 15}]})
    assert [s['difference'] for s in result['steps']] == [10, 5, 15]
    with pytest.raises(ValueError, match="Hess"):
        energy_cycle({'levels': [{'id': 'a'}, {'id': 'b'}, {'id': 'c'}],
                      'steps': [{'from': 'a', 'to': 'b', 'value': 10}, {'from': 'b', 'to': 'c', 'value': 5},
                                {'from': 'a', 'to': 'c', 'value': 20}]})


def test_jablonski_enforces_spin_and_direction():
    result = jablonski({'preset': 'standard'})
    kinds = {t['kind']: t['radiative'] for t in result['transitions']}
    assert kinds['fluorescence'] and kinds['phosphorescence'] and not kinds['isc'] and not kinds['ic']
    with pytest.raises(ValueError, match='multiplicity'):
        jablonski({'states': [{'id': 'S0', 'energy': 0}, {'id': 'T1', 'energy': 2}],
                   'transitions': [{'from': 'T1', 'to': 'S0', 'kind': 'fluorescence'}]})
    with pytest.raises(ValueError, match='higher'):
        jablonski({'states': [{'id': 'S0', 'energy': 0}, {'id': 'S1', 'energy': 2}],
                   'transitions': [{'from': 'S1', 'to': 'S0', 'kind': 'absorption'}]})


# ---------------------------------------------------------------- orbitals
@pytest.mark.parametrize('species,charge,order,unpaired', [
    ('H2', 0, 1, 0), ('He2', 0, 0, 0), ('Li2', 0, 1, 0), ('Be2', 0, 0, 0), ('B2', 0, 1, 2), ('C2', 0, 2, 0),
    ('N2', 0, 3, 0), ('O2', 0, 2, 2), ('F2', 0, 1, 0), ('Ne2', 0, 0, 0), ('O2', -1, 1.5, 1), ('N2', 1, 2.5, 1),
    ('O2', 1, 2.5, 1), ('O2', -2, 1, 0)])
def test_mo_presets_give_textbook_bond_orders_and_magnetism(species, charge, order, unpaired):
    result = mo_diagram({'preset': species, 'charge': charge})
    assert result['bond_order'] == order and result['unpaired_electrons'] == unpaired


def test_mo_diagram_conserves_electrons():
    with pytest.raises(ValueError, match='conserved'):
        mo_diagram({'left': {'levels': [{'id': '1s', 'energy': 0, 'electrons': 1}]},
                    'right': {'levels': [{'id': '1s', 'energy': 0, 'electrons': 1}]},
                    'center': {'levels': [{'label': 's', 'energy': -1, 'electrons': 1, 'from': ['1s']}]}})


def test_hund_and_pauli():
    assert hund_fill(4, 3) == ['ud', 'u', 'u'] and hund_fill(3, 3) == ['u', 'u', 'u']
    with pytest.raises(ValueError, match='Hund'):
        electron_configuration_diagram({'subshells': [{'label': '2p', 'boxes': ['ud', 'u', '']}]})
    with pytest.raises(ValueError, match='Pauli'):
        electron_configuration_diagram({'subshells': [{'label': '2p', 'boxes': ['uu', '', '']}]})
    fe = electron_configuration_diagram({'element': 'Fe'})
    assert fe['unpaired_electrons'] == 4 and fe['abbreviated'] == '[Ar] 3d6 4s2'


@pytest.mark.parametrize('spec,unpaired,cfse', [
    ({'ion': 'Fe2+', 'spin': 'low'}, 0, -2.4), ({'ion': 'Fe2+', 'spin': 'high'}, 4, -0.4),
    ({'ion': 'Fe3+', 'spin': 'high'}, 5, 0.0), ({'d_electrons': 3}, 3, -1.2),
    ({'geometry': 'tetrahedral', 'd_electrons': 7}, 3, -1.2)])
def test_crystal_field_fill_and_stabilisation(spec, unpaired, cfse):
    result = crystal_field(spec)
    assert result['unpaired_electrons'] == unpaired and result['cfse']['splitting_units'] == cfse


def test_square_planar_d8_is_diamagnetic():
    result = crystal_field({'geometry': 'square_planar', 'ion': 'Pt2+'})
    assert result['unpaired_electrons'] == 0 and result['spin_only_moment_bm'] == 0


# ---------------------------------------------------------------- Lewis structures
def test_lewis_counts_electrons_and_formal_charge():
    water = lewis_structure({'atoms': [{'id': 'H1', 'element': 'H', 'at': [0, 0]}, {'id': 'O', 'element': 'O', 'at': [1, 0]},
                                       {'id': 'H2', 'element': 'H', 'at': [2, 0]}], 'bonds': [['H1', 'O'], ['O', 'H2']]})
    oxygen = next(a for a in water['species'][0]['atoms'] if a['atom'] == 'O')
    assert oxygen['nonbonding_electrons'] == 4 and oxygen['shell'] == 8 and water['species'][0]['valence_electrons'] == 8
    ammonium = lewis_structure({'atoms': [{'id': 'N', 'element': 'N', 'at': [1, 1], 'formal_charge': 1}]
                                + [{'id': h, 'element': 'H', 'at': at} for h, at in
                                   (('a', [0, 1]), ('b', [2, 1]), ('c', [1, 0]), ('d', [1, 2]))],
                                'bonds': [['a', 'N'], ['N', 'b'], ['c', 'N'], ['N', 'd']]})
    assert ammonium['species'][0]['charge'] == 1
    with pytest.raises(ValueError, match='not enough'):
        lewis_structure({'atoms': [{'id': 'H', 'element': 'H', 'at': [0, 0]}, {'id': 'O', 'element': 'O', 'at': [1, 0]}],
                         'bonds': [['H', 'O', 2]]})
    flagged = lewis_structure({'atoms': [{'id': 'C', 'element': 'C', 'at': [0, 0]}, {'id': 'O', 'element': 'O', 'at': [1, 0]}],
                               'bonds': [['C', 'O', 2]]})
    assert any('fewer than 8' in w for w in flagged['warnings'])


# ---------------------------------------------------------------- unit cells
@pytest.mark.parametrize('preset,counts,cn', [
    ('sc', {'Po': 1}, {'Po': 6}), ('bcc', {'Fe': 2}, {'Fe': 8}), ('fcc', {'Cu': 4}, {'Cu': 12}),
    ('diamond', {'C': 8}, {'C': 4}), ('cscl', {'Cs': 1, 'Cl': 1}, {'Cs': 8, 'Cl': 8}),
    ('nacl', {'Na': 4, 'Cl': 4}, {'Na': 6, 'Cl': 6}), ('zinc_blende', {'Zn': 4, 'S': 4}, {'Zn': 4, 'S': 4}),
    ('fluorite', {'Ca': 4, 'F': 8}, {'Ca': 8, 'F': 4}), ('perovskite', {'Ca': 1, 'Ti': 1, 'O': 3}, {'Ti': 6})])
def test_unit_cell_counts_and_coordination(preset, counts, cn):
    result = unit_cell({'preset': preset})
    assert result['atoms_per_cell'] == counts
    for element, number in cn.items():
        assert result['coordination'][element]['coordination_number'] == number
    assert check_cdxml(result['cdxml'])['ok']


def test_unit_cell_presets_are_all_tested():
    assert set(PRESETS) == {'sc', 'bcc', 'fcc', 'diamond', 'cscl', 'nacl', 'zinc_blende', 'fluorite', 'perovskite'}


def test_unit_cell_expands_a_custom_site_to_its_boundary_copies():
    result = unit_cell({'sites': [{'element': 'X', 'frac': [0.5, 0, 0]}]})   # edge centres: 4 x 1/4
    assert result['atoms_per_cell'] == {'X': 1} and result['atoms_drawn'] == 4
    with pytest.raises(ValueError, match='twice'):
        unit_cell({'sites': [{'element': 'X', 'frac': [0, 0, 0]}, {'element': 'X', 'frac': [1, 0, 0]}]})


# ---------------------------------------------------------------- process diagrams
def test_flowchart_finds_return_loops_and_columns():
    result = flowchart({'nodes': [{'id': 'a'}, {'id': 'b'}, {'id': 'c'}],
                        'edges': [{'from': 'a', 'to': 'b'}, {'from': 'b', 'to': 'c'}, {'from': 'c', 'to': 'b'}]})
    assert result['return_loops'] == ['c -> b']
    assert [result['layout'][n]['col'] for n in 'abc'] == [0, 1, 2]


@pytest.mark.parametrize('kind,anode', [('galvanic', '-'), ('electrolytic', '+')])
def test_electrochemical_cell_polarity(kind, anode):
    result = electrochemical_cell({'type': kind, 'anode': {'electrode': 'Zn'}, 'cathode': {'electrode': 'Cu'}})
    assert result['anode']['process'] == 'oxidation' and result['anode']['polarity'] == anode
    assert result['cathode']['process'] == 'reduction'


def test_periodic_table_draws_every_element():
    result = periodic_table({})
    assert result['elements_drawn'] == 118
    assert periodic_table({'highlight': [{'group': 17}]})['highlighted'] == ['F', 'Cl', 'Br', 'I', 'At', 'Ts']


# ---------------------------------------------------------------- compose and circles
def _benzene():
    import math
    atoms = ''.join(f'<n id="{i + 1}" p="{50 + 14.4 * math.cos(math.pi / 3 * i):.2f} '
                    f'{50 + 14.4 * math.sin(math.pi / 3 * i):.2f}"/>' for i in range(6))
    bonds = ''.join(f'<b id="{10 + i}" B="{i + 1}" E="{(i + 1) % 6 + 1}"{" Order=\"2\"" if i % 2 == 0 else ""}/>'
                    for i in range(6))
    return ('<?xml version="1.0" encoding="UTF-8" ?><CDXML BondLength="14.4"><colortable><color r="0" g="0" b="1"/>'
            '</colortable><fonttable><font id="3" charset="iso-8859-1" name="Arial"/></fonttable><page id="1">'
            f'<fragment id="100">{atoms}{bonds}</fragment><t id="200" p="40 80" color="2"><s font="3" size="10">A</s></t>'
            '</page></CDXML>')


def test_add_ring_circles_turns_kekule_benzene_into_a_circle():
    result = add_ring_circles(_benzene())
    assert result['circles'] == 1 and result['bonds_made_single'] == 3
    root = ET.fromstring(result['cdxml'].split('\n', 2)[2])
    assert not [b for b in root.iter('b') if b.get('Order') == '2']
    assert [g for g in root.iter('graphic') if g.get('OvalType') == 'Circle']


def test_merge_renumbers_ids_and_keeps_colours():
    merged = merge_cdxml([{'cdxml': _benzene()}, {'cdxml': _benzene()}])
    root = ET.fromstring(merged['cdxml'].split('\n', 2)[2])
    ids = [n.get('id') for n in root.iter() if n.get('id') and n.tag != 'font']
    assert len(ids) == len(set(ids))
    texts = [t for t in root.iter('t')]
    colors = {c: tuple(col.attrib.values()) for c, col in enumerate(root.find('colortable').findall('color'), start=2)}
    assert all(colors[int(t.get('color'))] == ('0', '0', '1') for t in texts)
    assert not [e for e in check_cdxml(merged['cdxml'], palette=None)['errors'] if e['code'] in ('duplicate_id', 'dangling_reference', 'color_index')]
    for bond in root.iter('b'):          # every bond still joins atoms of its own copy
        fragment = next(f for f in root.iter('fragment') if bond in list(f))
        atom_ids = {n.get('id') for n in fragment.findall('n')}
        assert bond.get('B') in atom_ids and bond.get('E') in atom_ids


# ---------------------------------------------------------------- checker
def test_checker_finds_what_chemdraw_drops():
    bad = ('<?xml version="1.0" encoding="utf8" ?><CDXML BondLength="20"><fonttable><font id="3" name="Times"/></fonttable>'
           '<page id="1"><fragment id="1"><n id="2" p="0 0"/></fragment><curve id="3" CurvePoints="0 0 1 1 2 2"/>'
           '<t id="4" p="0 0"><s font="3" size="10">\u4e2d\u6587</s></t><b id="5" B="9" E="2"/></page></CDXML>')
    codes = {e['code'] for e in check_cdxml(bad)['errors']}
    assert {'xml_encoding', 'setting', 'blank_structure', 'curve_points', 'cjk_font', 'dangling_reference'} <= codes
    assert check_cdxml(bad, ignore=['setting'])['errors'] and 'setting' not in {
        e['code'] for e in check_cdxml(bad, ignore=['setting'])['errors']}


# ---------------------------------------------------------------- ProgID
def test_progid_override_and_fallback(monkeypatch):
    from cdxml_toolkit.chemdraw import progid
    monkeypatch.setenv('CHEMDRAW_PROGID', 'ChemDraw_x64.Application')
    assert progid.chemdraw_progid() == 'ChemDraw_x64.Application'
    monkeypatch.delenv('CHEMDRAW_PROGID')
    progid._detected.cache_clear()
    monkeypatch.setattr(progid, 'registered', lambda name: name == 'ChemDraw_x64.Application')
    assert progid.chemdraw_progid() == 'ChemDraw_x64.Application'
    progid._detected.cache_clear()


# ---------------------------------------------------------------- Hückel and phase diagrams
@pytest.mark.parametrize('system,energies,rule', [
    ('butadiene', [1.618, 0.618, -0.618, -1.618], None), ('benzene', [2, 1, 1, -1, -1, -2], 'aromatic (4n+2)'),
    ('cyclobutadiene', [2, 0, 0, -2], 'antiaromatic (4n)'), ('allyl_radical', [1.4142, 0, -1.4142], None)])
def test_huckel_energies_and_aromaticity(system, energies, rule):
    from cdxml_toolkit.scientific.huckel import huckel_mo
    result = huckel_mo({'system': system})
    assert [o['x'] for o in result['orbitals']] == pytest.approx(energies, abs=1e-3)
    assert result.get('huckel_rule') == rule
    assert check_cdxml(result['cdxml'])['ok']


def test_huckel_frontier_orbitals_and_nodes():
    from cdxml_toolkit.scientific.huckel import huckel_mo
    butadiene = huckel_mo({'system': 'butadiene'})
    assert (butadiene['homo'], butadiene['lumo'], butadiene['somo']) == (2, 3, [])
    assert [o['nodes_between_atoms'] for o in butadiene['orbitals']] == [0, 1, 2, 3]
    square = huckel_mo({'system': 'cyclobutadiene'})
    assert square['somo'] == [2, 3] and square['lumo'] == 4 and square['unpaired_electrons'] == 2
    assert huckel_mo({'system': 'benzene'})['delocalisation_energy_beta'] == pytest.approx(2.0)
    with pytest.raises(ValueError, match='do not fit'):
        huckel_mo({'atoms': 2, 'electrons': 5})


def test_phase_diagram_orders_its_points():
    from cdxml_toolkit.scientific.phase import phase_diagram
    result = phase_diagram({'triple': [216.6, 5.1], 'critical': [304, 73], 'marks': [{'P': 1, 'label': '1 atm'}]})
    assert result['warnings'] == [] and check_cdxml(result['cdxml'])['ok']
    with pytest.raises(ValueError, match='critical point'):
        phase_diagram({'triple': [300, 5], 'critical': [200, 70]})
