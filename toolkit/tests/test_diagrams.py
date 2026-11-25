"""Diagram generators: the chemistry each report states, the errors a self-contradicting spec
raises, and the house rules every bundled example meets. None of this needs ChemDraw or RDKit."""
import json
import re
import math
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
    """The operation an example file is for: the longest operation name its stem starts with."""
    matches = [op for op in DIAGRAM_OPERATIONS if name == op or name.startswith(op + '_')]
    return max(matches, key=len) if matches else name


@pytest.mark.parametrize('example', sorted(p.stem for p in EXAMPLES.glob('*.json')
                                           if _operation(p.stem) in DIAGRAM_OPERATIONS))
def test_every_example_meets_the_house_rules_without_layout_warnings(example):
    spec = json.loads((EXAMPLES / f'{example}.json').read_text(encoding='utf-8'))
    result = build(_operation(example), spec)
    check = check_cdxml(result['cdxml'], allow_sizes=[8, 6])
    assert check['ok'], check['errors']
    # No line, arrow, curve or other caption may run through a caption (the checker's estimate
    # matches the text boxes ChemDraw writes when it saves the file).
    assert [w for w in check['warnings'] if w['code'] == 'text_over_ink'] == []
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
    assert '\u2021' in result['cdxml']                     # transition states are marked
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
    # A backslash inside an f-string expression is a syntax error before Python 3.12,
    # so the optional double bond attribute is built outside the f-string.
    bonds = ''.join(
        '<b id="%d" B="%d" E="%d"%s/>' % (10 + i, i + 1, (i + 1) % 6 + 1,
                                          ' Order="2"' if i % 2 == 0 else '')
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
           '<page id="1"><fragment id="1"><n id="2" p="0 0"/></fragment><curve id="3" CurvePoints="0 0 1 1 2"/>'
           '<t id="4" p="0 0"><s font="3" size="10">\u4e2d\u6587</s></t><b id="5" B="9" E="2"/></page></CDXML>')
    codes = {e['code'] for e in check_cdxml(bad)['errors']}
    assert {'xml_encoding', 'setting', 'blank_structure', 'curve_points', 'cjk_font', 'dangling_reference'} <= codes
    assert check_cdxml(bad, ignore=['setting'])['errors'] and 'setting' not in {
        e['code'] for e in check_cdxml(bad, ignore=['setting'])['errors']}


# ---------------------------------------------------------------- Huckel and phase diagrams
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


def _e_axis_arrow_head_tail_y(cdxml):
    """The ``y`` of the energy-axis arrow's head and tail, from the finished CDXML text."""
    root = ET.fromstring(cdxml.split('\n', 2)[2])
    for node in root.iter('arrow'):
        if node.get('ArrowheadHead') == 'Full':
            return (float(node.get('Head3D').split()[1]), float(node.get('Tail3D').split()[1]))
    raise AssertionError('no E-axis arrow with a full head found')


def test_huckel_levels_ordered_by_energy():
    """Drawn levels must obey the energy axis: the most bonding level sits at the page bottom.

    The report exposes the drawn ``y`` (page coordinates, ``y`` down), so these assert physical
    ordering — not snapshot values — and lock the fix for the inverted-level bug.
    """
    from cdxml_toolkit.scientific.huckel import huckel_mo

    butadiene = huckel_mo({'system': 'butadiene'})['orbitals']
    assert [o['x'] for o in butadiene] == pytest.approx([1.618, 0.618, -0.618, -1.618], abs=1e-3)
    assert all(a['y'] > b['y'] for a, b in zip(butadiene, butadiene[1:]))

    benzene = huckel_mo({'system': 'benzene'})['orbitals']
    assert benzene[0]['y'] > benzene[1]['y'] == benzene[2]['y'] > benzene[3]['y'] == benzene[4]['y'] > benzene[5]['y']

    cyclobutadiene = huckel_mo({'system': 'cyclobutadiene'})['orbitals']
    assert cyclobutadiene[0]['y'] > cyclobutadiene[1]['y'] == cyclobutadiene[2]['y'] > cyclobutadiene[3]['y']

    pyridine = huckel_mo({'system': 'pyridine'})['orbitals']
    assert all(a['y'] >= b['y'] for a, b in zip(pyridine, pyridine[1:]))

    head_y, tail_y = _e_axis_arrow_head_tail_y(huckel_mo({'system': 'butadiene'})['cdxml'])
    assert head_y < tail_y                                    # E arrow points up the page


def _hetero_label_gaps(cdxml):
    """Size-8 captions (heteroatom symbols) and their min gaps to every orbital circle/curve box."""
    from cdxml_toolkit.scientific.canvas import caption_box
    root = ET.fromstring(cdxml.split('\n', 2)[2])
    labels, graphics = [], []
    for t in root.iter('t'):
        if not t.get('p'):
            continue
        size = max((float(s.get('size', 10) or 10) for s in t.iter('s')), default=10.0)
        if size != 8:
            continue
        x, y = (float(v) for v in t.get('p').split()[:2])
        markup = ''.join((s.text or '') for s in t.iter('s'))
        labels.append((markup, caption_box(x, y, markup, size, t.get('Justification', 'Left'))))
    for g in root.iter('graphic'):
        if g.get('BoundingBox') and g.get('OvalType') in ('Circle', 'Circle Filled'):
            v = [float(q) for q in g.get('BoundingBox').split()[:4]]
            graphics.append((min(v[0], v[2]), min(v[1], v[3]), max(v[0], v[2]), max(v[1], v[3])))
    for c in root.iter('curve'):
        v = [float(q) for q in c.get('CurvePoints', '').split()]
        if len(v) < 6:
            continue
        xs, ys = v[0::2], v[1::2]
        graphics.append((min(xs), min(ys), max(xs), max(ys)))
    gaps = []
    for _markup, box in labels:
        for gb in graphics:
            dx = max(box[0] - gb[2], gb[0] - box[2], 0.0)
            dy = max(box[1] - gb[3], gb[1] - box[3], 0.0)
            gaps.append(math.hypot(dx, dy))
    return labels, gaps


def test_huckel_hetero_labels_clear_orbital_graphics():
    """Heteroatom labels sit outside the orbital ink of their own and neighbouring levels.

    The report only exposes energy ordering, so this parses the drawn CDXML. It first checks the
    size-8 captions are exactly the declared heteroatom symbols (count included), then that each
    stays at least half a font size (4 pt) clear of every orbital circle and curve box.
    """
    from cdxml_toolkit.scientific.huckel import huckel_mo, HETERO_SYSTEMS, HETEROATOMS

    for name in HETERO_SYSTEMS:
        result = huckel_mo({'system': name})
        labels, gaps = _hetero_label_gaps(result['cdxml'])
        _, _, kinds = HETERO_SYSTEMS[name]
        expected = sorted(HETEROATOMS[kind]['symbol'] for kind in kinds.values())
        found = sorted(markup for markup, _ in labels)
        assert found == expected, f'{name}: heteroatom labels {found} != declared {expected}'
        for gap in gaps:
            assert gap >= 4.0, f'{name}: a heteroatom label is {gap:.2f} pt from orbital ink; need >= 4'


def test_huckel_custom_heteroatom_index_clears_graphics():
    """A public ``heteroatoms`` spec (not a preset name) must avoid the same overlap.

    Regresses the index-4 pyridine case, where text_over_ink reports zero while the label truly
    overlaps an orbital disc; the assertion therefore goes through the geometric boxes.
    """
    from cdxml_toolkit.scientific.huckel import huckel_mo

    result = huckel_mo({'atoms': 6, 'cyclic': True, 'electrons': 6,
                        'heteroatoms': [{'index': 4, 'type': 'N_pyridine'}]})
    labels, gaps = _hetero_label_gaps(result['cdxml'])
    assert [markup for markup, _ in labels] == ['N']
    assert gaps and min(gaps) >= 4.0, f'index-4 label gaps {[round(g, 2) for g in gaps]} need >= 4'


def test_phase_diagram_orders_its_points():
    from cdxml_toolkit.scientific.phase import phase_diagram
    result = phase_diagram({'triple': [216.6, 5.1], 'critical': [304, 73], 'marks': [{'P': 1, 'label': '1 atm'}]})
    assert result['warnings'] == [] and check_cdxml(result['cdxml'])['ok']
    with pytest.raises(ValueError, match='critical point'):
        phase_diagram({'triple': [300, 5], 'critical': [200, 70]})


# ---------------------------------------------------------------- hexagonal cells and binary diagrams
@pytest.mark.parametrize('preset,counts,cn,packing', [
    ('hcp', {'Mg': 6}, {'Mg': 12}, 0.7405), ('wurtzite', {'S': 6, 'Zn': 6}, {'S': 4, 'Zn': 4}, None),
    ('simple_hexagonal', {'X': 3}, {'X': 8}, 0.6046), ('fcc', {'Cu': 4}, {'Cu': 12}, 0.7405),
    ('bcc', {'Fe': 2}, {'Fe': 8}, 0.6802), ('diamond', {'C': 8}, {'C': 4}, 0.3401)])
def test_cells_count_coordinate_and_pack(preset, counts, cn, packing):
    result = unit_cell({'preset': preset})
    assert result['atoms_per_cell'] == counts
    assert {e: v['coordination_number'] for e, v in result['coordination'].items()} == cn
    assert result.get('packing_fraction') == (pytest.approx(packing, abs=2e-4) if packing else None)
    assert check_cdxml(result['cdxml'])['ok']


def test_hexagonal_atoms_do_not_hide_each_other():
    import itertools, math, re
    for preset in ('hcp', 'wurtzite', 'simple_hexagonal'):
        text = unit_cell({'preset': preset})['cdxml']
        centres = [tuple(map(float, m.split()[:2])) for m in re.findall(r'OvalType="Circle" [^>]*?Center3D="([^"]+)"', text)]
        assert min(math.dist(a, b) for a, b in itertools.combinations(centres, 2)) > 10.0


def test_binary_eutectic_fields_and_lever_rule():
    from cdxml_toolkit.scientific.binary import binary_phase_diagram
    result = binary_phase_diagram({'type': 'eutectic', 'components': ['Pb', 'Sn'], 'melting': [327, 232],
                                   'eutectic': [0.619, 183], 'solubility': [0.183, 0.977],
                                   'points': [{'x': 0.4, 'T': 230}, {'x': 0.4, 'T': 150}, {'x': 0.8, 'T': 300}]})
    assert set(result['fields']) == {'L', 'L + \u03b1', 'L + \u03b2', '\u03b1', '\u03b1 + \u03b2', '\u03b2'}
    fields = [p['field'] for p in result['points']]
    assert fields == ['L + \u03b1', '\u03b1 + \u03b2', 'L']
    for point in result['points']:
        assert sum(q['fraction'] for q in point['phases']) == pytest.approx(1.0)
        if len(point['phases']) == 2:      # the lever rule: the mixture's composition is recovered
            assert sum(q['fraction'] * q['composition'] for q in point['phases']) == pytest.approx(point['x'], abs=2e-3)
    assert result['warnings'] == [] and check_cdxml(result['cdxml'], allow_sizes=[8])['ok']
    with pytest.raises(ValueError, match='below both melting points'):
        binary_phase_diagram({'melting': [300, 200], 'eutectic': [0.5, 250]})
    with pytest.raises(ValueError, match='either side'):
        binary_phase_diagram({'melting': [300, 200], 'eutectic': [0.5, 150], 'solubility': [0.6, 0.9]})


def test_binary_isomorphous_uses_measured_curves():
    from cdxml_toolkit.scientific.binary import binary_phase_diagram
    result = binary_phase_diagram({'type': 'isomorphous', 'melting': [1085, 1455],
                                   'liquidus': [[0.4, 1290]], 'solidus': [[0.4, 1220]],
                                   'points': [{'x': 0.4, 'T': 1250}]})
    liquid, solid = result['points'][0]['phases']
    assert result['points'][0]['field'] == 'L + \u03b1' and liquid['composition'] < 0.4 < solid['composition']
    with pytest.raises(ValueError, match='above the solidus'):
        binary_phase_diagram({'type': 'isomorphous', 'melting': [1000, 1400],
                              'liquidus': [[0.5, 1100]], 'solidus': [[0.5, 1150]]})


# ---------------------------------------------------------------- redox, correlation, Walsh, surfaces
_FE = {'element': 'Fe', 'activity': 1e-6, 'species': [
    {'name': 'Fe2+', 'dGf': -78.90, 'charge': 2}, {'name': 'Fe3+', 'dGf': -4.7, 'charge': 3},
    {'name': 'Fe2O3', 'dGf': -742.2, 'atoms': 2, 'O': 3}, {'name': 'Fe3O4', 'dGf': -1015.4, 'atoms': 3, 'O': 4}]}


def test_pourbaix_fields_follow_from_gibbs_energies():
    from cdxml_toolkit.scientific.redox import pourbaix
    result = pourbaix({**_FE, 'points': [{'pH': 1, 'E': 0.3}, {'pH': 7, 'E': 0.0}, {'pH': 2, 'E': -0.8},
                                         {'pH': 0.5, 'E': 1.2}]})
    assert set(result['fields']) == {'Fe2+', 'Fe3+', 'Fe2O3', 'Fe3O4', 'Fe'}
    assert [p['stable'] for p in result['points']] == ['Fe2+', 'Fe2O3', 'Fe', 'Fe3+']
    assert result['warnings'] == [] and check_cdxml(result['cdxml'], allow_sizes=[8])['ok']


def test_pourbaix_fe2_fe_boundary_is_the_nernst_potential():
    from cdxml_toolkit.scientific.redox import pourbaix
    # Fe2+/Fe at a = 1e-6: E = dGf/(2F) + (RT/2F) ln a = -0.409 - 0.177 V, independent of pH.
    below = pourbaix({**_FE, 'points': [{'pH': 2, 'E': -0.60}]})['points'][0]['stable']
    above = pourbaix({**_FE, 'points': [{'pH': 2, 'E': -0.57}]})['points'][0]['stable']
    assert (below, above) == ('Fe', 'Fe2+')


def test_latimer_derives_couples_and_finds_disproportionation():
    from cdxml_toolkit.scientific.redox import frost, latimer
    spec = json.loads((EXAMPLES / 'latimer.json').read_text(encoding='utf-8'))
    result = latimer(spec)
    assert result['disproportionate'] == ['MnO{sub:4}{sup:2\u2212}', 'Mn{sup:3+}']
    assert [d['E'] for d in result['derived']] == pytest.approx([1.693, 1.508], abs=1e-3)
    points = frost({k: v for k, v in spec.items() if k != 'derived'})['series'][0]['points']
    assert [p['NE'] for p in points] == pytest.approx([5.18, 4.62, 0.10, -0.85, -2.36, 0.0], abs=1e-2)
    with pytest.raises(ValueError, match='highest oxidation state'):
        latimer({'species': [{'name': 'a', 'oxidation_state': 1}, {'name': 'b', 'oxidation_state': 2}], 'potentials': [0.1]})


def test_ellingham_crossings():
    from cdxml_toolkit.scientific.redox import ellingham
    result = ellingham(json.loads((EXAMPLES / 'ellingham.json').read_text(encoding='utf-8')))
    crossing = {tuple(sorted(x['lines'])): x for x in result['crossings']}
    co = crossing[tuple(sorted(['2C + O{sub:2} \u2192 2CO', 'C + O{sub:2} \u2192 CO{sub:2}']))]
    assert co['T'] == pytest.approx(981, abs=2) and co['lower_above'].endswith('2CO')
    assert not any('MgO' in name for x in result['crossings'] for name in x['lines'])   # Mg crosses above 2273 K
    with pytest.raises(ValueError, match='kink'):
        ellingham({'lines': [{'label': 'x', 'dH': -500, 'dS': -100, 'segments': [{'from_T': 1000, 'dH': -400, 'dS': -100}]}]})


@pytest.mark.parametrize('preset,allowed', [
    ('butadiene_conrotatory', True), ('butadiene_disrotatory', False), ('hexatriene_conrotatory', False),
    ('hexatriene_disrotatory', True), ('cycloaddition_2+2', False), ('cycloaddition_4+2', True)])
def test_woodward_hoffmann_correlations(preset, allowed):
    from cdxml_toolkit.scientific.correlation import correlation_diagram
    result = correlation_diagram({'preset': preset})
    assert result['allowed_thermally'] is allowed
    assert check_cdxml(result['cdxml'], allow_sizes=[8])['ok']


def test_correlation_conserves_symmetry():
    from cdxml_toolkit.scientific.correlation import correlation_diagram
    with pytest.raises(ValueError, match='conserved'):
        correlation_diagram({'left': [{'label': 'a', 'symmetry': 'S'}], 'right': [{'label': 'b', 'symmetry': 'A'}],
                             'electrons': 2})


@pytest.mark.parametrize('electrons,geometry', [(4, 'linear'), (5, 'bent'), (6, 'bent'), (8, 'bent')])
def test_walsh_ah2(electrons, geometry):
    from cdxml_toolkit.scientific.correlation import walsh_diagram
    assert walsh_diagram({'preset': 'AH2', 'electrons': electrons})['geometry'] == geometry


def test_mueller_brown_stationary_points():
    from cdxml_toolkit.scientific.surface import potential_surface
    result = potential_surface({'preset': 'mueller_brown'})
    minima = [p for p in result['stationary_points'] if p['kind'] == 'minimum']
    saddles = [p for p in result['stationary_points'] if p['kind'] != 'minimum']
    assert [round(p['energy'], 1) for p in minima] == [-146.7, -108.2, -80.8]
    assert [round(p['energy'], 1) for p in saddles] == [-72.2, -40.7]
    assert (minima[0]['x'], minima[0]['y']) == pytest.approx((-0.558, 1.442), abs=2e-3)
    assert check_cdxml(result['cdxml'], allow_sizes=[8])['ok']


def test_grid_surface_finds_a_saddle():
    from cdxml_toolkit.scientific.surface import potential_surface
    xs = [i / 10 for i in range(-20, 21)]
    z = [[(x * x - 1) ** 2 + 2 * y * y for x in xs] for y in xs]       # double well, saddle at (0, 0)
    result = potential_surface({'grid': {'x': xs, 'y': xs, 'z': z}, 'levels': 8})
    kinds = sorted((p['kind'], round(p['x'], 1), round(p['y'], 1)) for p in result['stationary_points'])
    assert kinds == [('minimum', -1.0, 0.0), ('minimum', 1.0, 0.0), ('transition state', 0.0, 0.0)]


def test_a_formal_charge_is_drawn_on_its_atom():
    # The charge used to be bookkeeping only: a nitro group came out neutral on paper.
    ammonium = lewis_structure({'style': 'lines', 'atoms': [{'id': 'N', 'element': 'N', 'at': [0, 0], 'formal_charge': 1}]
                                + [{'id': h, 'element': 'H', 'at': at} for h, at in
                                   (('a', [-1, 0]), ('b', [0, -1]), ('c', [1, 0]), ('d', [0, 1]))],
                                'bonds': [['a', 'N'], ['N', 'b'], ['c', 'N'], ['N', 'd']]})
    xml = ammonium['cdxml']
    assert re.search(r'face="64"[^>]*>\+<', xml), 'no superscript plus on the ammonium nitrogen'
    water = lewis_structure({'style': 'lines',
                             'atoms': [{'id': 'H1', 'element': 'H', 'at': [0, 0]},
                                       {'id': 'O', 'element': 'O', 'at': [1, 0]},
                                       {'id': 'H2', 'element': 'H', 'at': [2, 0]}],
                             'bonds': [['H1', 'O'], ['O', 'H2']]})
    assert not re.search(r'face="64"', water['cdxml']), 'a neutral molecule should carry no superscript charge'


def test_lewis_hex_grid_places_pairs_between_bonds():
    from cdxml_toolkit.scientific.lewis import _groups
    result = lewis_structure(json.loads((EXAMPLES / 'lewis_carbonate.json').read_text(encoding='utf-8')))
    atoms = {a['atom']: a for a in result['species'][0]['atoms']}
    assert result['species'][0]['charge'] == -2 and atoms['O2']['nonbonding_electrons'] == 6
    assert _groups({'id': 'O', 'nonbonding': 4, 'bonds': {0: None, 180: None}}) == {90: 2, 270: 2}
    with pytest.raises(ValueError, match='neighbours on the hex grid'):
        lewis_structure({'grid': 'hex', 'atoms': [{'id': 'a', 'element': 'H', 'at': [0, 0]},
                                                  {'id': 'b', 'element': 'H', 'at': [1, 1]}], 'bonds': [['a', 'b']]})


@pytest.mark.parametrize('system,energies', [
    ('formaldehyde', [1.618, -0.618]), ('pyridine', [2.107, 1.167, 1.0, -0.841, -1.0, -1.934])])
def test_huckel_heteroatoms(system, energies):
    from cdxml_toolkit.scientific.huckel import huckel_mo
    result = huckel_mo({'system': system})
    assert [o['x'] for o in result['orbitals']] == pytest.approx(energies, abs=1e-3)
    assert 'delocalisation_energy_beta' not in result and result['heteroatoms']


# ---------------------------------------------------------------- Tanabe-Sugano
def test_tanabe_sugano_free_ion_terms_are_exact():
    from cdxml_toolkit.scientific.tanabe import levels
    cb = 4.42
    free = levels(2, 0.0, cb)
    terms = {f'{m}{label}': round(e - free[0][0], 6) for e, m, label in free}
    assert terms == pytest.approx({'3F': 0, '3P': 15, '1D': 5 + 2 * cb, '1G': 12 + 2 * cb, '1S': 22 + 7 * cb})


@pytest.mark.parametrize('n,first', [(3, '4T{sub:2g}'), (8, '3T{sub:2g}')])
def test_first_spin_allowed_band_of_d3_and_d8_is_delta(n, first):
    from cdxml_toolkit.scientific.tanabe import levels
    for delta in (5.0, 12.0, 30.0):
        lv = levels(n, delta, 4.5)
        ground = lv[0]
        allowed = [(e - ground[0], f'{m}{label}') for e, m, label in lv[1:] if m == ground[1]]
        assert allowed[0] == (pytest.approx(delta, abs=1e-6), first)


@pytest.mark.parametrize('n,weak,strong,crossover', [
    (2, '3T{sub:1g}', '3T{sub:1g}', None), (3, '4A{sub:2g}', '4A{sub:2g}', None),
    (4, '5E{sub:g}', '3T{sub:1g}', 27.3), (5, '6A{sub:1g}', '2T{sub:2g}', 28.3),
    (6, '5T{sub:2g}', '1A{sub:1g}', 20.0), (7, '4T{sub:1g}', '2E{sub:g}', 21.7), (8, '3A{sub:2g}', '3A{sub:2g}', None)])
def test_ground_states_and_spin_crossover(n, weak, strong, crossover):
    from cdxml_toolkit.scientific.tanabe import tanabe_sugano
    result = tanabe_sugano({'d_electrons': n, 'samples': 80})
    assert (result['weak_field_ground'], result['strong_field_ground']) == (weak, strong)
    if crossover is None:
        assert result['spin_crossover'] == []
    else:
        assert result['spin_crossover'][0]['delta_over_b'] == pytest.approx(crossover, abs=0.15)
    assert check_cdxml(result['cdxml'], allow_sizes=[8, 6])['ok']


def test_band_fit_recovers_delta_for_d3():
    from cdxml_toolkit.scientific.tanabe import tanabe_sugano
    fit = tanabe_sugano({'d_electrons': 3, 'bands': [17400, 24500], 'samples': 60})['band_fit'][0]
    assert fit['Delta'] == pytest.approx(17400, abs=1) and fit['assignment'] == ['4T{sub:2g}', '4T{sub:1g}']
    assert 600 < fit['B'] < 900


# ---------------------------------------------------------------- species-list energy-profile format
_TWO_STEP = {'style': 'smooth', 'unit': 'kJ/mol', 'energy_symbol': 'E', 'y_label_style': 'slash', 'decimals': 0,
             'show_plus': True, 'reference': 'R',
             'species': [{'id': 'R', 'label': 'R', 'energy': 10}, {'id': 'TS1', 'energy': 85, 'ts': True},
                         {'id': 'I', 'energy': 40}, {'id': 'TS2', 'energy': 70, 'ts': True},
                         {'id': 'P', 'energy': -30, 'label_position': 'below'}],
             'annotations': [{'from': 'R', 'to': 'TS1', 'label': 'E{sub:a} = {value}'},
                             {'from': 'P', 'to': 'R', 'heads': 'from', 'x': 4.5, 'label': '\u0394H'}],
             'unit_width': 60}


def test_species_list_format_is_translated():
    from cdxml_toolkit.scientific.energy import from_species_spec
    spec, ignored = from_species_spec(_TWO_STEP)
    points = spec['pathways'][0]['points']
    assert [p['energy'] for p in points] == [0, 75, 30, 60, -40]                 # relative to the reference
    assert [p['kind'] for p in points] == ['minimum', 'ts', 'minimum', 'ts', 'minimum']
    assert [p['x'] for p in points] == [0, 0.25, 0.5, 0.75, 1.0]
    assert points[-1]['label_side'] == 'below'
    assert spec['arrows'][0]['label'] == 'E{sub:a} = +75'                         # {value}, signed, 0 decimals
    assert (spec['arrows'][1]['from'], spec['arrows'][1]['to'], spec['arrows'][1]['heads']) == ('R', 'P', 'to')
    assert spec['arrows'][1]['x'] == pytest.approx(1.125)
    assert spec['y_label'] == '{i:E} / (kJ mol{sup:\u22121})' and spec['label_bold'] and not spec['ts_marker']
    assert ignored == ['unit_width']


def test_species_list_profile_builds_and_reports_barriers():
    result = energy_profile(_TWO_STEP)
    barriers = result['pathways'][0]['barriers']
    assert [b['forward_barrier'] for b in barriers] == [75, 30]
    assert result['translated_from'] == 'species-list format' and result['layout_keys_ignored'] == ['unit_width']
    assert result['warnings'] == [] and check_cdxml(result['cdxml'], allow_sizes=[8, 6])['ok']


@pytest.mark.parametrize('end', ['first', 'last'])
def test_gap_arrow_beside_an_end_point_gets_room_for_its_label(end):
    # A torsion scan: the barrier arrow stands beside the first (or last) minimum, where the
    # page margin, not the curve, has to hold the label.
    species = [{'id': 'A', 'label': 'coplanar 0°', 'energy': 28},
               {'id': 'B', 'label': 'perpendicular ±90°', 'energy': 0},
               {'id': 'C', 'label': 'coplanar 180°', 'energy': 26}]
    target = 'A' if end == 'first' else 'C'
    result = energy_profile({'style': 'smooth', 'unit': 'kJ/mol', 'species': species,
                             'annotations': [{'from': 'B', 'to': target, 'label': 'rotation barrier'}]})
    assert result['warnings'] == []
    assert check_cdxml(result['cdxml'])['ok']


@pytest.mark.parametrize('name', ['energy-profile.spec.json', 'levels.json'])
def test_preview_energy_profiles_build_cleanly(name):
    preview = Path(__file__).resolve().parents[2] / 'skill' / 'chemdraw' / 'assets' / 'preview' / name
    if not preview.exists():
        pytest.skip('the Skill previews are not beside this toolkit checkout')
    result = energy_profile(json.loads(preview.read_text(encoding='utf-8')))
    assert result['warnings'] == []
    assert check_cdxml(result['cdxml'], ignore=['palette'])['ok']


def test_profile_accepts_a_conformer_template_as_structure():
    result = energy_profile({'points': [{'label': 'A', 'energy': 0, 'structure': {'template': 'chair'}},
                                        {'label': 'TS', 'energy': 40},
                                        {'label': 'B', 'energy': 5, 'structure': {'template': 'chair', 'flip': True}}]})
    assert result['structures'] == ['A', 'B']
    root = ET.fromstring(result['cdxml'].split('\n', 2)[2])
    assert len(list(root.iter('fragment'))) == 2 and check_cdxml(result['cdxml'], allow_sizes=[8, 6])['ok']


# ---------------------------------------------------------------- native-render hazards
def _doc(body, page='<page id="1" BoundingBox="0 0 200 100" Width="200" Height="100" WidthPages="2" HeightPages="2">'):
    return ('<?xml version="1.0" encoding="UTF-8" ?><CDXML BondLength="14.4" LineWidth="0.6" BoldWidth="2" '
            'HashSpacing="2.5" MarginWidth="1.6" LabelSize="10" CaptionSize="10" ChainAngle="120" LabelFont="3" '
            'CaptionFont="3"><fonttable><font id="3" charset="iso-8859-1" name="Arial"/>'
            '<font id="5" charset="Unicode" name="SimSun"/></fonttable>' + page + body + '</page></CDXML>')


def test_checker_flags_what_chemdraw_drops_or_wraps():
    codes = lambda text: {i['code'] for i in check_cdxml(text)['errors'] + check_cdxml(text)['warnings']}
    mixed = _doc('<t id="2" p="10 10"><s font="5" size="10">\u4e2dA</s></t>')
    assert 'mixed_script_run' in codes(mixed)
    markup = _doc('<t id="2" p="10 10"><s font="3" size="10">E{sub:a}</s></t>')
    assert 'unresolved_markup' in codes(markup)
    polyline = _doc('<curve id="2" CurvePoints="0 0 0 0 0 0 10 10 10 10 10 10" Closed="no" FillType="None"/>')
    assert 'curve_type_missing' in codes(polyline)
    assert 'curve_type_missing' not in codes(polyline.replace('<curve id="2"', '<curve id="2" CurveType="0"'))
    saved_by_chemdraw = _doc('<curve id="2" CurvePoints="0 0 0 0 0 0 10 10 10 10 10 10"/>')   # no Closed/FillType
    assert 'curve_type_missing' not in codes(saved_by_chemdraw)
    bezier = _doc('<curve id="2" CurvePoints="0 0 0 0 5 -5 10 -5 15 0 15 0"/>')      # an electron arrow
    assert 'curve_type_missing' not in codes(bezier)
    head = _doc('<arrow id="2" Head3D="50 10 0" Tail3D="10 10 0" HeadSize="1000" ArrowheadCenterSize="500" '
                'ArrowheadWidth="250"/>')
    assert 'arrowhead_ratio' in codes(head)
    outside = _doc('<t id="2" p="350 10"><s font="3" size="10">A</s></t>')
    assert 'outside_page' in codes(outside)
    one_page = _doc('<t id="2" p="650 10"><s font="3" size="10">A</s></t>',
                    '<page id="1" BoundingBox="0 0 700 100" Width="700" Height="100">')
    assert 'text_may_wrap' in codes(one_page)


def test_generated_pages_are_wide_enough_not_to_wrap():
    for result in (build('periodic_table', {}), build('flowchart', json.loads(
            (EXAMPLES / 'flowchart.json').read_text(encoding='utf-8')))):
        page = ET.fromstring(result['cdxml'].split('\n', 2)[2]).find('page')
        width = float(page.get('Width'))
        assert int(page.get('WidthPages')) >= max(2, math.ceil(width / 500) + 1)
        assert int(page.get('HeightPages')) >= 2
        assert check_cdxml(result['cdxml'], allow_sizes=[8, 6])['ok']


def test_new_checks_pass_the_bundled_native_files():
    """CDXML files already rendered natively for the README must not trip the render-hazard checks."""
    root = Path(__file__).resolve().parents[2]
    for name in ('assets/readme/scientific/apparatus.cdxml', 'assets/readme/scientific/kinetics.cdxml',
                 'assets/readme/scientific/nmr.cdxml', 'assets/readme/engines/nmr/caffeine/comparison.cdxml'):
        path = root / name
        if not path.exists():
            pytest.skip('README assets are not in this checkout')
        found = check_cdxml(str(path), palette=None)
        codes = {i['code'] for i in found['errors'] + found['warnings']}
        assert not codes & {'curve_type_missing', 'text_may_wrap', 'mixed_script_run', 'unresolved_markup',
                            'text_over_ink'}, (name, codes)


def _page(body):
    return ('<?xml version="1.0" encoding="UTF-8" ?><CDXML><fonttable><font id="3" charset="utf-8" name="Arial"/>'
            '</fonttable><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/></colortable>'
            f'<page id="1" Width="300" Height="200">{body}</page></CDXML>')


def test_text_over_ink_finds_a_line_through_a_caption_and_nothing_beside_it():
    from cdxml_toolkit.cdxml_check import text_over_ink
    caption = '<t id="2" p="50 50"><s font="3" size="10" face="0">Fe</s><s font="3" size="10" face="32">2</s></t>'
    through = '<graphic id="3" GraphicType="Line" BoundingBox="40 47 90 47"/>'
    beside = '<graphic id="3" GraphicType="Line" BoundingBox="40 56 90 56"/>'
    assert text_over_ink(_page(caption + through)) == ["a line runs through 'Fe{sub:2}'"]
    assert text_over_ink(_page(caption + beside)) == []
    assert check_cdxml(_page(caption + through))['warnings'][-1]['code'] == 'text_over_ink'


def test_text_over_ink_trusts_chemdraw_boxes_and_white_patches():
    from cdxml_toolkit.cdxml_check import text_over_ink
    line = '<graphic id="3" GraphicType="Line" BoundingBox="60 40 60 60"/>'
    # ChemDraw's own box (written when it saved the file) ends before the line.
    measured = '<t id="2" p="50 50" BoundingBox="50 42 58 52"><s font="3" size="10" face="0">MMMM</s></t>'
    assert text_over_ink(_page(line + measured)) == []
    # A white patch drawn after the line and under the caption hides the line there.
    caption = '<t id="4" p="55 50"><s font="3" size="10" face="0">TS</s></t>'
    patch = ('<curve id="5" CurvePoints="53 40 53 40 53 40 70 40 70 40 70 40 70 54 70 54 70 54 53 54 53 54 53 54"'
             ' Closed="yes" FillType="Solid" CurveType="129" color="2"/>')
    assert text_over_ink(_page(line + caption)) == ["a line runs through 'TS'"]
    assert text_over_ink(_page(line + patch + caption)) == []
    assert text_over_ink(_page(patch + line + caption)) == ["a line runs through 'TS'"]


def test_caption_widths_follow_arial():
    from cdxml_toolkit.scientific.canvas import text_width
    # ChemDraw measures "tautomerization" in 13 pt Arial as 88.9 pt wide.
    assert text_width('tautomerization', 13) == pytest.approx(88.9, abs=0.3)

def test_a_single_triple_curve_is_accepted():
    """Six numbers is one whole (handle, point, handle) triple — the smallest curve
    ChemDraw draws, and what Clipware stationery uses. It must not be reported."""
    good = ('<?xml version="1.0" encoding="utf-8" ?><CDXML BondLength="14.4">'
            '<fonttable><font id="3" name="Arial"/></fonttable>'
            '<page id="1"><curve id="3" CurvePoints="0 0 10 20 20 0"/></page></CDXML>')
    codes = {e['code'] for e in check_cdxml(good)['errors']}
    assert 'curve_points' not in codes

def test_a_setting_rounded_by_chemdraw_is_tolerated():
    """ChemDraw's own SaveAs writes 2.5 as 2.49 and 1.6 as 1.59.

    An exact comparison failed every file that had been through the native save, so the
    rule stops being a signal. 0.01 is rounding; a real deviation still reports.
    """
    import re as _re

    def codes(**overrides):
        text = _doc('')
        for key, value in overrides.items():
            text = _re.sub(r'%s="[^"]*"' % key, '%s="%s"' % (key, value), text, count=1)
        return {i['code'] for i in check_cdxml(text)['errors']}

    assert 'setting' not in codes(HashSpacing='2.49', MarginWidth='1.59')
    assert 'setting' in codes(HashSpacing='2.4')

def test_a_methyl_halide_is_not_a_bare_terminal_carbon():
    """A line ending at Br reads as CH3Br; that is the skeletal convention, not a defect.

    The warning belongs on a fragment where *both* ends are unlabelled: a bare line that
    identifies nothing.
    """
    halide = _doc('<fragment id="1"><n id="2" p="20 20"/><n id="3" p="40 20" Element="35"/>'
                  '<b id="4" B="2" E="3"/></fragment>')
    bare_line = _doc('<fragment id="1"><n id="2" p="20 20"/><n id="3" p="40 20"/>'
                     '<b id="4" B="2" E="3"/></fragment>')
    codes = lambda text: {i['code'] for i in check_cdxml(text)['errors'] + check_cdxml(text)['warnings']}
    assert 'bare_terminal_carbon' not in codes(halide)
    assert 'bare_terminal_carbon' in codes(bare_line)


@pytest.mark.parametrize('head, shape, flagged', [
    ('Full', ('1000', '875', '250'), False),
    ('Full', ('1000', '1000', '500'), True),
    ('HalfLeft', ('1000', '1000', '500'), True),    # a solid barb that reads as two electrons
    ('HalfLeft', ('600', '525', '150'), False),     # ChemDraw's own half head
    ('HalfLeft', ('1000', '700', '400'), True),
])
def test_a_half_head_keeps_the_house_proportions(head, shape, flagged):
    arrow = (f'<arrow id="5" ArrowheadHead="{head}" HeadSize="{shape[0]}" ArrowheadCenterSize="{shape[1]}" '
             f'ArrowheadWidth="{shape[2]}" Head3D="60 40 0" Tail3D="20 40 0"/>')
    codes = [w['code'] for w in check_cdxml(_page(arrow))['warnings']]
    assert ('arrowhead_ratio' in codes) is flagged


def _arc_arrow(ident, tail, head, centre, angle, end='Head'):
    r = math.dist(tail, centre)
    heads = 'ArrowheadHead="Full"' if end == 'Head' else 'ArrowheadHead="None" ArrowheadTail="Full"'
    return (f'<arrow id="{ident}" {heads} HeadSize="1000" ArrowheadCenterSize="875" ArrowheadWidth="250" '
            f'AngularSize="{angle}" Tail3D="{tail[0]} {tail[1]} 0" Head3D="{head[0]} {head[1]} 0" '
            f'Center3D="{centre[0]} {centre[1]} 0" MajorAxisEnd3D="{centre[0] + r} {centre[1]} 0" '
            f'MinorAxisEnd3D="{centre[0]} {centre[1] - r} 0"/>')


def test_check_reports_curved_arrows_that_run_into_each_other():
    # A quarter circle about (50, 50), r 30, from its east point anticlockwise to its north point
    # (50, 20); another about (50, -10) from its west point anticlockwise to its south point, the
    # same (50, 20): the two heads land together.
    first = _arc_arrow(10, (80, 50), (50, 20), (50, 50), 90)
    meeting = _arc_arrow(11, (20, -10), (50, 20), (50, -10), 90)
    elsewhere = _arc_arrow(12, (230, 50), (200, 20), (200, 50), 90)
    codes = [w['code'] for w in check_cdxml(_page(first + meeting))['warnings']]
    assert 'arrows_collide' in codes
    assert 'arrows_collide' not in [w['code'] for w in check_cdxml(_page(first + elsewhere))['warnings']]


def test_check_ignores_curved_arrows_that_only_share_a_tail():
    # Both leave the east point (80, 50): one turns up, the other (written from its target) down.
    up = _arc_arrow(10, (80, 50), (50, 20), (50, 50), 90)
    down = _arc_arrow(11, (50, 80), (80, 50), (50, 50), 90, end='Tail')
    assert 'arrows_collide' not in [w['code'] for w in check_cdxml(_page(up + down))['warnings']]
