"""Geometry and CDXML contract tests for native ring variable-attachment R groups."""
from __future__ import annotations

import json
import math
from xml.etree import ElementTree as ET

import pytest

from cdxml_toolkit.chemistry_semantics import read_molecules
from cdxml_toolkit.mcp_runtime.figure_tools import compose_chemical_figure


def _side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def _midpoint(a, b):
    return tuple((start + end) / 2 for start, end in zip(a, b))


def _measure_variable_attachments(cdxml_path, ring_atom_indices):
    root = ET.parse(str(cdxml_path)).getroot()
    frag = root.find('./page/fragment')
    assert frag is not None
    assert frag.find('n[@NodeType="MultiAttachment"]') is None

    nodes = frag.findall('n')
    node_by_id = {n.get('id'): n for n in nodes}
    bonds = frag.findall('b')
    bond_by_id = {b.get('id'): b for b in bonds}

    def pt(node):
        return tuple(map(float, node.get('p').split()))

    ring_nodes = [nodes[i] for i in ring_atom_indices]
    ring_edges = []
    n_ring = len(ring_nodes)
    for idx in range(n_ring):
        u_id = ring_nodes[idx].get('id')
        v_id = ring_nodes[(idx + 1) % n_ring].get('id')
        edge_bond = next(b for b in bonds if {b.get('B'), b.get('E')} == {u_id, v_id})
        ring_edges.append((edge_bond, pt(ring_nodes[idx]), pt(ring_nodes[(idx + 1) % n_ring])))

    anchors = [n for n in nodes if n.get('NodeType') == 'VariableAttachment']
    results = []
    for anchor in anchors:
        aid = anchor.get('id')
        vbond = next(b for b in bonds if aid in (b.get('B'), b.get('E')))
        lid = vbond.get('E') if vbond.get('B') == aid else vbond.get('B')
        label = node_by_id[lid]
        a_pt, l_pt = pt(anchor), pt(label)
        bond_len = math.dist(a_pt, l_pt)

        crossed_edges = [
            (eb, p1, p2) for eb, p1, p2 in ring_edges
            if _side(p1, p2, a_pt) * _side(p1, p2, l_pt) < 0
            and _side(a_pt, l_pt, p1) * _side(a_pt, l_pt, p2) < 0
        ]
        assert len(crossed_edges) == 1
        crossed, p1, p2 = crossed_edges[0]
        assert vbond.get('CrossingBonds') == crossed.get('id')
        assert vbond.get('id') in (crossed.get('CrossingBonds') or '').split()

        t_el = label.find('t')
        text = ''.join(s.text or '' for s in t_el.iter('s')) if t_el is not None else ''
        results.append({
            'anchor': anchor,
            'label': label,
            'bond': vbond,
            'crossed': crossed,
            'text': text,
            'edge_length': math.dist(p1, p2),
            'inside': math.dist(a_pt, _midpoint(p1, p2)),
            'outside': math.dist(l_pt, _midpoint(p1, p2)),
            'midpoint_error': abs(_side(a_pt, l_pt, _midpoint(p1, p2))) / bond_len,
            'normal_x': (l_pt[0] - a_pt[0]) / bond_len,
            'normal_y': (l_pt[1] - a_pt[1]) / bond_len,
            'crossed_double': crossed.get('Order') == '2',
        })
    return root, frag, results


@pytest.mark.parametrize('name, smiles, ring_atoms, groups, expected_texts', [
    (
        'mono',
        'c1ccccc1',
        [0, 1, 2, 3, 4, 5],
        [{'number': 1, 'atoms': [0, 1, 2, 3, 4, 5], 'color': '#C81E1E'}],
        ['R\u00b9'],
    ),
    (
        'di',
        'c1ccccc1',
        [0, 1, 2, 3, 4, 5],
        [
            {'number': 1, 'atoms': [0, 1, 2, 3, 4, 5], 'color': '#C81E1E'},
            {'number': 2, 'atoms': [0, 1, 2, 3, 4, 5], 'color': '#C81E1E'},
        ],
        ['R\u00b9', 'R\u00b2'],
    ),
    (
        'pyridine',
        'c1ccncc1',
        [0, 1, 2, 3, 4, 5],
        [{'number': 1, 'atoms': [0, 1, 2, 4, 5], 'color': '#C81E1E'}],
        ['R\u00b9'],
    ),
    (
        'thiophene',
        'c1ccsc1',
        [0, 1, 2, 3, 4],
        [{'number': 1, 'atoms': [0, 1, 2, 4], 'color': '#C81E1E'}],
        ['R\u00b9'],
    ),
])
def test_ring_r_group_geometry_and_crossing_rules(tmp_path, name, smiles, ring_atoms, groups, expected_texts):
    manifest = {
        'style': {'bond_length': 14.4, 'font_size': 9.0, 'font': 'Arial'},
        'objects': [{
            'type': 'molecule',
            'id': name,
            'smiles': smiles,
            'position': [100, 100],
            'rotation': 90 if name in ('mono', 'di', 'pyridine') else 0,
            'ring_r_groups': groups,
        }],
    }
    mpath = tmp_path / f'{name}.json'
    out = tmp_path / f'{name}.cdxml'
    mpath.write_text(json.dumps(manifest), encoding='utf-8')
    report = compose_chemical_figure(str(mpath), str(out))
    assert report['ok'] is True, report

    root, frag, results = _measure_variable_attachments(out, ring_atoms)
    assert len(results) == len(groups)
    assert [r['text'] for r in results] == expected_texts

    used_crossed = set()
    for r in results:
        assert r['midpoint_error'] < 0.02
        assert r['inside'] == pytest.approx(0.325 * r['edge_length'], abs=0.025)
        assert r['outside'] > 0.27 * r['edge_length']
        assert r['crossed_double'] is True, (name, r['crossed'].attrib)
        assert r['crossed'].get('id') not in used_crossed
        used_crossed.add(r['crossed'].get('id'))
        # Bond must stay scaffold black, not red like the R label
        assert r['bond'].get('color') is None
        assert r['label'].get('color') is not None
        assert r['label'].get('color') != frag.get('color')
        t_el = r['label'].find('t')
        assert t_el is not None and t_el.get('InterpretChemically') == 'no'

    if name == 'di':
        first, second = results
        assert first['normal_x'] == pytest.approx(-1)
        assert second['normal_x'] == pytest.approx(0.5, abs=0.01)
        assert second['normal_y'] == pytest.approx(math.sqrt(3) / 2, abs=0.01)
        assert second['outside'] > 0.48 * second['edge_length']
        assert 0 < float(second['crossed'].get('MarginWidth')) < float(root.get('MarginWidth'))
        crossing_bonds = [b for b in frag.findall('b') if b.get('CrossingBonds')]
        assert len(crossing_bonds) == 4

    # RDKit scaffold readback succeeds cleanly
    mols = read_molecules(out)
    assert len(mols) == 1


def test_variable_attachment_is_normal_to_a_skewed_edge():
    from rdkit import Chem

    from cdxml_toolkit.layout.ring_r_groups import ring_r_group_layout, ring_r_group_sites

    mol = Chem.MolFromSmiles('C1=CCCCC1')
    points = [(10, 10), (10, 24), (24, 33), (37, 27), (33, 8), (20, 3)]
    groups = ring_r_group_sites(mol, [{'number': 1, 'atoms': [0, 1]}])
    (placement,) = ring_r_group_layout(mol, points, groups, 14.4, 9)

    assert placement.crossed == mol.GetBondBetweenAtoms(0, 1).GetIdx()
    # The ring centroid is above this edge's midpoint. The connector stays horizontal.
    assert placement.anchor[1] == pytest.approx(17)
    assert placement.label[1] == pytest.approx(17)
