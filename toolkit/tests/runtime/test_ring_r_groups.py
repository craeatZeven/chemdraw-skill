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
    ring_pts = [pt(n) for n in ring_nodes]
    center = (sum(p[0] for p in ring_pts) / len(ring_pts),
              sum(p[1] for p in ring_pts) / len(ring_pts))
    ring_radius = sum(math.dist(p, center) for p in ring_pts) / len(ring_pts)

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
        anchor_dist = math.dist(a_pt, center)

        crossed_edges = [
            eb for eb, p1, p2 in ring_edges
            if _side(p1, p2, a_pt) * _side(p1, p2, l_pt) < 0
            and _side(a_pt, l_pt, p1) * _side(a_pt, l_pt, p2) < 0
        ]
        assert len(crossed_edges) == 1
        crossed = crossed_edges[0]
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
            'bond_ratio': bond_len / ring_radius,
            'anchor_ratio': anchor_dist / ring_radius,
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
        assert 0.30 <= r['bond_ratio'] <= 0.55, (name, r['bond_ratio'])
        assert r['anchor_ratio'] >= 0.45, (name, r['anchor_ratio'])
        assert r['crossed_double'] is True, (name, r['crossed'].attrib)
        assert r['crossed'].get('id') not in used_crossed
        used_crossed.add(r['crossed'].get('id'))
        # Bond must stay scaffold black, not red like the R label
        assert r['bond'].get('color') is None
        assert r['label'].get('color') is not None
        assert r['label'].get('color') != frag.get('color')
        t_el = r['label'].find('t')
        assert t_el is not None and t_el.get('InterpretChemically') == 'no'

    # RDKit scaffold readback succeeds cleanly
    mols = read_molecules(out)
    assert len(mols) == 1
