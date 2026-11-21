"""The ring-centre R attachment: its geometry, and that it still reads back as the brief's structure.

Witness: assets/readme/markush/markush.cdxml, whose MultiAttachment sits on a ring's centroid
and reaches its R label through the middle of a double bond. Documented in
docs/rgroup-attachment-geometry.md.
"""
from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET

import pytest
from rdkit import Chem

from cdxml_toolkit import overlay_geometry as G
from cdxml_toolkit.chemistry_semantics import document_inventory, semantic_key
from cdxml_toolkit.mcp_runtime.figure_scope import _validate
from cdxml_toolkit.mcp_runtime.figure_tools import compose_chemical_figure
from cdxml_toolkit.rdkit_utils import collapse_multi_attachments

BOND = 14.4


def hexagon(bond=BOND):
    """A regular hexagon's points in ring order, plus its bond orders (alternating)."""
    xy = {i: (bond * math.cos(math.pi / 3 * i), bond * math.sin(math.pi / 3 * i)) for i in range(6)}
    orders = {frozenset((i, (i + 1) % 6)): (2 if i % 2 == 0 else 1) for i in range(6)}
    return xy, orders


def test_the_label_sits_on_the_ring_centre_side_of_a_double_bond():
    xy, orders = hexagon()
    placement = G.ring_attachment_placement(list(range(6)), xy, orders, attach=0)
    centre = placement["centre"]
    assert centre == pytest.approx((0.0, 0.0), abs=1e-6)
    a, b = placement["crossed"]
    assert orders[frozenset((a, b))] == 2, "the convention crosses a double bond"
    midpoint = ((xy[a][0] + xy[b][0]) / 2, (xy[a][1] + xy[b][1]) / 2)
    reach = math.dist(centre, placement["label"])
    assert reach == pytest.approx(G.R_ATTACHMENT_BOND_RATIO * BOND, abs=1e-6)
    # The line from the centre to the label runs through the crossed bond's midpoint.
    direction = ((placement["label"][0] - centre[0]) / reach, (placement["label"][1] - centre[1]) / reach)
    along = (midpoint[0] - centre[0]) * direction[0] + (midpoint[1] - centre[1]) * direction[1]
    perp = abs((midpoint[0] - centre[0]) * direction[1] - (midpoint[1] - centre[1]) * direction[0])
    assert along > 0 and perp == pytest.approx(0.0, abs=1e-6)


def test_the_declared_atom_leads_the_attachment_list_so_readback_is_stable():
    xy, orders = hexagon()
    placement = G.ring_attachment_placement(list(range(6)), xy, orders, attach=3)
    assert placement["attachments"][0] == 3
    assert sorted(placement["attachments"]) == list(range(6))


def test_a_substituted_ring_atom_is_not_a_free_position():
    xy, orders = hexagon()
    blocked = {1, 2, 3, 4}
    placement = G.ring_attachment_placement(list(range(6)), xy, orders, attach=0, blocked=blocked)
    assert set(placement["attachments"]) == {0, 5}


def test_a_ring_with_no_open_edge_keeps_the_vertex_bond():
    xy, orders = hexagon()
    # Every ring atom but the declared one carries something, so no bond has two free atoms.
    blocked = {1, 2, 3, 4, 5}
    assert G.ring_attachment_placement(list(range(6)), xy, orders, attach=0, blocked=blocked) is None


def _compose(tmp_path, smiles, r_attachment=None):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "version": 1,
        "style": {"preset": "acs1996"},
        "objects": [{"type": "molecule", "id": "m1", "smiles": smiles, "position": [150, 110],
                     **({"r_attachment": r_attachment} if r_attachment else {})}],
    }), encoding="utf-8")
    out = tmp_path / "figure.cdxml"
    result = compose_chemical_figure(str(manifest), str(out))
    assert result["ok"], result
    return out


def _ring_centroid(fragment, mark):
    """The centroid of the ring the MultiAttachment belongs to, walked from its bonds.

    The node sits at the centre of the *whole* ring, not of the positions still free: the
    markush witness has three attachment atoms out of six ring atoms and its node is on the
    six-membered centroid.
    """
    nodes = {n.get("id"): n for n in fragment.findall("n")}
    skip = {mark.get("id")}
    adjacency = {}
    for b in fragment.findall("b"):
        a, e = b.get("B"), b.get("E")
        if a in skip or e in skip:
            continue
        adjacency.setdefault(a, []).append(e)
        adjacency.setdefault(e, []).append(a)
    start = (mark.get("Attachments") or "").split()[0]
    order, previous, current = [start], None, start
    for _ in range(12):
        following = [j for j in adjacency.get(current, []) if j != previous]
        if not following:
            break
        previous, current = current, following[0]
        if current == start:
            break
        order.append(current)
    pts = [[float(v) for v in nodes[i].get("p").split()] for i in order]
    return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))


@pytest.mark.parametrize("smiles", ["[*:1]c1ccccc1", "[*:1]c1ccc([*:2])cc1", "[*:1]c1cccnc1"])
def test_centre_mode_writes_a_multiattachment_on_the_ring_centroid(tmp_path, smiles):
    out = _compose(tmp_path, smiles, "centre")
    root = ET.fromstring(out.read_text(encoding="utf-8").split("?>", 1)[-1])
    marks = [n for n in root.iter("n") if n.get("NodeType") == "MultiAttachment"]
    assert len(marks) == smiles.count("[*:")
    for mark in marks:
        cx, cy = [float(v) for v in mark.get("p").split()]
        fragment = next(f for f in root.iter("fragment") if mark in list(f.findall("n")))
        expected_x, expected_y = _ring_centroid(fragment, mark)
        assert cx == pytest.approx(expected_x, abs=1e-6)
        assert cy == pytest.approx(expected_y, abs=1e-6)
        ring_ids = (mark.get("Attachments") or "").split()
        assert len(ring_ids) >= 3
        nodes = {n.get("id"): n for n in root.iter("n")}
        # It is bonded to the R label, and the bond runs from the centroid.
        link = next(b for b in root.iter("b") if mark.get("id") in (b.get("B"), b.get("E")))
        other = link.get("E") if link.get("B") == mark.get("id") else link.get("B")
        assert nodes[other].get("NodeType") == "GenericNickname"
        lx, ly = [float(v) for v in nodes[other].get("p").split()]
        assert math.hypot(lx - cx, ly - cy) > 0


@pytest.mark.parametrize("smiles", ["[*:1]c1ccccc1", "[*:1]c1ccc([*:2])cc1", "[*:1]c1cccnc1"])
def test_the_drawn_structure_reads_back_as_the_briefs_own(tmp_path, smiles):
    out = _compose(tmp_path, smiles, "centre")
    found = document_inventory(out)
    wanted = semantic_key(Chem.MolFromSmiles(smiles))
    assert dict(found) == {wanted: 1}, (dict(found), wanted)


def test_vertex_mode_writes_no_multiattachment(tmp_path):
    out = _compose(tmp_path, "[*:1]c1ccccc1")
    root = ET.fromstring(out.read_text(encoding="utf-8").split("?>", 1)[-1])
    assert not [n for n in root.iter("n") if n.get("NodeType") == "MultiAttachment"]
    assert dict(document_inventory(out)) == {semantic_key(Chem.MolFromSmiles("[*:1]c1ccccc1")): 1}


def test_a_metal_ligand_multiattachment_is_not_folded():
    page = ('<CDXML><page id="1"><fragment id="1">'
            '<n id="10" p="0 0" Element="26"/>'
            '<n id="11" p="10 10" NodeType="MultiAttachment" Attachments="12 13"/>'
            '<n id="12" p="20 0"/><n id="13" p="20 20"/>'
            '<b id="14" B="10" E="11" Order="dative"/></fragment></page></CDXML>')
    frag = ET.fromstring(page).find("page").find("fragment")
    folded = collapse_multi_attachments(frag)
    assert any(n.get("NodeType") == "MultiAttachment" for n in folded.findall("n"))


def test_an_r_marks_multiattachment_is_folded_onto_the_declared_atom():
    page = ('<CDXML><page id="1"><fragment id="1">'
            '<n id="10" p="0 0" NodeType="MultiAttachment" Attachments="11 12"/>'
            '<n id="11" p="10 0"/><n id="12" p="0 10"/>'
            '<n id="13" p="20 0" NodeType="GenericNickname" GenericNickname="R1"/>'
            '<b id="14" B="10" E="13"/></fragment></page></CDXML>')
    frag = ET.fromstring(page).find("page").find("fragment")
    folded = collapse_multi_attachments(frag)
    assert not [n for n in folded.findall("n") if n.get("NodeType") == "MultiAttachment"]
    bond = folded.find("b")
    assert {bond.get("B"), bond.get("E")} == {"11", "13"}


def test_the_brief_key_is_validated():
    brief = {"kind": "scope", "species": [{"id": "a", "smiles": "c1ccccc1"}],
             "r_attachment": "sideways"}
    assert any("r_attachment" in e for e in _validate(brief))
    brief["r_attachment"] = "centre"
    assert _validate(brief) == []
