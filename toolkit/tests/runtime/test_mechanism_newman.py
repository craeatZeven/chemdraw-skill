"""Newman projections placed among the structures of a mechanism scheme.

BJOC 22-125, scheme 4 and 22-130, scheme 25 carry their stereochemical argument with Newman
projections drawn inside the scheme, next to ordinary structures, linked by step arrows.
"""
from __future__ import annotations

import copy
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.mcp_runtime import figure_mechanism as M
from cdxml_toolkit.render.renderer import _IDGen

ANTI = {"front": ["Me", "H", "H"], "rear": ["Br", "H", "Me"], "conformation": "staggered"}
GAUCHE = {"front": ["Me", "H", "H"], "rear": ["H", "Me", "Br"], "conformation": "staggered"}


def _xml(frag, start=5000):
    ids = _IDGen(start)
    return [ET.fromstring(x) for x in M._newman_xml(frag, ids.next)], ids


def test_a_newman_projection_is_sized_like_the_structures_around_it():
    frag = M._newman_frag("n1", ANTI, "A")
    width = frag.bbox[2] - frag.bbox[0]
    # spokes of 1.4 bond lengths: wider than a C-C bond, narrower than a six-membered ring and its labels
    assert 30 < width < 70
    assert frag.ref.label == "A" and frag.newman["front"] == ["Me", "H", "H"]


def test_its_xml_takes_fresh_ids_and_keeps_the_rear_bonds_under_the_disc():
    frag = M._newman_frag("n1", ANTI, None)
    M.R._shift_atoms(frag.atoms, 200.0, 100.0)
    elements, ids = _xml(frag)
    fragment = next(e for e in elements if e.tag == "fragment")
    discs = [e for e in elements if e.tag == "graphic"]
    every = [e for e in elements if e is not fragment] + list(fragment.iter())
    numbers = [int(e.get("id")) for e in every if e.get("id")]
    assert len(numbers) == len(set(numbers)) and min(numbers) >= 5000
    nodes = {n.get("id") for n in fragment.iter("n")}
    assert all(b.get("B") in nodes and b.get("E") in nodes for b in fragment.iter("b"))
    filled = next(d for d in discs if d.get("OvalType") == "Circle Filled")
    cx, cy, _ = map(float, filled.get("Center3D").split())
    assert (round(cx, 2), round(cy, 2)) == (200.0, 100.0)
    z = {b.get("id"): int(b.get("Z")) for b in fragment.iter("b")}
    rear = [z[b.get("id")] for b in fragment.iter("b") if b.get("Z") and int(b.get("Z")) < int(filled.get("Z"))]
    front = [z[b.get("id")] for b in fragment.iter("b") if int(b.get("Z")) > int(filled.get("Z"))]
    assert len(rear) >= 3 and len(front) >= 3  # three rear spokes beneath the disc, the front ones above it


def test_a_newman_species_is_checked():
    brief = {"kind": "mechanism", "species": [{"id": "n", "newman": {"front": ["Me", "H"], "rear": ["H", "H", "H"]}}],
             "steps": [{"substrates": ["n"], "products": ["n"]}]}
    assert any("front" in e for e in M._validate(brief))
    brief["species"][0]["newman"] = {"front": ["Me", "H", "H"], "rear": ["H", "H", "H"], "spin": 3}
    assert any("newman" in e for e in M._validate(brief))


def test_an_arrow_cannot_point_at_a_newman_projection():
    brief = {"kind": "mechanism",
             "species": [{"id": "n", "newman": ANTI}, {"id": "p", "smiles": "C=C"}],
             "steps": [{"substrates": ["n"], "products": ["p"],
                        "moves": [{"from": {"species": "n", "atom": 0}, "to": {"species": "n", "atom": 1}}]}]}
    assert any("Newman" in e for e in M._validate(brief))


def test_a_step_with_a_newman_projection_may_have_no_moves():
    brief = {"kind": "mechanism", "species": [{"id": "a", "newman": ANTI}, {"id": "b", "newman": GAUCHE}],
             "steps": [{"substrates": ["a"], "products": ["b"], "arrow": "reversible"}]}
    assert M._validate(brief) == []


# --- whole briefs (these need RDKit) ----------------------------------------------------------

def test_two_projections_lay_out_among_the_structures_without_overlap(tmp_path, monkeypatch):
    monkeypatch.setenv("CHEMDRAW_OFFLINE", "1")
    from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure
    brief = {
        "kind": "mechanism",
        "species": [
            {"id": "sub", "smiles": "CC(Br)CC", "label": "1"},
            {"id": "anti", "newman": ANTI, "label": "1-anti"},
            {"id": "gauche", "newman": GAUCHE, "label": "1-gauche"},
        ],
        "steps": [
            {"substrates": ["sub"], "products": ["anti"]},
            {"substrates": ["anti"], "products": ["gauche"], "arrow": "reversible"},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "newman.cdxml"))
    assert report["ok"] is True, report.get("errors")
    assert report["lint"]["issues"] == []
    assert report["readback"]["status"] == "matched"
    root = ET.parse(tmp_path / "newman.cdxml").getroot()
    discs = [g for g in root.iter("graphic") if g.get("OvalType") == "Circle Filled"]
    assert len(discs) == 2


@pytest.mark.parametrize("offset", [0, 120, 240, -120, 365, 125])
def test_a_rear_set_on_top_of_the_front_set_is_refused(offset):
    # dihedral_offset is the rear set's angle from the front set (60 is staggered), not an
    # addition to it: 120 puts every rear spoke and label on a front one.
    spec = {"front": ["CH3", "H", "H"], "rear": ["CH3", "H", "H"], "dihedral_offset": offset}
    errors = M._newman_errors("g", spec)
    assert any("on top of the front" in e for e in errors), errors


@pytest.mark.parametrize("spec", [
    {"conformation": "staggered"}, {"conformation": "eclipsed"}, {"dihedral_offset": 60},
    {"dihedral_offset": 14}, {"dihedral_offset": -14}, {"dihedral_offset": 180},
])
def test_every_label_of_an_accepted_projection_has_a_place_of_its_own(spec):
    spec = {"front": ["CH3", "H", "H"], "rear": ["H", "CH3", "H"], **spec}
    assert M._newman_errors("p", spec) == []
    frag = M._newman_frag("p", spec, "p")
    places = [(round(a["x"], 1), round(a["y"], 1)) for a in frag.atoms if a.get("show_label")]
    assert len(places) == 6 and len(set(places)) == 6
