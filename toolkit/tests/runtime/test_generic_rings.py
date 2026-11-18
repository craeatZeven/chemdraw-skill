"""An unspecified ring drawn as a dashed arc through two of a structure's atoms (parity item P9).

Witness: BJOC 22-125 scheme 4, where X and the carbon beside it carry a dashed ring that stands
for any ring fused there; nothing is drawn inside it.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit import overlay_geometry as G
from cdxml_toolkit.mcp_runtime import figure_draw, figure_house, figure_mechanism


def _pairs(text: str) -> list[tuple[float, float]]:
    v = [float(x) for x in text.split()]
    return list(zip(v[0::2], v[1::2]))


A, B = (100.0, 100.0), (114.4, 100.0)   # a bond, 14.4 pt
INSIDE = (107.2, 112.0)                   # the rest of the molecule is below the bond


@pytest.mark.parametrize("members", [5, 6, 7])
def test_the_arc_runs_the_long_way_from_one_atom_to_the_other(members):
    pts = _pairs(G.generic_ring_points(A, B, INSIDE, members=members))
    anchors = pts[1::3]
    assert anchors[0] == pytest.approx(A, abs=0.01) and anchors[-1] == pytest.approx(B, abs=0.01)
    r = 14.4 / (2 * math.sin(math.pi / members))          # the ring's circumradius
    cy = 100.0 - math.sqrt(r * r - 7.2 * 7.2)               # the centre is above the bond, away from INSIDE
    for x, y in anchors:
        assert math.hypot(x - 107.2, y - cy) == pytest.approx(r, abs=0.05)
    # Every point of the arc but its ends lies on the far side of the bond: the short arc under
    # the bond, where the solid bond is, is not drawn.
    assert all(y < 100.0 - 0.5 for _, y in anchors[1:-1])
    turn = [(math.atan2(y2 - cy, x2 - 107.2) - math.atan2(y1 - cy, x1 - 107.2) + math.pi) % (2 * math.pi) - math.pi
            for (x1, y1), (x2, y2) in zip(anchors, anchors[1:])]
    assert all(t > 0 for t in turn) or all(t < 0 for t in turn)   # one way round
    sweep = abs(sum(turn))
    assert math.degrees(sweep) == pytest.approx(360 - 360 / members, abs=0.5)


def test_an_arc_stops_short_of_an_atom_label():
    pts = _pairs(G.generic_ring_points(A, B, INSIDE, members=6, clear_a=5.0, clear_b=0.0))
    first, last = pts[1], pts[-2]
    assert math.dist(first, A) == pytest.approx(5.0, abs=0.3)
    assert last == pytest.approx(B, abs=0.01)


def test_the_box_holds_the_arc():
    box = G.generic_ring_box(A, B, INSIDE, members=6)
    pts = _pairs(G.generic_ring_points(A, B, INSIDE, members=6))
    assert box[0] <= min(x for x, _ in pts) and box[2] >= max(x for x, _ in pts)
    assert box[1] <= min(y for _, y in pts) and box[3] >= max(y for _, y in pts)
    assert box[3] <= 100.0 + G.OVERLAY_LINE + 0.5   # nothing of it below the bond


@pytest.mark.parametrize("value, message", [
    ("1,2", "generic_rings must be a list"),
    ([{"atoms": [1]}], "two atoms"),
    ([{"atoms": [1, 1]}], "two atoms"),
    ([{"atoms": [1, -2]}], "two atoms"),
    ([{"atoms": [1, 2], "members": 2}], "members"),
    ([{"atoms": [1, 2], "colour": "red"}], "unknown keys"),
])
def test_generic_ring_entries_are_checked(value, message):
    errors = figure_house.generic_ring_errors({"id": "I", "generic_rings": value})
    assert any(message in e for e in errors), errors


def test_scheme_and_mechanism_briefs_accept_generic_rings():
    species = {"id": "I", "smiles": "OC(=O)C1CC#C", "generic_rings": [{"atoms": [3, 4]}]}
    scheme = {"species": [species, {"id": "P", "smiles": "CC"}],
              "steps": [{"substrates": ["I"], "products": ["P"]}]}
    assert figure_draw._validate(scheme) == []
    mechanism = {"kind": "mechanism", "species": [species], "steps": [{"substrates": ["I"], "products": ["I"],
                                                                      "moves": []}]}
    assert not [e for e in figure_mechanism._validate(mechanism) if "generic" in e or "unknown keys" in e]


def test_the_renderer_counts_a_generic_ring_in_a_structures_box():
    from cdxml_toolkit.render.renderer import _overlay_box
    atoms = [{"index": 1, "x": A[0], "y": A[1]}, {"index": 2, "x": B[0], "y": B[1]},
             {"index": 3, "x": INSIDE[0], "y": INSIDE[1]}]
    box = _overlay_box(atoms, [{"atoms": [0, 1], "shape": "ring", "away": [2], "members": 6}])
    assert box == pytest.approx(G.generic_ring_box(A, B, INSIDE, members=6))


def test_a_scheme_draws_one_dashed_arc_for_the_ring(tmp_path):
    # 4-pentynoic acid with a ring through C2-C3 standing for the cyclic substrate of the witness.
    brief = {"species": [{"id": "I", "smiles": "OC(=O)CCC#C", "label": "I", "generic_rings": [{"atoms": [3, 4]}]},
                         {"id": "P", "smiles": "OC(=O)CCC=C", "label": "P"}],
             "steps": [{"substrates": ["I"], "products": ["P"], "above": ["[Au]"]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "ring.cdxml"))
    assert report["ok"], report
    root = ET.parse(tmp_path / "ring.cdxml").getroot()
    arcs = [c for c in root.iter("curve") if c.get("LineType") == "Dashed"]
    assert len(arcs) == 1 and arcs[0].get("CurveType") == "2"   # an open arc, dashed
    assert "generic rings (1)" in report["style"]["applied"]


def test_a_ring_through_atoms_that_are_not_bonded_is_refused(tmp_path):
    brief = {"species": [{"id": "I", "smiles": "OC(=O)CCC#C", "generic_rings": [{"atoms": [0, 4]}]},
                         {"id": "P", "smiles": "CC"}],
             "steps": [{"substrates": ["I"], "products": ["P"]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "bad.cdxml"))
    assert not report["ok"] and any("bonded" in e for e in report["errors"])


def _segment_gap(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy or 1.0)))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


def test_an_open_chain_is_laid_out_with_room_for_its_ring(tmp_path):
    # The first render drew 4-pentynoic acid zigzag and the arc ran through the alkyne, with the
    # compound number inside it: a chain is now depicted as if the ring were there.
    brief = {"species": [{"id": "S", "smiles": "OC(=O)CCC#C", "label": "1", "generic_rings": [{"atoms": [3, 4]}]},
                         {"id": "P", "smiles": "C=C1CCC(=O)O1", "label": "2"}],
             "steps": [{"substrates": ["S"], "products": ["P"]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "chain.cdxml"))
    assert report["ok"], report
    root = ET.parse(tmp_path / "chain.cdxml").getroot()
    arc = next(c for c in root.iter("curve") if c.get("LineType") == "Dashed")
    pts = _pairs(arc.get("CurvePoints"))[1::3]
    frag = min(root.iter("fragment"), key=lambda f: min(float(n.get("p").split()[0]) for n in f.iter("n")))
    nodes = {n.get("id"): tuple(float(v) for v in n.get("p").split()[:2]) for n in frag.iter("n")}
    bonds = [(nodes[b.get("B")], nodes[b.get("E")]) for b in frag.iter("b")]
    inner = pts[2:-2]   # leave the ends, which sit on the two ring atoms
    assert min(_segment_gap(p, a, b) for p in inner for a, b in bonds) > 4.0
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    for t in root.iter("t"):
        if "".join(t.itertext()).strip() == "1":
            x, y = (float(v) for v in t.get("p").split()[:2])
            assert not (min(xs) < x < max(xs) and min(ys) < y < max(ys))
    assert not [w for w in report["warnings"] if "overlay" in w.lower()]


def test_the_arc_goes_to_the_side_clear_of_the_molecules_other_bonds():
    # A zigzag chain whose neighbours put the "away" side right across another bond: the arc
    # is drawn on whichever side keeps clear of the molecule, by distance, not by direction.
    a, b = (100.0, 100.0), (114.4, 100.0)
    blocking = [((100.0, 112.0), (114.4, 112.0)), ((96.0, 106.0), (118.0, 120.0))]   # bonds below the pair
    points, clearance, members = G.choose_generic_ring(a, b, away=(107.2, 90.0), members=6,
                                                       bonds=blocking + [(a, b)])
    anchors = _pairs(points)[1::3]
    assert all(y < 100.0 + 0.5 for _, y in anchors)       # it went up, away from the bonds
    assert clearance > 3.0 and members == 6


def test_a_ring_that_fits_nowhere_is_made_smaller_then_reported():
    a, b = (100.0, 100.0), (114.4, 100.0)
    cage = [((80.0, 85.0), (135.0, 85.0)), ((80.0, 115.0), (135.0, 115.0))]   # bonds 15 pt above and below
    points, clearance, members = G.choose_generic_ring(a, b, away=(107.2, 90.0), members=6,
                                                       bonds=cage + [(a, b)])
    assert members < 6 or clearance < 3.0


def _linear_page(overlays):
    from cdxml_toolkit.render import renderer as R
    from cdxml_toolkit.render.schema import ArrowContent, SchemeDescriptor, StepDescriptor, StructureRef
    structures = {"S": StructureRef(id="S", smiles="OC(=O)CCC#C", label="1", overlays=overlays),
                  "P": StructureRef(id="P", smiles="CCO", label="2")}
    steps = [StepDescriptor(substrates=["S"], products=["P"], above_arrow=ArrowContent(text=["[Au]"]))]
    xml, _ = R._layout_linear(SchemeDescriptor(structures=structures, steps=steps, layout="linear"), R._IDGen())
    return ET.fromstring(f"<page>{xml}</page>")


def test_a_linear_layout_keeps_the_arrow_and_the_number_off_a_generic_ring():
    # The open-chain render: the arc cleared the bonds, but the step arrow ran into it and the
    # number sat inside it, because only a cycle made room for it. Every layout does now: the
    # arc's box is part of the structure's.
    ring = ({"atoms": [1, 3], "shape": "ring", "away": [0, 2, 4], "members": 6},)
    page = _linear_page(ring)
    frag = page.findall("fragment")[0]
    where = {n.get("id"): tuple(float(v) for v in n.get("p").split()) for n in frag.findall("n")}
    nodes = list(where.values())
    rest = [nodes[0], nodes[2], nodes[4]]
    away = (sum(p[0] for p in rest) / 3, sum(p[1] for p in rest) / 3)
    bonds = [(where[b.get("B")], where[b.get("E")]) for b in frag.findall("b")]
    points, _, _ = G.choose_generic_ring(nodes[1], nodes[3], away, 6, bonds)
    xs = [p[0] for p in _pairs(points)]; ys = [p[1] for p in _pairs(points)]
    arrow = page.find("arrow")
    tail_x = min(float(arrow.get("Tail3D").split()[0]), float(arrow.get("Head3D").split()[0]))
    arrow_y = float(arrow.get("Tail3D").split()[1])
    assert tail_x > max(xs) or not (min(ys) <= arrow_y <= max(ys))   # the arrow starts past the arc
    label = next(t for t in page.findall("t") if "".join(t.itertext()).strip() == "1")
    lx, ly = (float(v) for v in label.get("p").split()[:2])
    # The number's glyphs (8.4 pt above its baseline) stay out of the arc's box.
    assert not (min(xs) - 4 < lx < max(xs) + 4 and ly - 8.4 < max(ys) and ly > min(ys))
