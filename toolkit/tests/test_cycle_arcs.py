"""A catalytic cycle drawn as a ring of arcs, each arrow on the circle the species sit on (P16).

Witness: BJOC 22-125 scheme 4 and 22-130 scheme 5 draw a cycle's steps as arcs of one circle;
straight chords make a three- or four-step cycle read as a polygon.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.cdxml_check import _curved_arrows
from cdxml_toolkit.render import renderer as R
from cdxml_toolkit.render.schema import ArrowContent, SchemeDescriptor, StepDescriptor, StructureRef


def _cycle(styles=("solid", "solid", "solid", "solid")):
    ids = ["A", "B", "C", "D"]
    smiles = ["CCO", "CCN", "CCCl", "CCBr"]
    structures = {i: StructureRef(id=i, smiles=s, label=i) for i, s in zip(ids, smiles)}
    steps = [StepDescriptor(substrates=[ids[k]], products=[ids[(k + 1) % 4]], arrow_style=styles[k],
                            above_arrow=ArrowContent(text=[f"step {k + 1}"])) for k in range(4)]
    xml, _ = R._layout_catalytic_cycle(SchemeDescriptor(structures=structures, steps=steps, layout="catalytic-cycle"),
                                       R._IDGen())
    return ET.fromstring(f"<page>{xml}</page>")


def _point(value):
    return tuple(float(v) for v in value.split()[:2])


def test_every_step_is_an_arc_of_one_circle():
    page = _cycle()
    arcs = [a for a in page.findall("arrow") if a.get("AngularSize")]
    assert len(arcs) == 4 and not [a for a in page.findall("arrow") if not a.get("AngularSize")]
    centres = {a.get("Center3D") for a in arcs}
    assert len(centres) == 1
    centre = _point(centres.pop())
    radii = {round(math.dist(centre, _point(a.get(end))), 1) for a in arcs for end in ("Tail3D", "Head3D")}
    assert len(radii) == 1
    frags = page.findall("fragment")
    mids = []
    for f in frags:   # each species is placed by the middle of its atoms' extent
        pts = [_point(n.get("p")) for n in f.iter("n")]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        mids.append(((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2))
    # The species sit on that circle too.
    # Each species is placed by the middle of its box, which a wide label (Br, OH) moves off the
    # middle of its atoms by a few points: 5.5 was measured on a native render.
    assert all(abs(math.dist(centre, m) - radii.copy().pop()) < 8.0 for m in mids)


def test_each_arc_points_at_the_next_species_and_starts_clear_of_both():
    page = _cycle()
    frags = page.findall("fragment")
    boxes = [tuple(float(v) for v in f.get("BoundingBox").split()) for f in frags]
    root = ET.fromstring(f"<CDXML>{ET.tostring(page, encoding='unicode')}</CDXML>")
    paths = [pts for _, pts in _curved_arrows(root)]
    assert len(paths) == 4
    for k, pts in enumerate(paths):
        here, there = boxes[k], boxes[(k + 1) % 4]
        centre = lambda b: ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
        assert math.dist(pts[-1], centre(there)) < math.dist(pts[0], centre(there))   # the head is at the next one
        for x, y in pts:
            for b in (here, there):
                assert not (b[0] < x < b[2] and b[1] < y < b[3])   # never inside a structure's box


def _radius(arc):
    return math.dist(_point(arc.get("Center3D")), _point(arc.get("Tail3D")))


def test_an_equilibrium_step_is_a_pair_of_concentric_half_arcs():
    # A reversible step of a cycle: two arcs on the cycle's circle, one each way, as a straight
    # equilibrium arrow is two half-headed shafts (BJOC 22-125 scheme 4, I <=> II <=> III).
    page = _cycle(styles=("solid", "reversible", "solid", "solid"))
    arrows = page.findall("arrow")
    assert all(a.get("AngularSize") for a in arrows) and len(arrows) == 5
    single = [a for a in arrows if "Full" in (a.get("ArrowheadHead"), a.get("ArrowheadTail"))]
    pair = [a for a in arrows if a not in single]
    assert len(single) == 3 and len(pair) == 2
    assert {a.get("Center3D") for a in arrows} == {single[0].get("Center3D")}
    radius = _radius(single[0])
    outer, inner = sorted(pair, key=_radius, reverse=True)
    assert _radius(outer) - radius == pytest.approx(R.CYCLE_EQUILIBRIUM_GAP / 2, abs=0.02)
    assert radius - _radius(inner) == pytest.approx(R.CYCLE_EQUILIBRIUM_GAP / 2, abs=0.02)
    # Same span, opposite ways: the outer runs on to C, the inner back to B.
    assert float(outer.get("AngularSize")) == pytest.approx(float(inner.get("AngularSize")), abs=0.01)
    assert (outer.get("ArrowheadHead"), outer.get("ArrowheadTail")) == ("None", R.CYCLE_EQUILIBRIUM_BARB)
    assert (inner.get("ArrowheadHead"), inner.get("ArrowheadTail")) == (R.CYCLE_EQUILIBRIUM_BARB, "None")


def test_the_pair_keeps_the_p1_gap_from_every_structure():
    page = _cycle(styles=("solid", "reversible", "solid", "solid"))
    root = ET.fromstring(f"<CDXML>{ET.tostring(page, encoding='unicode')}</CDXML>")
    atoms = [_point(n.get("p")) for n in page.iter("n")]
    for _, pts in _curved_arrows(root):
        assert min(math.dist(p, a) for p in pts for a in atoms) >= 3.0


def test_the_pair_arcs_are_written_from_one_span():
    ids = R._IDGen()
    parts, outer, inner, forward_id = R._build_cycle_equilibrium(100.0, 100.0, 50.0, -1.0, -0.2, ids)
    assert len(parts) == 2
    a, b = (ET.fromstring(x) for x in parts)
    # The outer arc runs clockwise on the page from -1.0 to -0.2: its head is the -0.2 end.
    assert math.dist(outer[-1], (100 + 51.5 * math.cos(-0.2), 100 + 51.5 * math.sin(-0.2))) < 0.05
    # The inner arc runs back: its head is the -1.0 end, on the smaller circle.
    assert math.dist(inner[-1], (100 + 48.5 * math.cos(-1.0), 100 + 48.5 * math.sin(-1.0))) < 0.05
    assert a.get("ArrowheadTail") == R.CYCLE_EQUILIBRIUM_BARB and a.get("ArrowheadHead") == "None"
    assert b.get("ArrowheadHead") == R.CYCLE_EQUILIBRIUM_BARB and b.get("ArrowheadTail") == "None"
    # ChemDraw draws an arc from Tail3D anticlockwise to Head3D; for the inner arc that path
    # starts at -0.2 and its head is where it ends.
    assert math.dist(_point(b.get("Head3D")), inner[-1]) < 0.05
    assert math.dist(_point(b.get("Tail3D")), inner[0]) < 0.05
    assert a.get("Center3D") == b.get("Center3D") == "100.00 100.00 0"
    assert a.get("id") == str(forward_id)


def test_the_step_labels_stand_outside_the_ring():
    page = _cycle()
    arc = next(a for a in page.findall("arrow") if a.get("AngularSize"))
    centre = _point(arc.get("Center3D"))
    radius = math.dist(centre, _point(arc.get("Tail3D")))
    for t in page.findall("t"):
        text = "".join(t.itertext())
        if text.startswith("step"):
            assert math.dist(centre, _point(t.get("p"))) > radius
