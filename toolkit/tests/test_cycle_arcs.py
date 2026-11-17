"""A catalytic cycle drawn as a ring of arcs, each arrow on the circle the species sit on (P16).

Witness: BJOC 22-125 scheme 4 and 22-130 scheme 5 draw a cycle's steps as arcs of one circle;
straight chords make a three- or four-step cycle read as a polygon.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET

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
    assert all(abs(math.dist(centre, m) - radii.copy().pop()) < 5.0 for m in mids)  # a label widens a box


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


def test_an_equilibrium_step_keeps_its_straight_double_arrow():
    page = _cycle(styles=("solid", "reversible", "solid", "solid"))
    straight = [a for a in page.findall("arrow") if not a.get("AngularSize")]
    assert len(straight) == 1 and straight[0].get("ArrowheadHead") == "HalfLeft"
    assert len([a for a in page.findall("arrow") if a.get("AngularSize")]) == 3


def test_the_step_labels_stand_outside_the_ring():
    page = _cycle()
    arc = next(a for a in page.findall("arrow") if a.get("AngularSize"))
    centre = _point(arc.get("Center3D"))
    radius = math.dist(centre, _point(arc.get("Tail3D")))
    for t in page.findall("t"):
        text = "".join(t.itertext())
        if text.startswith("step"):
            assert math.dist(centre, _point(t.get("p"))) > radius
