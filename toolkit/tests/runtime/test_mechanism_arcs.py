"""Curved electron arrows as ChemDraw writes them: an <arrow> with AngularSize on a circle.

The convention is read from toolkit/tests/data/native-arc-arrows.cdxml, saved by ChemDraw:
the arc runs from Tail3D anticlockwise on the page (y down) through AngularSize degrees to
Head3D, around Center3D, with MajorAxisEnd3D = centre + (r, 0) and MinorAxisEnd3D = centre +
(0, -r). An arrow that turns clockwise is the same arc with its ends swapped and the head on
the tail end.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from cdxml_toolkit.mcp_runtime import figure_mechanism as M

NATIVE = Path(__file__).resolve().parents[1] / "data" / "native-arc-arrows.cdxml"


def _xy(value):
    x, y, *_ = map(float, value.split())
    return x, y


def _walk(arrow, count=60):
    """Points of an <arrow> arc, following the native convention from Tail3D."""
    cx, cy = _xy(arrow.get("Center3D"))
    tx, ty = _xy(arrow.get("Tail3D"))
    r = math.hypot(tx - cx, ty - cy)
    start = math.atan2(-(ty - cy), tx - cx)  # angle with y up: anticlockwise as seen on the page
    sweep = math.radians(float(arrow.get("AngularSize")))
    return [(cx + r * math.cos(start + sweep * k / count), cy - r * math.sin(start + sweep * k / count))
            for k in range(count + 1)]


def test_the_convention_reproduces_every_arc_chemdraw_saved():
    text = re.sub(r"<!DOCTYPE[^>]*>", "", NATIVE.read_text(encoding="utf-8"))
    arrows = [a for a in ET.fromstring(text.encode()).iter("arrow") if a.get("AngularSize")]
    assert len(arrows) == 7
    for arrow in arrows:
        end = _walk(arrow)[-1]
        assert math.dist(end, _xy(arrow.get("Head3D"))) < 0.05
        cx, cy = _xy(arrow.get("Center3D"))
        r = math.dist((cx, cy), _xy(arrow.get("Tail3D")))
        assert math.dist(_xy(arrow.get("MajorAxisEnd3D")), (cx + r, cy)) < 0.05
        assert math.dist(_xy(arrow.get("MinorAxisEnd3D")), (cx, cy - r)) < 0.05


@pytest.mark.parametrize("side", [1.0, -1.0])
@pytest.mark.parametrize("dst", [(60.0, 0.0), (0.0, 40.0), (-30.0, -25.0)])
def test_an_arc_leaves_its_source_and_reaches_its_target_on_the_side_asked(side, dst):
    arc = M._arc((0.0, 0.0), 0.0, dst, 0.0, side, 1.0)
    pts = arc["pts"]
    assert math.dist(pts[0], (0.0, 0.0)) < 1e-6 and math.dist(pts[-1], dst) < 1e-6
    # the middle of the arc lies on the side of the chord that ``side`` names (as _bezier did)
    ux, uy = dst[0] / math.hypot(*dst), dst[1] / math.hypot(*dst)
    mid = pts[len(pts) // 2]
    assert (mid[0] * -uy + mid[1] * ux) * side > 0


def test_the_written_arrow_follows_the_native_convention_whichever_way_it_turns():
    for side in (1.0, -1.0):
        arc = M._arc((10.0, 50.0), 0.0, (70.0, 50.0), 0.0, side, 1.0)
        xml = M._arc_xml({"arc": arc, "electrons": 2, "pts": arc["pts"]}, iter(range(100, 200)).__next__)
        arrow = ET.fromstring(xml)
        walked = _walk(arrow, count=2)
        # the drawn arc is the routed one: same ends, same middle
        assert math.dist(walked[1], arc["pts"][len(arc["pts"]) // 2]) < 0.05
        assert math.dist(walked[-1], _xy(arrow.get("Head3D"))) < 0.05
        head_end = "Head3D" if arrow.get("ArrowheadHead", "None") != "None" else "Tail3D"
        assert math.dist(_xy(arrow.get(head_end)), (70.0, 50.0)) < 0.05  # the head is on the target


def test_a_one_electron_move_has_a_half_head():
    arc = M._arc((0.0, 0.0), 0.0, (50.0, 0.0), 0.0, 1.0, 1.0)
    arrow = ET.fromstring(M._arc_xml({"arc": arc, "electrons": 1, "pts": arc["pts"]}, iter(range(9)).__next__))
    heads = {arrow.get("ArrowheadHead", "None"), arrow.get("ArrowheadTail", "None")} - {"None"}
    assert heads == {"HalfLeft"}


def test_an_arc_is_scored_as_crossing_a_bond_it_cuts_and_not_one_it_clears():
    arc = M._arc((0.0, 0.0), 0.0, (60.0, 0.0), 0.0, 1.0, 1.0)
    apex = arc["pts"][len(arc["pts"]) // 2]
    across = {"key": "b", "label": "a bond", "a": (apex[0], apex[1] - 10), "b": (apex[0], apex[1] + 10),
              "atoms": {("x", 1), ("x", 2)}}
    clear = dict(across, a=(apex[0] + 120, apex[1] - 10), b=(apex[0] + 120, apex[1] + 10))
    assert M._bond_hits(arc["pts"], [across], set()) == [across]
    assert M._bond_hits(arc["pts"], [clear], set()) == []
    # a bond at one of the arrow's own atoms is never in its way
    assert M._bond_hits(arc["pts"], [across], {("x", 1)}) == []


def test_a_move_colour_is_checked():
    brief = {"kind": "mechanism", "species": [{"id": "a", "smiles": "C"}],
             "steps": [{"substrates": ["a"], "products": ["a"],
                        "moves": [{"from": {"species": "a", "atom": 0}, "to": {"species": "a", "atom": 0},
                                   "color": "green"}]}]}
    assert any("color" in e for e in M._validate(brief))
    brief["steps"][0]["moves"][0]["color"] = "#2E7D32"
    assert not any("color" in e for e in M._validate(brief))
