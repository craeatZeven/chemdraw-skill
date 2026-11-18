"""Steric clash marks: a short red double zigzag across the gap between two groups that crowd
each other (BJOC 22-125 scheme 4: between Ar and the ring in V-1, the ring and [Rh] in V-2/V-3)."""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit import clash_marks as C
from cdxml_toolkit.mcp_runtime import figure_mechanism as M


def _proj(point, origin, axis):
    return (point[0] - origin[0]) * axis[0] + (point[1] - origin[1]) * axis[1]


def test_a_mark_is_two_zigzags_across_the_line_between_the_groups():
    a, b = (0.0, 0.0), (30.0, 0.0)
    lines, gap = C.clash_mark(a, b, 6.5, 3.5)
    assert gap == pytest.approx(30.0 - 6.5 - 3.5)
    assert len(lines) == 2 and all(len(line) == C.TEETH + 1 for line in lines)
    middle = 6.5 + gap / 2   # the middle of the clear gap, not of the two atoms
    for line in lines:
        xs = [p[0] for p in line]
        ys = [p[1] for p in line]
        # Each runs across the a-b line, centred on it ...
        assert max(ys) - min(ys) == pytest.approx(C.LENGTH) and (max(ys) + min(ys)) / 2 == pytest.approx(0.0)
        # ... and zigzags along it, within the mark's thickness around the gap's middle.
        assert all(abs(x - middle) <= C.SEPARATION / 2 + C.AMPLITUDE + 1e-9 for x in xs)
        assert len({round(x, 6) for x in xs}) == 2
    centres = sorted(sum(p[0] for p in line) / len(line) for line in lines)
    assert centres[1] - centres[0] == pytest.approx(C.SEPARATION, abs=C.AMPLITUDE)


def test_a_mark_turns_with_the_groups():
    a, b = (0.0, 0.0), (20.0, 20.0)
    lines, _ = C.clash_mark(a, b, 3.5, 3.5)
    axis = (1 / math.sqrt(2), 1 / math.sqrt(2))
    across = (-axis[1], axis[0])
    for line in lines:
        spread = [_proj(p, a, across) for p in line]
        assert max(spread) - min(spread) == pytest.approx(C.LENGTH)


def test_a_gap_too_narrow_for_the_mark_is_reported():
    _, gap = C.clash_mark((0.0, 0.0), (12.0, 0.0), 6.5, 3.5)
    assert gap < C.MIN_GAP


def test_the_curve_points_draw_straight_segments():
    line = [(0.0, 0.0), (1.0, 2.0), (0.0, 4.0)]
    values = [float(v) for v in C.curve_points(line).split()]
    triples = [tuple(values[i:i + 6]) for i in range(0, len(values), 6)]
    # One (control in, anchor, control out) triple per point, all three on the point.
    assert triples == [(x, y, x, y, x, y) for x, y in line]


def _errors(species, steps=None):
    brief = {"kind": "mechanism", "species": species,
             "steps": steps or [{"substrates": [species[0]["id"]], "products": [species[-1]["id"]]}]}
    return M._validate(brief)


def test_clashes_are_checked():
    ok = [{"id": "a", "smiles": "CC(C)(C)c1ccccc1C", "clashes": [{"between": [1, 9]}]},
          {"id": "n", "newman": {"front": ["Me", "H", "Ph"], "rear": ["H", "Me", "H"]},
           "clashes": [{"between": ["front:0", "rear:1"], "color": "#FF0000"}]}]
    assert not any("clash" in e for e in _errors(ok))
    for bad in ([{"between": [1]}], [{"between": [1, 1]}], [{"between": [1, 2], "colour": "red"}],
                "1-9", [{"between": ["front:0", 2]}], [{"between": [1, 2], "color": "red"}]):
        species = [{**ok[0], "clashes": bad}]
        assert any("clash" in e for e in _errors(species)), bad
    for bad in ([{"between": ["front:3", "rear:0"]}], [{"between": [0, 1]}], [{"between": ["side:0", "rear:0"]}]):
        species = [{**ok[1], "clashes": bad}]
        assert any("clash" in e for e in _errors(species)), bad


def test_newman_positions_name_the_drawn_atoms():
    spec = {"front": ["Me", "H", "OH"], "rear": ["H", "Me", "H"]}
    frag = M._newman_frag("n", spec, None)
    front0 = M._clash_atom(frag, "front:0")
    rear1 = M._clash_atom(frag, "rear:1")
    assert front0["symbol"] == "Me" and front0["role"] == "front_substituent"
    assert rear1["symbol"] == "Me" and rear1["role"] == "rear_substituent"


def test_marks_are_drawn_in_their_colour(tmp_path):
    # Needs RDKit. tBu and the ortho methyl of 2-tert-butyltoluene crowd each other.
    brief = {"kind": "mechanism",
             "species": [{"id": "a", "smiles": "CC(C)(C)c1ccccc1C", "label": "1", "clashes": [{"between": [1, 10]}]},
                         {"id": "b", "smiles": "CC(C)(C)c1ccccc1C", "label": "1'"}],
             "steps": [{"substrates": ["a"], "products": ["b"], "arrow": "resonance"}]}
    from cdxml_toolkit.mcp_runtime import figure_draw
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "c.cdxml"))
    assert report["ok"] and report["clash_marks"] == 1, report
    root = ET.parse(tmp_path / "c.cdxml").getroot()
    table = [tuple(round(float(c.get(k)) * 255) for k in "rgb") for c in root.find("colortable").findall("color")]
    marks = [c for c in root.iter("curve") if c.get("color") and table[int(c.get("color")) - 2] == (255, 0, 0)]
    assert len(marks) == 2
