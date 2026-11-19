"""The page composer: rows of scientific drawings with labels, arrows and curved arrows."""
import json
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.scientific.lewis import lewis_structure
from cdxml_toolkit.scientific.page import compose_page, runs


def _water(tmp_path, name="water.cdxml"):
    spec = {"style": "lines", "atoms": [{"id": "H1", "element": "H", "at": [0, 0]},
                                       {"id": "O", "element": "O", "at": [1, 0]},
                                       {"id": "H2", "element": "H", "at": [2, 0]}],
            "bonds": [["H1", "O"], ["O", "H2"]]}
    target = tmp_path / name
    lewis_structure(spec, output=str(target))
    return str(target)


def test_a_label_subscripts_its_atom_counts():
    assert runs("CH3NO2") == [("CH", "0"), ("3", "32"), ("NO", "0"), ("2", "32")]


def test_a_caret_charge_becomes_a_superscript():
    assert ("2-", "64") in runs("CO3^2-")


def test_two_parts_sit_side_by_side_on_one_page(tmp_path):
    part = _water(tmp_path)
    result = compose_page({"rows": [{"parts": [part, part], "gap": 40.0, "labels": ["(a)", "(b)"]}],
                           "out": str(tmp_path / "page.cdxml")},
                          output=str(tmp_path / "page.cdxml"))
    boxes = result["rows"][0]["panels"]
    assert len(boxes) == 2 and boxes[1][0] > boxes[0][2], "the second part starts right of the first"
    page = ET.fromstring(result["cdxml"].split("?>", 1)[-1]).find("page")
    texts = ["".join(s.text or "" for s in t.iter("s")) for t in page.iter("t")]
    assert "(a)" in texts and "(b)" in texts


def test_a_resonance_arrow_sits_on_the_bond_line(tmp_path):
    part = _water(tmp_path)
    result = compose_page({"rows": [{"parts": [part, part], "gap": 60.0, "arrows": [[0, 1]]}]})
    page = ET.fromstring(result["cdxml"].split("?>", 1)[-1]).find("page")
    arrows = list(page.iter("arrow"))
    assert len(arrows) == 1 and arrows[0].get("ArrowheadHead") == "Full" and arrows[0].get("ArrowheadTail") == "Full"
    lines = [(float(g.get("BoundingBox").split()[1]), float(g.get("BoundingBox").split()[3]))
             for g in page.iter("graphic") if g.get("BoundingBox")
             and abs(float(g.get("BoundingBox").split()[3]) - float(g.get("BoundingBox").split()[1])) < 1.5]
    arrow_y = float(arrows[0].get("Tail3D").split()[1])
    assert any(abs((y0 + y1) / 2 - arrow_y) < 0.5 for y0, y1 in lines), "the arrow is not on a bond line"


def test_an_equilibrium_sign_is_two_half_arrows_eight_points_apart(tmp_path):
    part = _water(tmp_path)
    result = compose_page({"rows": [{"parts": [part, part], "gap": 60.0,
                                     "arrows": [["equilibrium", 0, 1]]}]})
    page = ET.fromstring(result["cdxml"].split("?>", 1)[-1]).find("page")
    ys = sorted(float(a.get("Tail3D").split()[1]) for a in page.iter("arrow"))
    assert len(ys) == 2 and abs((ys[1] - ys[0]) - 8.0) < 0.1
    heads = sorted(a.get("ArrowheadHead") for a in page.iter("arrow"))
    assert heads == ["HalfLeft", "None"]


def test_a_row_can_carry_a_margin_label_and_a_curved_arrow(tmp_path):
    part = _water(tmp_path)
    result = compose_page({"rows": [{"parts": [part], "margin_label": "(a)",
                                     "curves": [{"centre": [10.0, 10.0], "radius": 20.0, "from": 100.0, "to": 60.0}]}]})
    page = ET.fromstring(result["cdxml"].split("?>", 1)[-1]).find("page")
    texts = ["".join(s.text or "" for s in t.iter("s")) for t in page.iter("t")]
    assert "(a)" in texts
    curved = [a for a in page.iter("arrow") if a.get("AngularSize")]
    assert curved and curved[0].get("ArrowheadHead") == "Full"


def test_a_gap_too_small_for_an_arrow_is_reported(tmp_path):
    part = _water(tmp_path)
    result = compose_page({"rows": [{"parts": [part, part], "gap": 10.0, "arrows": [[0, 1]]}]})
    assert result["warnings"], "a gap that leaves no room must warn instead of drawing a backwards arrow"


def test_an_unknown_arrow_kind_is_refused(tmp_path):
    part = _water(tmp_path)
    with pytest.raises(ValueError, match="resonance"):
        compose_page({"rows": [{"parts": [part, part], "gap": 60.0, "arrows": [["dotted", 0, 1]]}]})
