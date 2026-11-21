"""CDXML geometry and phase contracts for hybrid-orbital and overlap drawings."""
from __future__ import annotations

import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.scientific.orbital_lobes import lobe, orbital_overlap, orbital_shapes, p_orbital
from cdxml_toolkit.scientific.page import compose_page


SKILL_SCRIPTS = Path(__file__).resolve().parents[2] / "skill" / "chemdraw" / "scripts"
if str(SKILL_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SKILL_SCRIPTS))
import acs_builder  # noqa: E402


def _root(result):
    return ET.fromstring(result["cdxml"])


def _point(node, attribute):
    return tuple(float(value) for value in node.get(attribute).split())


def _lobes(root):
    return [node for node in root.iter("graphic")
            if node.get("GraphicType") == "Oval" and node.get("OvalType") in (None, "Shaded")]


def _axis_angle(node):
    center = _point(node, "Center3D")
    major = _point(node, "MajorAxisEnd3D")
    return math.degrees(math.atan2(major[1] - center[1], major[0] - center[0]))


def _angle_between(a, b):
    ca, cb = _point(a, "Center3D"), _point(b, "Center3D")
    va = [m - c for m, c in zip(_point(a, "MajorAxisEnd3D"), ca)]
    vb = [m - c for m, c in zip(_point(b, "MajorAxisEnd3D"), cb)]
    cosine = sum(x * y for x, y in zip(va, vb)) / (math.sqrt(sum(x * x for x in va)) * math.sqrt(sum(x * x for x in vb)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _colors(root):
    table = {0: "#000000", 1: "#FFFFFF"}
    for index, color in enumerate(root.findall("./colortable/color"), start=2):
        rgb = [round(float(color.get(channel)) * 255) for channel in "rgb"]
        table[index] = "#" + "".join(f"{value:02X}" for value in rgb)
    return table


def test_ellipse_lobe_direction_coefficient_fill_and_phase_colors():
    from cdxml_toolkit.scientific.canvas import Canvas

    canvas = Canvas()
    lobe(canvas, (50, 60), angle_deg=30, size=40, width=16, coefficient=0.5, phase=1)
    lobe(canvas, (100, 60), angle_deg=-90, size=30, width=12, phase=-1, fill="hollow")
    p_orbital(canvas, (150, 60), size=30, width=12, phase=1)
    root = _root({"cdxml": canvas.to_string()})
    ellipses = _lobes(root)
    assert len(ellipses) == 4
    assert _axis_angle(ellipses[0]) == pytest.approx(30, abs=0.01)
    assert _point(ellipses[0], "Center3D") == pytest.approx((58.66, 65.0, 0))
    assert math.dist(_point(ellipses[0], "Center3D")[:2], _point(ellipses[0], "MajorAxisEnd3D")[:2]) == pytest.approx(10, abs=0.001)
    assert ellipses[0].get("OvalType") == "Shaded"
    assert ellipses[1].get("OvalType") is None
    colors = _colors(root)
    assert colors[int(ellipses[0].get("color"))] == "#C0392B"
    assert colors[int(ellipses[1].get("color"))] == "#0072B2"
    assert [colors[int(node.get("color"))] for node in ellipses[2:]] == ["#C0392B", "#0072B2"]
    assert [_axis_angle(node) for node in ellipses[2:]] == pytest.approx([-90, 90])


@pytest.mark.parametrize("hybridisation,count,expected_angle", [
    ("sp3", 4, 109.5),
    ("sp2", 5, 120.0),
    ("sp", 6, 180.0),
    ("unequal", 4, 109.5),
])
def test_orbital_shapes_panel_counts_angles_axes_and_arcs(hybridisation, count, expected_angle):
    result = orbital_shapes({"hybridisation": hybridisation})
    root = _root(result)
    ellipses = _lobes(root)
    assert len(ellipses) == count
    hybrids = ellipses[:4 if hybridisation in ("sp3", "unequal") else 3 if hybridisation == "sp2" else 2]
    angles = [_angle_between(a, b) for i, a in enumerate(hybrids) for b in hybrids[i + 1:]]
    assert all(angle == pytest.approx(expected_angle, abs=1.0) for angle in angles)
    arcs = [node for node in root.iter("arrow") if node.get("AngularSize") is not None]
    assert len(arcs) == 1
    assert abs(float(arcs[0].get("AngularSize"))) == pytest.approx(expected_angle, abs=1.0)
    assert any(str(expected_angle) in "".join(node.itertext()) for node in root.iter("t"))
    assert len([node for node in root.iter("graphic") if node.get("GraphicType") == "Line"]) >= len(hybrids)
    assert check_cdxml(result["cdxml"])["ok"]


def test_sp2_and_sp_include_unhybridised_p_pairs():
    sp2 = _lobes(_root(orbital_shapes({"hybridisation": "sp2"})))[3:]
    sp = _lobes(_root(orbital_shapes({"hybridisation": "sp"})))[2:]
    assert len(sp2) == 2 and _angle_between(*sp2) == pytest.approx(180, abs=1e-6)
    assert len(sp) == 4
    assert _angle_between(sp[0], sp[1]) == pytest.approx(180, abs=1e-6)
    assert _angle_between(sp[2], sp[3]) == pytest.approx(180, abs=1e-6)


def test_unequal_hybridisation_scales_each_lobe_and_labels_ch3cl():
    result = orbital_shapes({"hybridisation": "unequal", "lobe_sizes": [48, 32, 34, 36]})
    root = _root(result)
    lengths = [2 * math.dist(_point(node, "Center3D"), _point(node, "MajorAxisEnd3D")) for node in _lobes(root)]
    assert lengths == pytest.approx([48, 32, 34, 36], abs=0.03)
    text = " ".join("".join(node.itertext()) for node in root.iter("t"))
    assert "CH3Cl" in text and "109.5°" in text


def test_pi_overlap_has_four_lobes_with_matching_upper_and_opposite_lower_phases():
    result = orbital_overlap({"kind": "pi"})
    root = _root(result)
    ellipses = _lobes(root)
    assert len(ellipses) == 4
    axis = next(node for node in root.iter("graphic")
                if node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed")
    axis_box = _point_from_box(axis)
    assert axis_box[1] == axis_box[3]
    axis_y = axis_box[1]
    upper = [node for node in ellipses if _point(node, "Center3D")[1] < axis_y]
    lower = [node for node in ellipses if _point(node, "Center3D")[1] > axis_y]
    assert len(upper) == len(lower) == 2
    colors = _colors(root)
    upper_colors = [colors[int(node.get("color"))] for node in upper]
    lower_colors = [colors[int(node.get("color"))] for node in lower]
    assert upper_colors == ["#C0392B", "#C0392B"]
    assert lower_colors == ["#0072B2", "#0072B2"]
    assert set(upper_colors).isdisjoint(lower_colors)
    overlap = [node for node in root.iter("graphic") if node.get("OvalType") == "Filled"]
    assert len(overlap) == 2
    assert check_cdxml(result["cdxml"])["ok"]


def _point_from_box(node):
    values = [float(value) for value in node.get("BoundingBox").split()]
    return values


def test_sigma_overlap_is_head_on_and_marks_the_overlap_region():
    root = _root(orbital_overlap({"kind": "sigma"}))
    ellipses = _lobes(root)
    assert len(ellipses) == 2
    assert [_axis_angle(node) for node in ellipses] == pytest.approx([0, 180])
    assert len([node for node in root.iter("graphic") if node.get("OvalType") == "Filled"]) == 1
    assert any(node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed" for node in root.iter("graphic"))


def test_overlap_generator_defaults_to_sigma_and_pi_panels():
    result = orbital_overlap({})
    assert [panel["kind"] for panel in result["panels"]] == ["sigma", "pi"]
    assert len(_lobes(_root(result))) == 6


def test_acs_build_exposes_lobes_shapes_and_overlap():
    text = acs_builder.build_from_spec({"objects": [
        {"type": "lobe", "at": [30, 40], "angle": 20, "size": 30, "phase": -1},
        {"type": "orbital_shapes", "at": [100, 80], "hybridisation": "sp3"},
        {"type": "orbital_overlap", "at": [250, 80], "kind": "pi"},
    ]})
    root = ET.fromstring(text)
    assert len(_lobes(root)) == 9
    assert len([node for node in root.iter("graphic") if node.get("OvalType") == "Filled"]) == 2
    assert check_cdxml(text)["ok"]
    with pytest.raises(ValueError, match="unknown key"):
        acs_builder.build_from_spec({"objects": [{"type": "orbital_shapes", "hybridisation": "sp3", "angel": 30}]})


def test_orbital_panels_compose_without_losing_ellipse_geometry(tmp_path):
    shape_path = tmp_path / "shape.cdxml"
    overlap_path = tmp_path / "overlap.cdxml"
    orbital_shapes({"hybridisation": "sp3"}, output=shape_path)
    orbital_overlap({"kind": "pi"}, output=overlap_path)
    result = compose_page({"rows": [{"parts": [str(shape_path), str(overlap_path)]}]})
    root = _root(result)
    assert len(_lobes(root)) == 8
    assert len([node for node in root.iter("graphic") if node.get("OvalType") == "Filled"]) == 2
