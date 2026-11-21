"""CDXML geometry, phase, composition and compatibility contracts for orbital drawings."""
from __future__ import annotations

import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml, _curved_arrows
from cdxml_toolkit.scientific.canvas import caption_box
from cdxml_toolkit.scientific.orbital_lobes import lobe, orbital_overlap, orbital_shapes, p_orbital
from cdxml_toolkit.scientific.page import compose_page


ORBITAL_PALETTE = ["#000000", "#FFFFFF", "#C0392B", "#0072B2"]


SKILL_SCRIPTS = Path(__file__).resolve().parents[2] / "skill" / "chemdraw" / "scripts"
if str(SKILL_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SKILL_SCRIPTS))
import acs_builder  # noqa: E402


def _root(result):
    return ET.fromstring(result["cdxml"])


def _point(node, attribute):
    return tuple(float(value) for value in node.get(attribute).split())


def _box(node):
    return tuple(float(value) for value in node.get("BoundingBox").split())


def _lobes(root):
    return [node for node in root.iter("graphic")
            if node.get("GraphicType") == "Oval" and node.get("OvalType") in (None, "Shaded")]


def _overlaps(root):
    return [node for node in root.iter("graphic") if node.get("OvalType") == "Filled"]


def _axis_vector(node):
    return tuple(m - c for m, c in zip(_point(node, "MajorAxisEnd3D"), _point(node, "Center3D")))


def _axis_angle(node):
    x, y, _ = _axis_vector(node)
    return math.degrees(math.atan2(y, x))


def _angle(a, b, dimensions=3):
    va, vb = _axis_vector(a)[:dimensions], _axis_vector(b)[:dimensions]
    cosine = sum(x * y for x, y in zip(va, vb)) / (math.sqrt(sum(x * x for x in va)) * math.sqrt(sum(x * x for x in vb)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _intersection(a, b):
    aa, bb = _box(a), _box(b)
    return max(aa[0], bb[0]), max(aa[1], bb[1]), min(aa[2], bb[2]), min(aa[3], bb[3])


def _minor_vector(node):
    return tuple(m - c for m, c in zip(_point(node, "MinorAxisEnd3D"), _point(node, "Center3D")))


def _ellipse_value(node, xy):
    center = _point(node, "Center3D")
    major = _axis_vector(node)
    minor = _minor_vector(node)
    major_length = math.hypot(*major)
    minor_length = math.hypot(*minor)
    delta = tuple(xy[i] - center[i] for i in range(2))
    along = sum(delta[i] * major[i] for i in range(2)) / major_length
    across = sum(delta[i] * minor[i] for i in range(2)) / minor_length
    return (along / major_length) ** 2 + (across / minor_length) ** 2


def _assert_nonempty_intersection(a, b, overlap):
    expected = _intersection(a, b)
    assert expected[2] > expected[0] and expected[3] > expected[1]
    box = _box(overlap)
    assert expected[0] - 0.02 <= box[0] <= box[2] <= expected[2] + 0.02
    assert expected[1] - 0.02 <= box[1] <= box[3] <= expected[3] + 0.02
    center = _point(overlap, "Center3D")
    major, minor = _axis_vector(overlap), _minor_vector(overlap)
    for degrees in (0, 45, 90, 135, 180, 225, 270, 315):
        theta = math.radians(degrees)
        xy = tuple(center[i] + major[i] * math.cos(theta) + minor[i] * math.sin(theta)
                   for i in range(2))
        assert _ellipse_value(a, xy) <= 1.0 + 1e-6
        assert _ellipse_value(b, xy) <= 1.0 + 1e-6


def _assert_arc_endpoints(root):
    arcs = [node for node in root.iter("arrow") if node.get("AngularSize") is not None]
    sampled = _curved_arrows(root)
    assert len(sampled) == len(arcs)
    for (_, points), arc in zip(sampled, arcs):
        assert math.dist(points[-1], _point(arc, "Head3D")[:2]) < 0.5


def _markup(node):
    pieces = []
    for run in node.iter("s"):
        face = int(run.get("face", "0")) & 96
        prefix = "{sub:" if face == 32 else "{sup:" if face == 64 else ""
        pieces.append(f"{prefix}{run.text or ''}{'}' if prefix else ''}")
    return "".join(pieces)


def _caption_extent(node):
    x, y = _point(node, "p")[:2]
    size = max(float(run.get("size", 10)) for run in node.iter("s"))
    align = node.get("Justification", "Left")
    return caption_box(x, y, _markup(node), size, align if align in ("Left", "Center", "Right") else "Left")


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
    ("sp3", 4, 109.4712),
    ("sp2", 5, 120.0),
    ("sp", 6, 180.0),
])
def test_orbital_shapes_panel_counts_angles_axes_and_arcs(hybridisation, count, expected_angle):
    result = orbital_shapes({"hybridisation": hybridisation})
    root = _root(result)
    ellipses = _lobes(root)
    assert len(ellipses) == count
    hybrid_count = 4 if hybridisation == "sp3" else 3 if hybridisation == "sp2" else 2
    hybrids = ellipses[:hybrid_count]
    angles = [_angle(a, b) for i, a in enumerate(hybrids) for b in hybrids[i + 1:]]
    assert all(angle == pytest.approx(expected_angle, abs=1.0) for angle in angles)
    arcs = [node for node in root.iter("arrow") if node.get("AngularSize") is not None]
    assert len(arcs) == 1
    projected = _angle(hybrids[0], hybrids[1], dimensions=2)
    assert abs(float(arcs[0].get("AngularSize"))) == pytest.approx(projected, abs=0.02)
    text = " ".join("".join(node.itertext()) for node in root.iter("t"))
    if hybridisation == "sp3":
        assert projected == pytest.approx(90, abs=0.02)
        assert "3D angle 109.5° (projection 90.0°)" in text
    else:
        assert f"{expected_angle:.1f}°" in text
    assert len([node for node in root.iter("graphic") if node.get("GraphicType") == "Line"]) >= len(hybrids)
    _assert_arc_endpoints(root)
    assert check_cdxml(result["cdxml"], palette=ORBITAL_PALETTE)["ok"]


def test_sp2_p_axis_is_perpendicular_to_the_hybrid_plane():
    ellipses = _lobes(_root(orbital_shapes({"hybridisation": "sp2"})))
    hybrids, p_pair = ellipses[:3], ellipses[3:]
    assert len(p_pair) == 2 and _angle(*p_pair) == pytest.approx(180, abs=0.02)
    assert all(_angle(p, hybrid) == pytest.approx(90, abs=0.02) for p in p_pair for hybrid in hybrids)
    assert all(vector[0] == pytest.approx(0) and vector[1] == pytest.approx(0) and abs(vector[2]) > 1
               for vector in map(_axis_vector, p_pair))


def test_sp_p_axes_are_opposite_pairs_and_mutually_perpendicular():
    ellipses = _lobes(_root(orbital_shapes({"hybridisation": "sp"})))
    first, second = ellipses[2:4], ellipses[4:6]
    assert _angle(*first) == pytest.approx(180, abs=0.02)
    assert _angle(*second) == pytest.approx(180, abs=0.02)
    assert _angle(first[0], second[0]) == pytest.approx(90, abs=0.02)
    assert all(vector[2] == pytest.approx(0) for vector in map(_axis_vector, first))
    assert all(vector[0] == pytest.approx(0) and vector[1] == pytest.approx(0) and abs(vector[2]) > 1
               for vector in map(_axis_vector, second))


def test_unequal_hybridisation_has_distinct_hch_and_hccl_angles_and_matching_arcs():
    result = orbital_shapes({"hybridisation": "unequal", "lobe_sizes": [48, 32, 34, 36]})
    root = _root(result)
    ellipses = _lobes(root)
    cl, hydrogens = ellipses[0], ellipses[1:]
    lengths = [2 * math.dist(_point(node, "Center3D"), _point(node, "MajorAxisEnd3D")) for node in ellipses]
    assert lengths == pytest.approx([48, 32, 34, 36], abs=0.03)
    cl_angles = [_angle(cl, hydrogen) for hydrogen in hydrogens]
    assert cl_angles == pytest.approx([108.4] * 3, abs=0.25)
    assert result["angles_deg"]["H-C-Cl"] == pytest.approx(cl_angles[0], abs=0.01)
    assert [_angle(a, b) for i, a in enumerate(hydrogens) for b in hydrogens[i + 1:]] == pytest.approx([110.3] * 3, abs=0.1)
    arcs = [node for node in root.iter("arrow") if node.get("AngularSize") is not None]
    assert len(arcs) == 2
    for arc, a, b in zip(arcs, (cl, hydrogens[0]), (hydrogens[0], hydrogens[1])):
        assert abs(float(arc.get("AngularSize"))) == pytest.approx(_angle(a, b, dimensions=2), abs=0.02)
    text = " ".join("".join(node.itertext()) for node in root.iter("t"))
    assert "CH3Cl" in text and f"H-C-Cl 3D {cl_angles[0]:.1f}°" in text and "H-C-H 3D 110.3°" in text
    _assert_arc_endpoints(root)
    assert check_cdxml(result["cdxml"], palette=ORBITAL_PALETTE)["ok"]


def test_pi_overlap_has_four_lobes_with_matching_upper_and_opposite_lower_phases():
    result = orbital_overlap({"kind": "pi"})
    root = _root(result)
    ellipses = _lobes(root)
    assert len(ellipses) == 4
    axis = next(node for node in root.iter("graphic")
                if node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed")
    axis_y = _box(axis)[1]
    upper = [node for node in ellipses if _point(node, "Center3D")[1] < axis_y]
    lower = [node for node in ellipses if _point(node, "Center3D")[1] > axis_y]
    assert len(upper) == len(lower) == 2
    colors = _colors(root)
    upper_colors = [colors[int(node.get("color"))] for node in upper]
    lower_colors = [colors[int(node.get("color"))] for node in lower]
    assert upper_colors == ["#C0392B", "#C0392B"]
    assert lower_colors == ["#0072B2", "#0072B2"]
    assert set(upper_colors).isdisjoint(lower_colors)
    assert len(_overlaps(root)) == 2
    assert check_cdxml(result["cdxml"], palette=ORBITAL_PALETTE)["ok"]


@pytest.mark.parametrize("size,width,separation", [(64, 22, 16), (80, 30, 24), (50, 28, 14)])
def test_sigma_overlap_region_equals_the_lobe_bbox_intersection(size, width, separation):
    root = _root(orbital_overlap({"kind": "sigma", "size": size, "width": width, "separation": separation}))
    ellipses = _lobes(root)
    assert len(ellipses) == 2
    assert [_axis_angle(node) for node in ellipses] == pytest.approx([0, 180])
    assert _box(ellipses[0]) != _box(ellipses[1])
    _assert_nonempty_intersection(*ellipses, _overlaps(root)[0])
    assert any(node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed" for node in root.iter("graphic"))


@pytest.mark.parametrize("size,width,separation", [(64, 22, 16), (80, 30, 24), (50, 28, 14), (64, 2, 0.5)])
def test_pi_overlap_regions_equal_same_phase_bbox_intersections(size, width, separation):
    result = orbital_overlap({"kind": "pi", "size": size, "width": width, "separation": separation})
    root = _root(result)
    axis = next(node for node in root.iter("graphic")
                if node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed")
    axis_y = _box(axis)[1]
    upper = sorted((node for node in _lobes(root) if _point(node, "Center3D")[1] < axis_y),
                   key=lambda node: _point(node, "Center3D")[0])
    lower = sorted((node for node in _lobes(root) if _point(node, "Center3D")[1] > axis_y),
                   key=lambda node: _point(node, "Center3D")[0])
    overlaps = sorted(_overlaps(root), key=lambda node: _point(node, "Center3D")[1])
    assert len(upper) == len(lower) == len(overlaps) == 2
    _assert_nonempty_intersection(*upper, overlaps[0])
    _assert_nonempty_intersection(*lower, overlaps[1])


@pytest.mark.parametrize("separation", [63.99, 63.99999])
def test_sigma_near_coincident_overlap_is_nondegenerate(separation):
    root = _root(orbital_overlap({"kind": "sigma", "separation": separation}))
    ellipses = _lobes(root)
    overlap = _overlaps(root)[0]
    assert min(math.hypot(*_axis_vector(overlap)), math.hypot(*_minor_vector(overlap))) >= 0.6
    _assert_nonempty_intersection(*ellipses, overlap)


@pytest.mark.parametrize("separation", [0.01, 0.00001])
def test_pi_near_coincident_overlap_is_nondegenerate(separation):
    result = orbital_overlap({"kind": "pi", "separation": separation})
    root = _root(result)
    axis = next(node for node in root.iter("graphic")
                if node.get("GraphicType") == "Line" and node.get("LineType") == "Dashed")
    axis_y = _box(axis)[1]
    upper = sorted((node for node in _lobes(root) if _point(node, "Center3D")[1] < axis_y),
                   key=lambda node: _point(node, "Center3D")[0])
    lower = sorted((node for node in _lobes(root) if _point(node, "Center3D")[1] > axis_y),
                   key=lambda node: _point(node, "Center3D")[0])
    overlaps = sorted(_overlaps(root), key=lambda node: _point(node, "Center3D")[1])
    for sources, overlap in ((upper, overlaps[0]), (lower, overlaps[1])):
        assert min(math.hypot(*_axis_vector(overlap)), math.hypot(*_minor_vector(overlap))) >= 0.6
        _assert_nonempty_intersection(*sources, overlap)


def test_impossible_overlap_separations_are_rejected():
    with pytest.raises(ValueError, match="separation"):
        orbital_overlap({"kind": "sigma", "size": 64, "separation": 64})
    with pytest.raises(ValueError, match="separation"):
        orbital_overlap({"kind": "pi", "width": 22, "separation": 22})


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
    assert len(_overlaps(root)) == 2
    assert check_cdxml(text, palette=ORBITAL_PALETTE)["ok"]
    with pytest.raises(ValueError, match="unknown key"):
        acs_builder.build_from_spec({"objects": [{"type": "orbital_shapes", "hybridisation": "sp3", "angel": 30}]})


def test_orbital_panels_compose_with_valid_colors_and_page_frame(tmp_path):
    shape_path = tmp_path / "shape.cdxml"
    overlap_path = tmp_path / "overlap.cdxml"
    orbital_shapes({"hybridisation": "sp3"}, output=shape_path)
    orbital_overlap({"kind": "pi"}, output=overlap_path)
    result = compose_page({"rows": [{"parts": [str(shape_path), str(overlap_path)]}]})
    root = _root(result)
    assert len(_lobes(root)) == 8
    assert len(_overlaps(root)) == 2
    colors = _colors(root)
    assert {"#C0392B", "#0072B2"} <= set(colors.values())
    for node in root.iter("graphic"):
        if node.get("color"):
            assert int(node.get("color")) in colors
    checked = check_cdxml(result["cdxml"], palette=ORBITAL_PALETTE)
    assert checked["errors"] == []
    assert not [warning for warning in checked["warnings"] if warning["code"] == "outside_page"]
    page = root.find("page")
    width, height = float(page.get("Width")), float(page.get("Height"))
    for node in root.iter():
        for attribute in ("p", "Head3D", "Tail3D", "Center3D"):
            if node.get(attribute):
                x, y = _point(node, attribute)[:2]
                assert 0 <= x <= width and 0 <= y <= height
    for node in root.iter("t"):
        x0, y0, x1, y1 = _caption_extent(node)
        assert 0 <= x0 <= x1 <= width and 0 <= y0 <= y1 <= height


def test_phase_colors_do_not_relax_the_default_text_palette():
    text = acs_builder.build_from_spec({"objects": [
        {"type": "text", "at": [20, 20], "text": "ordinary red text", "color": "#C0392B"}
    ]})
    assert "palette" in {error["code"] for error in check_cdxml(text)["errors"]}
    third = acs_builder.build_from_spec({"objects": [
        {"type": "text", "at": [20, 20], "text": "unlisted colour", "color": "#123456"}
    ]})
    assert "palette" in {error["code"] for error in check_cdxml(third, palette=ORBITAL_PALETTE)["errors"]}


@pytest.mark.parametrize("spec,expected", [
    ({"objects": [{"type": "orbital", "at": [205.93, 69.76], "kind": "p", "shaded": "shaded",
                   "angle_deg": -90, "length": 6.8, "width": 2.3, "color": "#0072B2"}]},
     '<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">\n<CDXML BoundingBox="197.13 58.66 214.73 78.56" BondLength="14.4" LabelFont="3" LabelSize="10" LabelFace="96" CaptionFont="3" CaptionSize="10" LineWidth="0.6" BoldWidth="2" HashSpacing="2.5" MarginWidth="1.6" BondSpacing="18" ChainAngle="120" LabelJustification="Auto" CaptionJustification="Left" color="0" bgcolor="1"><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/><color r="0" g="0.4471" b="0.698"/></colortable><fonttable><font id="3" charset="iso-8859-1" name="Arial"/><font id="5" charset="Unicode" name="SimSun"/></fonttable><page id="1">\n<graphic id="1001" GraphicType="Orbital" OrbitalType="p" BoundingBox="205.93 69.76 205.93 62.96" Center3D="205.93 69.76 0" MajorAxisEnd3D="205.93 62.96 0" MinorAxisEnd3D="208.23 69.76 0" color="4"/>\n</page></CDXML>\n'),
    ({"objects": [{"type": "orbital", "at": [213.13, 82.24], "kind": "p", "shaded": "shaded",
                   "angle_deg": -90, "length": 6.0, "width": 2.0, "color": "#336699"}]},
     '<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">\n<CDXML BoundingBox="205.13 72.24 221.13 90.24" BondLength="14.4" LabelFont="3" LabelSize="10" LabelFace="96" CaptionFont="3" CaptionSize="10" LineWidth="0.6" BoldWidth="2" HashSpacing="2.5" MarginWidth="1.6" BondSpacing="18" ChainAngle="120" LabelJustification="Auto" CaptionJustification="Left" color="0" bgcolor="1"><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/><color r="0.2" g="0.4" b="0.6"/></colortable><fonttable><font id="3" charset="iso-8859-1" name="Arial"/><font id="5" charset="Unicode" name="SimSun"/></fonttable><page id="1">\n<graphic id="1001" GraphicType="Orbital" OrbitalType="p" BoundingBox="213.13 82.24 213.13 76.24" Center3D="213.13 82.24 0" MajorAxisEnd3D="213.13 76.24 0" MinorAxisEnd3D="215.13 82.24 0" color="4"/>\n</page></CDXML>\n'),
])
def test_native_p_orbital_build_output_is_byte_compatible(spec, expected):
    assert acs_builder.build_from_spec(spec) == expected
