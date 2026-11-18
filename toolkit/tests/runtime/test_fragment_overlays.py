"""Dashed overlays that mark one fragment in several structures at once (parity item P9)."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdxml_toolkit.mcp_runtime import figure_draw, figure_house

DATA = Path(__file__).resolve().parents[1] / "data"


def _floats(text: str) -> list[float]:
    return [float(v) for v in text.split()]


def _pairs(text: str) -> list[tuple[float, float]]:
    v = _floats(text)
    return list(zip(v[0::2], v[1::2]))


def test_a_circle_overlay_is_written_the_way_chemdraw_writes_a_closed_one():
    # native-dashed-overlay.cdxml: ChemDraw's own closed dashed circle, CurveType="3" with
    # Closed="yes" and the first (in, anchor, out) repeated last. (Its first version was
    # CurveType="2", which ChemDraw draws open, with a gap at the upper right.)
    sample = (DATA / "native-dashed-overlay.cdxml").read_text(encoding="utf-8")
    curve = re.search(r'<curve\s+id="9001"[^>]*>', sample, re.S).group(0)
    native = _floats(re.search(r'CurvePoints="([^"]+)"', curve).group(1))
    assert 'CurveType="3"' in curve and 'Closed="yes"' in curve
    ours = _floats(figure_house.overlay_curve_points([(70.0, 95.0)], shape="circle", radius=36.0))
    assert len(native) == len(ours) == 54
    for a, b in zip(ours, native):
        assert a == pytest.approx(b, abs=0.02)
    assert ours[48:] == ours[:6]  # the last segment returns to the start


def test_an_outline_ends_where_it_starts():
    chain = [(100.0, 100.0), (112.47, 92.8), (124.94, 100.0), (137.41, 92.8)]
    pts = _pairs(figure_house.overlay_curve_points(chain, shape="outline"))
    assert pts[-3:] == pts[:3]


def test_a_circle_overlay_clears_every_marked_atom():
    hexagon = [(100 + 14.4 * math.cos(math.pi / 3 * k), 100 + 14.4 * math.sin(math.pi / 3 * k)) for k in range(6)]
    pts = _pairs(figure_house.overlay_curve_points(hexagon, shape="circle"))
    anchors = pts[1::3][:-1]
    radii = {round(math.hypot(x - 100, y - 100), 2) for x, y in anchors}
    assert len(radii) == 1 and radii.pop() >= 14.4 + figure_house.OVERLAY_PAD - 0.01


def test_an_outline_overlay_follows_the_atoms_at_a_distance():
    chain = [(100.0, 100.0), (112.47, 92.8), (124.94, 100.0), (137.41, 92.8)]
    pts = _pairs(figure_house.overlay_curve_points(chain, shape="outline"))
    assert len(pts) % 3 == 0 and len(pts) >= 9
    anchors = pts[1::3][:-1]
    for ax, ay in anchors:
        nearest = min(math.hypot(ax - x, ay - y) for x, y in chain)
        assert nearest == pytest.approx(figure_house.OVERLAY_PAD, abs=0.5)
    # Every atom lies inside the outline (ray casting on the anchors).
    for x, y in chain:
        inside = False
        for (x1, y1), (x2, y2) in zip(anchors, anchors[1:] + anchors[:1]):
            if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
                inside = not inside
        assert inside


def test_an_outline_of_one_or_two_atoms_falls_back_to_a_circle():
    one = figure_house.overlay_curve_points([(50.0, 50.0)], shape="outline")
    assert len(_floats(one)) == 54


@pytest.mark.parametrize("overlays, message", [
    ("c1ccccc1", "overlays must be a list"),
    ([{}], "overlays[1] needs a fragment"),
    ([{"fragment": "c1ccccc1", "shape": "star"}], "shape must be one of"),
    ([{"fragment": "c1ccccc1", "species": ["X"]}], "names 'X'"),
    ([{"fragment": "c1ccccc1", "style": "dotted"}], "unknown keys"),
])
def test_overlay_entries_are_checked(overlays, message):
    errors = figure_house.overlay_errors({"overlays": overlays}, {"A", "B"})
    assert any(message in e for e in errors), errors


def test_the_brief_accepts_overlays():
    brief = {"species": [{"id": "A", "smiles": "c1ccccc1C"}, {"id": "B", "smiles": "c1ccccc1CO"}],
             "steps": [{"substrates": ["A"], "products": ["B"], "above": ["[O]"]}],
             "overlays": [{"fragment": "c1ccccc1"}]}
    assert figure_draw._validate(brief) == []


# Toluene twice, the ring drawn at different places: the overlay must sit on the ring in both.
_TWO = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML LineWidth="0.6"><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/></colortable>
<page id="1">{frags}</page></CDXML>"""


def _ring_fragment(fid: int, cx: float, cy: float, extra: str) -> str:
    nodes, bonds = [], []
    for k in range(6):
        a = math.pi / 3 * k
        nodes.append(f'<n id="{fid + 1 + k}" p="{cx + 14.4 * math.cos(a):.2f} {cy + 14.4 * math.sin(a):.2f}"/>')
        order = ' Order="2"' if k % 2 == 0 else ""
        bonds.append(f'<b id="{fid + 20 + k}" B="{fid + 1 + k}" E="{fid + 1 + (k + 1) % 6}"{order}/>')
    nodes.append(f'<n id="{fid + 7}" p="{cx + 28.8:.2f} {cy:.2f}"{extra}/>')
    bonds.append(f'<b id="{fid + 30}" B="{fid + 1}" E="{fid + 7}"/>')
    return f'<fragment id="{fid}">' + "".join(nodes + bonds) + "</fragment>"


def test_every_named_structure_gets_the_overlay_in_the_same_place():
    cdxml = _TWO.format(frags=_ring_fragment(100, 100, 100, "") + _ring_fragment(200, 260, 140, ' Element="8"'))
    brief = {"species": [{"id": "A", "smiles": "Cc1ccccc1"}, {"id": "B", "smiles": "Oc1ccccc1"}],
             "overlays": [{"fragment": "c1ccccc1"}]}
    grounded = [{"id": "A", "smiles": "Cc1ccccc1"}, {"id": "B", "smiles": "Oc1ccccc1"}]
    out, applied, warnings = figure_house.apply_fragment_overlays(cdxml, brief, grounded)
    page = ET.fromstring(str(out).split("?>", 1)[-1]).find("page")
    curves = page.findall("curve")
    assert applied == ["fragment overlays (2)"] and warnings == []
    assert len(curves) == 2
    offsets = []
    for curve, (cx, cy) in zip(curves, [(100, 100), (260, 140)]):
        assert curve.get("LineType") == "Dashed" and curve.get("CurveType") == "3" and curve.get("Closed") == "yes"
        assert curve.get("FillType") is None
        anchors = _pairs(curve.get("CurvePoints"))[1::3][:-1]
        mx = sum(x for x, _ in anchors) / len(anchors)
        my = sum(y for _, y in anchors) / len(anchors)
        offsets.append((round(mx - cx, 2), round(my - cy, 2)))
    assert offsets[0] == offsets[1] == (0.0, 0.0)


def test_a_named_structure_without_the_fragment_is_refused():
    cdxml = _TWO.format(frags=_ring_fragment(100, 100, 100, ""))
    brief = {"species": [{"id": "A", "smiles": "Cc1ccccc1"}],
             "overlays": [{"fragment": "C1CCCCC1", "species": ["A"]}]}
    with pytest.raises(ValueError, match="does not contain"):
        figure_house.apply_fragment_overlays(cdxml, brief, [{"id": "A", "smiles": "Cc1ccccc1"}])


CYCLE = {
    "layout": "catalytic-cycle",
    "species": [
        {"id": "I", "smiles": "c1ccc(cc1)[Pd]Br", "label": "I"},
        {"id": "II", "smiles": "c1ccc(cc1)[Pd]C", "label": "II"},
        {"id": "III", "smiles": "Cc1ccccc1", "label": "III"},
        {"id": "MeB", "smiles": "CB(O)O"},
        {"id": "Pd", "smiles": "[Pd]"},
        {"id": "PhBr", "smiles": "Brc1ccccc1"},
    ],
    "steps": [
        {"substrates": ["I", "MeB"], "products": ["II"], "above": ["transmetalation"]},
        {"substrates": ["II"], "products": ["III", "Pd"], "above": ["reductive elimination"]},
        {"substrates": ["Pd", "PhBr"], "products": ["I"], "above": ["oxidative addition"]},
    ],
    "overlays": [{"fragment": "c1ccccc1", "species": ["I", "II", "III"]}],
}


def test_a_cycle_marks_the_ring_in_every_intermediate(tmp_path):
    report = figure_draw.draw_figure(CYCLE, output_path=str(tmp_path / "cycle.cdxml"))
    assert report["ok"], report
    root = ET.parse(tmp_path / "cycle.cdxml").getroot()
    dashed = [c for c in root.iter("curve") if c.get("LineType") == "Dashed"]
    assert len(dashed) == 3
    radii = set()
    for curve in dashed:
        anchors = _pairs(curve.get("CurvePoints"))[1::3][:-1]
        mx = sum(x for x, _ in anchors) / len(anchors)
        my = sum(y for _, y in anchors) / len(anchors)
        radii.add(round(math.hypot(anchors[0][0] - mx, anchors[0][1] - my), 1))
    assert len(radii) == 1


def _overlay_page(*extra: str) -> str:
    ring = figure_house.overlay_curve_points([(100.0, 100.0)], shape="circle", radius=20.0)
    return ('<?xml version="1.0" encoding="UTF-8" ?><CDXML><page id="1">'
            f'<curve id="7" LineType="Dashed" CurveType="3" Closed="yes" CurvePoints="{ring}"/>' + "".join(extra) + "</page></CDXML>")


def test_an_arrow_or_caption_on_an_overlay_is_reported():
    from cdxml_toolkit.cdxml_check import check_cdxml, overlays_collide
    crossing = '<arrow id="20" Head3D="140 100 0" Tail3D="60 100 0" ArrowheadHead="Full"/>'
    clear = '<arrow id="21" Head3D="260 100 0" Tail3D="200 100 0" ArrowheadHead="Full"/>'
    on_rim = '<t id="30" p="96 83"><s size="10">reductive elimination</s></t>'
    far = '<t id="31" p="200 200"><s size="10">transmetalation</s></t>'
    found = overlays_collide(ET.fromstring(_overlay_page(crossing, clear, on_rim, far).split("?>", 1)[1]))
    assert len(found) == 2, found
    assert any('<arrow id="20">' in f and "overlay" in f for f in found)
    assert any("reductive elimination" in f for f in found)
    codes = {w["code"] for w in check_cdxml(_overlay_page(crossing, on_rim), palette=None)["warnings"]}
    assert "overlay_collides" in codes


def test_an_overlay_with_nothing_on_it_is_not_reported():
    from cdxml_toolkit.cdxml_check import overlays_collide
    clear = '<arrow id="21" Head3D="260 100 0" Tail3D="200 100 0" ArrowheadHead="Full"/>'
    assert overlays_collide(ET.fromstring(_overlay_page(clear).split("?>", 1)[1])) == []


def test_the_overlay_box_covers_the_drawn_outline():
    from cdxml_toolkit.render.renderer import _overlay_box
    atoms = [{"index": i + 1, "x": 100 + 14.4 * math.cos(math.pi / 3 * i), "y": 100 + 14.4 * math.sin(math.pi / 3 * i)}
             for i in range(6)] + [{"index": 7, "x": 300.0, "y": 300.0}]
    for shape in ("circle", "outline"):
        box = _overlay_box(atoms, [{"atoms": list(range(6)), "shape": shape}])
        pts = _pairs(figure_house.overlay_curve_points([(a["x"], a["y"]) for a in atoms[:6]], shape=shape))
        assert box[0] <= min(x for x, _ in pts) and box[2] >= max(x for x, _ in pts)
        assert box[1] <= min(y for _, y in pts) and box[3] >= max(y for _, y in pts)
        assert box[2] < 150  # the unmarked atom is not part of it
    assert _overlay_box(atoms, []) is None


def test_cycle_arrows_and_labels_keep_off_an_intermediates_overlay():
    from cdxml_toolkit.cdxml_check import overlays_collide
    from cdxml_toolkit.render import renderer as R
    from cdxml_toolkit.render.schema import ArrowContent, SchemeDescriptor, StepDescriptor, StructureRef
    ring = "COc1ccc(cc1)"
    structures = {
        "I": StructureRef(id="I", smiles=ring + "[Pd]I", label="I", overlays=({"atoms": list(range(8)), "shape": "circle"},)),
        "II": StructureRef(id="II", smiles=ring + "[Pd]c1ccccc1", label="II",
                           overlays=({"atoms": list(range(8)), "shape": "circle"},)),
        "Pd": StructureRef(id="Pd", smiles="[Pd]", label="Pd(0)"),
    }
    steps = [StepDescriptor(substrates=["Pd"], products=["I"], above_arrow=ArrowContent(text=["oxidative addition"])),
             StepDescriptor(substrates=["I"], products=["II"], above_arrow=ArrowContent(text=["transmetalation"])),
             StepDescriptor(substrates=["II"], products=["Pd"], above_arrow=ArrowContent(text=["reductive elimination"]))]
    scheme = SchemeDescriptor(structures=structures, steps=steps, layout="catalytic-cycle")
    xml, _ = R._layout_catalytic_cycle(scheme, R._IDGen())
    page = ET.fromstring(f"<page>{xml}</page>")
    frags = {f.get("id"): f for f in page.findall("fragment")}
    curves = []
    for frag in frags.values():
        nodes = frag.findall("n")
        if len(nodes) < 8:  # the bare Pd(0) has one atom; I and II carry the anisyl ring
            continue
        pts = [tuple(float(v) for v in n.get("p").split()) for n in nodes[:8]]
        curves.append(f'<curve id="9{len(curves)}" LineType="Dashed" CurveType="3" Closed="yes" '
                      f'CurvePoints="{figure_house.overlay_curve_points(pts)}"/>')
    assert len(curves) == 2
    for curve in curves:
        page.append(ET.fromstring(curve))
    assert overlays_collide(ET.fromstring(f"<CDXML>{ET.tostring(page, encoding='unicode')}</CDXML>")) == []


def test_the_check_reads_an_open_dashed_curve_as_drawn():
    # A curve without the closed bit is drawn open, so nothing meets its missing last segment.
    from cdxml_toolkit.cdxml_check import overlays_collide
    native = re.search(r'CurvePoints="([^"]+)"', _overlay_page()).group(1).rsplit(" ", 6)[0]  # drop the repeat
    through_gap = '<arrow id="22" Head3D="135 85.5 0" Tail3D="108 96.7 0" ArrowheadHead="Full"/>'
    open_page = ('<CDXML><page id="1"><curve id="7" LineType="Dashed" CurveType="2" '
                 f'CurvePoints="{native}"/>{through_gap}</page></CDXML>')
    closed_page = open_page.replace('CurveType="2"', 'CurveType="3"')
    assert overlays_collide(open_page) == []
    assert len(overlays_collide(closed_page)) == 1
