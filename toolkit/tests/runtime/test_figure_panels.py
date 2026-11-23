"""Several schemes on one page, each under a boxed title, with a colour legend (P20).

Witness: BJOC 22-103 figure 1 — panels (a), (b), (c), each under a boxed one-line title,
substrates coloured by class, and a footer legend (R1 = H, MeO, Cl in one colour, R2 = ... in
another).
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.mcp_runtime import figure_draw, figure_panels

SCHEME = {"species": [{"id": "A", "smiles": "CCO"}, {"id": "B", "smiles": "CC=O"}],
          "steps": [{"substrates": ["A"], "products": ["B"], "above": ["PCC"]}]}


def _brief(**extra):
    return {"kind": "panels", "panels": [{"title": "(a) oxidation", "brief": SCHEME},
                                         {"title": "(b) the same again", "brief": SCHEME}], **extra}


def test_a_panels_brief_is_checked():
    assert figure_panels._validate(_brief(legend=["R1 = H, MeO", [{"text": "R2", "color": "#1F4E9E"}, {"text": " = CF3"}]])) == []
    for brief, message in ((_brief(panels=[]), "panels"),
                           ({"kind": "panels", "panels": [{"brief": SCHEME}]}, "title"),
                           ({"kind": "panels", "panels": [{"title": "x", "brief": {"kind": "scope"}}]}, "scheme, mechanism, table or energy_profile"),
                           (_brief(legend="R1 = H"), "legend"),
                           (_brief(colour="red"), "Unknown brief keys")):
        assert any(message in e for e in figure_panels._validate(brief)), brief


def test_panels_stack_under_titles_of_one_width():
    sizes = [(300.0, 80.0), (220.0, 120.0)]
    titles, panels, width, bottom = figure_panels._stack(sizes, title_height=16.0, gap=8.0, panel_gap=16.0)
    assert width == 300.0
    assert titles[0] == (0.0, 0.0) and panels[0] == (0.0, 24.0)
    assert titles[1] == (0.0, 24.0 + 80.0 + 16.0) and panels[1] == (0.0, titles[1][1] + 24.0)
    assert bottom == panels[1][1] + 120.0


def test_draw_figure_routes_panels():
    report = figure_draw.draw_figure({"kind": "panels", "panels": []})
    assert report["ok"] is False and any("panels" in e for e in report["errors"])


def test_three_panels_are_drawn_on_one_page_with_titles_and_a_legend(tmp_path):
    brief = _brief(legend=[[{"text": "R1", "color": "#1F4E9E"}, {"text": " = H, MeO, Cl"}]])
    brief["panels"].append({"title": "(c) this work", "brief": SCHEME})
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "panels.cdxml"))
    assert report["ok"], report
    root = ET.parse(tmp_path / "panels.cdxml").getroot()
    rects = [g for g in root.iter("graphic") if g.get("GraphicType") == "Rectangle"]
    assert len(rects) == 3
    widths = {round(float(r.get("BoundingBox").split()[2]) - float(r.get("BoundingBox").split()[0]), 1) for r in rects}
    assert len(widths) == 1
    texts = ["".join(t.itertext()).strip() for t in root.iter("t")]
    assert "(a) oxidation" in texts and "(c) this work" in texts and "R1 = H, MeO, Cl" in texts
    assert len(root.findall(".//fragment")) == 6
    assert [p["readback"] for p in report["panels"]] == ["matched"] * 3


# Panel (a) of the witness is a gallery: structures with their names, no arrows. Panel (b) is
# split into sections, each opened by a boxed author line on a dashed rule.

GALLERY = [{"id": "g1", "smiles": "c1ccccc1O", "caption": ["phenol", "antiseptic"]},
           {"id": "g2", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "caption": "aspirin\nanalgesic agent",
            "color": "#1F4E9E"}]


def test_gallery_and_sections_panels_are_checked():
    ok = {"kind": "panels", "panels": [
        {"title": "(a) examples", "gallery": GALLERY},
        {"title": "(b) earlier work", "sections": [{"title": "A. Author 2009 [35]", "brief": SCHEME},
                                                   {"title": "B. Author 2014 [36]", "brief": SCHEME}]}]}
    assert figure_panels._validate(ok) == []
    for panel, message in (({"title": "x"}, "one of brief, gallery or sections"),
                           ({"title": "x", "brief": SCHEME, "gallery": GALLERY}, "one of brief, gallery or sections"),
                           ({"title": "x", "gallery": []}, "gallery"),
                           ({"title": "x", "gallery": [{"id": "g", "smiles": "C", "colour": "red"}]}, "unknown keys"),
                           ({"title": "x", "gallery": [{"id": "g"}]}, "smiles or a name"),
                           ({"title": "x", "sections": [{"brief": SCHEME}]}, "title"),
                           ({"title": "x", "sections": [{"title": "s", "brief": {"kind": "scope"}}]}, "scheme, mechanism, table or energy_profile")):
        errors = figure_panels._validate({"kind": "panels", "panels": [panel]})
        assert any(message in e for e in errors), (panel, errors)


def test_gallery_rows_wrap_at_the_width_and_spread_to_the_widest_row():
    sizes = [(100.0, 60.0), (120.0, 80.0), (90.0, 50.0), (150.0, 70.0)]
    at, width, height = figure_panels._gallery_rows(sizes, max_width=300.0, gap=20.0)
    # 100 + 20 + 120 = 240 fits; adding 90 would make 330 > 300, so the third starts a new row.
    assert [round(y, 1) for _, y in at] == [0.0, 0.0, 80.0 + 20.0, 80.0 + 20.0]
    assert width == 90.0 + 20.0 + 150.0   # the widest row's width
    # Each row runs from 0 to the widest row's right edge.
    assert at[0][0] == 0.0 and at[1][0] + 120.0 == pytest.approx(width)
    assert at[2][0] == 0.0 and at[3][0] + 150.0 == pytest.approx(width)
    assert height == 80.0 + 20.0 + 70.0


def test_a_lone_gallery_cell_is_not_stretched():
    at, width, height = figure_panels._gallery_rows([(100.0, 60.0)], max_width=300.0, gap=20.0)
    assert at == [(0.0, 0.0)] and width == 100.0 and height == 60.0


def test_a_section_title_is_boxed_on_a_dashed_rule_across_the_figure():
    objects = figure_panels._section_rule("L. Grimaud 2009 [35]", width=500.0, text_width=90.0)
    dashed = [o for o in objects if o["type"] == "line"]
    boxes = [o for o in objects if o["type"] == "rectangle"]
    (text,) = [o for o in objects if o["type"] == "text"]
    assert len(dashed) == 2 and all(o["style"] == "dashed" for o in dashed) and len(boxes) == 1
    lead, tail = sorted(dashed, key=lambda o: o["start"][0])
    box = boxes[0]
    assert lead["start"][0] == 0.0 and lead["end"][0] < box["start"][0]
    assert tail["start"][0] > box["end"][0] and tail["end"][0] == 500.0
    assert box["end"][0] - box["start"][0] >= 0.6 * 500.0 - 1   # the witness's box spans most of the figure
    assert box["start"][0] < text["position"][0] < box["end"][0]
    assert lead["start"][1] == pytest.approx((box["start"][1] + box["end"][1]) / 2)


def test_a_gallery_and_sections_are_drawn(tmp_path):
    brief = {"kind": "panels", "panels": [
        {"title": "(a) examples", "gallery": GALLERY},
        {"title": "(b) earlier work", "sections": [{"title": "A. Author 2009 [35]", "brief": SCHEME},
                                                   {"title": "B. Author 2014 [36]", "brief": SCHEME}]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "panels.cdxml"))
    assert report["ok"], report
    assert report["panels"][0]["layout"] == "gallery" and report["panels"][0]["readback"] == "matched"
    root = ET.parse(tmp_path / "panels.cdxml").getroot()
    texts = ["".join(t.itertext()).strip() for t in root.iter("t")]
    assert "phenol" in texts and "analgesic agent" in texts and "A. Author 2009 [35]" in texts
    rects = [g for g in root.iter("graphic") if g.get("GraphicType") == "Rectangle"]
    assert len(rects) == 4   # two panel titles, two section titles
    assert len([a for a in root.iter("arrow") if a.get("ArrowheadHead") == "Full"]) == 2   # none in the gallery
