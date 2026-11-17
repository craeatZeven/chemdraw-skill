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
                           ({"kind": "panels", "panels": [{"title": "x", "brief": {"kind": "scope"}}]}, "scheme or a mechanism"),
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
