"""Where ChemDraw sets a label's glyphs, and charge symbols kept 3 pt from them (P42).

Measured on ChemDraw's own saves (assets/coverage/*.native.cdxml, 10 pt Arial): the glyphs run
from 4.6 pt above the atom to 4.0 pt below it; a left-justified label centres its first
character on the atom, a right-justified one its last."""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdxml_toolkit import label_geometry as G
from cdxml_toolkit.mcp_runtime import figure_house

NATIVE = Path(__file__).resolve().parents[3] / "assets" / "coverage"


@pytest.mark.parametrize("name, text, justification", [
    ("stereo_model.native.cdxml", "OH", "Left"), ("stereo_model.native.cdxml", "Me", "Right"),
    ("stereo_model.native.cdxml", "Br", "Left"), ("process_route.native.cdxml", "Br", "Right"),
    ("stereo_model.native.cdxml", "CHO", "Left"), ("preview_scheme.native.cdxml", "O-", "Right"),
])
def test_the_model_matches_chemdraws_own_label_boxes(name, text, justification):
    root = ET.parse(NATIVE / name).getroot()
    for node in root.iter("n"):
        t = node.find("t")
        if t is None or t.get("LabelJustification") != justification or "".join(t.itertext()).strip() != text:
            continue
        ax, ay = map(float, node.get("p").split()[:2])
        native = tuple(map(float, t.get("BoundingBox").split()))
        model = G.label_box(text, ax, ay, justification)
        # Horizontally to a tenth of a point; vertically the model is the taller of ChemDraw's boxes.
        assert model[0] == pytest.approx(native[0], abs=0.15) and model[2] == pytest.approx(native[2], abs=0.15)
        assert model[1] <= native[1] + 0.01 and model[3] >= native[3] - 0.01
        return
    pytest.skip(f"{text} {justification} not in {name}")


def _page(symbol_box, label_at=(100.0, 100.0), text="[Rh]", justification="Right"):
    x1, y1, x2, y2 = symbol_box
    return f"""<?xml version="1.0" encoding="UTF-8" ?>
<CDXML><page id="1"><fragment id="2">
<n id="3" p="{label_at[0]} {label_at[1]}" NodeType="GenericNickname" GenericNickname="{text}" Charge="-1">
<t id="4" p="0 0" LabelJustification="{justification}"><s>{text}</s></t></n>
<n id="5" p="{label_at[0] + 14.4} {label_at[1]}"/><b id="6" B="3" E="5"/>
<graphic id="7" BoundingBox="{x1} {y1} {x2} {y2}" GraphicType="Symbol" SymbolType="CircleMinus">
<represent attribute="Charge" object="3"/></graphic></fragment></page></CDXML>"""


def _clearance(text_xml):
    root = ET.fromstring(str(text_xml).split("?>", 1)[-1])
    return figure_house.charge_symbol_clearance(root)


def test_a_symbol_on_its_label_is_moved_clear():
    # A right-justified [Rh] runs left of its atom: a symbol written there sits on the text.
    before = _page((86.0, 92.0, 94.0, 100.0))
    assert min(c for _, c in _clearance(before)) < 0
    after, moved, crowded = figure_house.clear_charge_symbols(before)
    assert moved == 1 and crowded == []
    assert min(c for _, c in _clearance(after)) >= 3.0 - 1e-6


def test_a_symbol_already_clear_is_left_where_it_is():
    page = _page((112.0, 80.0, 120.0, 88.0))
    assert min(c for _, c in _clearance(page)) >= 3.0
    after, moved, _ = figure_house.clear_charge_symbols(page)
    assert moved == 0 and str(after) == page


def test_a_newman_projection_keeps_every_charge_3_pt_from_every_label_and_bond():
    from cdxml_toolkit.templates.conformers import newman_projection
    ring = {"root": 1, "atoms": [{"index": 1, "symbol": "C", "x": 0.0, "y": 0.0},
                                 {"index": 2, "symbol": "O", "x": 7.2, "y": -12.47, "charge": 1},
                                 {"index": 3, "symbol": "C", "x": 14.4, "y": 0.0}],
            "bonds": [{"atom1": 1, "atom2": 2, "order": 2}, {"atom1": 1, "atom2": 3, "order": 1},
                      {"atom1": 2, "atom2": 3, "order": 1}]}
    for rear, charges in ((("H", "[Rh]", "Me"), [0, -1, 0]), (("[Rh]", "Me", "H"), [-1, 0, 0]),
                          (("Me", "H", "[Rh]"), [0, 0, -1])):
        drawn = newman_projection(front=("*", "H", "Me"), rear=rear, cx=100.0, cy=100.0,
                                  bond_length=20.16, attach=ring, charges={"rear": charges})
        root = ET.fromstring(drawn["cdxml"].split("?>", 1)[-1])
        clearances = figure_house.charge_symbol_clearance(root)
        assert len(clearances) == 2
        assert min(c for _, c in clearances) >= 3.0 - 1e-6, (rear, clearances)
