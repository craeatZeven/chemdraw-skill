"""P24: charges as ChemDraw's circled signs on ordinary structures, when a figure asks for it.

Most corpus figures that show charges circle them (BJOC 21-208 s2, 6-73 f1, 6-84 s4, 7-37 f2);
some write a plain sign (6-118 s1, 7-45 s22), so the default stays the superscript sign."""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.render import renderer as R

ALKOXIDE = ([{"index": 1, "symbol": "O", "x": 0.0, "y": 0.0, "charge": -1, "num_hydrogens": 0},
             {"index": 2, "symbol": "C", "x": 12.47, "y": 7.2},
             {"index": 3, "symbol": "C", "x": 24.94, "y": 0.0}],
            [{"index": 1, "atom1": 1, "atom2": 2, "order": 1}, {"index": 2, "atom1": 2, "atom2": 3, "order": 1}])
CATION = ([{"index": 1, "symbol": "C", "x": 0.0, "y": 0.0, "charge": 1},
           {"index": 2, "symbol": "C", "x": 12.47, "y": 7.2},
           {"index": 3, "symbol": "C", "x": -12.47, "y": 7.2},
           {"index": 4, "symbol": "C", "x": 0.0, "y": -14.4}],
          [{"index": i, "atom1": 1, "atom2": i + 1, "order": 1} for i in (1, 2, 3)])


def _fragment(drawing, style=None):
    atoms, bonds = drawing
    if style:
        with R.charge_style(style):
            xml, _, _ = R._build_fragment([dict(a) for a in atoms], [dict(b) for b in bonds], R._IDGen())
    else:
        xml, _, _ = R._build_fragment([dict(a) for a in atoms], [dict(b) for b in bonds], R._IDGen())
    return ET.fromstring(xml)


def _symbols(fragment):
    return [g for g in fragment.iter("graphic") if g.get("GraphicType") == "Symbol"]


def test_the_default_is_the_superscript_sign():
    fragment = _fragment(ALKOXIDE)
    (oxygen,) = [n for n in fragment.iter("n") if n.get("Charge")]
    assert [s.text for s in oxygen.iter("s")] == ["O", "-"]
    assert _symbols(fragment) == []


def test_circled_puts_the_sign_in_a_symbol_that_represents_the_charge():
    fragment = _fragment(ALKOXIDE, "circled")
    (oxygen,) = [n for n in fragment.iter("n") if n.get("Charge")]
    assert oxygen.get("Charge") == "-1" and "".join(oxygen.itertext()).strip() == "O"
    (symbol,) = _symbols(fragment)
    assert symbol.get("SymbolType") == "CircleMinus"
    assert symbol.find("represent").attrib == {"attribute": "Charge", "object": oxygen.get("id")}


def test_a_charged_carbon_stays_unlabelled_with_its_sign_beside_it():
    fragment = _fragment(CATION, "circled")
    (carbon,) = [n for n in fragment.iter("n") if n.get("Charge")]
    assert carbon.get("Charge") == "1" and carbon.find("t") is None and carbon.get("Element") in (None, "6")
    (symbol,) = _symbols(fragment)
    assert symbol.get("SymbolType") == "CirclePlus"
    x1, y1, x2, y2 = map(float, symbol.get("BoundingBox").split())
    centre = ((x1 + x2) / 2, (y1 + y2) / 2)
    at = {n.get("id"): tuple(map(float, n.get("p").split())) for n in fragment.iter("n")}
    def gap(p, a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
        return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)
    assert min(gap(centre, at[b.get("B")], at[b.get("E")]) for b in fragment.iter("b")) - (x2 - x1) / 2 >= 3.0


def test_a_charge_beyond_one_keeps_its_superscript():
    atoms, bonds = ALKOXIDE
    fragment = _fragment(([{**atoms[0], "charge": -2}] + atoms[1:], bonds), "circled")
    (oxygen,) = [n for n in fragment.iter("n") if n.get("Charge")]
    assert [s.text for s in oxygen.iter("s")] == ["O", "2-"] and _symbols(fragment) == []


def test_the_layout_box_holds_the_sign():
    atoms, _ = ALKOXIDE
    (symbol,) = _symbols(_fragment(ALKOXIDE, "circled"))
    x1, y1, x2, y2 = map(float, symbol.get("BoundingBox").split())
    with R.charge_style("circled"):
        box = R._fragment_bbox([dict(a) for a in atoms])
    assert box[0] <= x1 and box[1] <= y1 and x2 <= box[2] and y2 <= box[3]


def test_the_style_is_a_known_word():
    with pytest.raises(ValueError):
        with R.charge_style("boxed"):
            pass

def test_a_lone_charged_atom_gets_a_symbol(tmp_path):
    # No bonds, no neighbours: every clearance term is empty, which used to be min() of nothing.
    from cdxml_toolkit.mcp_runtime import figure_draw
    brief = {"kind": "scheme",
             "species": [{"id": "a", "smiles": "[CH3+]", "label": "1"}],
             "steps": [{"substrates": ["a"], "products": ["a"]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "one.cdxml"))
    assert report["ok"], report
