"""Species labels that are formulas get subscripts in a scheme; compound numbers do not (P17)."""

from __future__ import annotations

import xml.etree.ElementTree as ET

from cdxml_toolkit.render import renderer as R


def _runs(label: str) -> list[tuple[str, str]]:
    xml, _ = R._build_label_element(label, 100.0, 100.0, R._IDGen(10))
    t = ET.fromstring(xml)
    return [(s.text, s.get("face")) for s in t.iter("s")]


def test_a_formula_label_is_written_with_subscripts():
    runs = _runs("PhB(OH)2")
    assert ("2", "32") in runs and "".join(text for text, _ in runs) == "PhB(OH)2"


def test_a_compound_number_is_written_plain():
    for label in ("TS2", "INT3", "IM1", "3a", "L1"):
        runs = _runs(label)
        assert len(runs) == 1 and runs[0] == (label, "0" if runs[0][1] == "0" else runs[0][1])
        assert runs[0][1] not in ("32", "96")


def test_a_bold_house_label_keeps_its_subscripts():
    from cdxml_toolkit.mcp_runtime.figure_house import _with_weight
    assert _with_weight("32", bold=True) == "33" and _with_weight("96", bold=True) == "97"
    assert _with_weight("33", bold=False) == "32" and _with_weight("0", bold=True) == "1"
    assert _with_weight(None, bold=True) == "1"


def test_a_numbered_r_group_is_drawn_as_r_with_its_number_raised():
    atoms = [{"index": 1, "symbol": "*", "x": 0.0, "y": 0.0, "atomic_number": 0, "generic_number": 2},
             {"index": 2, "symbol": "C", "x": 14.4, "y": 0.0, "atomic_number": 6}]
    bonds = [{"index": 1, "atom1": 1, "atom2": 2, "order": 1}]
    xml, _, _ = R._build_fragment(atoms, bonds, R._IDGen(10))
    node = ET.fromstring(xml).find("n")
    assert node.get("GenericNickname") == "R2"
    assert [(s.text, s.get("face")) for s in node.iter("s")][-1] == ("2", "64")
