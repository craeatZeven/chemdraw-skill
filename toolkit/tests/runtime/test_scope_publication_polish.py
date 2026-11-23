"""TDD tests for publication-grade substrate scope polish derived from cross-journal corpus analysis."""
from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest
from PIL import Image

from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setenv("CHEMDRAW_OFFLINE", "1")


def _root(report) -> ET.Element:
    text = open(report["pages"][0]["path"], encoding="utf-8").read()
    return ET.fromstring(text[text.index("<CDXML"):])


def _texts(root: ET.Element) -> list[str]:
    return ["".join(t.itertext()) for t in root.iter("t")]


def test_empty_title_omits_title_element(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "species": [
            {"id": "3a", "smiles": "Cc1ccccc1", "label": "3a", "values": {"yield": "85%"}},
            {"id": "3b", "smiles": "COc1ccccc1", "label": "3b", "values": {"yield": "78%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "no-title"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    texts = _texts(_root(report))
    assert "Substrate scope" not in texts


def test_centered_header_and_colored_header_substrates(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "width_mm": 180,
        "columns": 4,
        "row_align": "left",
        "header": {
            "substrates": [
                {"smiles": "c1ccc([*:1])cc1", "label": "1", "color": "#0B4F8A"},
                {"smiles": "BrC1CCCCC1", "label": "2", "color": "#8B1E24"},
            ],
            "products": [
                {"smiles": "c1ccc(C2CCCCC2)cc1", "label": "3"},
            ],
            "above": ["Fe(acac)3 (5 mol%)"],
            "below": ["MCH/TMEDA, 0 °C, 1 h"],
        },
        "species": [
            {"id": f"3{c}", "smiles": smi, "label": f"3{c}", "values": {"yield": f"{70+i}%"}}
            for i, (c, smi) in enumerate([
                ("a", "Cc1ccc(C2CCCCC2)cc1"),
                ("b", "CCc1ccc(C2CCCCC2)cc1"),
                ("c", "COc1ccc(C2CCCCC2)cc1"),
                ("d", "Fc1ccc(C2CCCCC2)cc1"),
            ])
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "header-color"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    root = _root(report)
    page_w = float(root.get("BoundingBox").split()[2])
    # Header arrow should be roughly centered over the page (within 25% of center, not flush left)
    arrows = [a for a in root.iter("arrow") if a.get("ArrowheadHead") == "Full"]
    assert len(arrows) == 1
    head_x = float(arrows[0].get("Head3D").split()[0])
    tail_x = float(arrows[0].get("Tail3D").split()[0])
    arrow_mid = (head_x + tail_x) / 2.0
    assert abs(arrow_mid - page_w / 2.0) < page_w * 0.25
    # Both header colors should appear in the colortable
    colors = root.find("colortable").findall("color")
    rgbs = {(round(float(c.get("r", 0)) * 255), round(float(c.get("g", 0)) * 255), round(float(c.get("b", 0)) * 255))
            for c in colors}
    assert (11, 79, 138) in rgbs  # #0B4F8A
    assert (139, 30, 36) in rgbs  # #8B1E24


def test_row_baseline_alignment_columns_and_centered_last_row(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "columns": 3,
        "row_align": "center",
        "species": [
            # Row 1: markedly different vertical heights (tall biphenyl-naphthyl vs short benzene)
            {"id": "3a", "smiles": "c1ccccc1", "label": "3a", "values": {"yield": "90%"}},
            {"id": "3b", "smiles": "c1ccc2cc3ccccc3cc2c1", "label": "3b", "values": {"yield": "85%"}},
            {"id": "3c", "smiles": "Cc1ccccc1", "label": "3c", "values": {"yield": "82%"}},
            # Row 2: 1 compound only -> should be horizontally centered!
            {"id": "3d", "smiles": "COc1ccccc1", "label": "3d", "values": {"yield": "77%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "aligned-rows"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    root = _root(report)
    page_w = float(root.get("BoundingBox").split()[2])
    label_pos = {}
    for t in root.iter("t"):
        txt = "".join(t.itertext()).strip()
        if txt in {"3a", "3b", "3c", "3d"}:
            x, y = map(float, t.get("p").split()[:2])
            label_pos[txt] = (x, y)
    assert set(label_pos) == {"3a", "3b", "3c", "3d"}
    # Row 1 labels share the exact same Y baseline!
    assert abs(label_pos["3a"][1] - label_pos["3b"][1]) < 0.5
    assert abs(label_pos["3b"][1] - label_pos["3c"][1]) < 0.5
    # Row 2 single compound is on a lower row AND horizontally centered (x far right of 3a's x)
    assert label_pos["3d"][1] > label_pos["3a"][1] + 20
    assert label_pos["3d"][0] > label_pos["3a"][0] + 40
    assert abs(label_pos["3d"][0] - page_w / 2.0) < 40


def test_scope_header_combines_text_rotation_color_and_ring_marks(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "columns": 3,
        "rmark_in_ring": True,
        "header": {
            "substrates": [
                {"smiles": "c1ccc([*:1])cc1", "rotation": 180, "color": "#0B4F8A"},
                {"text": "S8"},
            ],
            "products": [{"smiles": "Cc1ccccc1"}],
        },
        "species": [{"id": "3a", "smiles": "Cc1ccccc1", "values": {"yield": "88%"}}],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "combined-header"))
    assert report["ok"] is True, report.get("errors")
    assert report["readback"]["status"] == "matched"
    assert report["pages"][0]["grid"]["columns"] == 3
    root = _root(report)
    assert "S8" in _texts(root)
    assert any(n.get("GenericNickname") == "R1" and n.get("Visible") == "no" for n in root.iter("n"))
    assert any(
        tuple(round(float(c.get(key, 0)) * 255) for key in ("r", "g", "b")) == (11, 79, 138)
        for c in root.findall("colortable/color")
    )


def test_fragment_colors_multi_synthon_and_secondary_stays_grey(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "highlight_varying": "#D55E00",
        "fragment_colors": [
            {"smarts": "C1CCCCC1", "color": "#8B1E24"},
            {"smarts": "c1ccccc1", "color": "#0B4F8A"},
        ],
        "species": [
            {"id": "3a", "smiles": "Cc1ccc(C2CCCCC2)cc1", "label": "3a", "values": {"yield": "82%"}},
            {"id": "3b", "smiles": "COc1ccc(C2CCCCC2)cc1", "label": "3b", "values": {"yield": "76%"}},
            {"id": "3b'", "smiles": "COc1cccc(C2CCCCC2)c1", "label": "3b'", "secondary_of": "3b",
             "values": {"yield": "8%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "synthon-colors"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    root = _root(report)
    colors = root.find("colortable").findall("color")
    rgbs = {(round(float(c.get("r", 0)) * 255), round(float(c.get("g", 0)) * 255), round(float(c.get("b", 0)) * 255))
            for c in colors}
    assert (139, 30, 36) in rgbs   # #8B1E24 cyclohexyl
    assert (11, 79, 138) in rgbs   # #0B4F8A phenyl
    # Secondary product 3b' must stay grey (#808080) and NOT carry #D55E00 or #8B1E24
    color_idx_to_rgb = {
        str(k + 2): (round(float(c.get("r", 0)) * 255), round(float(c.get("g", 0)) * 255), round(float(c.get("b", 0)) * 255))
        for k, c in enumerate(colors)
    }
    groups = list(root.find("page").findall("group"))
    sec_group = next(g for g in groups if any("3b'" in "".join(t.itertext()) for t in g.findall("t")))
    sec_frag = sec_group.find("fragment")
    sec_colors = {color_idx_to_rgb.get(el.get("color")) for el in list(sec_frag.findall("n")) + list(sec_frag.findall("b"))
                  if el.get("color")}
    assert sec_colors <= {(128, 128, 128)}


def test_footnote_superscripts_subscripts_and_italic_notes(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "footnotes": [
            "a Isolated yield in H2O at 80 °C.",
            "b Determined by chiral HPLC.",
        ],
        "species": [
            {"id": "3a", "smiles": "Cc1ccccc1", "label": "3a^a", "note": "*from thymol*",
             "values": {"yield": "84%"}},
            {"id": "3b", "smiles": "COc1ccccc1", "label": "3b^b", "note": "*from carvacrol*",
             "values": {"yield": "79%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "typography"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    root = _root(report)
    # Footnote element should have italic superscript 'a' (face bit 64|2 = 66) and subscript '2' (face bit 32)
    foot_t = next(t for t in root.iter("t") if "Isolated yield" in "".join(t.itertext()))
    runs = foot_t.findall("s")
    assert any((r.text or "") == "a" and (int(r.get("face", "0")) & 64) for r in runs)
    assert any((r.text or "") == "2" and (int(r.get("face", "0")) & 32) for r in runs)
    # Italic note "*from thymol*" renders as "from thymol" with italic bit (face & 2)
    note_t = next(t for t in root.iter("t") if "from thymol" in "".join(t.itertext()))
    assert "*" not in "".join(note_t.itertext())
    assert any("from thymol" in (r.text or "") and (int(r.get("face", "0")) & 2) for r in note_t.findall("s"))


def test_scope_preserves_charge_runs_and_italics_at_cell_width(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "footnotes": ["a With TBA^+ Nu^\u2212 in MeCN."],
        "species": [{"id": "20", "smiles": "CS", "group": "With TBA^+ Nu^\u2212",
                     "note": "*From N-Ac-L-cysteine*"}],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "charges"))
    assert report["checks_passed"], report
    root = _root(report)
    headings = [t for t in root.iter("t") if "TBA" in "".join(t.itertext())]
    assert len(headings) == 2
    for heading in headings:
        charge_runs = [s for s in heading.findall("s") if s.text in ("+", "\u2212")]
        assert [s.text for s in charge_runs] == ["+", "\u2212"]
        assert all(int(s.get("face", "0")) & 64 for s in charge_runs)
        assert all(not int(s.get("face", "0")) & 2 for s in charge_runs)
    notes = [t for t in root.iter("t") if "cysteine" in "".join(t.itertext())]
    assert len(notes) == 1
    assert "".join(notes[0].itertext()) == "From N-Ac-L-cysteine"
    assert any(int(s.get("face", "0")) & 2 for s in notes[0].findall("s"))


def test_scope_header_labels_and_footer_reaction_preserve_inventory(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "header": {
            "substrates": [{"smiles": "CC(=O)Oc1ccccc1", "attached_label": {"smarts": "Oc1ccccc1", "text": "OPh"},
                            "counterion": "^\u2212OTf"}],
            "products": [{"smiles": "CC[*:1]", "variable_label": "Nu", "variable_label_color": "#1B4D89"}],
            "above": [[{"text": "1.0 equiv. "}, {"text": "H\u2013Nu", "color": "#1B4D89"}]],
        },
        "species": [{"id": "1", "smiles": "CCO"}],
        "footer_scheme": {"width": 220, "height": 65, "objects": [
            {"type": "molecule", "smiles": "CCCl", "position": [35, 30]},
            {"type": "arrow", "start": [75, 30], "end": [145, 30]},
            {"type": "molecule", "smiles": "CCBr", "position": [185, 30]},
        ]},
    }
    report = draw_figure(brief, output_path=str(tmp_path / "reaction"))
    assert report["ok"], report
    assert report["checks_passed"], report
    assert report["readback"]["structures_drawn"] == 5
    root = _root(report)
    assert root.find(".//n[@NodeType='Fragment']/fragment") is not None
    label = next(node.find("t/s") for node in root.findall(".//n[@NodeType='GenericNickname']")
                 if node.findtext("t/s") == "Nu")
    condition = next(run for run in root.iter("s") if run.text == "H\u2013Nu")
    assert label.get("color") == condition.get("color") not in (None, "0")
    assert {"OPh", "Nu", "\u2212OTf"} <= set(_texts(root))
    assert len(root.findall(".//arrow[@ArrowheadHead='Full']")) == 2


def test_scope_composed_header_keeps_individual_bond_colors(tmp_path):
    spec = {"kind": "scope", "species": [{"id": "1", "smiles": "CCC"}],
            "header_scheme": {"width": 90, "height": 75, "objects": [
                {"type": "molecule", "smiles": "CCOC", "position": [45, 40],
                 "bond_colors": {"0": "#1B4D89", "2": "#800000"}},
            ]}}
    report = draw_figure(spec, output_path=str(tmp_path / "header"))
    assert report["ok"] and report["checks_passed"], report
    assert report["readback"]["structures_drawn"] == 2
    bonds = _root(report).findall("./page/fragment/b")
    assert len(bonds) == 3
    assert len({bond.get("color") for bond in bonds}) == 3


def test_scope_counts_disconnected_salt_components(tmp_path):
    brief = {"kind": "scope", "species": [{"id": "1", "smiles": "C[N+](C)(C)C.[Cl-]"}]}
    report = draw_figure(brief, output_path=str(tmp_path / "salt"))
    assert report["ok"], report
    assert report["checks_passed"], report
    assert report["readback"]["structures_drawn"] == 2


def test_scope_preserves_signed_exponent_and_embeds_a_source_crop(tmp_path):
    source = tmp_path / "source.png"
    Image.new("RGB", (50, 40), "white").save(source)
    brief = {"kind": "scope", "species": [{"id": "1", "smiles": "CC", "note": "35.3 kcal mol^\u22121^a"}],
             "footer_scheme": {"width": 100, "height": 90, "objects": [
                 {"type": "image", "file": str(source), "crop": [5, 5, 35, 25],
                  "start": [20, 20], "end": [80, 60]},
             ]}}
    report = draw_figure(brief, output_path=str(tmp_path / "inset"))
    assert report["ok"] and report["checks_passed"], report
    assert report["readback"]["structures_drawn"] == 1
    root = _root(report)
    exponent = next(run for run in root.iter("s") if run.text == "\u22121")
    assert int(exponent.get("face", "0")) & 64
    assert not int(exponent.get("face", "0")) & 2
    embedded = root.find(".//embeddedobject")
    from io import BytesIO
    with Image.open(BytesIO(bytes.fromhex(embedded.get("PNG")))) as crop:
        assert crop.size == (30, 20)


def test_band_style_pill_and_panel_border(tmp_path):
    brief = {
        "kind": "scope",
        "title": "",
        "band_style": "pill",
        "border": True,
        "group_colors": {"Aryl halides": "#1B4D89", "Heteroaryl halides": "#8B1E24"},
        "species": [
            {"id": "3a", "smiles": "Cc1ccccc1", "label": "3a", "group": "Aryl halides", "values": {"yield": "88%"}},
            {"id": "3b", "smiles": "c1ccncc1", "label": "3b", "group": "Heteroaryl halides", "values": {"yield": "75%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "pill-border"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    root = _root(report)
    graphics = list(root.iter("graphic"))
    rects = [g for g in graphics if g.get("GraphicType") == "Rectangle"]
    # Expect at least 2 pill rectangles (one per group) + 1 outer border rectangle = 3
    assert len(rects) >= 3


def test_definition_contracted_group_with_highlight_varying_has_no_bond_warnings(tmp_path):
    from cdxml_toolkit.cdxml_check import text_over_ink

    brief = {
        "kind": "scope",
        "title": "",
        "highlight_varying": "#C81E1E",
        "species": [
            {"id": "3a", "smiles": "Cc1ccc(NC(=O)c2nc3ccccc3s2)cc1", "label": "3a", "values": {"yield": "64%"}},
            {"id": "3p", "smiles": "CCCCCCCNC(=O)c1nc2ccccc2s1", "label": "3p", "definition": "Alk = C7H15",
             "values": {"yield": "63%"}},
            {"id": "3q", "smiles": "CCCCCCCCCCCCCCCCCCNC(=O)c1nc2ccccc2s1", "label": "3q",
             "definition": "Alk = C18H37", "drawn_with": "3p", "values": {"yield": "60%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "def-highlight"))
    assert report["ok"] is True, report.get("errors")
    assert report["checks_passed"] is True
    assert not any("not bonded" in w for w in report["warnings"]), report["warnings"]
    assert text_over_ink(_root(report)) == []


def test_check_cdxml_accepts_scope_caption_size_when_in_allow_sizes(tmp_path):
    from cdxml_toolkit.cdxml_check import check_cdxml

    brief = {
        "kind": "scope",
        "title": "",
        "style": "nature_portfolio",
        "species": [
            {"id": "3a", "smiles": "Cc1ccccc1", "label": "3a", "values": {"yield": "85%"}},
            {"id": "3b", "smiles": "COc1ccccc1", "label": "3b", "values": {"yield": "78%"}},
        ],
    }
    report = draw_figure(brief, output_path=str(tmp_path / "nature-scope"))
    assert report["ok"] is True
    chk = check_cdxml(report["pages"][0]["path"], profile="nature_portfolio", palette=None, allow_sizes=(7.0, 8.0, 10.0))
    assert chk["ok"] is True, chk["errors"]
