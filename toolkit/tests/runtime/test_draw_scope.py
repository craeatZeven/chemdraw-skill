from __future__ import annotations

import math
import copy
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdxml_toolkit.mcp_runtime import figure_scope
from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure

SCOPE = {
    "kind": "scope",
    "title": "Substrate scope",
    "conditions": "Conditions as reported by the user.",
    "footnotes": ["a Isolated yield."],
    "species": [
        {"id": "3a", "smiles": "Cc1ccccc1", "label": "3a", "values": {"yield": "85%", "ee": "92%"}},
        {"id": "3b", "smiles": "COc1ccccc1", "label": "3b", "values": {"yield": "78%", "ee": "90%"}},
        {"id": "3c", "smiles": "Clc1ccccc1", "label": "3c", "values": {"yield": "81%"}},
        {"id": "3d", "smiles": "Fc1ccccc1", "label": "3d", "values": {"yield": "74%", "ee": "88%"}},
    ],
}


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setenv("CHEMDRAW_OFFLINE", "1")


def _brief(**changes):
    brief = copy.deepcopy(SCOPE)
    brief.update(changes)
    return brief


def _pages(report):
    return [ET.parse(page["path"]).getroot() for page in report["pages"]]


def _texts(report):
    return [" ".join("".join(t.itertext()).split()) for root in _pages(report) for t in root.iter("t")]


def test_a_scope_brief_is_drawn_into_a_folder_of_editable_pages(tmp_path):
    out = tmp_path / "scope"

    report = draw_figure(SCOPE, output_path=str(out))

    assert report["ok"] is True
    assert report["kind"] == "scope"
    assert report["output_dir"] == str(out)
    assert (out / "page-001.cdxml").is_file()
    assert [Path(p["path"]).name for p in report["pages"]] == ["page-001.cdxml"]
    assert len(report["species"]) == 4
    assert {s["source"] for s in report["species"]} == {"user"}
    assert report["readback"]["status"] == "matched"
    assert report["checks_passed"] is True


def test_labels_and_values_stay_with_their_structures(tmp_path):
    report = draw_figure(SCOPE, output_path=str(tmp_path / "scope"))
    texts = _texts(report)

    for expected in ("3a", "yield: 85%", "ee: 92%", "Substrate scope", "Conditions as reported by the user."):
        assert expected in texts


def test_a_missing_value_is_shown_as_a_dash_and_named_in_the_notes(tmp_path):
    report = draw_figure(SCOPE, output_path=str(tmp_path / "scope"))

    assert "ee: —" in _texts(report)
    assert any("3c" in note and "ee" in note for note in report["notes"])


def test_the_columns_follow_the_fields_the_brief_names(tmp_path):
    brief = _brief(fields=[{"key": "ee", "label": "ee", "unit": "%"}, "yield"])
    for entry in brief["species"]:
        entry["values"] = {k: v.rstrip("%") for k, v in entry["values"].items()}

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    texts = _texts(report)

    assert [f["key"] for f in report["fields"]] == ["ee", "yield"]
    assert "ee: 92 %" in texts
    assert "yield: 85" in texts


def test_groups_become_headings_in_first_seen_order(tmp_path):
    brief = _brief()
    for entry, group in zip(brief["species"], ("para", "ortho", "para", "ortho")):
        entry["group"] = group

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    texts = _texts(report)

    assert texts.index("para") < texts.index("ortho")


def test_many_compounds_continue_on_further_pages(tmp_path):
    brief = _brief(width_mm=90, height_mm=90)
    brief["species"] = [
        {"id": f"c{i}", "smiles": "c1ccccc1" + "C" * (i % 5), "label": str(i), "values": {"yield": f"{50 + i}%"}}
        for i in range(1, 25)
    ]

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))

    assert report["ok"] is True
    assert len(report["pages"]) > 1
    assert report["readback"]["status"] == "matched"
    assert sum(p["compounds"] for p in report["pages"]) == 24


def test_native_rendering_is_not_claimed_when_it_did_not_run(tmp_path):
    report = draw_figure(SCOPE, output_path=str(tmp_path / "scope"))

    assert report["native"]["status"] == "not_run"
    assert any("not been rendered" in n for n in report["notes"])


def test_a_native_render_lists_the_pngs(tmp_path, monkeypatch):
    def fake(stage, pages):
        (stage / "png").mkdir()
        for page in pages:
            (stage / "png" / page["file"].replace(".cdxml", ".png")).write_bytes(b"png")
        return {"status": "rendered_pending_visual_review", "formats": ["png", "svg"]}

    monkeypatch.setattr("cdxml_toolkit.mcp_runtime.publication_figure._native_previews", fake)

    report = draw_figure(SCOPE, output_path=str(tmp_path / "scope"), render=True)

    assert report["native"]["status"] == "rendered"
    assert [Path(p).name for p in report["native"]["png_paths"]] == ["page-001.png"]


def test_without_an_output_path_a_fresh_folder_is_made(tmp_path, monkeypatch):
    monkeypatch.setattr("tempfile.tempdir", str(tmp_path))

    first = draw_figure(SCOPE)
    second = draw_figure(SCOPE)

    assert first["output_dir"] != second["output_dir"]
    assert Path(first["output_dir"]).is_dir() and Path(second["output_dir"]).is_dir()


def test_an_existing_folder_is_never_overwritten(tmp_path):
    out = tmp_path / "scope"
    out.mkdir()
    (out / "keep.txt").write_text("mine", encoding="utf-8")

    report = draw_figure(SCOPE, output_path=str(out))

    assert report["ok"] is False
    assert (out / "keep.txt").read_text(encoding="utf-8") == "mine"
    assert not (out / "page-001.cdxml").exists()


def test_a_name_that_cannot_be_resolved_stops_the_run(tmp_path):
    brief = _brief()
    brief["species"][0] = {"id": "3a", "name": "no-such-compound-zzz", "values": {"yield": "85%"}}

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))

    assert report["ok"] is False
    assert "no-such-compound-zzz" in " ".join(report["errors"])
    assert not (tmp_path / "scope").exists()


@pytest.mark.parametrize(
    "change, message",
    [
        ({"steps": [{}]}, "steps"),
        ({"surprise": 1}, "Unknown brief keys"),
        ({"fields": [{"label": "x"}]}, "key"),
        ({"fields": ["yield", "yield"]}, "duplicate"),
        ({"footnotes": "a single string"}, "footnotes"),
        ({"width_mm": 5}, "width_mm"),
    ],
)
def test_a_malformed_scope_brief_is_refused_with_a_reason(tmp_path, change, message):
    report = draw_figure(_brief(**change), output_path=str(tmp_path / "scope"))

    assert report["ok"] is False
    assert message in " ".join(report["errors"])


def test_a_value_for_an_unlisted_field_is_refused(tmp_path):
    brief = _brief(fields=["yield"])

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))

    assert report["ok"] is False
    assert "ee" in " ".join(report["errors"])


def test_duplicate_ids_are_refused(tmp_path):
    brief = _brief()
    brief["species"][1]["id"] = "3a"

    report = draw_figure(brief, output_path=str(tmp_path / "scope"))

    assert report["ok"] is False
    assert "duplicate" in " ".join(report["errors"])


def test_the_scope_checks_are_the_publication_checks(tmp_path):
    report = draw_figure(SCOPE, output_path=str(tmp_path / "scope"))

    assert report["checks"]["data_binding"]["status"] == "passed"
    assert report["lint"]["issues"] == []


def test_the_brief_validator_names_the_scope_keys():
    assert {"fields", "conditions", "footnotes", "width_mm", "height_mm"} <= figure_scope.SCOPE_KEYS


# --- the varying group ----------------------------------------------------------------------

def _coloured_atoms(root, colour):
    """Element symbols (or 'C' for a bare carbon) of the atoms and bond ends drawn in ``colour``."""
    table = [tuple(round(float(c.get(k)) * 255) for k in "rgb") for c in root.find("colortable").findall("color")]
    want = tuple(int(colour[i:i + 2], 16) for i in (1, 3, 5))
    index = str(table.index(want) + 2) if want in table else None
    nodes = {n.get("id"): n for n in root.iter("n")}
    atoms = [n for n in nodes.values() if n.get("color") == index]
    bonds = [b for b in root.iter("b") if b.get("color") == index]
    return atoms, bonds, nodes


def test_highlight_varying_colours_what_differs_from_the_common_core(tmp_path):
    # BJOC 22-127, scheme 2: the substituent that changes is drawn in a second colour in every
    # product, so the reader need not compare structures to find it.
    report = draw_figure(_brief(highlight_varying=True), output_path=str(tmp_path / "scope"))
    assert report["ok"] is True
    assert report["highlighted"] == {"3a": 1, "3b": 2, "3c": 1, "3d": 1}  # CH3, OCH3, Cl, F
    (root,) = _pages(report)
    atoms, bonds, nodes = _coloured_atoms(root, "#FF6600")
    assert len(atoms) == 5
    # Each substituent's bond to the ring is coloured too, so a methyl drawn as a line shows.
    assert len(bonds) == 5
    ring = [n for n in nodes.values() if n.get("Element") in (None, "6") and n.get("color") is None]
    assert len(ring) >= 24  # the four benzene rings stay black


def test_a_species_own_colours_win_over_the_highlight(tmp_path):
    brief = _brief(highlight_varying="#C62828")
    brief["species"][2]["atom_colors"] = [{"atoms": [0], "color": "#336699"}]
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is True and "3c" not in report["highlighted"]
    (root,) = _pages(report)
    assert len(_coloured_atoms(root, "#C62828")[0]) == 4
    assert len(_coloured_atoms(root, "#336699")[0]) == 1


def test_highlight_varying_is_checked():
    for value in ("yes", "#FF66", 3):
        report = draw_figure(_brief(highlight_varying=value))
        assert report["ok"] is False and any("highlight_varying" in e for e in report["errors"])


def test_highlight_varying_needs_something_in_common(tmp_path):
    brief = _brief(highlight_varying=True)
    brief["species"] = [{"id": "x", "smiles": "CCO", "label": "1", "values": {"yield": "1%"}},
                        {"id": "y", "smiles": "c1ccncc1", "label": "2", "values": {"yield": "2%"}}]
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is True and report["highlighted"] == {}
    assert any("highlight_varying" in w for w in report["warnings"])


def test_a_secondary_product_is_drawn_beside_its_main_product(tmp_path):
    brief = _brief()
    brief["species"].insert(1, {"id": "3a-minor", "smiles": "Cc1cccc(C)c1", "label": "3a'",
                                "values": {"yield": "9%"}, "secondary_of": "3a"})
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is True and report["readback"]["status"] == "matched"
    assert report["checks_passed"] is True
    (root,) = _pages(report)
    table = [tuple(round(float(c.get(k)) * 255) for k in "rgb") for c in root.find("colortable").findall("color")]
    grey = str(table.index((0x80, 0x80, 0x80)) + 2)
    assert sum(f.get("color") == grey for f in root.iter("fragment")) == 1


def test_secondary_of_must_name_another_species():
    brief = _brief()
    brief["species"][0]["secondary_of"] = "nobody"
    report = draw_figure(brief)
    assert report["ok"] is False and any("secondary_of" in e for e in report["errors"])


def test_a_secondary_product_does_not_make_later_ids_look_duplicated():
    brief = {"kind": "scope", "fields": [{"key": "yield"}], "species": [
        {"id": "3a", "smiles": "CCO", "values": {"yield": "90%"}},
        {"id": "3a'", "smiles": "CCCO", "secondary_of": "3a", "values": {"yield": "5%"}},
        {"id": "3b", "smiles": "CCN", "values": {"yield": "80%"}},
    ]}
    assert figure_scope._validate(brief) == []


def test_a_definition_is_written_to_its_own_unnamed_line(tmp_path):
    brief = _brief()
    brief["species"][0]["definition"] = "Alk = C7H15"
    report = draw_figure(brief, output_path=str(tmp_path / "family"))
    assert report["ok"], report
    root = ET.parse(report["pages"][0]["path"]).getroot()   # a scope brief writes a directory, not a file
    lines = ["".join(s.text or "" for s in t.iter("s")) for t in root.iter("t")]
    assert "Alk = C7H15" in lines
    assert not any(line.startswith("definition") or line == "—" for line in lines)


def test_a_defined_group_is_drawn_as_its_label_and_reads_back_whole(tmp_path):
    # BJOC 22-127 scheme 2: 3p and 3q draw "NH-Alk" and say under it what Alk is, instead of a
    # heptyl or octadecyl chain that would widen the cell.
    brief = {"kind": "scope", "fields": ["yield"], "species": [
        {"id": "3p", "smiles": "CCCCCCCNC(=O)c1nc2ccccc2s1", "label": "3p", "definition": "Alk = C7H15",
         "values": {"yield": "63%"}},
        {"id": "3q", "smiles": "CCCCCCCCCCCCCCCCCCNC(=O)c1nc2ccccc2s1", "label": "3q",
         "definition": "Alk = C18H37", "values": {"yield": "60%"}},
        {"id": "3a", "smiles": "Cc1ccc(NC(=O)c2nc3ccccc3s2)cc1", "label": "3a", "values": {"yield": "64%"}}]}
    report = draw_figure(brief, output_path=str(tmp_path / "alk"))
    assert report["ok"], report
    assert report["readback"]["status"] == "matched"
    root = ET.parse(report["pages"][0]["path"]).getroot()
    labelled = [n for n in root.iter("n") if n.get("NodeType") == "Fragment"
                and "".join(n.find("t").itertext()).strip() == "Alk"]
    assert len(labelled) == 2
    project = json.loads((tmp_path / "alk" / "project.json").read_text(encoding="utf-8"))
    widths = {c: project["compounds"][c]["bbox"][2] - project["compounds"][c]["bbox"][0] for c in ("3p", "3q", "3a")}
    # The shared structure no longer widens the cell; the caption does — "3q, Alk = C18H37, 60%"
    # is longer than "3a, 88%". Measured 124.55 against 104.05, so 20 pt was 0.5 pt too tight.
    assert widths["3q"] <= widths["3a"] + 24


def test_a_species_drawn_with_another_lists_its_line_under_that_structure(tmp_path):
    # BJOC 22-127, fig. 3: 3q (Alk = C18H37) shares 3p's NH-Alk structure and adds one line.
    brief = _brief()
    brief["species"][0]["definition"] = "Alk = C7H15"
    brief["species"].insert(1, {"id": "3a-long", "smiles": "CCCCCCCCCCc1ccccc1", "label": "3e",
                                "definition": "Alk = C18H37", "values": {"yield": "60%"}, "drawn_with": "3a"})
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is True, report
    assert report["readback"]["status"] == "matched" and report["checks_passed"] is True
    (root,) = _pages(report)
    assert len(list(root.iter("fragment"))) == 4   # 3a, 3b, 3c, 3d; 3e has none of its own
    lines = ["".join(s.text or "" for s in t.iter("s")) for t in root.iter("t")]
    assert "3e" in lines and "Alk = C18H37" in lines


@pytest.mark.parametrize("drawn_with, message", [
    ("nobody", "another species"), ("3a", "drawn with its own"), ])
def test_drawn_with_must_name_a_species_drawn_on_its_own(drawn_with, message):
    brief = _brief()
    brief["species"][1]["drawn_with"] = "3a"
    brief["species"][2]["drawn_with"] = drawn_with if drawn_with == "nobody" else "3b"
    report = draw_figure(brief)
    assert report["ok"] is False and any(message in e for e in report["errors"]), report["errors"]


def test_header_r_marks_share_the_highlight_colour(tmp_path):
    # BJOC 22-127 scheme 2: R1 and R2 over the grid in the colour of what varies in it.
    brief = _brief(highlight_varying=True, header={
        "substrates": [{"smiles": "c1ccc([*:1])cc1", "label": "1"}],
        "products": [{"smiles": "Cc1ccc([*:1])cc1", "label": "3"}]})
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is True, report
    (root,) = _pages(report)
    table = [tuple(round(float(c.get(k)) * 255) for k in "rgb") for c in root.find("colortable").findall("color")]
    orange = str(table.index((0xFF, 0x66, 0x00)) + 2)
    marks = [n for n in root.iter("n") if n.get("NodeType") == "GenericNickname"]
    assert marks and all(n.get("color") == orange for n in marks)


@pytest.mark.parametrize("columns", [0, 9, 2.5, True, "4"])
def test_columns_is_a_small_whole_number(columns):
    report = draw_figure(_brief(columns=columns))
    assert report["ok"] is False and any("columns" in e for e in report["errors"]), report["errors"]


def test_four_columns_that_do_not_fit_report_the_arithmetic(tmp_path):
    # The cells are made too wide by their captions, which do not depend on how the structures
    # are depicted: a long value is wider than the 110.6 pt a column has at 180 mm.
    brief = _brief(columns=4)
    for entry in brief["species"]:
        entry["values"] = {key: "not determined after aqueous workup" for key in entry["values"]}
    report = draw_figure(brief, output_path=str(tmp_path / "scope"))
    assert report["ok"] is False
    (message,) = [e for e in report["errors"] if "columns do not fit" in e]
    for part in ("4 columns do not fit", "pt a column", "the page has 478.2 pt",
                 "width_mm to at least", "the type would be"):
        assert part in message, (part, message)
    assert not (tmp_path / "scope").exists()

def test_a_minority_compound_does_not_shrink_the_core_for_the_rest(tmp_path):
    # BJOC 22-127 scheme 2: a dozen N-aryl amides and three alkyl ones. The aryl ring is shared by
    # 75% of the compounds, so it is scaffold and the aryl rows colour their substituent only;
    # before the majority rule one alkyl row shrank the core and every aryl ring was coloured.
    brief = {"kind": "scope", "highlight_varying": True, "fields": [{"key": "yield"}],
             "species": [
                 {"id": "a", "smiles": "Cc1ccc(NC(=O)c2nc3ccccc3s2)cc1", "values": {"yield": "90%"}},
                 {"id": "b", "smiles": "COc1ccc(NC(=O)c2nc3ccccc3s2)cc1", "values": {"yield": "80%"}},
                 {"id": "c", "smiles": "Fc1ccc(NC(=O)c2nc3ccccc3s2)cc1", "values": {"yield": "70%"}},
                 {"id": "d", "smiles": "CCCCNC(=O)c2nc3ccccc3s2", "values": {"yield": "60%"}}]}
    report = draw_figure(brief, output_path=str(tmp_path / "mixed"))
    assert report["ok"], report
    assert report["highlighted"]["a"] == 1     # the methyl
    assert report["highlighted"]["b"] == 2     # O and C of the methoxy
    assert report["highlighted"]["c"] == 1     # the fluorine
    assert report["highlighted"]["d"] >= 4     # the whole butyl chain is what varies there

def test_an_empty_title_asks_for_no_title_at_all(tmp_path):
    # A journal figure carries no title; the caption does that job. "" used to fall back to the
    # default "Substrate scope" because the default was applied with `or`.
    brief = _brief()
    brief["title"] = ""
    report = draw_figure(brief, output_path=str(tmp_path / "untitled"))
    assert report["ok"], report
    root = ET.parse(report["pages"][0]["path"]).getroot()
    texts = ["".join(s.text or "" for s in t.iter("s")).strip() for t in root.iter("t")]
    assert "Substrate scope" not in texts
    assert all("scope" not in text.lower() for text in texts)

def test_depiction_rotation_turns_the_whole_table(tmp_path):
    # The first compound sets the orientation of every cell, so turning it turns the table: a source
    # figure with the benzothiazole on the left needs 180 where RDKit draws it on the right. The
    # vector from the ring's nitrogen to its sulfur says which way round the heterocycle sits.
    def brief(rotation=None):
        out = {"kind": "scope", "title": "", "fields": [{"key": "yield"}],
               "species": [{"id": "a", "smiles": "O=C(Nc1ccccc1)c1nc2ccccc2s1", "values": {"yield": "64%"}},
                           {"id": "b", "smiles": "Cc1ccccc1NC(=O)c1nc2ccccc2s1", "values": {"yield": "61%"}}]}
        if rotation is not None:
            out["depiction_rotation"] = rotation
        return out

    straight = draw_figure(brief(), output_path=str(tmp_path / "straight"))
    turned = draw_figure(brief(180), output_path=str(tmp_path / "turned"))
    assert turned["ok"], turned

    def n_to_s_angle(report):
        root = ET.parse(report["pages"][0]["path"]).getroot()
        where = {}
        for node in root.find(".//fragment").iter("n"):
            if node.get("p") and node.get("Element") in ("7", "16"):
                x, y = node.get("p").split()[:2]
                where[node.get("Element")] = (float(x), float(y))
        assert set(where) == {"7", "16"}, where
        (nx, ny), (sx, sy) = where["7"], where["16"]
        return math.degrees(math.atan2(sy - ny, sx - nx))

    change = abs((n_to_s_angle(straight) - n_to_s_angle(turned) + 180) % 360 - 180)
    assert 150 < change < 210, change

def test_a_general_reaction_keeps_a_short_real_r_bond(tmp_path):
    brief = _brief(rmark_in_ring=True, header={
        "substrates": [{"smiles": "c1ccc([*:1])cc1", "label": "1"}],
        "products": [{"smiles": "Cc1ccc([*:1])cc1", "label": "3"}]})
    report = draw_figure(brief, output_path=str(tmp_path / "inside"))
    assert report["ok"] is True, report
    (root,) = _pages(report)
    where = {n.get("id"): [float(v) for v in n.get("p").split()[:2]]
             for n in root.iter("n") if n.get("p")}
    marks = [n for n in root.iter("n") if n.get("NodeType") == "GenericNickname"]
    assert marks, "no R mark was drawn"
    for mark in marks:
        bonds = [b for b in root.iter("b") if mark.get("id") in (b.get("B"), b.get("E"))]
        assert len(bonds) == 1, "the R mark keeps exactly one bond"
        bond = bonds[0]
        other = bond.get("E") if bond.get("B") == mark.get("id") else bond.get("B")
        mx, my = where[mark.get("id")]
        ox, oy = where[other]
        # Header structures are drawn larger than the grid, so the stub scales with them.
        assert math.hypot(mx - ox, my - oy) < 22.0, "the mark is still a whole bond length away"

def test_an_r_mark_reads_r_then_a_raised_number(tmp_path):
    # "R1" written as R plus the superscript character. A face-64 run cannot be used: ChemDraw draws
    # one as a raised *prefix*, so the mark read "1R".
    brief = _brief(header={
        "substrates": [{"smiles": "c1ccc([*:1])cc1", "label": "1"}],
        "products": [{"smiles": "Cc1ccc([*:1])cc1", "label": "3"}]})
    report = draw_figure(brief, output_path=str(tmp_path / "runs"))
    assert report["ok"] is True, report
    (root,) = _pages(report)
    marks = [n for n in root.iter("n") if n.get("NodeType") == "GenericNickname"]
    assert marks, "no R mark was drawn"
    runs = [s for s in marks[0].iter("s")]
    assert runs, "no label runs"
    written = "".join(s.text or "" for s in runs)
    assert written.startswith("R"), written
    assert written[1:] in ("\u00b9", "\u00b2", "\u00b3"), written
    assert all(s.get("face") != "64" for s in runs), "a face-64 run would be drawn as a prefix"


def test_scope_ring_edge_marks_keep_the_bonded_r_in_chemical_readback(tmp_path):
    from cdxml_toolkit.chemistry_semantics import validate_fragment

    fixture = Path(__file__).resolve().parents[3] / "assets/parity/p22-four-columns-wide.json"
    source = json.loads(fixture.read_text(encoding="utf-8"))
    report = draw_figure(source, output_path=str(tmp_path / "ring-edge"))
    assert report["ok"] and report["checks"]["chemistry"]["status"] == "passed", report
    assert report["readback"]["status"] == "matched", report["readback"]
    (root,) = _pages(report)
    page = root.find("page")
    labels = [t for t in page.findall("t") if t.findtext("s") == "R\u00b9"]
    assert len(labels) == 2
    table = root.find("colortable").findall("color")
    fragments = [f for f in page.findall("fragment") if f.find('n[@GenericNickname="R1"]') is not None]
    assert len(fragments) == 2
    smiles = [source["header"]["substrates"][0]["smiles"], source["header"]["products"][0]["smiles"]]
    for fragment, label, structure in zip(fragments, labels, smiles):
        run = label.find("s")
        assert label.get("InterpretChemically") == "no"
        red = table[int(run.get("color")) - 2]
        assert tuple(float(red.get(k)) for k in "rgb") == (1.0, 0.0, 0.0)
        mark = fragment.find('n[@GenericNickname="R1"]')
        assert mark.get("Visible") == "no"
        assert mark.findtext("t/s") == "R1"
        nodes = {n.get("id"): tuple(map(float, n.get("p").split()[:2])) for n in fragment.findall("n")}
        bonds = [(b, b.get("B"), b.get("E")) for b in fragment.findall("b")]
        [(attachment, begin, end)] = [item for item in bonds if mark.get("id") in item[1:]]
        assert attachment.get("Visible") == "no"
        anchor = end if begin == mark.get("id") else begin
        (ax, ay) = nodes[anchor]
        ring_neighbours = [nodes[e if b == anchor else b] for bond, b, e in bonds
                           if anchor in (b, e) and bond.get("Order") == "2"]
        # The stroke is the normal of the edge it crosses: its midpoint sits on that edge and its
        # direction is square to it, at whatever angle the ring was depicted.
        strokes = []
        for line in page.findall('arrow[@ArrowheadHead="None"]'):
            x0, y0 = map(float, line.get("Tail3D").split()[:2])
            x1, y1 = map(float, line.get("Head3D").split()[:2])
            strokes.append(((x0 + x1) / 2, (y0 + y1) / 2, math.degrees(math.atan2(y1 - y0, x1 - x0))))
        assert strokes, "no ring-edge stroke was drawn"
        crossings = []
        for nx, ny in ring_neighbours:
            mid = ((ax + nx) / 2, (ay + ny) / 2)
            edge = math.degrees(math.atan2(ny - ay, nx - ax))
            for mx, my, stroke in strokes:
                # The stroke reaches further into the ring than out of it, so its midpoint sits
                # a point or so inside the edge midpoint.
                if math.dist(mid, (mx, my)) < 2.5:
                    crossings.append(abs(90 - abs((stroke - edge + 90) % 180 - 90)))
        assert crossings and min(crossings) < 5, (crossings, ring_neighbours)
        assert label.get("InterpretChemically") == "no"
        assert validate_fragment(structure, ET.tostring(fragment, encoding="unicode"))["status"] == "preserved"
        fragment.remove(attachment)
        with pytest.raises(ValueError, match="chemical_semantics_changed"):
            validate_fragment(structure, ET.tostring(fragment, encoding="unicode"))
    r2_marks = page.findall('.//n[@GenericNickname="R2"]')
    assert len(r2_marks) == 2 and all(n.get("Visible") != "no" for n in r2_marks)

def test_a_ring_edge_mark_crosses_its_edge_squarely():
    # The stroke is a decoration, but it still has to look drawn across one ring edge: its direction
    # must be that edge's normal, whatever angle the ring sits at.
    import math as _math

    from rdkit import Chem
    from rdkit.Chem import AllChem

    from cdxml_toolkit.publication.layout import _ring_edge_marks, style_settings

    mol = Chem.MolFromSmiles("c1ccc([*:1])cc1")
    AllChem.Compute2DCoords(mol)
    conformer = mol.GetConformer()
    angle = _math.radians(20)          # a ring at 20 degrees: an axis-aligned stroke would be askew
    c, s = _math.cos(angle), _math.sin(angle)
    xy = []
    for i in range(mol.GetNumAtoms()):
        x, y = conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y
        xy.append((x * c - y * s, x * s + y * c))

    style = style_settings({})
    marks = _ring_edge_marks(xy, mol, style["bond_length"], style)
    assert marks, "no ring mark was placed"
    _nickname, start, end, _position, _align, _label, edge = next(iter(marks.values()))
    stroke = _math.atan2(end[1] - start[1], end[0] - start[0])
    edge_angle = _math.atan2(edge[3] - edge[1], edge[2] - edge[0])
    difference = abs(((_math.degrees(stroke - edge_angle) + 90) % 180) - 90)
    assert abs(90 - difference) < 5, (edge_angle, stroke, difference)
