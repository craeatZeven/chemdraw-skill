"""Tests for publication palettes, profile-aware checking, retrosynthesis arrows,
first-class optimization/SAR tables, per-band scope colouring, and phantom-fragment guards."""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure
from cdxml_toolkit.mcp_runtime.figure_validation import compare_inventories, inspect_document
from cdxml_toolkit.publication.house_styles import APPLIED_FIGURE_CONVENTIONS
from cdxml_toolkit.publication.palettes import PALETTES, resolve_palette


class TestPublicationPalettes:
    def test_builtin_palettes_have_required_roles(self):
        for name in ("okabe_ito", "nature", "acs_jacs", "wiley_angew", "beilstein", "rsc"):
            assert name in PALETTES
            p = resolve_palette(name)
            for key in ("primary", "secondary", "accent", "new_bond", "muted", "band_colors", "allowed"):
                assert key in p, f"{name} missing {key}"
            assert "#000000" in p["allowed"]
            assert "#FFFFFF" in p["allowed"]
            assert "#808080" in p["allowed"]
            assert len(p["band_colors"]) >= 3

    def test_resolve_palette_from_house_style_collects_card_colors(self):
        p = resolve_palette("nature_portfolio")
        assert "#0078B2" in p["allowed"]
        assert "#FF0000" in p["allowed"]
        p_jacs = resolve_palette("acs_jacs")
        assert "#800000" in p_jacs["allowed"]
        assert "#004080" in p_jacs["allowed"]


class TestProfileAwareChecker:
    def test_check_accepts_house_style_profile_and_palette(self, tmp_path):
        out = tmp_path / "jacs.cdxml"
        brief = {
            "style": "acs_jacs",
            "species": [
                {"id": "A", "smiles": "c1ccccc1Br", "label": "1"},
                {"id": "B", "smiles": "OB(O)c1ccc(C)cc1", "label": "2"},
                {"id": "P", "smiles": "Cc1ccc(-c2ccccc2)cc1", "label": "3"},
            ],
            "steps": [
                {"substrates": ["A", "B"], "products": ["P"],
                 "above": ["Pd(PPh3)4 (5 mol %)"], "below": ["K2CO3, DMF, 80 °C"], "yield": "88%"}
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        report = check_cdxml(out, profile="acs_jacs")
        assert report["ok"], report["errors"]
        assert report["profile"] == "acs_jacs"

    def test_check_accepts_measurement_preset_profile(self, tmp_path):
        out = tmp_path / "wiley.cdxml"
        brief = {
            "style": "wiley",
            "species": [
                {"id": "A", "smiles": "CCO", "label": "1"},
                {"id": "B", "smiles": "CC=O", "label": "2"},
            ],
            "steps": [{"substrates": ["A"], "products": ["B"], "above": ["DMP"], "below": ["CH2Cl2"]}],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        report = check_cdxml(out, profile="wiley")
        assert report["ok"], report["errors"]

    def test_check_accepts_named_palette_string(self, tmp_path):
        out = tmp_path / "okabe.cdxml"
        brief = {
            "species": [
                {"id": "A", "smiles": "c1ccccc1O", "label": "1", "color": "#0072B2"},
                {"id": "B", "smiles": "c1ccccc1OC", "label": "2", "color": "#D55E00"},
            ],
            "steps": [{"substrates": ["A"], "products": ["B"], "above": ["MeI, K2CO3"]}],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        report = check_cdxml(out, profile="acs1996", palette="okabe_ito")
        assert report["ok"], report["errors"]


class TestPhantomFragmentGuard:
    def test_condition_text_marks_interpret_chemically_no(self, tmp_path):
        out = tmp_path / "mech.cdxml"
        brief = {
            "kind": "mechanism",
            "species": [
                {"id": "A", "smiles": "CC=O", "label": "1"},
                {"id": "OH", "smiles": "[OH-]"},
                {"id": "P", "smiles": "CC([O-])O", "label": "2"},
            ],
            "steps": [{
                "substrates": ["A", "OH"],
                "products": ["P"],
                "above": ["H2O"],
                "moves": [
                    {"from": {"species": "OH", "atom": 0}, "to": {"species": "A", "atom": 1}},
                    {"from": {"species": "A", "bond": [1, 2]}, "to": {"species": "A", "atom": 2}},
                ],
            }],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        root = ET.parse(out).getroot()
        atom_texts = {id(t) for f in root.iter("fragment") for t in f.iter("t")}
        free_texts = [t for t in root.iter("t") if id(t) not in atom_texts]
        assert free_texts
        for t in free_texts:
            assert t.get("InterpretChemically") == "no", ET.tostring(t, encoding="unicode")

    def test_compare_inventories_flags_added_phantom_molecule(self):
        before = {
            "molecules": [
                {"smiles": "CC=O", "connectivity": "CC=O", "semantic_key": "CC=O",
                 "formula": "C2H4O", "formal_charge": 0, "assigned_centers": 0},
            ],
        }
        after = {
            "molecules": [
                {"smiles": "CC=O", "connectivity": "CC=O", "semantic_key": "CC=O",
                 "formula": "C2H4O", "formal_charge": 0, "assigned_centers": 0},
                {"smiles": "O", "connectivity": "O", "semantic_key": "O",
                 "formula": "H2O", "formal_charge": 0, "assigned_centers": 0},
            ],
        }
        diff = compare_inventories(before, after)
        assert diff["status"] == "changed"
        assert diff["molecule_count_delta"] == 1
        assert diff["added_molecules"] == ["O"]
        assert diff["removed_molecules"] == []


class TestRetrosynthesisArrows:
    def test_draw_figure_supports_retro_arrow(self, tmp_path):
        out = tmp_path / "retro.cdxml"
        brief = {
            "species": [
                {"id": "T", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "1"},
                {"id": "S1", "smiles": "Oc1ccccc1C(=O)O", "label": "2"},
                {"id": "S2", "smiles": "CC(=O)OC(C)=O", "label": "3"},
            ],
            "steps": [
                {"substrates": ["T"], "products": ["S1", "S2"],
                 "arrow": "retro", "above": ["acetylation"]},
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["checks_passed"], res["warnings"]
        root = ET.parse(out).getroot()
        arrows = list(root.iter("arrow"))
        assert len(arrows) == 1
        assert arrows[0].get("LineType") == "Double"
        assert arrows[0].get("ArrowheadType") == "Angle"
        assert arrows[0].get("ArrowShaftSpacing") == "300"


class TestOptimizationTableFigure:
    def test_draw_optimization_table_clean_booktabs_with_header_and_highlight(self, tmp_path):
        out = tmp_path / "opt_table.cdxml"
        brief = {
            "kind": "table",
            "title": "Optimization of the coupling conditions",
            "header": {
                "substrates": [
                    {"smiles": "Brc1ccccc1", "label": "1a"},
                    {"smiles": "OB(O)c1ccc(C)cc1", "label": "2a"},
                ],
                "products": [{"smiles": "Cc1ccc(-c2ccccc2)cc1", "label": "3a"}],
                "above": ["Pd(OAc)2 (2 mol %)"],
                "below": ["base, solvent, T, 12 h"],
            },
            "columns": [
                {"key": "entry", "label": "Entry"},
                {"key": "base", "label": "Base"},
                {"key": "solvent", "label": "Solvent"},
                {"key": "temp", "label": "T (°C)"},
                {"key": "yield", "label": "Yield (%)^a"},
            ],
            "rows": [
                {"entry": "1", "base": "K2CO3", "solvent": "DMF", "temp": "80", "yield": "62"},
                {"entry": "2", "base": "Cs2CO3", "solvent": "DMF", "temp": "80", "yield": "88"},
                {"entry": "3", "base": "Cs2CO3", "solvent": "1,4-dioxane", "temp": "100",
                 "yield": "94^b", "highlight": True},
                {"entry": "4", "base": "DBU", "solvent": "MeCN", "temp": "80", "yield": "trace"},
            ],
            "footnotes": [
                "^a Determined by 1H NMR using 1,3,5-trimethoxybenzene as internal standard.",
                "^b Isolated yield.",
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["kind"] == "table"
        assert res["checks_passed"], res["warnings"]
        assert res["readback"]["status"] == "matched"
        # Check with house-rule checker: zero text_over_ink, zero errors!
        chk = check_cdxml(out, profile="acs1996", allow_sizes=(8,))
        assert chk["ok"], chk["errors"]
        assert not any(w["code"] == "text_over_ink" for w in chk["warnings"]), chk["warnings"]
        # Verify chemical subscripts (K2CO3 -> subscript 2 and 3) and superscript footnote (^a, ^b)
        root = ET.parse(out).getroot()
        sub_runs = [s.text for s in root.iter("s") if int(s.get("face", "0") or 0) & 32]
        sup_runs = [s.text for s in root.iter("s") if int(s.get("face", "0") or 0) & 64]
        assert "2" in sub_runs and "3" in sub_runs
        assert "a" in sup_runs and "b" in sup_runs
        # All table caption <t> elements must carry InterpretChemically="no"
        atom_texts = {id(t) for f in root.iter("fragment") for t in f.iter("t")}
        for t in root.iter("t"):
            if id(t) not in atom_texts:
                assert t.get("InterpretChemically") == "no"

    def test_optimization_table_applies_house_style_conventions(self, tmp_path):
        out = tmp_path / "opt_nature.cdxml"
        brief = {
            "kind": "table",
            "style": "nature_portfolio",
            "columns": ["Entry", "Catalyst", "Yield (%)"],
            "rows": [
                {"Entry": "1", "Catalyst": "L1", "Yield (%)": "45"},
                {"Entry": "2", "Catalyst": "L2", "Yield (%)": "92", "highlight": True},
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["style"]["house"] == "nature_portfolio"
        assert "optimization_table" in APPLIED_FIGURE_CONVENTIONS


class TestPerBandScopeTheming:
    def test_scope_group_colors_tint_bands_and_varying_groups(self, tmp_path):
        out_dir = tmp_path / "banded_scope"
        brief = {
            "kind": "scope",
            "title": "Banded substrate scope",
            "highlight_varying": True,
            "group_colors": {
                "para substituents": "#0078B2",
                "C-nucleophiles": "#D55E00",
            },
            "fields": [{"key": "yield", "label": "Yield"}],
            "species": [
                {"id": "1a", "smiles": "Cc1ccc(C(=O)N)cc1", "label": "1a",
                 "group": "para substituents", "values": {"yield": "88%"}},
                {"id": "1b", "smiles": "COc1ccc(C(=O)N)cc1", "label": "1b",
                 "group": "para substituents", "values": {"yield": "91%"}},
                {"id": "2a", "smiles": "FC(F)(F)c1ccc(C(=O)N)cc1", "label": "2a",
                 "group": "C-nucleophiles", "values": {"yield": "76%"}},
            ],
        }
        res = draw_figure(brief, str(out_dir), render=False)
        assert res["ok"], res
        assert res["checks_passed"], res["warnings"]
        page_path = Path(res["pages"][0]["path"])
        chk = check_cdxml(page_path, profile="acs1996", palette="nature")
        assert chk["ok"], chk["errors"]
        used_colors = set(chk["summary"]["colors"].keys())
        assert "#0078B2" in used_colors
        assert "#D55E00" in used_colors

    def test_scope_house_style_variation_color_applied(self):
        assert "variation_color" in APPLIED_FIGURE_CONVENTIONS.get("substrate_scope", ())
