"""TDD tests for Phase 2 multi-type publication figures:
- Multi-panel side-by-side layout, rule-flanked titles, accent colour, embedded tables & energy profiles, TOC check
- First-class kind: "energy_profile" in draw_figure with house-style / palette integration
- Retrosynthesis wavy/dashed bond disconnection marks without breaking chemistry readback
- Italic stereochemical prefixes on compound labels (cis-3a, rac-2a, E-1a)
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure
from cdxml_toolkit.publication.house_styles import APPLIED_FIGURE_CONVENTIONS
from cdxml_toolkit.publication.layout import caption_runs
from cdxml_toolkit.scientific.energy import energy_profile


class TestEnergyProfileIntegration:
    def test_scientific_energy_profile_uses_named_palette(self, tmp_path):
        out = tmp_path / "ep_pal.cdxml"
        spec = {
            "palette": "okabe_ito",
            "style": "levels",
            "unit": "kcal/mol",
            "pathways": [
                {"label": "Path A", "points": [
                    {"id": "R", "label": "SM", "energy": 0.0},
                    {"id": "TS1", "label": "TS1", "energy": 18.4, "kind": "ts"},
                    {"id": "P", "label": "Prod", "energy": -6.2},
                ]},
                {"label": "Path B", "points": [
                    {"id": "R", "label": "SM", "energy": 0.0},
                    {"id": "TS2", "label": "TS2", "energy": 24.1, "kind": "ts"},
                    {"id": "P", "label": "Prod", "energy": -6.2},
                ]},
            ],
        }
        res = energy_profile(spec, output=str(out))
        assert Path(res["output_path"]).is_file()
        chk = check_cdxml(out, profile="acs1996", palette="okabe_ito", allow_sizes=(8, 9))
        assert chk["ok"], chk["errors"]
        used = set(chk["summary"]["colors"].keys())
        assert "#0072B2" in used
        assert "#D55E00" in used

    def test_draw_figure_supports_kind_energy_profile_with_house_style(self, tmp_path):
        out = tmp_path / "draw_ep.cdxml"
        brief = {
            "kind": "energy_profile",
            "style": "houk_group",
            "profile_style": "levels",
            "unit": "kcal/mol",
            "pathways": [
                {"label": "major", "points": [
                    {"id": "SM", "label": "1", "energy": 0.0},
                    {"id": "TS", "label": "TS1", "energy": 15.2, "kind": "ts"},
                    {"id": "PR", "label": "2", "energy": -10.5},
                ]},
                {"label": "minor", "points": [
                    {"id": "SM", "label": "1", "energy": 0.0},
                    {"id": "TS2", "label": "TS2", "energy": 21.0, "kind": "ts"},
                    {"id": "PR", "label": "2'", "energy": -8.1},
                ]},
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["kind"] == "energy_profile"
        assert res["checks_passed"], res["warnings"]
        assert Path(res["output_path"]).is_file()
        assert "energy_profile" in APPLIED_FIGURE_CONVENTIONS


class TestMultiPanelUpgrades:
    def test_panels_side_by_side_with_scheme_table_and_energy_profile(self, tmp_path):
        out = tmp_path / "multi_fig.cdxml"
        brief = {
            "kind": "panels",
            "layout": "side_by_side",
            "columns": 2,
            "title_style": "rule",
            "panels": [
                {
                    "title": "(a) Model reaction",
                    "brief": {
                        "species": [
                            {"id": "A", "smiles": "Brc1ccccc1", "label": "1"},
                            {"id": "P", "smiles": "Cc1ccccc1", "label": "2"},
                        ],
                        "steps": [{"substrates": ["A"], "products": ["P"], "above": ["MeMgBr, [Ni]"]}],
                    },
                },
                {
                    "title": "(b) Optimization",
                    "accent": True,
                    "brief": {
                        "kind": "table",
                        "columns": ["Entry", "Ligand", "Yield (%)"],
                        "rows": [
                            {"Entry": "1", "Ligand": "L1", "Yield (%)": "42"},
                            {"Entry": "2", "Ligand": "L2", "Yield (%)": "93", "highlight": True},
                        ],
                    },
                },
                {
                    "title": "(c) Free-energy profile",
                    "brief": {
                        "kind": "energy_profile",
                        "profile_style": "levels",
                        "unit": "kcal/mol",
                        "pathways": [
                            {"label": "Ni", "points": [
                                {"id": "R", "label": "1", "energy": 0.0},
                                {"id": "TS", "label": "TS", "energy": 14.0, "kind": "ts"},
                                {"id": "P", "label": "2", "energy": -8.0},
                            ]},
                        ],
                    },
                },
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["layout"] == "side_by_side"
        assert len(res["panels"]) == 3
        assert res["checks_passed"], res["warnings"]
        # Verify panel (b) sits to the right of panel (a) (x > 0)
        assert res["panel_positions"][1][0] > res["panel_positions"][0][0]
        assert "concept_panel" in APPLIED_FIGURE_CONVENTIONS

    def test_panels_toc_check_reports_dimensions(self, tmp_path):
        out = tmp_path / "toc.cdxml"
        brief = {
            "kind": "panels",
            "title_style": "plain",
            "toc": True,
            "panels": [
                {
                    "title": "TOC summary",
                    "brief": {
                        "species": [
                            {"id": "A", "smiles": "CCO", "label": "1"},
                            {"id": "B", "smiles": "CC=O", "label": "2"},
                        ],
                        "steps": [{"substrates": ["A"], "products": ["B"], "above": ["[O]"]}],
                    },
                },
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert "toc" in res
        assert "width_pt" in res["toc"] and "height_pt" in res["toc"]
        assert isinstance(res["toc"]["fits_acs_toc"], bool)


class TestRetrosynthesisDisconnectionsAndStereoLabels:
    def test_species_disconnection_draws_wavy_cut_preserving_readback(self, tmp_path):
        out = tmp_path / "retro_disc.cdxml"
        brief = {
            "species": [
                {
                    "id": "T",
                    "smiles": "CC(=O)Oc1ccccc1",
                    "label": "1",
                    "disconnections": [{"bond": [1, 3], "color": "#DC3220", "style": "wavy"}],
                },
                {"id": "S1", "smiles": "CC(=O)Cl", "label": "2"},
                {"id": "S2", "smiles": "Oc1ccccc1", "label": "3"},
            ],
            "steps": [
                {"substrates": ["T"], "products": ["S1", "S2"], "arrow": "retro", "above": ["esterification"]},
            ],
        }
        res = draw_figure(brief, str(out), render=False)
        assert res["ok"], res
        assert res["readback"]["status"] == "matched", res["readback"]
        assert res["checks_passed"], res["warnings"]
        root = ET.parse(out).getroot()
        curves = list(root.iter("curve"))
        assert len(curves) >= 1, "Expected wavy disconnection <curve>"
        # CurvePoints must have a valid multiple of 6 numbers so ChemDraw never drops it
        pts = (curves[0].get("CurvePoints") or "").split()
        assert len(pts) >= 6 and len(pts) % 6 == 0
        assert "retrosynthesis" in APPLIED_FIGURE_CONVENTIONS

    def test_caption_runs_italicizes_stereochemical_prefix_on_bold_compound_label(self):
        plain, runs = caption_runs("cis-3a", bold=True)
        assert plain == "cis-3a"
        assert runs is not None
        assert runs[0] == {"text": "cis", "italic": True}
        assert any(r.get("bold") and "3a" in r["text"] for r in runs)

        plain2, runs2 = caption_runs("rac-2a", bold=True)
        assert plain2 == "rac-2a"
        assert runs2 is not None
        assert runs2[0] == {"text": "rac", "italic": True}
