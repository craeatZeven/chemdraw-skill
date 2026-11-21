from __future__ import annotations

import argparse
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from scripts import run_scope_benchmark as benchmark

try:
    from PIL import Image
except ImportError:
    Image = None


class ScopeBenchmarkTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)

    def _run_journal(self, native_by_case, *, render):
        cases = []
        for slug in native_by_case:
            (self.root / f"{slug}.json").write_text(
                json.dumps({"kind": "scope", "slug": slug}), encoding="utf-8"
            )
            cases.append({
                "slug": slug,
                "journal": "Test Journal",
                "figure": "Scheme 1",
                "archetype": "classic_header_grid",
                "check_profile": "acs1996",
                "brief": f"{slug}.json",
                "corpus_png": "reference.png",
            })

        def draw_scope(brief, output_path, render):
            return {
                "ok": True,
                "checks_passed": True,
                "readback": {"status": "matched"},
                "native": native_by_case[brief["slug"]],
            }

        modules = {
            "cdxml_toolkit.cdxml_check": SimpleNamespace(check_cdxml=None),
            "cdxml_toolkit.mcp_runtime.figure_draw": SimpleNamespace(draw_figure=None),
            "cdxml_toolkit.mcp_runtime.figure_scope": SimpleNamespace(draw_scope=draw_scope),
        }
        out = self.root / "results"
        args = argparse.Namespace(only=None, limit=None, output_dir=str(out),
                                  render=render, side_by_side=False)
        with patch.dict("sys.modules", modules), patch.dict("os.environ"), \
                patch.object(benchmark, "ROOT", self.root), \
                patch.object(benchmark, "_load_manifest", return_value={"cases": cases}), \
                redirect_stdout(io.StringIO()):
            self.assertEqual(benchmark._run_in_toolkit(args), 0)
        return json.loads((out / "benchmark-report.json").read_text(encoding="utf-8"))

    def test_journal_no_render_leaves_metrics_null(self):
        report = self._run_journal({"audit": {"status": "not_run"}}, render=False)
        case, = report["cases"]
        self.assertTrue(case["passed"])
        self.assertIsNone(case["pixel_metrics"])
        self.assertIsNone(case["pixel_metrics_normalized"])
        self.assertIsNone(case["visual_metrics_error"])
        self.assertEqual(report["pixel_metrics_normalized"]["cases_measured"], 0)
        self.assertIsNone(report["pixel_metrics_normalized"]["rmse_median"])

    @unittest.skipUnless(Image is not None, "Pillow is required to compare PNGs")
    def test_journal_measures_without_side_by_side_and_reports_unavailable_pairs(self):
        reference = self.root / "reference.png"
        Image.new("RGB", (10, 10), "black").save(reference)
        narrow = self.root / "narrow.png"
        Image.new("RGB", (10, 20), "black").save(narrow)
        corrupt = self.root / "corrupt.png"
        corrupt.write_bytes(b"not a PNG")
        report = self._run_journal({
            "narrow": {"status": "rendered", "png_paths": [str(narrow), str(reference)]},
            "corrupt": {"status": "rendered", "png_paths": [str(corrupt)]},
            "no-png": {"status": "not_run", "reason": "ChemDraw unavailable"},
            "identical": {"status": "rendered", "png_path": str(reference)},
        }, render=True)

        cases = {case["slug"]: case for case in report["cases"]}
        self.assertEqual(report["passed"], 4)
        self.assertEqual(report["diff_threshold_grey"], 16)
        self.assertIsNone(cases["narrow"]["side_by_side"])
        raw = cases["narrow"]["pixel_metrics"]
        self.assertFalse(raw["size_match"])
        self.assertIsNone(raw["rmse"])
        normalized = cases["narrow"]["pixel_metrics_normalized"]
        self.assertEqual(normalized["canvas_size"], [900, 900])
        self.assertEqual(normalized["render_scaled_size"], [450, 900])
        self.assertEqual(normalized["aspect_ratio_delta"], 0.5)
        self.assertEqual(normalized["diff_pixel_ratio"], 0.5)
        self.assertAlmostEqual(normalized["rmse"], 180.312, places=3)
        self.assertEqual(cases["identical"]["pixel_metrics"]["rmse"], 0.0)
        self.assertEqual(cases["identical"]["pixel_metrics_normalized"]["rmse"], 0.0)

        self.assertIsNone(cases["corrupt"]["pixel_metrics_normalized"])
        self.assertIsNotNone(cases["corrupt"]["visual_metrics_error"])
        self.assertIsNone(cases["no-png"]["pixel_metrics_normalized"])
        self.assertEqual(cases["no-png"]["visual_metrics_error"], "ChemDraw unavailable")
        metrics = report["pixel_metrics_normalized"]
        self.assertEqual(metrics["cases_measured"], 2)
        self.assertEqual(metrics["rmse_median"], 90.156)
        self.assertEqual(metrics["diff_pixel_ratio_median"], 0.25)
        self.assertEqual(metrics["rmse_worst_case"], "narrow")
        self.assertEqual(metrics["diff_pixel_ratio_worst_case"], "narrow")

    @unittest.skipUnless(Image is not None, "Pillow is required to compare PNGs")
    def test_journal_composites_transparency_and_applies_threshold_strictly(self):
        reference = self.root / "reference.png"
        render = self.root / "render.png"
        Image.new("RGBA", (2, 1), (0, 0, 0, 0)).save(reference)
        image = Image.new("RGB", (2, 1))
        image.putdata([(239, 239, 239), (238, 238, 238)])
        image.save(render)

        raw, _ = benchmark._journal_pixel_metrics(reference, render)
        self.assertEqual(raw["rmse"], 16.508)
        self.assertEqual(raw["diff_pixel_ratio"], 0.5)
