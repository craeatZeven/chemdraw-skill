"""Regression checks on the source content of the journal fixtures."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


BRIEFS = Path(__file__).resolve().parents[1] / "assets" / "scope-benchmark" / "briefs"


def brief(slug):
    return json.loads((BRIEFS / f"{slug}.json").read_text(encoding="utf-8"))


class JournalFidelityTests(unittest.TestCase):
    def test_beilstein_requests_a_fixed_four_column_grid_and_left_scaffold(self):
        spec = brief("bjoc-articles-22-127-scheme2")
        self.assertEqual(spec["columns"], 4)
        self.assertNotIn("row_align", spec)
        self.assertIn("depiction_rotation", spec)
        self.assertEqual([entry["id"] for entry in spec["species"]], [f"3{letter}" for letter in "abcdefghijklmnopqrst"])

    def test_helicenes_include_the_six_ring_series_and_two_source_insets(self):
        spec = brief("nature-chemistry-s41557-023-01174-5-fig3")
        series = [entry["id"] for entry in spec["species"] if entry["group"] == "[6]Helicenes"]
        self.assertEqual(series, [f"2{suffix}" for suffix in "rstuvwx"])
        images = [item for item in spec["footer_scheme"]["objects"] if item["type"] == "image"]
        self.assertEqual([item["id"] for item in images], ["xray-2p", "xray-2q"])

    def test_compare_includes_the_last_page(self):
        try:
            from PIL import Image
        except ImportError:
            self.skipTest("Pillow is required for native image comparisons")
        from scripts.run_scope_benchmark import _stitch_side_by_side

        with TemporaryDirectory() as directory:
            original, first, last, output = (Path(directory) / f"{name}.png" for name in ("orig", "1", "2", "out"))
            Image.new("RGB", (40, 40), "white").save(original)
            Image.new("RGB", (40, 40), "blue").save(first)
            Image.new("RGB", (60, 40), "red").save(last)
            self.assertTrue(_stitch_side_by_side(original, first, output, "Figure", additional_drawn_pngs=(last,)))
            with Image.open(output) as comparison:
                self.assertEqual(comparison.getpixel((comparison.width - 60, comparison.height - 80)), (255, 0, 0))
