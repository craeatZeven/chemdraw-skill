"""Regression checks for the published journal benchmark fixtures."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


BRIEFS = Path(__file__).resolve().parents[1] / "assets" / "scope-benchmark" / "briefs"
BLUE = "#1B4D89"


def brief(slug):
    return json.loads((BRIEFS / f"{slug}.json").read_text(encoding="utf-8"))


class JournalFidelityTests(unittest.TestCase):
    def test_thianthrenium_scope_preserves_chemical_labels_and_provenance(self):
        spec = brief("nature-synthesis-s44160-025-00821-8-fig4")
        entries = {entry["id"]: entry for entry in spec["species"]}
        self.assertEqual(len(entries), 20)
        self.assertEqual(set(entries), {"2'", *map(str, range(3, 22))})
        self.assertEqual(entries["9"]["abbreviate"], ["NCS"])
        self.assertEqual(entries["11"]["abbreviate"], ["N3"])
        self.assertEqual(entries["20"]["note"], "From *N*-Ac-L-cysteine")
        self.assertEqual(spec["header"]["products"][0]["variable_label_color"], BLUE)
        runs = spec["header"]["above"][0]
        self.assertIn({"text": "H\u2013Nu", "color": BLUE}, runs)

    def test_iron_scope_preserves_synthon_boundaries_and_primary_section(self):
        spec = brief("nature-synthesis-s44160-025-00771-1-fig4")
        entries = spec["species"]
        self.assertEqual(len(entries), 24)
        groups = {entry["id"]: entry["group"] for entry in entries}
        self.assertIn("secondary", groups["48"])
        self.assertTrue(all("primary" in groups[str(i)] for i in range(49, 54)))
        header = spec["header_scheme"]["objects"]
        self.assertEqual([item["text"] for item in header if item.get("id", "").startswith("synthon-")],
                         ["Ar", "alkyl", "Ar", "alkyl"])
        self.assertIn("fe-marker", {item.get("id") for item in header})
        molecules = {item["id"]: item for item in spec["footer_scheme"]["objects"] if item["type"] == "molecule"}
        self.assertEqual(molecules["55"]["color"], "#800000")
        self.assertEqual(molecules["55"]["highlight_color"], "#0B417C")
        for label in ("54", "55"):
            coords = molecules[label]["coordinates"]
            self.assertGreater(coords[16][1], coords[6][1])

    def test_beilstein_has_four_columns_with_source_scaffold_orientation(self):
        spec = brief("bjoc-articles-22-127-scheme2")
        self.assertEqual(spec["columns"], 4)
        self.assertNotIn("row_align", spec)
        reference = spec["species"][0]["coordinates"]
        self.assertLess(reference[16][1], reference[9][1])
        self.assertLess(reference[16][0], reference[5][0])
        self.assertLess(reference[7][1], reference[6][1])
        self.assertEqual([entry["id"] for entry in spec["species"]],
                         [f"3{letter}" for letter in "abcdefghijklmnopqrst"])

    def test_helicenes_include_the_six_ring_series_and_two_source_insets(self):
        spec = brief("nature-chemistry-s41557-023-01174-5-fig3")
        series = [entry["id"] for entry in spec["species"] if entry["group"] == "[6]Helicenes"]
        self.assertEqual(series, [f"2{suffix}" for suffix in "rstuvwx"])
        images = [item for item in spec["footer_scheme"]["objects"] if item["type"] == "image"]
        self.assertEqual([item["id"] for item in images], ["xray-2p", "xray-2q"])
        product = next(item for item in spec["footer_scheme"]["objects"] if item.get("id") == "2z")
        self.assertIn("coordinates", product)
        self.assertTrue(all("\u2091\u2099" in item["values"]["barrier"] for item in spec["species"]
                            if "barrier" in item["values"]))

    def test_hexahelicene_2r_uses_the_supporting_information_graph(self):
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is required to verify the supporting-information structure")

        spec = brief("nature-chemistry-s41557-023-01174-5-fig3")
        entry = next(item for item in spec["species"] if item["id"] == "2r")
        molecule = Chem.MolFromSmiles(entry["smiles"])
        expected = Chem.MolFromSmiles("Fc1c(F)c2ccc3ccc4ccccc4c3c2c2c1ccc1ccccc12")
        self.assertEqual(Chem.MolToSmiles(molecule), Chem.MolToSmiles(expected))
        self.assertEqual(len(entry["coordinates"]), molecule.GetNumAtoms())

    def test_comparison_includes_every_redraw_page(self):
        try:
            from PIL import Image
        except ImportError:
            self.skipTest("Pillow is required for image comparisons")
        from scripts.run_scope_benchmark import _stitch_side_by_side

        with TemporaryDirectory() as directory:
            original, first, last, output = (Path(directory) / f"{name}.png"
                                              for name in ("original", "first", "last", "comparison"))
            Image.new("RGB", (40, 40), "white").save(original)
            Image.new("RGB", (40, 40), "blue").save(first)
            Image.new("RGB", (60, 40), "red").save(last)
            self.assertTrue(_stitch_side_by_side(original, first, output, "Figure", additional_drawn_pngs=(last,)))
            with Image.open(output) as comparison:
                self.assertEqual(comparison.getpixel((comparison.width - 60, comparison.height - 80)), (255, 0, 0))
