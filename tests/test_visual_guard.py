from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from tests import check_visual


class VisualGuardTests(unittest.TestCase):
    def test_size_threshold_is_inclusive_for_each_dimension(self):
        golden = {"size_px": [100, 200], "perceptual_hash": "0", "marks": []}
        within_tolerance = {**golden, "size_px": [102, 196]}
        self.assertFalse(check_visual.compare(golden, within_tolerance)[0])
        for size, dimension in (([103, 200], "width"), ([100, 195], "height")):
            with self.subTest(size=size):
                outside_tolerance = {**golden, "size_px": size}
                errors, _ = check_visual.compare(golden, outside_tolerance)
                self.assertIn(dimension, " ".join(errors))

    def test_hash_distance_counts_differing_bits(self):
        golden = {"size_px": [100, 200], "perceptual_hash": "0", "marks": []}
        actual = {**golden, "perceptual_hash": "f"}
        self.assertFalse(check_visual.compare(golden, actual, hash_limit=4)[0])
        self.assertIn("hash", " ".join(check_visual.compare(golden, actual, hash_limit=3)[0]))

    def test_hash_detects_recoloring_at_fixed_size(self):
        try:
            from PIL import Image, ImageDraw
        except ImportError:
            self.skipTest("Pillow is required for visual fingerprint tests")
        red = Image.new("RGB", (128, 128), "white")
        blue = red.copy()
        ImageDraw.Draw(red).line((16, 64, 112, 64), fill="#c81e1e", width=4)
        ImageDraw.Draw(blue).line((16, 64, 112, 64), fill="#1e1ec8", width=4)
        digest = check_visual.perceptual_hash(red)
        self.assertEqual(digest, check_visual.perceptual_hash(red.resize((256, 256))))
        self.assertGreater(check_visual.hamming_distance(digest, check_visual.perceptual_hash(blue)),
                           check_visual.HASH_LIMIT)

    def test_paper_header_hash_flags_a_change_that_the_page_hash_misses(self):
        try:
            from PIL import Image, ImageDraw
        except ImportError:
            self.skipTest("Pillow is required for visual fingerprint tests")
        before = Image.new("RGB", (704, 896), "white")
        after = before.copy()
        ImageDraw.Draw(before).rectangle((30, 40, 39, 49), fill="#c81e1e")
        ImageDraw.Draw(after).rectangle((30, 40, 39, 49), fill="#1e1ec8")
        golden = {**check_visual.image_metrics("paper", before), "marks": []}
        actual = {**check_visual.image_metrics("paper", after), "marks": []}
        errors, _ = check_visual.compare(golden, actual)
        self.assertEqual(len(errors), 1)
        self.assertIn("header hash", errors[0])

    def test_coordinates_come_from_the_crossed_ring_and_connector(self):
        cdxml = '''<CDXML><page><fragment>
            <n id="1" p="0 0"/><n id="2" p="6 0"/><n id="3" p="3 6"/>
            <n id="4" NodeType="VariableAttachment" Attachments="1 2" p="3 1"/>
            <n id="5" p="3 -2"><t><s>R1</s></t></n>
            <b id="6" B="1" E="2" Order="2" CrossingBonds="9"/>
            <b id="7" B="2" E="3"/><b id="8" B="3" E="1"/>
            <b id="9" B="5" E="4" CrossingBonds="6"/>
        </fragment></page></CDXML>'''
        with TemporaryDirectory() as directory:
            path = Path(directory) / "molecule.cdxml"
            path.write_text(cdxml, encoding="utf-8")
            mark, = check_visual.measure_marks(path)
        self.assertEqual(mark["ring_center_pt"], [3.0, 2.0])
        self.assertEqual(mark["anchor_pt"], [3.0, 1.0])
        self.assertEqual(mark["label_endpoint_pt"], [3.0, -2.0])
        golden = {"size_px": [100, 200], "perceptual_hash": "0", "marks": [mark]}
        moved = {**golden, "marks": [{**mark, "anchor_pt": [3.0, 1.2]}]}
        errors, _ = check_visual.compare(golden, moved)
        self.assertIn("anchor_pt", " ".join(errors))
