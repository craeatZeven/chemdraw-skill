"""Crops must reach the engine even when the temp directory has a non-ASCII name.

OpenCV's imwrite takes a narrow path string and returns False (no exception) for a
path such as ``C:\\Users\\<non-ASCII>\\AppData\\Local\\Temp``; the crop never appears
and the DECIMER engine is then handed paths that do not exist.
"""
import os

import numpy as np

from cdxml_toolkit.image.structure_from_image import _write_crop


def test_a_crop_is_written_under_a_non_ascii_directory(tmp_path):
    target = tmp_path / "\u4f7f\u7528\u8005" / "crop_000.png"
    target.parent.mkdir()
    image = np.full((12, 16, 3), 255, dtype=np.uint8)
    image[2:8, 3:9] = 0

    _write_crop(str(target), image)

    assert target.is_file() and target.stat().st_size > 0


def test_a_crop_written_to_an_ascii_path_is_still_a_png(tmp_path):
    target = tmp_path / "crop_ascii.png"

    _write_crop(str(target), np.zeros((8, 8, 3), dtype=np.uint8))

    assert target.is_file()
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
