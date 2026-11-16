"""Free text captions must stay text when ChemDraw opens the file on Windows.

ChemDraw reads a caption such as "H2O" or "MeOH" as a nickname fragment when the
document interprets text chemically, so a saved drawing gains a chemical entity
nobody drew.  The macOS automation path already guards against this
(``mac_automation.protect_captions``); the COM path did not.
"""
import xml.etree.ElementTree as ElementTree

from cdxml_toolkit.chemdraw.cdx_converter import mark_free_captions


def _document(caption: str) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" ?><CDXML BondLength="14.4">'
        '<page id="1">'
        '<fragment id="100"><n id="101" p="10 20"><t p="10 20"><s font="3" size="10">O</s></t></n>'
        '<n id="102" p="30 20"/><b id="103" B="101" E="102"/></fragment>'
        f'<t id="200" p="40 60"><s font="3" size="10">{caption}</s></t>'
        '</page></CDXML>'
    )


def _by_id(text: str) -> dict:
    root = ElementTree.fromstring(text[text.index("<CDXML"):])
    return {node.get("id"): node for node in root.iter()}


def test_a_free_caption_is_marked_as_plain_text():
    marked = mark_free_captions(_document("H2O"))
    assert _by_id(marked)["200"].get("InterpretChemically") == "no"


def test_an_atom_label_is_left_alone():
    marked = mark_free_captions(_document("OH"))
    label = _by_id(marked)["101"]
    assert label.get("InterpretChemically") is None
    assert label.find("t").get("InterpretChemically") is None


def test_a_caption_that_already_says_something_is_not_overwritten():
    source = _document("H2O").replace('<t id="200"', '<t id="200" InterpretChemically="yes"')
    marked = mark_free_captions(source)
    assert _by_id(marked)["200"].get("InterpretChemically") == "yes"
