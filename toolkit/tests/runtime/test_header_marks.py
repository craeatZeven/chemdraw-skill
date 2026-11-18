"""The R marks of a scope's general reaction take the colour that marks what varies below
(BJOC 22-127 scheme 2: R¹ and R² in red over a grid whose varying groups are red)."""
import xml.etree.ElementTree as ET

from cdxml_toolkit.mcp_runtime import figure_house

PAGE = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML><colortable><color r="1" g="1" b="1"/><color r="0" g="0" b="0"/></colortable><page id="1">
<fragment id="10">
 <n id="11" p="0 0"/>
 <n id="12" p="10 0" NodeType="GenericNickname" GenericNickname="R1"><t id="13" p="10 0">
  <s font="3" size="10" color="0" face="96">R</s><s font="3" size="10" color="0" face="64">1</s></t></n>
 <b id="14" B="11" E="12"/>
</fragment>
<group id="20"><fragment id="21">
 <n id="22" p="0 50"/>
 <n id="23" p="10 50" NodeType="GenericNickname" GenericNickname="R"><t id="24" p="10 50">
  <s font="3" size="10" color="0" face="96">R</s></t></n>
 <b id="25" B="22" E="23"/>
</fragment><t id="26" p="0 70"><s color="0">3a</s></t></group>
</page></CDXML>"""


def _rgb(root, index):
    entry = root.find("colortable").findall("color")[int(index) - 2]
    return tuple(round(float(entry.get(c)) * 255) for c in "rgb")


def test_header_r_marks_take_the_highlight_colour():
    text, painted = figure_house.paint_header_marks(PAGE, "#FF6600")
    assert painted == 1
    root = ET.fromstring(str(text).split("?>", 1)[-1])
    mark = root.find(".//n[@id='12']")
    assert _rgb(root, mark.get("color")) == (255, 102, 0)
    assert {_rgb(root, s.get("color")) for s in mark.iter("s")} == {(255, 102, 0)}   # R and its number


def test_marks_inside_a_compound_cell_are_left_alone():
    text, _ = figure_house.paint_header_marks(PAGE, "#FF6600")
    root = ET.fromstring(str(text).split("?>", 1)[-1])
    inside = root.find(".//n[@id='23']")
    assert inside.get("color") is None and {s.get("color") for s in inside.iter("s")} == {"0"}
    assert root.find(".//n[@id='11']").get("color") is None   # only the marks, not the rest


def test_a_page_without_marks_is_unchanged():
    text, painted = figure_house.paint_header_marks(PAGE.replace('NodeType="GenericNickname"', ''), "#FF6600")
    assert painted == 0
