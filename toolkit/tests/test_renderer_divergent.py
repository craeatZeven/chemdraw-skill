from __future__ import annotations

import xml.etree.ElementTree as ET

from cdxml_toolkit.mcp_runtime import figure_lint
from cdxml_toolkit.render import renderer
from cdxml_toolkit.render.parser import parse_yaml

DIVERGENT = """
layout: divergent
structures:
  S: {smiles: "O=Cc1ccccc1", label: "1"}
  P: {smiles: "OCc1ccccc1", label: "2"}
  Q: {smiles: "OC(=O)c1ccccc1", label: "3"}
  R: {smiles: "C=Cc1ccccc1", label: "4"}
steps:
  - {substrates: [S], products: [P], above_arrow: {text: ["NaBH4"]}}
  - {substrates: [S], products: [Q], above_arrow: {text: ["KMnO4"]}, below_arrow: {text: ["H2O, 80 C"]}}
  - {substrates: [S], products: [R], above_arrow: {text: ["Ph3P=CH2"]}, arrow_style: failed}
"""


def _render(tmp_path):
    source = tmp_path / "scheme.yaml"
    source.write_text(DIVERGENT, encoding="utf-8")
    cdxml = renderer.render(parse_yaml(str(source)))
    return cdxml, ET.fromstring(cdxml[cdxml.index("<CDXML"):])


def _boxes(root):
    return {f.get("id"): [float(v) for v in f.get("BoundingBox").split()] for f in root.iter("fragment")}


def _point(arrow, key):
    x, y, _ = (float(v) for v in arrow.get(key).split())
    return x, y


def test_every_branch_arrow_leaves_from_the_starting_material(tmp_path):
    _, root = _render(tmp_path)
    boxes = _boxes(root)
    arrows = {a.get("id"): a for a in root.iter("arrow")}
    steps = list(root.iter("step"))
    assert len(steps) == 3
    start = steps[0].get("ReactionStepReactants")
    assert all(step.get("ReactionStepReactants") == start for step in steps)
    tails = set()
    for step in steps:
        arrow = arrows[step.get("ReactionStepArrows")]
        (tail_x, tail_y), (head_x, head_y) = _point(arrow, "Tail3D"), _point(arrow, "Head3D")
        product = boxes[step.get("ReactionStepProducts")]
        # A horizontal arrow on the product's own row, pointing at it from the left.
        assert tail_y == head_y and product[1] <= head_y <= product[3]
        assert head_x < product[0]
        # No arrow starts under or beside another product: all tails share one column
        # right of the starting material.
        assert tail_x > boxes[start][2]
        tails.add(round(tail_x, 2))
    assert len(tails) == 1


def test_a_connector_joins_the_branch_tails(tmp_path):
    _, root = _render(tmp_path)
    tails = sorted(_point(a, "Tail3D")[1] for a in root.iter("arrow"))
    lines = [g for g in root.iter("graphic") if g.get("GraphicType") == "Line"]
    assert len(lines) == 1
    x1, y1, x2, y2 = (float(v) for v in lines[0].get("BoundingBox").split())
    assert x1 == x2 and (min(y1, y2), max(y1, y2)) == (tails[0], tails[-1])


def test_branch_conditions_are_tidied_like_other_layouts(tmp_path):
    _, root = _render(tmp_path)
    texts = ["".join(t.itertext()).strip() for t in root.iter("t")]
    assert "H2O, 80 \u00b0C" in texts


def test_the_branches_pass_the_layout_checks(tmp_path):
    cdxml, _ = _render(tmp_path)
    assert figure_lint.lint_cdxml(cdxml) == []


THREE_WAY = """
layout: divergent
structures:
  E: {smiles: "C=CC(=O)OC", label: "1"}
  N1: {smiles: "CC(=O)CC(=O)C", label: "2a"}
  N2: {smiles: "C(C#N)C#N", label: "2b"}
  N3: {smiles: "O=[N+]([O-])CC", label: "2c"}
  P1: {smiles: "CC(=O)C(CCC(=O)OC)C(=O)C", label: "3a"}
  P2: {smiles: "N#CC(CCC(=O)OC)C#N", label: "3b"}
  P3: {smiles: "CCC(CCC(=O)OC)[N+](=O)[O-]", label: "3c"}
  B: {smiles: "C1CCN2CCCN=C2CC1", label: "DBU"}
  X: {smiles: "OC(=O)C(O)c1ccccc1", label: "4"}
  Y: {smiles: "O=C(O)C(=O)c1ccccc1", label: "5"}
steps:
  - {substrates: [E, N1], products: [P1], below_arrow: {text: ["THF, rt"]}, yield: "91%"}
  - {substrates: [E, N2], products: [P2], above_arrow: {structures: [B], text: ["(10 mol%)"]}}
  - {substrates: [E, N3], products: [P3], above_arrow: {text: ["K2CO3"]}}
  - {substrates: [X], products: [Y], above_arrow: {text: ["PCC"]}}
"""


def _render_three_way(tmp_path):
    source = tmp_path / "three_way.yaml"
    source.write_text(THREE_WAY, encoding="utf-8")
    cdxml = renderer.render(parse_yaml(str(source)))
    return cdxml, ET.fromstring(cdxml[cdxml.index("<CDXML"):])


def test_each_branch_draws_its_co_substrate_as_a_reactant(tmp_path):
    _, root = _render_three_way(tmp_path)
    boxes = _boxes(root)
    arrows = {a.get("id"): a for a in root.iter("arrow")}
    steps = list(root.iter("step"))
    branches = [s for s in steps if len(s.get("ReactionStepReactants", "").split()) == 2]
    assert len(branches) == 3
    shared = {s.get("ReactionStepReactants").split()[0] for s in branches}
    assert len(shared) == 1
    for step in branches:
        co = boxes[step.get("ReactionStepReactants").split()[1]]
        arrow = arrows[step.get("ReactionStepArrows")]
        (tail_x, tail_y), (head_x, _) = _point(arrow, "Tail3D"), _point(arrow, "Head3D")
        # The co-substrate stands above its own arrow, within the arrow's length.
        assert co[3] < tail_y and tail_x <= co[0] and co[2] <= head_x


def test_reagent_structures_above_the_arrow_are_drawn(tmp_path):
    _, root = _render_three_way(tmp_path)
    above = [s.get("ReactionStepObjectsAboveArrow", "").split() for s in root.iter("step")]
    fragment_ids = {f.get("id") for f in root.iter("fragment")}
    assert any(any(i in fragment_ids for i in ids) for ids in above)


def test_a_step_without_the_shared_substrate_is_still_drawn(tmp_path):
    _, root = _render_three_way(tmp_path)
    labels = {"".join(t.itertext()).strip() for t in root.iter("t")}
    assert {"4", "5", "PCC"} <= labels
    assert len(list(root.iter("step"))) == 4


def test_the_three_way_scheme_passes_the_layout_checks(tmp_path):
    cdxml, _ = _render_three_way(tmp_path)
    assert figure_lint.lint_cdxml(cdxml) == []
