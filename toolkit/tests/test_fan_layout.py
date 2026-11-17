"""One starting material with arrows leaving it in several directions (P19).

Witness: BJOC 22-103 figure 1, panel (c): the isocyanide stands in the middle and one arrow
leaves it to the left and one to the right, each to its own product, a co-substrate above each.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.render import renderer as R
from cdxml_toolkit.render.schema import ArrowContent, SchemeDescriptor, StepDescriptor, StructureRef


def _fan(products: int, co_substrates: bool = False):
    smiles = {"SM": "CCN", "P1": "CCO", "P2": "CCCl", "P3": "CCBr", "P4": "CCI", "X1": "CS", "X2": "CF"}
    structures = {k: StructureRef(id=k, smiles=v, label=k) for k, v in smiles.items()}
    steps = []
    for k in range(products):
        subs = ["SM"] + ([f"X{k + 1}"] if co_substrates and k < 2 else [])
        steps.append(StepDescriptor(substrates=subs, products=[f"P{k + 1}"],
                                    above_arrow=ArrowContent(text=["base"]), below_arrow=ArrowContent(text=["catalyst"])))
    xml, _ = R._layout_fan(SchemeDescriptor(structures=structures, steps=steps, layout="fan"), R._IDGen())
    return ET.fromstring(f"<page>{xml}</page>")


def _box(frag):
    pts = [tuple(float(v) for v in n.get("p").split()) for n in frag.iter("n")]
    return min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts)


def _by_label(page):
    out = {}
    frags = page.findall("fragment")
    labels = [t for t in page.findall("t")]
    for f in frags:
        b = _box(f)
        below = [t for t in labels if b[0] - 5 <= float(t.get("p").split()[0]) <= b[2] + 5
                 and 0 < float(t.get("p").split()[1]) - b[3] < 30]
        name = "".join(below[0].itertext()).strip() if below else None
        out[name] = b
    return out


def test_two_products_leave_the_starting_material_to_the_left_and_to_the_right():
    page = _fan(2)
    boxes = _by_label(page)
    sm, p1, p2 = boxes["SM"], boxes["P1"], boxes["P2"]
    assert p1[2] < sm[0] and p2[0] > sm[2]                      # one each side
    arrows = page.findall("arrow")
    assert len(arrows) == 2
    heads = sorted(float(a.get("Head3D").split()[0]) - float(a.get("Tail3D").split()[0]) for a in arrows)
    assert heads[0] < 0 < heads[1]                               # one points left, one right
    for a in arrows:
        tail_x = float(a.get("Tail3D").split()[0])
        assert not (sm[0] <= tail_x <= sm[2])                    # it leaves from the edge


def test_a_third_product_goes_below_and_a_fourth_above():
    page = _fan(4)
    boxes = _by_label(page)
    sm = boxes["SM"]
    assert boxes["P3"][1] > sm[3] and boxes["P4"][3] < sm[1]
    vertical = [a for a in page.findall("arrow") if a.get("Head3D").split()[0] == a.get("Tail3D").split()[0]]
    assert len(vertical) == 2


def test_co_substrates_stand_over_their_arrows_and_count_as_reactants():
    page = _fan(2, co_substrates=True)
    steps = page.findall(".//scheme/step")
    assert len(steps) == 2 and all(len(s.get("ReactionStepReactants").split()) == 2 for s in steps)
    assert len(page.findall("fragment")) == 5


def test_a_brief_can_ask_for_a_fan():
    from cdxml_toolkit.mcp_runtime import figure_draw
    brief = {"species": [{"id": "SM", "smiles": "CCN"}, {"id": "A", "smiles": "CCO"}, {"id": "B", "smiles": "CCCl"}],
             "steps": [{"substrates": ["SM"], "products": ["A"]}, {"substrates": ["SM"], "products": ["B"]}],
             "layout": "fan"}
    assert figure_draw._validate(brief) == []
    assert figure_draw._route_layout(figure_draw._step_descriptors(brief), "fan")[0] == "fan"
    chain = {"species": brief["species"], "layout": "fan",
             "steps": [{"substrates": ["SM"], "products": ["A"]}, {"substrates": ["A"], "products": ["B"]}]}
    assert any("fan" in e for e in figure_draw._validate(chain))


def test_two_products_from_one_compound_default_to_a_fan_and_more_stay_stacked():
    from cdxml_toolkit.mcp_runtime import figure_draw
    two = {"species": [{"id": "SM", "smiles": "CCN"}, {"id": "A", "smiles": "CCO"}, {"id": "B", "smiles": "CCCl"}],
           "steps": [{"substrates": ["SM"], "products": ["A"]}, {"substrates": ["SM"], "products": ["B"]}]}
    steps = figure_draw._step_descriptors(two)
    assert figure_draw._route_layout(steps, None)[0] == "fan"
    assert figure_draw._route_layout(steps, "divergent")[0] == "divergent"
    three = dict(two, species=two["species"] + [{"id": "C", "smiles": "CCBr"}],
                 steps=two["steps"] + [{"substrates": ["SM"], "products": ["C"]}])
    assert figure_draw._route_layout(figure_draw._step_descriptors(three), None)[0] == "divergent"


def test_a_fan_too_wide_for_the_row_budget_is_drawn_stacked():
    page = _fan(2)
    xs = [float(n.get("p").split()[0]) for n in page.iter("n")]
    smiles = {"SM": "CCN", "P1": "CCO", "P2": "CCCl"}
    structures = {k: StructureRef(id=k, smiles=v, label=k) for k, v in smiles.items()}
    steps = [StepDescriptor(substrates=["SM"], products=[f"P{k}"], above_arrow=ArrowContent(text=["base"]))
             for k in (1, 2)]
    narrow = SchemeDescriptor(structures=structures, steps=steps, layout="fan",
                              max_row_width=(max(xs) - min(xs)) / 2)
    xml, _ = R._layout_fan(narrow, R._IDGen())
    arrows = ET.fromstring(f"<page>{xml}</page>").findall("arrow")
    assert all(float(a.get("Head3D").split()[0]) > float(a.get("Tail3D").split()[0]) for a in arrows)  # all point right


def test_the_scheme_records_when_a_fan_is_drawn_stacked():
    smiles = {"SM": "CCN", "P1": "CCO", "P2": "CCCl"}
    structures = {k: StructureRef(id=k, smiles=v, label=k) for k, v in smiles.items()}
    steps = [StepDescriptor(substrates=["SM"], products=[f"P{k}"]) for k in (1, 2)]
    scheme = SchemeDescriptor(structures=structures, steps=steps, layout="fan", max_row_width=10.0)
    R._layout_fan(scheme, R._IDGen())
    assert scheme.drawn_layout == "divergent"
