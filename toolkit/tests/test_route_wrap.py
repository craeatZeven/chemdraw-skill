"""A long route wraps by continuing its arrow on the next row, each structure drawn once (P18).

Witness: journal routes run left to right and wrap; the intermediate at the end of a row serves
the arrow that starts the next one, instead of being drawn again.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

from cdxml_toolkit.render import renderer as R
from cdxml_toolkit.render.schema import ArrowContent, SchemeDescriptor, StepDescriptor, StructureRef

SMILES = ["CCO", "CCN", "CCCl", "CCBr", "CCI", "CCS"]


def _route(wrap: str, width: float = 380.0):
    ids = [f"S{k}" for k in range(6)]
    structures = {i: StructureRef(id=i, smiles=s, label=str(k + 1)) for k, (i, s) in enumerate(zip(ids, SMILES))}
    steps = [StepDescriptor(substrates=[ids[k]], products=[ids[k + 1]],
                            above_arrow=ArrowContent(text=["reagent", "solvent"])) for k in range(5)]
    scheme = SchemeDescriptor(structures=structures, steps=steps, layout="sequential", wrap=wrap,
                              max_row_width=width)
    xml, _ = R._layout_sequential(scheme, R._IDGen())
    return ET.fromstring(f"<page>{xml}</page>")


def _rows(page):
    return sorted({round(float(a.get("Head3D").split()[1]), 1) for a in page.findall("arrow")})


def test_a_wrapped_route_draws_every_structure_once():
    page = _route("continue")
    assert len(page.findall("fragment")) == 6
    assert len(_rows(page)) >= 2


def test_the_next_row_starts_with_the_arrow_that_leaves_the_last_structure():
    page = _route("continue")
    rows = _rows(page)
    second = [a for a in page.findall("arrow") if round(float(a.get("Head3D").split()[1]), 1) == rows[1]]
    first_arrow = min(second, key=lambda a: float(a.get("Tail3D").split()[0]))
    tail_x = float(first_arrow.get("Tail3D").split()[0])
    row_frags = [f for f in page.findall("fragment")
                 if abs(sum(float(n.get("p").split()[1]) for n in f.iter("n")) / len(f.findall(".//n")) - rows[1]) < 30]
    assert row_frags and tail_x < min(float(n.get("p").split()[0]) for f in row_frags for n in f.iter("n"))


def test_continuing_fits_more_steps_on_a_row_than_repeating():
    assert len(_rows(_route("continue"))) <= len(_rows(_route("repeat")))
    assert len(_route("repeat").findall("fragment")) > 6   # the old way draws the shared ones again


def _brief(**extra):
    species = [{"id": f"S{k}", "smiles": s} for k, s in enumerate(SMILES)]
    steps = [{"substrates": [f"S{k}"], "products": [f"S{k + 1}"]} for k in range(5)]
    return {"species": species, "steps": steps, **extra}


def test_a_brief_route_continues_by_default_and_can_choose_another_wrap():
    from cdxml_toolkit.mcp_runtime import figure_draw
    grounded = [{"id": f"S{k}", "smiles": s, "label": None} for k, s in enumerate(SMILES)]
    assert figure_draw._scheme(_brief(), grounded, 504.0, 1.0).wrap == "continue"
    assert figure_draw._scheme(_brief(wrap="serpentine"), grounded, 504.0, 1.0).wrap == "serpentine"
    assert figure_draw._validate(_brief(wrap="repeat")) == []
    assert any("wrap" in e for e in figure_draw._validate(_brief(wrap="zigzag")))
