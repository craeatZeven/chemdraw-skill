"""A double bond whose geometry the source leaves open reads back with the one it was drawn in.

RDKit reads E/Z off the 2D coordinates of a drawn C=N or C=C, so acetophenone oxime given as
CC(=NO)c1ccccc1 comes back as C/C(=N\\O)c1ccccc1. That is the compound asked for, drawn one way;
the guards accepted it only when the source had fixed the geometry, so every oxime, imine and
hydrazone without one was refused (found on the P19 and P20 parity renders).
"""

from __future__ import annotations

from collections import Counter

import pytest
from rdkit import Chem

from cdxml_toolkit import chemistry_semantics as S
from cdxml_toolkit.cdxml_builder import build_molecule_cdxml


@pytest.mark.parametrize("smiles", ["CC(=NO)c1ccccc1", "CC=NNc1ccccc1", "CC=CC(=O)O"])
def test_an_open_double_bond_is_drawn_and_read_back_as_the_same_compound(smiles):
    source = Chem.MolFromSmiles(smiles)
    cdxml = build_molecule_cdxml(source)             # validates the fragment it writes
    key = S.semantic_key(source)
    assert S.document_inventory(cdxml, like=[key]) == Counter({key: 1})


def test_a_fixed_geometry_must_still_come_back_the_same():
    source = Chem.MolFromSmiles("C/C=C/C(=O)O")      # E fixed in the source
    other = Chem.MolFromSmiles("C/C=C\\C(=O)O")       # Z drawn
    assert S.same_compound(source, other) is None
    assert S.same_compound(source, source) == "exact"
    assert S.same_compound(Chem.MolFromSmiles("CC=CC(=O)O"), other) == "drawn_geometry"


def test_a_different_compound_is_not_excused():
    assert S.same_compound(Chem.MolFromSmiles("CC(=NO)c1ccccc1"), Chem.MolFromSmiles("CC(=NO)c1ccccn1")) is None


def test_the_fan_with_an_oxime_draws(tmp_path):
    from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure
    brief = {"species": [{"id": "SM", "smiles": "CC(=O)c1ccccc1"}, {"id": "NH2OH", "smiles": "NO"},
                         {"id": "B", "smiles": "CC(=NO)c1ccccc1"}, {"id": "A", "smiles": "CC(O)c1ccccc1"}],
             "steps": [{"substrates": ["SM"], "products": ["A"], "above": ["NaBH4"]},
                       {"substrates": ["SM", "NH2OH"], "products": ["B"], "above": ["NaOAc"]}]}
    report = draw_figure(brief, output_path=str(tmp_path / "oxime.cdxml"))
    assert report["ok"], report
    assert report["readback"]["status"] == "matched"
