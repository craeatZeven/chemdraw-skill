"""Comprehensive tests for conformer/projection templates, acs_builder power-user tricks,
and draw_figure ring_fills / atom_halos / ts_bracket decorations.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.mcp_runtime import figure_draw
from cdxml_toolkit.templates import (
    bicyclic_cage,
    build_conformer_cdxml,
    chair_cyclohexane,
    fischer_projection,
    haworth_pyranose,
    newman_projection,
)

SKILL_SCRIPTS = Path(__file__).resolve().parents[2] / "skill" / "chemdraw" / "scripts"
if str(SKILL_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SKILL_SCRIPTS))

import acs_builder  # noqa: E402


# ---------------------------------------------------------------------------
# 1. Conformer & Perspective Templates (conformers.py)
# ---------------------------------------------------------------------------

def test_chair_cyclohexane_chair1_and_chair2_geometry():
    c1 = chair_cyclohexane(
        flip=False,
        substituents={1: {"ax": "OH", "eq": "Me"}, 4: {"eq": "tBu"}},
        cx=180.0,
        cy=160.0,
        bond_length=14.4,
    )
    c2 = chair_cyclohexane(
        flip=True,
        substituents={1: {"ax": "OH", "eq": "Me"}, 4: {"ax": "tBu"}},
        cx=180.0,
        cy=160.0,
        bond_length=14.4,
    )

    assert c1["conformation"] == "chair1" and c2["conformation"] == "chair2"
    for chair in (c1, c2):
        pts = chair["ring_coords"]
        assert len(pts) == 6
        # Opposite ring bonds (0-1 vs 3-4, 1-2 vs 4-5, 2-3 vs 5-0) are strictly parallel
        for k in range(3):
            v1 = (pts[(k + 1) % 6][0] - pts[k][0], pts[(k + 1) % 6][1] - pts[k][1])
            v2 = (pts[(k + 4) % 6][0] - pts[k + 3][0], pts[(k + 4) % 6][1] - pts[k + 3][1])
            assert math.isclose(math.hypot(*v1), 14.4, abs_tol=1e-6)
            assert math.isclose(math.hypot(*v2), 14.4, abs_tol=1e-6)
            cross = v1[0] * v2[1] - v1[1] * v2[0]
            assert abs(cross) < 1e-6

        # Equatorial vectors are parallel to the (i +/- 2) ring bonds and point outward
        for pos in range(1, 7):
            eq_x, eq_y = chair["equatorial_vectors"][pos]
            assert math.isclose(math.hypot(eq_x, eq_y), 14.4, abs_tol=1e-6)
            ref_a = pts[pos % 6]
            ref_b = pts[(pos + 1) % 6]
            rv_x, rv_y = ref_b[0] - ref_a[0], ref_b[1] - ref_a[1]
            assert abs(eq_x * rv_y - eq_y * rv_x) < 1e-6
            outward = (pts[pos - 1][0] - 180.0) * eq_x + (pts[pos - 1][1] - 160.0) * eq_y
            assert outward > 0.0

    # Axial vectors are strictly vertical and invert on ring flip
    for pos in (1, 3, 5):
        assert c1["axial_vectors"][pos] == (0.0, 14.4)
        assert c2["axial_vectors"][pos] == (0.0, -14.4)
    for pos in (2, 4, 6):
        assert c1["axial_vectors"][pos] == (0.0, -14.4)
        assert c2["axial_vectors"][pos] == (0.0, 14.4)

    root = ET.fromstring(c1["cdxml"])
    labels = ["".join(n.itertext()) for n in root.findall(".//fragment/n") if list(n)]
    assert "OH" in labels and "Me" in labels and "tBu" in labels

    # Corner angle at C2 (between C2->C1 and C2->C3) is 40 deg (u=-60 deg, v=160 deg)
    pts2 = c2["ring_coords"]
    e21 = (pts2[0][0] - pts2[1][0], pts2[0][1] - pts2[1][1])
    e23 = (pts2[2][0] - pts2[1][0], pts2[2][1] - pts2[1][1])
    cos_corner = (e21[0] * e23[0] + e21[1] * e23[1]) / (math.hypot(*e21) * math.hypot(*e23))
    assert math.isclose(math.degrees(math.acos(cos_corner)), 40.0, abs_tol=1e-6)

    # Rear axial C6-H in chair2 passes behind front ring bond C3-C4 (front bond carries CrossingBonds)
    # and keeps >= 3.0 pt clearance from C3-C4 at default ACS bond length 14.4 pt.
    from cdxml_toolkit.label_geometry import label_box

    c2_h6 = chair_cyclohexane(
        flip=True,
        substituents={
            1: {"ax": {"label": "Br", "display": "Bold"}},
            4: {"eq": "Me"},
            6: {"ax": {"label": "H", "display": "WedgedHashBegin"}},
        },
        bond_length=14.4,
    )
    c2_root = ET.fromstring(c2_h6["cdxml"])
    c2_bonds = {b.get("id"): b for b in c2_root.findall(".//fragment/b")}
    c2_nodes = {n.get("id"): n for n in c2_root.findall(".//fragment/n")}
    h6_bond = next(b for b in c2_bonds.values() if b.get("Display") == "WedgedHashBegin")
    crossing_ring_bonds = [b for b in c2_bonds.values() if h6_bond.get("id") in (b.get("CrossingBonds") or "").split()]
    assert len(crossing_ring_bonds) == 1
    assert int(crossing_ring_bonds[0].get("Z")) > int(h6_bond.get("Z"))

    # Left-pointing equatorial substituent (C4 eq) uses Right justification
    me_node = next(n for n in c2_nodes.values() if "".join(n.itertext()) == "Me")
    assert me_node.find("t").get("LabelJustification") == "Right"

    h6_atom = next(a for a in c2_h6["atoms"] if a.get("symbol") == "H" and a.get("ring_position") == 6)
    hx1, hy1, hx2, hy2 = label_box("H", h6_atom["x"], h6_atom["y"], "Left")
    p3, p4 = c2_h6["ring_coords"][2], c2_h6["ring_coords"][3]
    min_d = min(
        math.hypot(
            max(hx1 - (p3[0] + t * (p4[0] - p3[0])), 0.0, (p3[0] + t * (p4[0] - p3[0])) - hx2),
            max(hy1 - (p3[1] + t * (p4[1] - p3[1])), 0.0, (p3[1] + t * (p4[1] - p3[1])) - hy2),
        )
        for t in (i / 40.0 for i in range(41))
    )
    assert min_d >= 3.0, f"rear axial H too close to front ring bond: {min_d:.2f} pt"


def test_newman_projection_staggered_and_eclipsed():
    stag = newman_projection(
        conformation="staggered",
        front=("Me", "H", "H"),
        rear=("Br", "H", "Me"),
        cx=180.0,
        cy=160.0,
        bond_length=24.0,
    )
    eclip = newman_projection(
        conformation="eclipsed",
        front=("Me", "H", "H"),
        rear=("Br", "H", "Me"),
        cx=180.0,
        cy=160.0,
        bond_length=24.0,
    )

    assert stag["front_angles"] == (-90.0, 30.0, 150.0)
    assert stag["rear_angles"] == (-30.0, 90.0, 210.0)
    assert eclip["rear_angles"] == (-76.0, 44.0, 164.0)

    expected_r = 0.48 * 24.0
    assert math.isclose(stag["circle_radius"], expected_r, abs_tol=1e-6)
    for sx, sy in stag["rear_spoke_starts"]:
        assert math.isclose(math.hypot(sx - 180.0, sy - 160.0), expected_r, abs_tol=1e-5)

    root = ET.fromstring(stag["cdxml"])
    oval = root.find('.//page/graphic[@GraphicType="Oval"]')
    assert oval is not None
    assert oval.get("OvalType") == "Circle Filled"
    assert oval.get("color") == "2"


def test_haworth_fischer_and_bicyclic_cages(tmp_path):
    haw = haworth_pyranose(anomer="beta", bond_length=18.0)
    h_root = ET.fromstring(haw["cdxml"])
    displays = [b.get("Display") for b in h_root.findall(".//fragment/b") if b.get("Display")]
    assert "Bold" in displays
    assert displays.count("WedgeBegin") == 2

    fisch = fischer_projection(
        centers=[{"left": "H", "right": "OH"}, {"left": "OH", "right": "H"}],
        top="CHO",
        bottom="CH2OH",
        bond_length=18.0,
    )
    f_root = ET.fromstring(fisch["cdxml"])
    f_texts = ["".join(n.itertext()) for n in f_root.findall(".//fragment/n") if list(n)]
    assert "CHO" in f_texts and "CH2OH" in f_texts and f_texts.count("OH") == 2

    from rdkit import Chem
    from cdxml_toolkit.chemistry_semantics import read_fragments
    from cdxml_toolkit.templates.conformers import _segment_crossing_params

    expected_cage_smiles = {
        "norbornane": "C1CC2CCC1C2",
        "bicyclo_222": "C1CC2CCC1CC2",
        "adamantane": "C1C2CC3CC1CC(C2)C3",
    }
    for scaf, expected_smi in expected_cage_smiles.items():
        bare = bicyclic_cage(scaffold=scaf, bond_length=16.0)
        mols = [m for _, m in read_fragments(bare["cdxml"])]
        assert len(mols) == 1
        assert Chem.MolToSmiles(mols[0]) == Chem.CanonSmiles(expected_smi), f"{scaf} graph mismatch"

        cage = bicyclic_cage(scaffold=scaf, substituents={1: "Br"}, bond_length=16.0)
        c_root = ET.fromstring(cage["cdxml"])
        node_xy = {
            n.get("id"): tuple(map(float, n.get("p").split()))
            for n in c_root.findall(".//fragment/n")
        }
        bond_map = {b.get("id"): b for b in c_root.findall(".//fragment/b")}
        crossing = [b for b in bond_map.values() if b.get("CrossingBonds")]
        assert len(crossing) >= 1, f"{scaf} missing CrossingBonds"
        for fb in crossing:
            a1, e1 = node_xy[fb.get("B")], node_xy[fb.get("E")]
            l1 = math.hypot(e1[0] - a1[0], e1[1] - a1[1])
            for rid in fb.get("CrossingBonds").split():
                rb = bond_map[rid]
                assert int(fb.get("Z")) > int(rb.get("Z")), f"{scaf} front bond Z must exceed rear bond Z"
                a2, e2 = node_xy[rb.get("B")], node_xy[rb.get("E")]
                l2 = math.hypot(e2[0] - a2[0], e2[1] - a2[1])
                params = _segment_crossing_params(a1, e1, a2, e2)
                assert params is not None, f"{scaf} CrossingBonds must be a real geometric intersection"
                t, u = params
                margin = min(t, 1.0 - t, u, 1.0 - u) * min(l1, l2)
                assert margin >= 0.25 * 16.0, f"{scaf} crossing too close to vertex ({margin:.2f} pt)"
                cos_a = abs(((e1[0] - a1[0]) * (e2[0] - a2[0]) + (e1[1] - a1[1]) * (e2[1] - a2[1])) / (l1 * l2))
                ang = math.degrees(math.acos(min(1.0, cos_a)))
                assert ang >= 40.0, f"{scaf} crossing angle too shallow ({ang:.1f} deg)"

    out_file = tmp_path / "chair_out.cdxml"
    cdxml_str = build_conformer_cdxml({"template": "chair", "flip": True}, output_path=out_file)
    assert out_file.is_file() and "<CDXML" in cdxml_str


# ---------------------------------------------------------------------------
# 2. acs_builder.py Power-User Tricks
# ---------------------------------------------------------------------------

def test_acs_builder_ring_fills_atom_halos_and_z_ordering(tmp_path):
    spec = {
        "objects": [
            {
                "type": "molecule",
                "smiles": "c1ccccc1O",
                "at": [150, 120],
                "ring_fills": [{"ring_index": 0, "color": "pastel_blue", "style": "polygon", "inset": 1.0}],
                "atom_halos": [{"atoms": [0, 6], "color": "pastel_yellow", "radius": 8.0}],
            }
        ]
    }
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps(spec), encoding="utf-8")
    out_path = acs_builder.build_from_spec(spec_path, tmp_path / "out.cdxml")
    root = ET.parse(out_path).getroot()
    page = root.find("page")
    assert page is not None

    children_tags = [el.tag for el in list(page)]
    frag_idx = children_tags.index("fragment")
    # Background curve and oval graphics must precede <fragment> in <page>
    assert "curve" in children_tags[:frag_idx]
    assert "graphic" in children_tags[:frag_idx]

    curve = page.find("curve")
    assert curve is not None
    assert curve.get("Closed") == "yes"
    assert curve.get("FillType") == "Solid"
    assert curve.get("LineType") == "None"
    assert curve.get("Z") == "1"

    halos = [g for g in page.findall("graphic") if g.get("OvalType") == "Circle Filled"]
    assert len(halos) == 2


def test_acs_builder_ferrocene_multi_attachment_ts_bracket_orbital_and_interaction():
    doc = acs_builder.Doc()
    # Ferrocene: Fe (atom 0) + two Cp rings (atoms 1..5 and 6..10)
    doc.add_smiles(
        "[Fe].c1cc[cH-]c1.c1cc[cH-]c1",
        at=(180, 140),
        multi_attachments=[
            {"metal_atom": 0, "ligand_atoms": [1, 2, 3, 4, 5], "bond_display": "WedgeBegin", "bond_order": "dative"},
            {"metal_atom": 0, "ligand_atoms": [6, 7, 8, 9, 10], "bond_display": "Dash", "bond_order": "dative"},
        ],
    )
    doc.add_ts_bracket((120, 90, 240, 190), superscript="\u2021")
    doc.add_orbital(180, 140, orbital_type="pShaded", angle_deg=-90.0, length=18.0, width=7.5, color="jacs_blue")
    doc.add_interaction((140, 140), (220, 140), style="dashed", color="jacs_red", label="2.05 A", clip_margin=3.5)
    doc.add_conformer({"template": "norbornane", "at": [320, 140]})

    cdxml = doc.to_cdxml()
    root = ET.fromstring(cdxml)

    ma_nodes = root.findall('.//fragment/n[@NodeType="MultiAttachment"]')
    assert len(ma_nodes) == 2
    assert all(len(n.get("Attachments", "").split()) == 5 for n in ma_nodes)

    dative_bonds = [b for b in root.findall(".//fragment/b") if b.get("Order") == "dative"]
    assert len(dative_bonds) == 2

    brackets = root.findall('.//page/graphic[@GraphicType="Bracket"]')
    assert len(brackets) == 2
    assert root.find(".//page/bracketedgroup") is not None
    all_text = " ".join("".join(t.itertext()) for t in root.findall(".//page/t"))
    assert "\u2021" in all_text
    assert "2.05 A" in all_text

    orbital = root.find('.//page/graphic[@GraphicType="Orbital"]')
    assert orbital is not None and orbital.get("OrbitalType") == "pShaded"

    dashed_lines = [
        g for g in root.findall('.//page/graphic[@GraphicType="Line"]') if g.get("LineType") == "Dashed"
    ]
    assert len(dashed_lines) == 1


# ---------------------------------------------------------------------------
# 3. draw_figure with ring_fills, atom_halos, and ts_bracket
# ---------------------------------------------------------------------------

def test_draw_figure_with_ring_fills_atom_halos_and_ts_bracket(tmp_path, monkeypatch):
    monkeypatch.setenv("CHEMDRAW_OFFLINE", "1")
    brief = {
        "species": [
            {
                "id": "SM",
                "smiles": "c1ccccc1Br",
                "label": "1",
                "ring_fills": [{"ring_index": 0, "color": "pastel_blue", "style": "polygon", "inset": 1.0}],
                "atom_halos": [{"atoms": [6], "color": "pastel_yellow", "radius": 7.5}],
            },
            {
                "id": "TS",
                "smiles": "c1ccccc1O",
                "label": "TS-1",
                "ts_bracket": True,
            },
        ],
        "steps": [
            {
                "substrates": ["SM"],
                "products": ["TS"],
                "above": ["Pd(PPh3)4, KOH"],
                "below": ["80 C, 4 h"],
            }
        ],
    }
    out_path = tmp_path / "decorated_scheme.cdxml"
    report = figure_draw.draw_figure(brief, output_path=str(out_path))
    assert report["ok"] is True
    assert report["lint"]["issues"] == []
    assert report["readback"]["status"] == "matched"
    assert report["checks_passed"] is True

    root = ET.parse(out_path).getroot()
    page = root.find("page")
    assert page is not None
    tags = [el.tag for el in list(page)]
    first_frag = tags.index("fragment")
    # Ring fill <curve> and atom halo <graphic> are placed before <fragment>
    assert "curve" in tags[:first_frag]
    assert "graphic" in tags[:first_frag]

    # Transition-state brackets and double dagger (\u2021) are present
    brackets = page.findall('graphic[@GraphicType="Bracket"]')
    assert len(brackets) == 2
    assert page.find("bracketedgroup") is not None
    page_texts = ["".join(t.itertext()) for t in page.findall("t")]
    assert "\u2021" in page_texts


def test_newman_projection_is_one_connected_substituted_ethane():
    from rdkit import Chem
    from cdxml_toolkit.chemistry_semantics import read_fragments

    res = newman_projection(front=("Me", "H", "H"), rear=("Br", "H", "Me"))
    mols = [m for _, m in read_fragments(res["cdxml"])]
    assert len(mols) == 1
    assert len(Chem.GetMolFrags(mols[0])) == 1
    assert Chem.MolToSmiles(mols[0]) == Chem.CanonSmiles("CCC(C)Br")
    root = ET.fromstring(res["cdxml"])
    assert root.find('.//page/graphic[@OvalType="Circle Filled"]') is not None
    assert len(root.findall(".//fragment/b")) == 7  # 6 spokes + the front-rear C-C bond


def test_conformer_labels_are_real_chemistry():
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from cdxml_toolkit.chemistry_semantics import read_fragments

    haw = [m for _, m in read_fragments(haworth_pyranose()["cdxml"])]
    assert len(haw) == 1 and rdMolDescriptors.CalcMolFormula(haw[0]) == "C6H12O6"
    fis = [m for _, m in read_fragments(fischer_projection(centers=[{"left": "H", "right": "OH"}])["cdxml"])]
    assert len(fis) == 1 and len(Chem.GetMolFrags(fis[0])) == 1
    assert Chem.MolToSmiles(fis[0], isomericSmiles=False) == Chem.CanonSmiles("O=CC(O)CO")


def test_haworth_lower_up_labels_sit_outside_the_ring():
    haw = haworth_pyranose()
    ring = [haw["ring_coords"][i] for i in (1, 2, 3, 4, 5, 6)]

    def inside(x, y):
        hit = False
        for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
            if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
                hit = not hit
        return hit

    ups = [a for a in haw["atoms"] if a.get("direction") == "up" and a["ring_position"] in (2, 3)]
    assert ups
    for a in ups:
        assert not inside(a["x"], a["y"])
        assert a["y"] < haw["ring_coords"][a["ring_position"]][1]


CONFORMER_SPECS = {
    "chair": {"substituents": {"1": {"ax": {"label": "Br", "display": "Bold"}, "eq": "H"}}},
    "newman": {"front": [{"label": "F", "display": "WedgeBegin"}, "H", "H"],
               "rear": [{"label": "F", "display": "WedgedHashBegin"}, "H", "H"]},
    "haworth": {"substituents": {"1": {"up": {"label": "OMe", "display": "Bold"}}}},
    "fischer": {"top": {"label": "CHO", "display": "Bold"},
                "centers": [{"left": {"label": "H", "display": "WedgedHashBegin"},
                             "right": {"label": "OH", "display": "WedgeBegin"}}],
                "bottom": "CH2OH"},
    "norbornane": {"substituents": {"7": {"label": "Br", "display": "WedgeBegin"}}},
    "bicyclo_222": {"substituents": {"1": {"label": "Br", "display": "WedgedHashBegin"}}},
    "adamantane": {"substituents": {"1": {"label": "Br", "display": "Bold"}}},
}


@pytest.mark.parametrize("template", list(CONFORMER_SPECS))
def test_every_conformer_template_draws_and_reads_back_in_a_scheme(tmp_path, template):
    spec = {"template": template, **CONFORMER_SPECS[template]}
    brief = {
        "species": [{"id": "a", "conformer": spec}, {"id": "b", "conformer": spec}],
        "steps": [{"substrates": ["a"], "products": ["b"], "arrow": "reversible"}],
    }
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / f"{template}.cdxml"))
    assert report["ok"] is True, report.get("errors")
    assert report["lint"]["issues"] == []
    assert report["readback"]["status"] == "matched", report["readback"]
    assert report["readback"]["fragment_roundtrip"] == "preserved"
    assert report["checks_passed"] is True, report["warnings"]
    root = ET.parse(report["output_path"]).getroot()
    assert float(root.get("BondLength")) == pytest.approx(14.4)
    assert len(root.findall(".//page/fragment")) == 2


@pytest.mark.parametrize("template, expected", [
    ("chair", {"Bold"}),
    ("newman", {"WedgeBegin", "WedgedHashBegin"}),
    ("haworth", {"Bold", "WedgeBegin"}),
    ("fischer", {"Bold", "WedgeBegin", "WedgedHashBegin"}),
    ("norbornane", {"WedgeBegin"}),
    ("bicyclo_222", {"WedgedHashBegin"}),
    ("adamantane", {"Bold"}),
])
def test_structured_substituents_write_bond_displays(template, expected):
    result = build_conformer_cdxml({"template": template, **CONFORMER_SPECS[template]})
    root = ET.fromstring(result)
    displays = {b.get("Display") for b in root.findall(".//fragment/b") if b.get("Display")}
    assert expected <= displays


def test_wedged_hash_alias_normalizes_to_canonical_cdxml_display():
    result = build_conformer_cdxml({
        "template": "chair",
        "flip": True,
        "substituents": {"2": {"ax": {"label": "H", "display": "WedgedHash"}}},
    })
    root = ET.fromstring(result)
    displays = {b.get("Display") for b in root.findall(".//fragment/b") if b.get("Display")}
    assert "WedgedHashBegin" in displays
    assert "WedgedHash" not in displays


@pytest.mark.parametrize("spec, message", [
    ({"template": "twist-boat"}, "template"),
    ({"template": "chair", "substituents": {"1": {"ax": {"label": "Br", "display": "Dash"}}}}, "display"),
    ({"template": "chair", "substituents": {"9": {"ax": "Br"}}}, "1..6"),
    ({"template": "newman", "front": [{"label": "F", "display": "Dash"}, "H", "H"]}, "display"),
])
def test_invalid_conformer_specs_are_refused_before_drawing(tmp_path, spec, message):
    brief = {"species": [{"id": "bad", "conformer": spec}],
             "steps": [{"substrates": ["bad"], "products": ["bad"]}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "bad.cdxml"))
    assert report["ok"] is False
    assert message.lower() in " ".join(report["errors"]).lower()
    assert not (tmp_path / "bad.cdxml").exists()


def test_conformer_is_exclusive_with_smiles_and_name():
    brief = {"species": [{"id": "bad", "smiles": "C", "conformer": {"template": "chair"}}],
             "steps": [{"substrates": ["bad"], "products": ["bad"]}]}
    errors = figure_draw._validate(brief)
    assert any("exactly one" in error for error in errors)


@pytest.mark.parametrize("template", ["norbornane", "bicyclo_222", "adamantane"])
def test_drawn_cage_crossing_bond_ids_are_remapped_to_real_rear_bonds(tmp_path, template):
    brief = {"species": [{"id": "cage", "conformer": {"template": template}},
                          {"id": "copy", "conformer": {"template": template}}],
             "steps": [{"substrates": ["cage"], "products": ["copy"], "arrow": "reversible"}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / f"{template}-ids.cdxml"))
    assert report["ok"] is True, report.get("errors")
    root = ET.parse(report["output_path"]).getroot()
    for fragment in root.findall(".//page/fragment"):
        bonds = {b.get("id"): b for b in fragment.findall("b")}
        crossing = [b for b in bonds.values() if b.get("CrossingBonds")]
        assert crossing
        for bond in crossing:
            for crossed_id in bond.get("CrossingBonds").split():
                assert crossed_id in bonds
                assert crossed_id != bond.get("id")
                assert bonds[crossed_id].get("B") in {n.get("id") for n in fragment.findall("n")}


def test_conformer_supports_ts_brackets_and_panels(tmp_path):
    chair = {"template": "chair", "substituents": {"1": {"ax": "Br"}}}
    brief = {"species": [{"id": "ts", "conformer": chair, "ts_bracket": True},
                          {"id": "copy", "conformer": chair}],
             "steps": [{"substrates": ["ts"], "products": ["copy"], "arrow": "reversible"}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "ts.cdxml"))
    assert report["ok"] and report["checks_passed"], report
    root = ET.parse(report["output_path"]).getroot()
    assert root.find(".//bracketedgroup") is not None
    bracketed = set(root.find(".//bracketedgroup").get("BracketedObjectIDs").split())
    chair_nodes = {n.get("id") for n in root.findall(".//page/fragment")[0].findall("n")}
    assert bracketed == chair_nodes

    cage_brief = {"species": [{"id": "cage", "conformer": {"template": "norbornane"}},
                               {"id": "copy", "conformer": {"template": "norbornane"}}],
                  "steps": [{"substrates": ["cage"], "products": ["copy"], "arrow": "reversible"}]}
    panels = {"kind": "panels", "panels": [{"title": "Conformer", "brief": cage_brief}]}
    panel_report = figure_draw.draw_figure(panels, output_path=str(tmp_path / "panel.cdxml"))
    assert panel_report["ok"] and panel_report["panels"][0]["checks_passed"], panel_report
    panel_root = ET.parse(panel_report["output_path"]).getroot()
    bond_ids = {b.get("id") for b in panel_root.iter("b")}
    assert all(ref in bond_ids for b in panel_root.iter("b") if b.get("CrossingBonds")
               for ref in b.get("CrossingBonds").split())


def test_unified_conformer_newman_works_in_mechanism(tmp_path):
    newman = {"template": "newman", "front": ["Me", "H", "H"], "rear": ["Br", "H", "Me"]}
    brief = {"kind": "mechanism",
             "species": [{"id": "anti", "conformer": newman}, {"id": "gauche", "conformer": newman}],
             "steps": [{"substrates": ["anti"], "products": ["gauche"], "arrow": "reversible"}]}
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "mechanism.cdxml"))
    assert report["ok"] and report["checks_passed"], report
    root = ET.parse(report["output_path"]).getroot()
    assert len([g for g in root.iter("graphic") if g.get("OvalType") == "Circle Filled"]) == 2


def test_figure_2_10_brief_draws_axial_bromide_beta_hydrogens_and_products(tmp_path):
    source = Path(__file__).resolve().parents[2] / "assets" / "parity" / "conformer-chair-fig2-10.json"
    brief = json.loads(source.read_text(encoding="utf-8"))
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "fig2-10.cdxml"))
    assert report["ok"] and report["checks_passed"], report
    root = ET.parse(report["output_path"]).getroot()
    displays = {b.get("Display") for b in root.iter("b") if b.get("Display")}
    assert {"WedgeBegin", "WedgedHashBegin"} <= displays
    assert "WedgedHash" not in displays
    bond_ids = {b.get("id") for b in root.iter("b")}
    crossing_bonds = [b for b in root.iter("b") if b.get("CrossingBonds")]
    assert len(crossing_bonds) == 2  # both cis and trans chairs cross C3-C4 over rear C6-H
    for b in crossing_bonds:
        assert all(ref in bond_ids for ref in b.get("CrossingBonds").split())
    for fragment in root.findall(".//page/fragment"):
        nodes = {n.get("id") for n in fragment.findall("n")}
        assert all(b.get("B") in nodes and b.get("E") in nodes for b in fragment.findall("b"))
    text = " ".join("".join(t.itertext()) for t in root.iter("t"))
    assert "cis (Br a, up)" in text
    assert "trans (Br a, up; Me a, down)" in text
    assert "1-methyl isomer" in text and "3-methyl isomer" in text


def test_figure_2_10_saved_products_have_expected_connectivity():
    from rdkit import Chem

    source = Path(__file__).resolve().parents[2] / "assets" / "parity" / "conformer-chair-fig2-10.cdxml"
    mols = {m.GetIntProp("CDX_FRAG_ID"): m for m in Chem.MolsFromCDXMLFile(str(source))}
    expected = {
        32: ("CC1=CCCCC1", "CTMHWPIWNRWQEG-UHFFFAOYSA-N"),
        50: ("CC1C=CCCC1", "UZPWKTCMUADILM-UHFFFAOYSA-N"),
        97: ("CC1C=CCCC1", "UZPWKTCMUADILM-UHFFFAOYSA-N"),
    }
    assert set(expected) <= set(mols)
    for fragment_id, (smiles, inchikey) in expected.items():
        assert Chem.MolToSmiles(mols[fragment_id]) == smiles
        assert Chem.MolToInchiKey(mols[fragment_id]) == inchikey


def test_figure_2_10_saved_substrates_preserve_geometry_and_read_back_cis_trans():
    from rdkit import Chem
    from cdxml_toolkit.chemistry_semantics import semantic_key

    source = Path(__file__).resolve().parents[2] / "assets" / "parity" / "conformer-chair-fig2-10.cdxml"
    root = ET.parse(source).getroot()
    fragments = {int(f.get("id")): f for f in root.findall(".//page/fragment")}

    def vectors(fragment_id, label):
        fragment = fragments[fragment_id]
        nodes = {n.get("id"): n for n in fragment.findall("n")}
        labelled = [n for n in nodes.values() if "".join(n.itertext()) == label]
        out = []
        for node in labelled:
            bond = next(b for b in fragment.findall("b") if node.get("id") in (b.get("B"), b.get("E")))
            ring_id = bond.get("E") if bond.get("B") == node.get("id") else bond.get("B")
            rx, ry = map(float, nodes[ring_id].get("p").split())
            x, y = map(float, node.get("p").split())
            out.append((round(x - rx, 2), round(y - ry, 2)))
        return sorted(out)

    assert vectors(4, "Br") == [(0.0, -14.4)]
    assert vectors(4, "H") == [(0.0, 14.4), (0.0, 14.4)]
    assert vectors(72, "Br") == [(0.0, -14.4)]
    assert vectors(72, "Me") == [(0.0, 14.4)]
    assert vectors(72, "H") == [(0.0, 14.4)]

    mols = {m.GetIntProp("CDX_FRAG_ID"): m for m in Chem.MolsFromCDXMLFile(str(source))}

    def configurations(mol):
        centres = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
        br = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Br").GetNeighbors()[0].GetIdx()
        methyl = next(a for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetDegree() == 1)
        me = methyl.GetNeighbors()[0].GetIdx()
        return centres[br], centres[me]

    cis = configurations(mols[4])
    trans = configurations(mols[72])
    assert "?" not in cis + trans
    assert cis[0] != cis[1]
    assert trans[0] == trans[1]
    assert semantic_key(mols[4]) != semantic_key(mols[72])


@pytest.mark.parametrize("kind", ["scheme", "mechanism"])
def test_smiles_wins_over_name_in_scheme_and_mechanism(tmp_path, monkeypatch, kind):
    from cdxml_toolkit.naming import mol_builder

    def unexpected_name_resolution(*args, **kwargs):
        raise AssertionError("a valid SMILES must not be replaced by name resolution")

    monkeypatch.setattr(mol_builder, "resolve_compound", unexpected_name_resolution)
    brief = {
        "kind": kind,
        "species": [{"id": "a", "smiles": "C", "name": "unresolvable-name"},
                    {"id": "b", "smiles": "CC", "name": "also-unresolvable"}],
        "steps": [{"substrates": ["a"], "products": ["b"], "arrow": "reversible"}],
    }
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / f"{kind}.cdxml"))
    assert report["ok"] is True, report.get("errors")


def test_legacy_newman_rejects_a_conformer_template_key(tmp_path):
    brief = {
        "kind": "mechanism",
        "species": [{"id": "bad", "newman": {"template": "chair"}}],
        "steps": [{"substrates": ["bad"], "products": ["bad"]}],
    }
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "bad.cdxml"))
    assert report["ok"] is False
    assert any("newman" in error and "template" in error for error in report["errors"])
    assert not (tmp_path / "bad.cdxml").exists()


@pytest.mark.parametrize("substituent", [17, "17"])
def test_numeric_conformer_substituents_are_rejected(tmp_path, substituent):
    brief = {
        "species": [{"id": "bad", "conformer": {"template": "chair", "substituents": {"1": {"ax": substituent}}}}],
        "steps": [{"substrates": ["bad"], "products": ["bad"]}],
    }
    report = figure_draw.draw_figure(brief, output_path=str(tmp_path / "bad.cdxml"))
    assert report["ok"] is False
    assert any("substituent" in error and "label" in error for error in report["errors"])
    assert not (tmp_path / "bad.cdxml").exists()


def test_structure_ref_eighth_positional_argument_remains_overlays():
    from cdxml_toolkit.render.schema import StructureRef

    overlays = ({"atoms": [0], "shape": "circle"},)
    ref = StructureRef("a", "C", None, None, None, "1", None, overlays)
    assert ref.overlays is overlays
    assert ref.conformer is None
