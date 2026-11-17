"""
renderer.py — Render a SchemeDescriptor to a CDXML document.

Supports:
  - linear: single step (substrates → products)
  - sequential: multi-step in a single row
  - wrap: repeat — multi-row L→R with repeated structures
  - wrap: serpentine — zigzag layout (L→R, R→L, L→R, ...) with vertical arrows

Uses RDKit for SMILES → 2D coords (no ChemDraw COM dependency).
Uses cdxml_builder infrastructure for fragment/arrow/text XML generation.
Uses reaction_cleanup-style layout logic for positioning.
"""

from __future__ import annotations

import json
import math
import os
from copy import deepcopy
from functools import partial
from dataclasses import dataclass, field, replace as dataclass_replace
from typing import Any, Dict, List, Optional, Tuple
from xml.sax.saxutils import escape as xml_escape

from ..constants import (
    ACS_BOND_LENGTH,
    ACS_BOND_LENGTH_STR,
    ACS_BOND_SPACING,
    ACS_BOLD_WIDTH,
    ACS_CAPTION_FACE,
    ACS_CAPTION_SIZE,
    ACS_CHAIN_ANGLE_STR,
    ACS_HASH_SPACING,
    ACS_LABEL_FACE,
    ACS_LABEL_FONT,
    ACS_LABEL_SIZE,
    ACS_LINE_WIDTH,
    ACS_MARGIN_WIDTH,
    CDXML_FOOTER,
    CDXML_HEADER,
    LAYOUT_ABOVE_GAP,
    LAYOUT_BELOW_GAP,
    LAYOUT_FRAG_GAP_BONDS,
    LAYOUT_INTER_GAP_BONDS,
    LAYOUT_PLUS_GAP,
    LAYOUT_ROW_GAP,
)
from ..cdxml_utils import LABEL_H_WIDTH, drawn_fragment_bbox, hydrogens_lead
from ..overlay_geometry import generic_ring_box, overlay_box
from ..text_formatting import build_formatted_s_xml, is_formula_label, normalize_condition_text

from .route_graph import convergent_sections
from .schema import (
    ArrowContent,
    RunArrowEntry,
    SchemeDescriptor,
    SectionDescriptor,
    StepDescriptor,
    StepRunArrows,
    StructureRef,
)


# ---------------------------------------------------------------------------
# Text metrics for Arial 10pt Bold in ChemDraw
# (measured via bbox investigation: 135 fragments, 237 texts, 83 arrows)
# ---------------------------------------------------------------------------

_CHAR_WIDTH = 4.7       # average character width (proportional Arial)
RUN_LABEL_GAP = 6.0     # white space between run-arrow labels on one line
_LINE_ADVANCE = 11.5    # line-to-line distance (ChemDraw ~1.15× multiplier)
# Shortest reaction arrow, in bond lengths; a house style's measured arrow length raises it.
MIN_ARROW_BONDS = 5.0
_CAP_HEIGHT = 9.1       # baseline to top of uppercase letters
_DESCENT = 2.0          # baseline to bottom of descenders

# Fragment bbox padding beyond atom center positions (half of measured excess)
_FRAG_PAD_W = 3.3       # ±3.3 pt per side (6.6 pt total width excess)
_FRAG_PAD_H = 1.45      # ±1.45 pt per side (2.9 pt total height excess)
_LABEL_CHAR_W = 6.0     # Arial 10 pt label character beyond the first
_LABEL_DIGIT_W = 4.0    # subscript count or charge sign


# ---------------------------------------------------------------------------
# Multi-row constants
# ---------------------------------------------------------------------------

# Vertical gap between bottom of row N (including run arrows) and top of row N+1.
# ~55 pts matches real Report-scheme-extr-2 inter-row spacing.
ROW_GAP = 55.0


# ---------------------------------------------------------------------------
# Source JSON loader (reaction_parser output)
# ---------------------------------------------------------------------------

def _load_source_json(path: str) -> Dict[str, Dict]:
    """
    Load reaction_parser JSON, build species lookup by ID and by role.

    Returns a dict mapping keys to species dicts. Keys include:
      - species ID ("sp_0", "sp_1", ...)
      - "SM" / "DP" shortcuts
      - lowercase species name
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    lookup: Dict[str, Dict] = {}
    for sp in data.get("species", []):
        sp_id = sp.get("id", "")
        if sp_id:
            lookup[sp_id] = sp
        if sp.get("is_sm"):
            lookup["SM"] = sp
        if sp.get("is_dp"):
            lookup["DP"] = sp
        # Register by lowercase name for flexible matching
        name = sp.get("name", "")
        if name:
            lookup[name.lower()] = sp
        # Register by CSV name too
        csv_name = sp.get("csv_name", "")
        if csv_name:
            lookup[csv_name.lower()] = sp
    return lookup


# ---------------------------------------------------------------------------
# ID generator (same pattern as cdxml_builder)
# ---------------------------------------------------------------------------

class _IDGen:
    """Simple incrementing integer ID generator."""
    def __init__(self, start: int = 1000):
        self._n = start

    def next(self) -> int:
        v = self._n
        self._n += 1
        return v


# ---------------------------------------------------------------------------
# Resolved structure: SMILES → atoms/bonds in CDXML points
# ---------------------------------------------------------------------------

@dataclass
class ResolvedFragment:
    """A structure that has been resolved to atom/bond data + XML."""
    ref: StructureRef
    atoms: List[Dict]
    bonds: List[Dict]
    xml: str = ""                # <fragment> XML string
    frag_id: int = 0             # XML element ID
    # Bounding box (CDXML points): min_x, min_y, max_x, max_y
    bbox: Tuple[float, float, float, float] = (0, 0, 0, 0)
    # Current center position
    cx: float = 0
    cy: float = 0


@dataclass
class ResolvedStep:
    """A step with all structures resolved and laid out."""
    descriptor: StepDescriptor
    substrates: List[ResolvedFragment]
    products: List[ResolvedFragment]
    above_structures: List[ResolvedFragment]
    above_text: List[str]
    below_text: List[str]
    below_structures: List[ResolvedFragment]
    # Arrow geometry (set during layout)
    arrow_tail_x: float = 0
    arrow_tail_y: float = 0
    arrow_head_x: float = 0
    arrow_head_y: float = 0


# ---------------------------------------------------------------------------
# SMILES → atom/bond dicts (using structure_from_image + normalize)
# ---------------------------------------------------------------------------

_DEPICTION_CACHE_MAX = 512
_depiction_cache: Dict[Tuple[Any, ...], Tuple[List[Dict], List[Dict]]] = {}


def _cache_get(key: Tuple[Any, ...]) -> Optional[Tuple[List[Dict], List[Dict]]]:
    cached = _depiction_cache.get(key)
    if cached is None:
        return None
    return deepcopy(cached[0]), deepcopy(cached[1])


def _cache_put(key: Tuple[Any, ...], atoms: List[Dict], bonds: List[Dict]) -> None:
    if len(_depiction_cache) >= _DEPICTION_CACHE_MAX:
        _depiction_cache.clear()
    _depiction_cache[key] = (deepcopy(atoms), deepcopy(bonds))


def _smiles_to_fragment_data(
    smiles: str,
    center_x: float = 200.0,
    center_y: float = 300.0,
    ring_hints: Tuple[Tuple[int, int, int], ...] = (),
) -> Optional[Tuple[List[Dict], List[Dict]]]:
    """
    Convert SMILES to atom/bond dicts in CDXML point coordinates.

    Returns (atoms, bonds) or None on failure.
    """
    cache_key = ("smiles", smiles, round(center_x, 4), round(center_y, 4), tuple(ring_hints))
    cached = _cache_get(cache_key)
    if cached is not None:
        return cached

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise RuntimeError("RDKit is required for SMILES→structure conversion")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    if not ring_hints or not _depict_with_rings(mol, ring_hints):
        AllChem.Compute2DCoords(mol)

    source_mol = Chem.Mol(mol)

    # Kekulize for explicit single/double bonds.
    # Use a copy because clearAromaticFlags=True corrupts the mol on failure.
    mol_kek = Chem.RWMol(mol)
    try:
        Chem.Kekulize(mol_kek, clearAromaticFlags=True)
        mol = mol_kek
    except Exception:
        # Kekulization failed -- use original mol with aromatic bonds.
        # _rdkit_mol_to_atom_bond_dicts handles AROMATIC -> order 2 fallback.
        pass

    # Extract atom/bond dicts
    atoms, bonds = _extract_mol_atom_bond_dicts(mol, source_mol=source_mol)

    # Normalize to CDXML coordinate space (scale + flip y + center)
    from ..image.structure_from_image import normalize_for_cdxml
    atoms, bonds = normalize_for_cdxml(atoms, bonds, center_x, center_y)
    _separate_components(atoms, bonds)

    _cache_put(cache_key, atoms, bonds)
    return atoms, bonds


def _depict_with_rings(mol: Any, ring_hints) -> bool:
    """2D coordinates for a structure drawn with generic rings, laid out as if each ring were
    there: a ring of ``members`` carbons is added through each bonded pair (i, j), the whole
    is depicted, and the added atoms are dropped. The arc then takes the ring's place instead
    of crossing a chain that RDKit drew zigzag across it. False when that cannot be done."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    try:
        work = Chem.RWMol(mol)
        for i, j, members in ring_hints:
            previous = i
            for _ in range(int(members) - 2):
                added = work.AddAtom(Chem.Atom(6))
                work.AddBond(previous, added, Chem.BondType.SINGLE)
                previous = added
            work.AddBond(previous, j, Chem.BondType.SINGLE)
        Chem.SanitizeMol(work)
        rdDepictor.Compute2DCoords(work)
        placed = work.GetConformer()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for k in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(k, placed.GetAtomPosition(k))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return True
    except Exception:
        return False


def _extract_mol_atom_bond_dicts(mol: Any, source_mol: Optional[Any] = None) -> Tuple[List[Dict], List[Dict]]:
    """Extract atom/bond dicts from *mol*, preserving dative bonds, radicals, and pre-Kekulize semantic_key."""
    from rdkit import Chem
    from ..chemistry_semantics import _elemental_metal, semantic_key
    from ..image.structure_from_image import _rdkit_mol_to_atom_bond_dicts

    source_cx = semantic_key(source_mol if source_mol is not None else mol)
    try:
        atoms, bonds = _rdkit_mol_to_atom_bond_dicts(mol)
    except ValueError as exc:
        if "unsupported_radical" in str(exc):
            from ..cdxml_builder import mol_to_atom_bond_dicts
            atoms, bonds = mol_to_atom_bond_dicts(mol, normalize=False)
        else:
            raise

    if atoms:
        atoms[0]["_source_cxsmiles"] = source_cx
    for i, atom in enumerate(mol.GetAtoms()):
        if i < len(atoms):
            n_rad = atom.GetNumRadicalElectrons()
            if n_rad > 0 and not _elemental_metal(atom):
                atoms[i]["num_radical_electrons"] = n_rad
                atoms[i]["radical"] = "Doublet" if n_rad == 1 else "Triplet"
                atoms[i]["num_hydrogens"] = atom.GetTotalNumHs(includeNeighbors=False)

    dative_right = {
        Chem.BondType.DATIVE,
        Chem.BondType.DATIVER,
        Chem.BondType.DATIVEONE,
    }
    for bi, bond in enumerate(mol.GetBonds()):
        if bi < len(bonds):
            bt = bond.GetBondType()
            if bt in dative_right:
                bonds[bi]["order"] = "dative"
            elif bt == Chem.BondType.DATIVEL:
                bonds[bi]["order"] = "dative"
                bonds[bi]["atom1"], bonds[bi]["atom2"] = bonds[bi]["atom2"], bonds[bi]["atom1"]
    return atoms, bonds


_ION_GAP = 3.0  # white space between the labels of two ions of a salt (pt)


def _separate_components(atoms: List[Dict], bonds: List[Dict]) -> None:
    """Move the parts of a salt apart where their labels would run together.

    RDKit places disconnected parts one bond apart whatever their labels, so
    Na+ and OH- printed as "OHNa+". Parts that already clear each other keep
    RDKit's placement; the drawing stays centred where it was.
    """
    owner = {a["index"]: a["index"] for a in atoms}

    def root(i):
        while owner[i] != i:
            i = owner[i]
        return i

    for b in bonds:
        owner[root(b["atom1"])] = root(b["atom2"])
    parts: Dict[int, List[Dict]] = {}
    for a in atoms:
        parts.setdefault(root(a["index"]), []).append(a)
    if len(parts) < 2:
        return
    before = _fragment_bbox(atoms)
    def extent(part):
        # A label is about 10 pt tall, taller than a bare carbon's box.
        x1, _, x2, _ = _fragment_bbox(part)
        reach = [5.0 if _label_overhang(a) or a.get("symbol", "C") != "C" else _FRAG_PAD_H for a in part]
        return (x1, min(a["y"] - r for a, r in zip(part, reach)),
                x2, max(a["y"] + r for a, r in zip(part, reach)))

    placed: List[List[Dict]] = []
    for part in sorted(parts.values(), key=lambda p: _fragment_bbox(p)[0]):
        for other in placed:
            box, other_box = extent(part), extent(other)
            level = box[1] < other_box[3] and other_box[1] < box[3]
            if level and box[0] < other_box[2] + _ION_GAP:
                _shift_atoms(part, other_box[2] + _ION_GAP - box[0], 0.0)
        placed.append(part)
    after = _fragment_bbox(atoms)
    _shift_atoms(atoms, (before[0] + before[2] - after[0] - after[2]) / 2.0, 0.0)


def _align_mol_to_reference(
    target_mol,
    ref_mol,
    center_x: float = 200.0,
    center_y: float = 300.0,
) -> Optional[Tuple[List[Dict], List[Dict]]]:
    """
    Align target_mol's 2D layout to ref_mol using MCS, then extract
    atom/bond dicts in CDXML coordinates.

    Uses RDKit's GenerateDepictionMatching2DStructure to orient the target
    so its shared scaffold matches the reference orientation.

    Returns (atoms, bonds) or None if alignment fails or MCS is too small.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdFMCS
    except ImportError:
        return None

    if target_mol is None or ref_mol is None:
        return None

    cache_key: Optional[Tuple[Any, ...]] = None
    try:
        cache_key = (
            "align",
            Chem.MolToCXSmiles(target_mol),
            Chem.MolToCXSmiles(ref_mol),
            round(center_x, 4),
            round(center_y, 4),
        )
        cached = _cache_get(cache_key)
        if cached is not None:
            return cached
    except Exception:
        cache_key = None

    # Find MCS between target and reference
    try:
        mcs_result = rdFMCS.FindMCS(
            [target_mol, ref_mol],
            timeout=5,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
        )
    except Exception:
        return None

    if mcs_result.numAtoms < 3:
        return None  # MCS too small for meaningful alignment

    # Build the MCS query mol
    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
    if mcs_mol is None:
        return None

    # Get substructure matches
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)
    target_match = target_mol.GetSubstructMatch(mcs_mol)

    if not ref_match or not target_match:
        return None

    # Build atom map: (ref_idx, target_idx) for MCS atoms
    atom_map = list(zip(ref_match, target_match))

    # Ensure reference has 2D coordinates
    if ref_mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(ref_mol)

    # Generate aligned 2D coordinates for target
    try:
        AllChem.Compute2DCoords(target_mol)
        AllChem.GenerateDepictionMatching2DStructure(
            target_mol, ref_mol, atomMap=atom_map
        )
    except Exception:
        # Alignment failed — fall back to standard coords
        AllChem.Compute2DCoords(target_mol)

    source_mol = Chem.Mol(target_mol)

    # Kekulize for explicit bonds
    mol_kek = Chem.RWMol(target_mol)
    try:
        Chem.Kekulize(mol_kek, clearAromaticFlags=True)
        target_mol = mol_kek
    except Exception:
        pass

    # Extract atom/bond dicts
    atoms, bonds = _extract_mol_atom_bond_dicts(target_mol, source_mol=source_mol)

    # Normalize to CDXML coordinate space
    from ..image.structure_from_image import normalize_for_cdxml
    atoms, bonds = normalize_for_cdxml(atoms, bonds, center_x, center_y)

    if cache_key is not None:
        _cache_put(cache_key, atoms, bonds)

    return atoms, bonds


# ---------------------------------------------------------------------------
# Fragment XML builder (adapted from cdxml_builder._build_fragment)
# ---------------------------------------------------------------------------

# Element number lookup
ELEMENT_NUMBERS = {
    "H":  1,  "He": 2,  "Li": 3,  "Be": 4,  "B":  5,  "C":  6,  "N":  7,
    "O":  8,  "F":  9,  "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14,
    "P":  15, "S":  16, "Cl": 17, "Ar": 18, "K":  19, "Ca": 20, "Ti": 22,
    "V":  23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Zr": 40, "Mo": 42, "Ru": 44, "Rh": 45, "Pd": 46,
    "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I":  53,
    "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "W":  74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83,
}

# Bond stereo
BOND_STEREO_ATTR = {
    1: "WedgeBegin",
    4: "Wavy",
    6: "WedgedHashBegin",
}


from contextvars import ContextVar

_fragment_receipts = ContextVar('fragment_receipts', default=None)


def _build_fragment(
    atoms: List[Dict],
    bonds: List[Dict],
    ids: _IDGen,
) -> Tuple[str, Dict[int, int], int]:
    """
    Build a <fragment> XML string from atom/bond dicts.

    Returns (xml_string, atom_id_map, fragment_xml_id).
    """
    atom_id_map: Dict[int, int] = {}
    frag_id = ids.next()

    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    bb_x1, bb_y1 = min(xs), min(ys)
    bb_x2, bb_y2 = max(xs), max(ys)

    lines: List[str] = []
    lines.append(
        f'<fragment id="{frag_id}" '
        f'BoundingBox="{bb_x1:.2f} {bb_y1:.2f} {bb_x2:.2f} {bb_y2:.2f}" '
        f'Z="{ids.next()}">'
    )

    for a in atoms:
        atom_xml_id = ids.next()
        atom_id_map[a["index"]] = atom_xml_id
        ax, ay = a["x"], a["y"]
        z = ids.next()

        sym = a.get("symbol", "C")
        elem_num = a.get('atomic_number', ELEMENT_NUMBERS.get(sym, 6))
        nh = a.get("num_hydrogens")
        charge = a.get("charge", 0)
        isotope = a.get("isotope")
        rad = a.get("radical")
        rad_attr_val: Optional[str] = None
        if isinstance(rad, str) and rad in ("Doublet", "Triplet", "Singlet", "Diradical"):
            rad_attr_val = "Triplet" if rad == "Diradical" else rad
        else:
            n_rad = int(rad) if isinstance(rad, (int, float)) else int(a.get("num_radical_electrons", 0))
            if n_rad == 1:
                rad_attr_val = "Doublet"
            elif n_rad >= 2:
                rad_attr_val = "Triplet"

        attrs = [
            f'id="{atom_xml_id}"',
            f'p="{ax:.2f} {ay:.2f}"',
            f'Z="{z}"',
        ]

        is_carbon = (sym == "C" and not charge and not isotope)
        is_variable = elem_num == 0
        if is_variable:
            # Element=0 is accepted by some parsers but cannot be opened by
            # ChemDraw. Use its native generic substituent node instead.
            if charge or isotope:
                raise ValueError('Charged or isotopically labelled variable atoms are unsupported')
            # A numbered R group ([*:2] in the SMILES) is R2, drawn R with a raised 2.
            number = a.get('generic_number')
            attrs.extend(['NodeType="GenericNickname"', f'GenericNickname="R{number or ""}"'])
        elif not is_carbon:
            attrs.append(f'Element="{elem_num}"')
            if nh is not None:
                attrs.append(f'NumHydrogens="{nh}"')
            if isotope:
                attrs.append(f'Isotope="{isotope}"')
            attrs.append('NeedsClean="yes"')
            if charge:
                attrs.append(f'Charge="{charge}"')
        elif rad_attr_val and nh is not None:
            attrs.append(f'NumHydrogens="{nh}"')

        if rad_attr_val:
            attrs.append(f'Radical="{rad_attr_val}"')

        if a.get('enhanced_stereo_type'):
            attrs.append(f'EnhancedStereoType="{a["enhanced_stereo_type"]}"')
            attrs.append(f'EnhancedStereoGroupNum="{a["enhanced_stereo_group"]}"')

        if is_carbon:
            lines.append(f'<n {" ".join(attrs)}/>')
        else:
            # Heteroatom needs a text label
            lines.append(f'<n {" ".join(attrs)}>')
            lx = ax - 3.25
            ly = ay + 3.52
            # The first letter is centred on the atom and the rest ("r" of Br) follows it.
            label_w = max(len(sym) * 5.5, 6.0)
            lbx1 = ax - 3.0
            lby1 = ay - 7.52
            lbx2 = lbx1 + label_w
            lby2 = ay
            tid = ids.next()
            lines.append(
                f'<t id="{tid}" p="{lx:.2f} {ly:.2f}" '
                f'BoundingBox="{lbx1:.2f} {lby1:.2f} {lbx2:.2f} {lby2:.2f}" '
                f'LabelJustification="Left">'
            )
            # Use isotope-specific symbol for display (e.g. D for deuterium)
            if is_variable:
                display_text = "R"
            elif sym == "H" and isotope == 2:
                display_text = "D"
            elif sym == "H" and isotope == 3:
                display_text = "T"
            else:
                display_text = sym
                if nh is not None and nh > 0:
                    display_text += "H" if nh == 1 else f"H{nh}"
            if isotope and not (sym == 'H' and isotope in (2, 3)):
                lines.append(
                    f'<s font="{ACS_LABEL_FONT}" size="{ACS_LABEL_SIZE}" '
                    f'color="0" face="64">{int(isotope)}</s>'
                )
            lines.append(
                f'<s font="{ACS_LABEL_FONT}" size="{ACS_LABEL_SIZE}" '
                f'color="0" face="{ACS_LABEL_FACE}">{xml_escape(display_text)}</s>'
            )
            if is_variable and a.get('generic_number'):
                lines.append(
                    f'<s font="{ACS_LABEL_FONT}" size="{ACS_LABEL_SIZE}" '
                    f'color="0" face="64">{int(a["generic_number"])}</s>'
                )
            if charge:
                charge_text = (str(abs(charge)) if abs(charge) > 1 else '') + ('+' if charge > 0 else '-')
                lines.append(
                    f'<s font="{ACS_LABEL_FONT}" size="{ACS_LABEL_SIZE}" '
                    f'color="0" face="64">{charge_text}</s>'
                )
            lines.append('</t>')
            lines.append('</n>')

    # Bonds
    for b in bonds:
        bid = ids.next()
        z = ids.next()
        a1 = atom_id_map.get(b["atom1"], 0)
        a2 = atom_id_map.get(b["atom2"], 0)
        order = b.get("order", 1)
        order_str_upper = str(order).upper()
        if order_str_upper.endswith("DATIVEL"):
            a1, a2 = a2, a1

        attrs = [
            f'id="{bid}"',
            f'Z="{z}"',
            f'B="{a1}"',
            f'E="{a2}"',
        ]
        if order == 2 or order == "2":
            attrs.append('Order="2"')
        elif order == 3 or order == "3":
            attrs.append('Order="3"')
        elif b.get("dative") or order in ("dative", 9) or "DATIVE" in order_str_upper:
            attrs.append('Order="dative"')

        cfg = b.get("cfg", 0)
        if cfg in BOND_STEREO_ATTR:
            attrs.append(f'Display="{BOND_STEREO_ATTR[cfg]}"')

        dp = b.get("double_pos")
        if dp:
            attrs.append(f'DoublePosition="{dp}"')

        lines.append(f'<b {" ".join(attrs)}/>')

    lines.append('</fragment>')
    fragment_xml = "\n".join(lines)
    source_cx = atoms[0].get('_source_cxsmiles') if atoms else None
    if source_cx:
        from ..chemistry_semantics import validate_fragment
        receipt = validate_fragment(source_cx, fragment_xml)
        collector = _fragment_receipts.get()
        if collector is not None:
            collector.append({'fragment_id': frag_id, **receipt})
    return fragment_xml, atom_id_map, frag_id


# ---------------------------------------------------------------------------
# Text builder
# ---------------------------------------------------------------------------

def _build_text_element(
    text_lines: List[str],
    x: float,
    y: float,
    ids: _IDGen,
    justification: str = "Center",
    use_formatting: bool = True,
    font_size: Optional[float] = None,
) -> Tuple[str, int]:
    """
    Build a standalone <t> element for condition text.

    Parameters
    ----------
    font_size : float, optional
        Override font size (default: ACS_CAPTION_SIZE, typically 10pt).

    Returns (xml_string, text_xml_id).
    """
    size = float(font_size if font_size is not None else ACS_CAPTION_SIZE)
    scale = size / float(ACS_CAPTION_SIZE)
    char_w = _CHAR_WIDTH * scale
    cap_h = _CAP_HEIGHT * scale
    descent = _DESCENT * scale
    line_adv = _LINE_ADVANCE * scale

    tid = ids.next()
    z = ids.next()

    max_chars = max((len(ln) for ln in text_lines), default=1)
    n = len(text_lines)
    w = max_chars * char_w

    if justification == "Left":
        bx1, bx2 = x, x + w
    elif justification == "Right":
        bx1, bx2 = x - w, x
    else:
        bx1, bx2 = x - w / 2.0, x + w / 2.0
    by1 = y - cap_h
    by2 = y + max(0, n - 1) * line_adv + descent

    parts = [
        f'<t id="{tid}" p="{x:.2f} {y:.2f}" '
        f'BoundingBox="{bx1:.2f} {by1:.2f} {bx2:.2f} {by2:.2f}" '
        f'Z="{z}" '
        f'CaptionJustification="{justification}" '
        f'Justification="{justification}" '
        f'LineHeight="auto">'
    ]

    if use_formatting:
        # Use chemistry-aware text formatting for each line.
        # ChemDraw requires \n INSIDE <s> text content for line breaks —
        # newlines between XML elements are treated as whitespace.
        formatted_runs = []
        for i, line in enumerate(text_lines):
            run = build_formatted_s_xml(
                line,
                font=ACS_LABEL_FONT,
                size=size,
                color="0",
            )
            if i < len(text_lines) - 1:
                # Inject \n before the closing </s> of this line's last run
                # so ChemDraw renders a line break
                last_close = run.rfind("</s>")
                if last_close >= 0:
                    run = run[:last_close] + "\n" + run[last_close:]
            formatted_runs.append(run)
        parts.append("".join(formatted_runs))
    else:
        text = "\n".join(xml_escape(ln) for ln in text_lines)
        parts.append(
            f'<s font="{ACS_LABEL_FONT}" size="{size}" '
            f'color="0" face="{ACS_CAPTION_FACE}">{text}</s>'
        )

    parts.append("</t>")
    return "\n".join(parts), tid


def _build_label_element(
    label: str,
    x: float,
    y: float,
    ids: _IDGen,
) -> Tuple[str, int]:
    """Build a compound label centered below a structure: a number ("1", "TS2") as written, a
    formula ("PhB(OH)2") with its subscripts."""
    return _build_text_element(
        [label], x, y, ids, justification="Center", use_formatting=is_formula_label(label),
    )


# ---------------------------------------------------------------------------
# Arrow builder
# ---------------------------------------------------------------------------

def _arrow_kind(style: Optional[str]) -> Optional[str]:
    """The ``kind`` for _build_arrow: only reversible and resonance change the shape."""
    return style if style in ("reversible", "resonance") else None

def _build_arrow(
    tail_x: float,
    tail_y: float,
    head_x: float,
    head_y: float,
    ids: _IDGen,
    dashed: bool = False,
    nogo: bool = False,
    kind: Optional[str] = None,
) -> Tuple[str, int]:
    """Build an <arrow> element. Returns (xml_string, arrow_xml_id).

    Parameters
    ----------
    nogo : bool
        If True, adds ``NoGo="Cross"`` — ChemDraw's native failed-arrow
        rendering (bold X through the arrow shaft).
    kind : "reversible" | "resonance" | None
        ``reversible`` is ChemDraw's equilibrium arrow (two half-headed
        shafts); ``resonance`` is one shaft with a full head at both ends.
    """
    aid = ids.next()
    z = ids.next()

    bx1 = min(tail_x, head_x)
    by1 = min(tail_y, head_y) - 4.0
    bx2 = max(tail_x, head_x)
    by2 = max(tail_y, head_y) + 4.0

    cx3 = (tail_x + head_x) / 2.0
    cy3 = tail_y + 100.0

    attrs = [
        f'id="{aid}"',
        f'BoundingBox="{bx1:.2f} {by1:.2f} {bx2:.2f} {by2:.2f}"',
        f'Z="{z}"',
        f'FillType="None"',
        f'ArrowheadHead="Full"',
        f'ArrowheadType="Solid"',
        f'HeadSize="1000"',
        f'ArrowheadCenterSize="875"',
        f'ArrowheadWidth="250"',
        f'Head3D="{head_x:.2f} {head_y:.2f} 0"',
        f'Tail3D="{tail_x:.2f} {tail_y:.2f} 0"',
        f'LineWidth="{ACS_LINE_WIDTH}"',
    ]

    if kind == "reversible":
        attrs = [a for a in attrs if not a.startswith('ArrowheadHead=')]
        attrs += ['ArrowheadHead="HalfLeft"', 'ArrowheadTail="HalfLeft"',
                  'ArrowShaftSpacing="300"']
    elif kind == "resonance":
        attrs.append('ArrowheadTail="Full"')
    if dashed:
        attrs.append('LineType="Dashed"')
    if nogo:
        attrs.append('NoGo="Cross"')

    xml = f'<arrow {" ".join(attrs)}/>'
    return xml, aid


def _build_failed_x(
    cx: float,
    cy: float,
    ids: _IDGen,
) -> str:
    """Build a bold 'X' text element centered on an arrow midpoint.

    .. deprecated::
        Use ``NoGo="Cross"`` on the ``<arrow>`` element instead.
        This function is kept for backwards compatibility.
    """
    tid = ids.next()
    z = ids.next()
    # Position slightly above arrow line so X sits centered on it
    py = cy + 3.5  # baseline offset — text anchor is at baseline
    return (
        f'<t id="{tid}" p="{cx:.2f} {py:.2f}" Z="{z}" '
        f'Justification="Center" InterpretChemically="no" '
        f'CaptionJustification="Center">\n'
        f'<s font="{ACS_LABEL_FONT}" size="12" color="0" face="1">X</s>\n'
        f'</t>'
    )


def _is_letter_condition(text_lines: List[str]) -> bool:
    """Check if below-arrow text is a letter condition like 'a' or 'b,c'."""
    if len(text_lines) != 1:
        return False
    import re
    return bool(re.match(r'^[a-z](,\s*[a-z])*$', text_lines[0].strip()))


def _build_letter_label(
    letter: str,
    cx: float,
    y: float,
    ids: _IDGen,
) -> Tuple[str, int]:
    """Build a small letter label centered above an arrow."""
    tid = ids.next()
    z = ids.next()
    # Position slightly above arrow line
    py = y - 6.0
    xml = (
        f'<t id="{tid}" p="{cx:.2f} {py:.2f}" Z="{z}" '
        f'Justification="Center" InterpretChemically="no" '
        f'CaptionJustification="Center">\n'
        f'<s font="{ACS_LABEL_FONT}" size="8" color="0" face="1">'
        f'{xml_escape(letter)}</s>\n'
        f'</t>'
    )
    return xml, tid


def _build_condition_key(
    condition_key: Dict[str, str],
    left_x: float,
    top_y: float,
    ids: _IDGen,
) -> Tuple[str, float]:
    """Build the condition key block below the scheme.

    Returns (xml_string, bottom_y).
    """
    xml_parts: List[str] = []
    y = top_y
    for letter in sorted(condition_key.keys()):
        text = condition_key[letter]
        tid = ids.next()
        z = ids.next()
        # Format: "(a) conditions text" — letter in italic, rest in regular
        # Use build_formatted_s_xml for chemistry-aware subscript formatting
        letter_run = (
            f'<s font="{ACS_LABEL_FONT}" size="{ACS_CAPTION_SIZE}" '
            f'color="0" face="2">({xml_escape(letter)})</s>'
        )
        text_run = build_formatted_s_xml(
            f" {text}",
            font=ACS_LABEL_FONT,
            size=ACS_CAPTION_SIZE,
            color="0",
        )
        xml_parts.append(
            f'<t id="{tid}" p="{left_x:.2f} {y:.2f}" Z="{z}" '
            f'InterpretChemically="no">\n'
            f'{letter_run}{text_run}\n'
            f'</t>'
        )
        y += 14.0  # line spacing
    return "\n".join(xml_parts), y


def _build_vertical_arrow(
    x: float,
    top_y: float,
    bottom_y: float,
    ids: _IDGen,
    condition_lines: Optional[List[str]] = None,
    condition_side: str = "right",
    dashed: bool = False,
    nogo: bool = False,
    kind: Optional[str] = None,
) -> Tuple[str, int]:
    """
    Build a vertical down arrow with optional condition text beside it.

    Parameters
    ----------
    x : float
        Horizontal position of the arrow.
    top_y : float
        Y coordinate of the arrow tail (top, in CDXML y-down coords).
    bottom_y : float
        Y coordinate of the arrow head (bottom).
    condition_lines : list of str, optional
        Condition text lines placed beside the arrow.
    condition_side : "right" or "left"
        Which side of the arrow to place condition text.
    nogo : bool
        If True, adds ``NoGo="Cross"`` for failed arrow.

    Returns (xml_string, arrow_xml_id).
    """
    arrow_xml, aid = _build_arrow(x, top_y, x, bottom_y, ids, dashed=dashed, nogo=nogo, kind=kind)

    if not condition_lines:
        return arrow_xml, aid

    # Place condition text beside the arrow
    mid_y = (top_y + bottom_y) / 2.0
    text_offset = 12.0  # horizontal gap from arrow shaft
    if condition_side == "right":
        text_x = x + text_offset
        justification = "Left"
    else:
        text_x = x - text_offset
        justification = "Right"

    # p.y is the first line baseline. To vertically center the visual
    # text block on the arrow midpoint:
    #   visual_top = p.y - _CAP_HEIGHT
    #   visual_bottom = p.y + (n-1)*_LINE_ADVANCE + _DESCENT
    #   visual_center = visual_top + text_h/2 = mid_y
    #   → p.y = mid_y - text_h/2 + _CAP_HEIGHT
    text_h = _estimate_text_height(condition_lines)
    text_y = mid_y - text_h / 2.0 + _CAP_HEIGHT

    txt_xml, _ = _build_text_element(
        condition_lines, text_x, text_y, ids, justification=justification,
    )

    return arrow_xml + "\n" + txt_xml, aid


def _build_run_arrow(
    tail_x: float,
    head_x: float,
    y: float,
    input_label: str,
    output_label: str,
    ids: _IDGen,
) -> Tuple[str, List[int]]:
    """
    Build a simple run arrow with input/output labels.

    Returns (xml_string, [element_ids]).
    """
    parts = []
    all_ids = []

    # Arrow
    arrow_xml, arrow_id = _build_arrow(tail_x, y, head_x, y, ids)
    parts.append(arrow_xml)
    all_ids.append(arrow_id)

    # Input label (left of arrow) — center text vertically on arrow.
    # p.y is baseline; visual center = p.y - (_CAP_HEIGHT - _DESCENT)/2
    # To center on arrow at y: p.y = y + (_CAP_HEIGHT - _DESCENT)/2
    label_baseline_y = y + (_CAP_HEIGHT - _DESCENT) / 2.0
    inp_xml, inp_id = _build_text_element(
        [input_label], tail_x - 5.0, label_baseline_y, ids,
        justification="Right", use_formatting=False,
    )
    parts.append(inp_xml)
    all_ids.append(inp_id)

    # Output label (right of arrow) — center text vertically on arrow
    if output_label:
        out_xml, out_id = _build_text_element(
            [output_label], head_x + 5.0, label_baseline_y, ids,
            justification="Left", use_formatting=False,
        )
        parts.append(out_xml)
        all_ids.append(out_id)

    return "\n".join(parts), all_ids


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def _label_overhang(atom: Dict) -> float:
    """Horizontal reach of an atom label beyond its centred first character.

    Labels start on the atom ("Br", "OH") and ChemDraw may flip hydrogens to
    the left ("HO"), so the extra characters can reach either side.
    """
    symbol = atom.get("symbol", "C")
    hydrogens = atom.get("num_hydrogens") or 0
    if symbol == "C" and not atom.get("charge") and not atom.get("isotope"):
        return 0.0
    extra = _LABEL_CHAR_W * (len(symbol) - 1)
    if hydrogens:
        extra += LABEL_H_WIDTH + (_LABEL_DIGIT_W if hydrogens > 1 else 0.0)
    if atom.get("charge"):
        extra += _LABEL_DIGIT_W
    return extra


def _fragment_bbox(atoms: List[Dict]) -> Tuple[float, float, float, float]:
    """Compute (min_x, min_y, max_x, max_y) from atom positions.

    Includes padding for bond lines and the first label character
    (measured: +6.6pt width, +2.9pt height) plus the width of longer
    labels such as "Br" or "OH".  A lone atom's label is drawn on one side,
    as ChemDraw writes it (OH-, Br- rightwards; H2O, HCl leftwards).
    """
    if len(atoms) == 1:
        atom = atoms[0]
        overhang = _label_overhang(atom)
        left = hydrogens_lead(atom.get("symbol", "C"), atom.get("num_hydrogens"))
        return (
            atom["x"] - _FRAG_PAD_W - (overhang if left else 0.0),
            atom["y"] - _FRAG_PAD_H,
            atom["x"] + _FRAG_PAD_W + (0.0 if left else overhang),
            atom["y"] + _FRAG_PAD_H,
        )
    return (
        min(a["x"] - _FRAG_PAD_W - _label_overhang(a) for a in atoms),
        min(a["y"] for a in atoms) - _FRAG_PAD_H,
        max(a["x"] + _FRAG_PAD_W + _label_overhang(a) for a in atoms),
        max(a["y"] for a in atoms) + _FRAG_PAD_H,
    )


def _bbox_width(bbox: Tuple[float, float, float, float]) -> float:
    return bbox[2] - bbox[0]


def _bbox_height(bbox: Tuple[float, float, float, float]) -> float:
    return bbox[3] - bbox[1]


def _bbox_center(bbox: Tuple[float, float, float, float]) -> Tuple[float, float]:
    return ((bbox[0] + bbox[2]) / 2.0, (bbox[1] + bbox[3]) / 2.0)


def _shift_atoms(atoms: List[Dict], dx: float, dy: float) -> None:
    """Translate all atoms in place."""
    for a in atoms:
        a["x"] += dx
        a["y"] += dy


def _estimate_text_width(lines: List[str]) -> float:
    """Estimate text block width in CDXML points."""
    if not lines:
        return 0
    return max(len(ln) for ln in lines) * _CHAR_WIDTH


def _estimate_text_height(lines: List[str]) -> float:
    """Estimate visual text block height (cap-height + line advances + descent)."""
    n = len(lines)
    if n == 0:
        return 0.0
    return _CAP_HEIGHT + max(0, n - 1) * _LINE_ADVANCE + _DESCENT


# ---------------------------------------------------------------------------
# Structure resolver
# ---------------------------------------------------------------------------

ORIENTATIONS = ("flip_horizontal", "flip_vertical", "rotate_180")


def _orient_fragment(frag: ResolvedFragment, orientation: Optional[str]) -> ResolvedFragment:
    """Mirror or turn a resolved structure about its centre, keeping its stereochemistry.

    A mirror image swaps wedges and hashes and the side of ring double bonds, so the
    drawn molecule is the same compound; a half turn needs neither.
    """
    if not orientation:
        return frag
    if orientation not in ORIENTATIONS:
        raise ValueError(f"orientation must be one of {list(ORIENTATIONS)}, not {orientation!r}")
    cx, cy = _bbox_center(_fragment_bbox(frag.atoms))
    flip_x = orientation in ("flip_horizontal", "rotate_180")
    flip_y = orientation in ("flip_vertical", "rotate_180")
    for atom in frag.atoms:
        if flip_x:
            atom["x"] = 2 * cx - atom["x"]
        if flip_y:
            atom["y"] = 2 * cy - atom["y"]
    if orientation != "rotate_180":
        swap_cfg = {1: 6, 6: 1}
        swap_side = {"Left": "Right", "Right": "Left"}
        for bond in frag.bonds:
            bond["cfg"] = swap_cfg.get(bond.get("cfg"), bond.get("cfg"))
            if bond.get("double_pos") in swap_side:
                bond["double_pos"] = swap_side[bond["double_pos"]]
    frag.bbox = _fragment_bbox(frag.atoms)
    frag.cx, frag.cy = _bbox_center(frag.bbox)
    return frag


def _resolve_structure(
    ref: StructureRef,
    center_x: float = 200.0,
    center_y: float = 300.0,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> ResolvedFragment:
    """Resolve a StructureRef, then mirror or turn it as ``ref.orientation`` asks."""
    frag = _resolve_structure_unoriented(ref, center_x, center_y, source_data, reference_mol)
    if not ref.orientation:
        return frag
    frag = _orient_fragment(frag, ref.orientation)
    frag.ref = dataclass_replace(frag.ref, orientation=ref.orientation)
    return frag


def _resolve_structure_unoriented(
    ref: StructureRef,
    center_x: float = 200.0,
    center_y: float = 300.0,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> ResolvedFragment:
    """
    Resolve a StructureRef to atoms/bonds in CDXML coordinates.

    Resolution order (three tiers):
      1. Source JSON — species ID or shorthand lookup
      2. Declared SMILES — from StructureRef.smiles
      3. Name resolution — reagent_db → PubChem cascade

    Parameters
    ----------
    reference_mol : RDKit Mol, optional
        When provided, the structure is aligned to this reference via MCS.
        This orients shared scaffolds to match the reference (product) layout.
    """
    smiles = ref.smiles
    label = ref.label

    # Tier 1: Source JSON lookup (by ref.id)
    if source_data and not smiles:
        # Try exact id, then lowercase
        sp = source_data.get(ref.id) or source_data.get(ref.id.lower())
        if sp:
            smiles = sp.get("smiles")
            # JSON label overrides ref label when ref has none
            if label is None and sp.get("label"):
                label = sp["label"]

    # Tier 2: Declared SMILES
    if smiles:
        # A structure with generic rings is depicted with room for them (see _depict_with_rings).
        ring_hints = tuple((o["atoms"][0], o["atoms"][1], int(o.get("members", 6)))
                           for o in getattr(ref, "overlays", ()) or () if o.get("shape") == "ring")
        # Try MCS alignment to reference if available
        if reference_mol is not None and not ring_hints:
            try:
                from rdkit import Chem
                target_mol = Chem.MolFromSmiles(smiles)
                if target_mol is not None:
                    aligned = _align_mol_to_reference(
                        target_mol, reference_mol, center_x, center_y,
                    )
                    if aligned is not None:
                        atoms, bonds = aligned
                        bbox = _fragment_bbox(atoms)
                        cx, cy = _bbox_center(bbox)
                        resolved_ref = StructureRef(
                            id=ref.id, smiles=smiles, name=ref.name, file=ref.file,
                            cdxml_id=ref.cdxml_id, label=label,
                        )
                        return ResolvedFragment(
                            ref=resolved_ref, atoms=atoms, bonds=bonds,
                            bbox=bbox, cx=cx, cy=cy,
                        )
            except Exception:
                pass  # Fall through to standard coords

        result = _smiles_to_fragment_data(smiles, center_x, center_y, ring_hints)
        if result is None:
            raise ValueError(f"Failed to generate 2D coords for '{ref.id}' (SMILES: {smiles})")
        atoms, bonds = result
        bbox = _fragment_bbox(atoms)
        cx, cy = _bbox_center(bbox)
        resolved_ref = StructureRef(
            id=ref.id, smiles=smiles, name=ref.name, file=ref.file,
            cdxml_id=ref.cdxml_id, label=label,
        )
        return ResolvedFragment(
            ref=resolved_ref, atoms=atoms, bonds=bonds, bbox=bbox, cx=cx, cy=cy,
        )

    # Tier 3: Name resolution (reagent_db cascade)
    name = ref.name
    # If source_data had a species with a name, use it for resolution
    if not name and source_data:
        sp = source_data.get(ref.id) or source_data.get(ref.id.lower())
        if sp:
            name = sp.get("name")

    if name:
        try:
            from ..resolve.reagent_db import get_reagent_db
            db = get_reagent_db()
            entry = db.entry_for_name(name.lower())
            if entry and entry.get("smiles"):
                smi = entry["smiles"]
                if isinstance(smi, list):
                    smi = smi[0]
                ref_copy = StructureRef(
                    id=ref.id, smiles=smi, label=label,
                )
                return _resolve_structure(
                    ref_copy, center_x, center_y,
                    reference_mol=reference_mol,
                )
        except Exception:
            pass
        raise ValueError(
            f"Cannot resolve structure '{ref.id}' by name '{name}'. "
            f"Provide a SMILES string instead."
        )

    if ref.file:
        raise NotImplementedError(
            f"File-based structure loading not yet implemented for '{ref.id}'"
        )

    # Tier 4: Try the ID itself as a compound name in reagent_db
    try:
        from ..resolve.reagent_db import get_reagent_db
        db = get_reagent_db()
        entry = db.entry_for_name(ref.id.lower())
        if entry and entry.get("smiles"):
            smi = entry["smiles"]
            if isinstance(smi, list):
                smi = smi[0]
            ref_copy = StructureRef(
                id=ref.id, smiles=smi, label=label,
            )
            return _resolve_structure(
                ref_copy, center_x, center_y,
                reference_mol=reference_mol,
            )
    except ImportError:
        pass

    raise ValueError(
        f"Structure '{ref.id}' has no smiles, name, or file — cannot resolve. "
        f"Provide explicit smiles: or name: in the structures block."
    )


# ---------------------------------------------------------------------------
# Layout engine: linear (single step)
# ---------------------------------------------------------------------------

def _layout_linear(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
    arrow_y: float = 300.0,
) -> Tuple[str, float]:
    """
    Layout a single-step scheme: substrates → products.

    Returns (inner_xml, lowest_y).
    """
    xml, lowest_y, _, _ = _layout_steps_row(scheme, scheme.steps, ids, start_x=100.0,
                                            arrow_y=arrow_y, source_data=source_data,
                                            reference_mol=reference_mol)
    return xml, lowest_y


# ---------------------------------------------------------------------------
# Layout engine: sequential (multi-step, with wrap:repeat support)
# ---------------------------------------------------------------------------

def _split_into_rows(
    steps: List[StepDescriptor],
    steps_per_row: Optional[int],
) -> List[List[StepDescriptor]]:
    """Split steps into row groups for wrap:repeat."""
    if steps_per_row is None or steps_per_row >= len(steps):
        return [steps]
    rows = []
    for i in range(0, len(steps), steps_per_row):
        rows.append(steps[i:i + steps_per_row])
    return rows


def _split_serpentine_rows(
    steps: List[StepDescriptor],
    steps_per_row: int,
) -> Tuple[List[List[StepDescriptor]], List[Optional[StepDescriptor]]]:
    """
    Split steps into horizontal rows + transition steps for serpentine layout.

    Returns (rows, transitions) where:
      - rows[i] = list of steps rendered horizontally in row i
      - transitions[i] = the step between row i and row i+1 (rendered as
        a vertical arrow), or None if row i is the last row.
    """
    rows: List[List[StepDescriptor]] = []
    transitions: List[Optional[StepDescriptor]] = []
    i = 0
    while i < len(steps):
        row = steps[i:i + steps_per_row]
        rows.append(row)
        i += steps_per_row
        if i < len(steps):
            # The next step becomes the vertical transition arrow
            transitions.append(steps[i])
            i += 1
        else:
            transitions.append(None)
    return rows, transitions


_ROW_START_X = 80.0


def _auto_steps_per_row(
    scheme: SchemeDescriptor,
    budget: float,
    source_data: Optional[Dict[str, Dict]],
    reference_mol: Optional[Any],
) -> int:
    """Largest steps-per-row whose widest row still fits in ``budget`` pt.

    Rows are measured by laying them out on a throwaway id counter, so the
    result matches what the real layout will produce. When even one step per
    row is wider than the budget, one step per row is returned and the caller's
    geometry check has to flag the overflow.
    """
    steps = scheme.steps
    token = _fragment_receipts.set(None)  # trial rows are not part of the drawing
    try:
        for per_row in range(len(steps), 0, -1):
            rows = (_split_serpentine_rows(steps, per_row)[0]
                    if scheme.wrap == "serpentine"
                    else _split_into_rows(steps, per_row))
            widest = 0.0
            for index, row in enumerate(rows):
                carried = scheme.wrap == "continue" and index > 0 and _continues(rows[index - 1], row)
                _, _, edge_x, _ = _layout_steps_row(
                    scheme, row, _IDGen(1), start_x=_ROW_START_X, arrow_y=300.0,
                    source_data=source_data, reference_mol=reference_mol,
                    skip_first_substrate_id=-1 if carried else None,
                    skip_first_substrate_cursor_x=_ROW_START_X if carried else None)
                widest = max(widest, edge_x - _ROW_START_X)
            if widest <= budget:
                return per_row
        return 1
    finally:
        _fragment_receipts.reset(token)


def _continues(previous: List[StepDescriptor], row: List[StepDescriptor]) -> bool:
    """Whether a row's first step starts from the structure the previous row ended with."""
    return bool(previous and row and previous[-1].products and row[0].substrates
                and row[0].substrates[0] == previous[-1].products[-1])


def _layout_sequential(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """
    Layout a multi-step sequential scheme.

    Supports wrap:continue (multi-row L->R; a row after the first starts with
    the arrow that leaves the structure the row above ended with, so every
    structure is drawn once), wrap:repeat (the same, but each row starts by
    drawing that structure again), wrap:serpentine (zigzag with vertical
    arrows between rows), and single-row (no wrapping).

    Returns (inner_xml, lowest_y).
    """
    wrap = scheme.wrap
    steps_per_row = scheme.steps_per_row
    if (steps_per_row is None and scheme.max_row_width
            and wrap in ("repeat", "continue", "serpentine") and len(scheme.steps) > 1):
        steps_per_row = _auto_steps_per_row(
            scheme, scheme.max_row_width, source_data, reference_mol)

    # Single row: no wrapping needed
    if (wrap not in ("repeat", "continue", "serpentine")
            or steps_per_row is None
            or steps_per_row >= len(scheme.steps)):
        xml, lowest_y, _, _ = _layout_steps_row(
            scheme, scheme.steps, ids, start_x=80.0, arrow_y=300.0,
            source_data=source_data, reference_mol=reference_mol)
        return xml, lowest_y

    if wrap == "serpentine":
        return _layout_serpentine(scheme, ids, source_data=source_data,
                                  reference_mol=reference_mol)

    # --- Multi-row wrap:repeat ---
    rows = _split_into_rows(scheme.steps, steps_per_row)

    # Assign run arrows to rows with step numbers adjusted to be 1-indexed
    # within each row.  Step N (1-indexed overall) -> row (N-1)//spr,
    # adjusted step = N - row_idx * spr.
    row_run_arrows: Dict[int, List[StepRunArrows]] = {}
    for sra in scheme.run_arrows:
        row_idx = (sra.step - 1) // steps_per_row
        adjusted_step = sra.step - row_idx * steps_per_row
        row_run_arrows.setdefault(row_idx, []).append(
            StepRunArrows(step=adjusted_step, runs=sra.runs))

    xml_parts: List[str] = []
    arrow_y = 300.0
    lowest_y = arrow_y
    previous_bottom: Optional[float] = None

    last_product_id: Optional[int] = None
    for row_idx, row_steps in enumerate(rows):
        row_ra = row_run_arrows.get(row_idx, [])
        # wrap:continue: the structure that ended the row above is not drawn again; this row
        # starts at the margin with the arrow that leaves it.
        carried = (wrap == "continue" and row_idx > 0 and last_product_id is not None
                   and _continues(rows[row_idx - 1], row_steps))

        def layout(at_y: float, row_ids: _IDGen, row_steps=row_steps, row_ra=row_ra, carried=carried,
                   carried_id=last_product_id):
            return _layout_steps_row(
                scheme, row_steps, row_ids,
                start_x=80.0, arrow_y=at_y,
                source_data=source_data,
                run_arrows=row_ra,
                reference_mol=reference_mol,
                skip_first_substrate_id=carried_id if carried else None,
                skip_first_substrate_cursor_x=80.0 if carried else None,
            )

        arrow_y = _next_arrow_y(layout, arrow_y, previous_bottom, LAYOUT_ROW_GAP)
        row_xml, row_lowest, _, row_info = layout(arrow_y, ids)
        last_product_id = row_info.get("last_product_id")
        xml_parts.append(row_xml)
        previous_bottom = _drawn_extent(row_xml)[1]
        lowest_y = max(row_lowest, previous_bottom)

    return "\n".join(xml_parts), lowest_y


def _next_arrow_y(layout, arrow_y: float, previous_bottom: Optional[float], gap: float) -> float:
    """Arrow height that puts a row's top ``gap`` below ``previous_bottom``.

    The row is measured by laying it out on a throwaway id counter.
    """
    if previous_bottom is None:
        return arrow_y
    token = _fragment_receipts.set(None)  # trial rows are not part of the drawing
    try:
        top, _ = _drawn_extent(layout(0.0, _IDGen(1))[0])
    finally:
        _fragment_receipts.reset(token)
    return previous_bottom + gap - top


def _drawn_extent(row_xml: str) -> Tuple[float, float]:
    """Top and bottom of what a laid-out row draws: structures with their labels, text, arrows."""
    from xml.etree import ElementTree as ET

    tops, bottoms = [], []
    for el in ET.fromstring(f"<row>{row_xml}</row>"):
        if el.tag == "fragment":
            box = drawn_fragment_bbox(el)
        elif el.get("BoundingBox"):
            box = tuple(float(v) for v in el.get("BoundingBox").split()[:4])
        else:
            continue
        if box:
            tops.append(box[1])
            bottoms.append(box[3])
    return (min(tops), max(bottoms)) if tops else (0.0, 0.0)


def _layout_serpentine(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """
    Layout a serpentine (zigzag) multi-row scheme.

    Row 1 flows L→R, vertical down arrow at right edge, Row 2 flows R→L,
    vertical down arrow at left edge, Row 3 flows L→R, etc.

    Transition steps between rows are rendered as vertical arrows with
    their conditions placed beside the arrow.

    Returns (inner_xml, lowest_y).
    """
    steps_per_row = scheme.steps_per_row
    left_margin = 80.0

    rows, transitions = _split_serpentine_rows(scheme.steps, steps_per_row)

    bond_len = ACS_BOND_LENGTH
    frag_gap = bond_len * LAYOUT_FRAG_GAP_BONDS * scheme.gap_scale
    inter_gap = max(bond_len * LAYOUT_INTER_GAP_BONDS, LAYOUT_PLUS_GAP) * scheme.gap_scale

    xml_parts: List[str] = []
    arrow_y = 300.0
    right_edge = 0.0  # will be set from first row's width
    last_product_frag_id: Optional[int] = None
    last_product_cursor_x: Optional[float] = None  # cursor pos after pre-placed substrate

    def _get_ref(sid: str) -> StructureRef:
        if sid in scheme.structures:
            return scheme.structures[sid]
        return StructureRef(id=sid)

    for row_idx, row_steps in enumerate(rows):
        direction = "ltr" if row_idx % 2 == 0 else "rtl"

        if direction == "ltr":
            start_x = left_margin
        else:
            start_x = right_edge

        # For rows after the first, the first substrate was placed by
        # the transition step's product. Pass skip_first_substrate_id
        # so _layout_steps_row doesn't redraw it.
        skip_id = last_product_frag_id if row_idx > 0 else None
        skip_cursor = last_product_cursor_x if row_idx > 0 else None

        row_xml, lowest_y, row_edge, row_info = _layout_steps_row(
            scheme, row_steps, ids,
            start_x=start_x, arrow_y=arrow_y,
            source_data=source_data,
            direction=direction,
            skip_first_substrate_id=skip_id,
            skip_first_substrate_cursor_x=skip_cursor,
            reference_mol=reference_mol,
        )
        xml_parts.append(row_xml)

        # After row 1, record the right edge for RTL alignment
        if row_idx == 0:
            right_edge = row_edge

        # --- Transition: vertical arrow + product of transition step ---
        transition = transitions[row_idx]
        if transition is not None:
            # Vertical arrow position: below the CENTER of the last product
            # in the current row (not at the row edge).
            if direction == "ltr":
                vert_x = row_info["last_product_cx"]
            else:
                vert_x = row_info["last_product_cx"]

            # Vertical arrow geometry:
            # - tail starts BELOW the last product (below its compound label)
            # - head ends with enough room for the next row's product
            last_prod_bottom = row_info["last_product_bottom"]
            vert_top_y = last_prod_bottom + 8.0  # below the last product + label
            vert_bottom_y = vert_top_y + ROW_GAP  # adequate vertical span

            # Condition text from the transition step
            cond_lines: List[str] = []
            if transition.below_arrow and transition.below_arrow.text:
                cond_lines = transition.below_arrow.text[:]
            if transition.above_arrow and transition.above_arrow.text:
                cond_lines = transition.above_arrow.text + cond_lines
            cond_lines = [normalize_condition_text(t) for t in cond_lines]
            if transition.yield_:
                cond_lines.append(transition.yield_)

            # Condition text placement: right side for right-edge arrows,
            # left side for left-edge arrows
            if direction == "ltr":
                cond_side = "right"
            else:
                cond_side = "left"

            vert_xml, _ = _build_vertical_arrow(
                vert_x, vert_top_y, vert_bottom_y, ids,
                condition_lines=cond_lines if cond_lines else None,
                condition_side=cond_side,
            )
            xml_parts.append(vert_xml)

            # Place the transition step's product below the vertical
            # arrow. This becomes the first substrate of the next row.
            prod_ref = _get_ref(transition.products[0])
            prod_frag = _resolve_structure(
                prod_ref,
                center_x=vert_x,
                center_y=vert_bottom_y + 25.0,
                source_data=source_data,
                reference_mol=reference_mol,
            )
            # Shift so product top is below arrowhead with clearance
            prod_h = _bbox_height(prod_frag.bbox)
            target_cy = vert_bottom_y + frag_gap * 1.0 + prod_h / 2.0
            dx = vert_x - prod_frag.cx
            dy = target_cy - prod_frag.cy
            _shift_atoms(prod_frag.atoms, dx, dy)
            prod_frag.bbox = _fragment_bbox(prod_frag.atoms)
            prod_frag.cx, prod_frag.cy = _bbox_center(prod_frag.bbox)

            frag_xml, _, frag_id = _build_fragment(
                prod_frag.atoms, prod_frag.bonds, ids,
            )
            xml_parts.append(frag_xml)

            # Compound label
            if prod_frag.ref.label:
                lbl_xml, _ = _build_label_element(
                    prod_frag.ref.label, prod_frag.cx,
                    prod_frag.bbox[3] + 14.0, ids,
                )
                xml_parts.append(lbl_xml)

            last_product_frag_id = frag_id

            # Compute cursor position for the next row's first step
            # (after the pre-placed substrate).
            # Next row direction is opposite of current.
            next_direction = "rtl" if direction == "ltr" else "ltr"
            if next_direction == "rtl":
                # RTL: cursor is at the LEFT edge of the substrate
                last_product_cursor_x = prod_frag.bbox[0] - inter_gap
            else:
                # LTR: cursor is at the RIGHT edge of the substrate
                last_product_cursor_x = prod_frag.bbox[2] + inter_gap

            # Next row's arrow_y = below the transition product
            arrow_y = prod_frag.cy
        else:
            last_product_frag_id = None
            last_product_cursor_x = None

    return "\n".join(xml_parts), lowest_y


# ---------------------------------------------------------------------------
# Layout engine: divergent (one SM → multiple products on branch rows)
# ---------------------------------------------------------------------------

# Gap between the bottom of one branch row and the top of the next (halved
# between rows: branch rows have no run arrows, so they can sit closer).
_DIVERGENT_BRANCH_GAP = 2.5 * ACS_BOND_LENGTH
_LABEL_ROOM_ABOVE = 20.0  # a compound number under a structure above an arrow


def _layout_divergent(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """
    Layout a divergent scheme: one SM gives multiple products.

    Detects the shared substrate (appears in multiple steps) and renders it
    once on the left.  Each step is a row with its own horizontal arrow; the
    arrow tails share one column joined by a vertical line.

    Returns (inner_xml, lowest_y).
    """
    steps = scheme.steps
    if not steps:
        return "", 300.0

    # --- Identify shared substrate ---
    # Count how many steps each substrate ID appears in
    from collections import Counter
    sub_counts: Counter = Counter()
    for step in steps:
        for sid in step.substrates:
            sub_counts[sid] += 1

    # The shared substrate is the one appearing in the most steps
    shared_id = sub_counts.most_common(1)[0][0] if sub_counts else None

    # If no substrate appears more than once, just stack steps vertically
    if shared_id is None or sub_counts[shared_id] < 2:
        # Fall back to stacking as independent rows
        return _layout_divergent_stacked(scheme, ids, source_data,
                                        reference_mol=reference_mol)

    # Separate steps: those sharing the substrate vs others
    shared_steps = [s for s in steps if shared_id in s.substrates]
    other_steps = [s for s in steps if shared_id not in s.substrates]

    # --- Resolve shared substrate once ---
    def _get_ref(sid: str) -> StructureRef:
        if sid in scheme.structures:
            return scheme.structures[sid]
        return StructureRef(id=sid)

    bond_len = ACS_BOND_LENGTH
    frag_gap = bond_len * LAYOUT_FRAG_GAP_BONDS * scheme.gap_scale
    inter_gap = max(bond_len * LAYOUT_INTER_GAP_BONDS, LAYOUT_PLUS_GAP) * scheme.gap_scale
    min_arrow_len = MIN_ARROW_BONDS * bond_len

    start_x = 100.0
    arrow_y = 300.0

    sm_ref = _get_ref(shared_id)
    sm_frag = _resolve_structure(sm_ref, center_x=200.0, center_y=arrow_y,
                                 source_data=source_data,
                                 reference_mol=reference_mol)
    # Position SM at start
    sm_w = _bbox_width(sm_frag.bbox)
    cx_target = start_x + sm_w / 2.0
    dx = cx_target - sm_frag.cx
    dy = arrow_y - sm_frag.cy
    _shift_atoms(sm_frag.atoms, dx, dy)
    sm_frag.bbox = _fragment_bbox(sm_frag.atoms)
    sm_frag.cx, sm_frag.cy = _bbox_center(sm_frag.bbox)

    xml_parts: List[str] = []
    all_frag_ids: List[int] = []
    step_metadata: List[Dict] = []

    # Build SM fragment
    sm_xml, _, sm_frag_id = _build_fragment(sm_frag.atoms, sm_frag.bonds, ids)
    xml_parts.append(sm_xml)
    all_frag_ids.append(sm_frag_id)

    # SM label
    if sm_frag.ref.label:
        lbl_xml, _ = _build_label_element(
            sm_frag.ref.label, sm_frag.cx, sm_frag.bbox[3] + 14.0, ids)
        xml_parts.append(lbl_xml)

    sm_right_x = sm_frag.bbox[2]
    lowest_y = arrow_y + 60.0

    # Every branch is a row with its own horizontal arrow.  The arrows share one
    # tail column right of the starting material, joined by a vertical connector,
    # so each product reads as made from the starting material (a column of
    # stacked vertical arrows reads as a chain, 1 -> 3 -> 4).
    # A branch's other substrates (the co-substrate of a Michael addition, say) and
    # any reagent structures stand above its arrow, in one row over the text.  The
    # co-substrates are recorded as reactants of the step, the reagents above it.
    branches = []
    for step in shared_steps:
        above_text = [normalize_condition_text(t) for t in step.above_arrow.text] if step.above_arrow else []
        below_text = [normalize_condition_text(t) for t in step.below_arrow.text] if step.below_arrow else []
        if step.yield_:
            below_text.append(step.yield_)
        products = [_resolve_structure(_get_ref(sid), source_data=source_data,
                                       reference_mol=reference_mol)
                    for sid in step.products]
        over = [(_resolve_structure(_get_ref(sid), source_data=source_data,
                                    reference_mol=reference_mol), "reactant")
                for sid in step.substrates if sid != shared_id]
        if step.above_arrow:
            over += [(_resolve_structure(_get_ref(sid), source_data=source_data,
                                         reference_mol=reference_mol), "above")
                     for sid in step.above_arrow.structures]
        branches.append((step, above_text, below_text, products, over))

    def _over_width(over) -> float:
        if not over:
            return 0.0
        return sum(_bbox_width(f.bbox) for f, _ in over) + inter_gap * (len(over) - 1)

    def _over_height(over) -> float:
        return max((_bbox_height(f.bbox) + (_LABEL_ROOM_ABOVE if f.ref.label else 0.0)
                    for f, _ in over), default=0.0)

    content_w = max((max(_estimate_text_width(above), _estimate_text_width(below), _over_width(over), 0)
                     for _, above, below, _, over in branches), default=0)
    arrow_len = max(content_w + 10.0, min_arrow_len)
    tail_x = sm_right_x + frag_gap * 0.3
    head_x = tail_x + arrow_len
    arrow_mid_x = (tail_x + head_x) / 2.0
    tail_ys: List[float] = []
    previous_bottom: Optional[float] = None

    for step, above_text, below_text, products, over in branches:
        # Extent of the row above and below its arrow line.
        above_h = 0.0
        text_h = 0.0
        if above_text:
            text_h = (max(0, len(above_text) - 1) * _LINE_ADVANCE + _DESCENT + _CAP_HEIGHT)
            above_h = LAYOUT_BELOW_GAP + text_h
        # The structures above the arrow sit over the text (4 pt apart), or
        # LAYOUT_ABOVE_GAP over the arrow when there is no text.
        over_base = (LAYOUT_BELOW_GAP + text_h + 4.0) if above_text else LAYOUT_ABOVE_GAP
        if over:
            above_h = over_base + _over_height(over)
        below_h = 0.0
        if below_text:
            below_h = (LAYOUT_BELOW_GAP + _CAP_HEIGHT
                       + max(0, len(below_text) - 1) * _LINE_ADVANCE + _DESCENT)
        for frag in products:
            above_h = max(above_h, frag.cy - frag.bbox[1])
            # The compound label sits 14 pt below the structure.
            below_h = max(below_h, frag.bbox[3] - frag.cy + (14.0 + _DESCENT if frag.ref.label else 0.0))
        if previous_bottom is None:
            branch_y = arrow_y
        else:
            branch_y = previous_bottom + _DIVERGENT_BRANCH_GAP / 2.0 + above_h
        tail_ys.append(branch_y)

        step_meta: Dict[str, Any] = {
            "reactant_ids": [sm_frag_id],
            "product_ids": [],
            "arrow_id": 0,
            "above_ids": [],
            "below_ids": [],
        }
        arrow_xml, arrow_id = _build_arrow(
            tail_x, branch_y, head_x, branch_y, ids,
            dashed=(step.arrow_style == "dashed"),
            nogo=(step.arrow_style == "failed"),
            kind=_arrow_kind(step.arrow_style))
        xml_parts.append(arrow_xml)
        step_meta["arrow_id"] = arrow_id

        # Above-arrow text (p.y = baseline of the first line), its visual bottom
        # LAYOUT_BELOW_GAP above the arrow like the text below it.
        if above_text:
            text_below_baseline = max(0, len(above_text) - 1) * _LINE_ADVANCE + _DESCENT
            text_y = branch_y - LAYOUT_BELOW_GAP - text_below_baseline
            txt_xml, txt_id = _build_text_element(above_text, arrow_mid_x, text_y, ids,
                                                  use_formatting=True)
            xml_parts.append(txt_xml)
            step_meta["above_ids"].append(txt_id)
        if below_text:
            text_y = branch_y + LAYOUT_BELOW_GAP + _CAP_HEIGHT
            txt_xml, txt_id = _build_text_element(below_text, arrow_mid_x, text_y, ids)
            xml_parts.append(txt_xml)
            step_meta["below_ids"].append(txt_id)

        over_x = arrow_mid_x - _over_width(over) / 2.0
        for frag, role in over:
            fw, fh = _bbox_width(frag.bbox), _bbox_height(frag.bbox)
            label_h = _LABEL_ROOM_ABOVE if frag.ref.label else 0.0
            _shift_atoms(frag.atoms, over_x + fw / 2.0 - frag.cx,
                         branch_y - over_base - label_h - fh / 2.0 - frag.cy)
            frag.bbox = _fragment_bbox(frag.atoms)
            frag.cx, frag.cy = _bbox_center(frag.bbox)
            ofrag_xml, _, ofrag_id = _build_fragment(frag.atoms, frag.bonds, ids)
            xml_parts.append(ofrag_xml)
            all_frag_ids.append(ofrag_id)
            step_meta["reactant_ids" if role == "reactant" else "above_ids"].append(ofrag_id)
            if frag.ref.label:
                lbl_xml, _ = _build_label_element(frag.ref.label, frag.cx, frag.bbox[3] + 14.0, ids)
                xml_parts.append(lbl_xml)
            over_x += fw + inter_gap

        cursor_x = head_x + frag_gap * 1.0
        for prod_frag in products:
            pw = _bbox_width(prod_frag.bbox)
            _shift_atoms(prod_frag.atoms, cursor_x + pw / 2.0 - prod_frag.cx, branch_y - prod_frag.cy)
            prod_frag.bbox = _fragment_bbox(prod_frag.atoms)
            prod_frag.cx, prod_frag.cy = _bbox_center(prod_frag.bbox)
            pfrag_xml, _, pfrag_id = _build_fragment(prod_frag.atoms, prod_frag.bonds, ids)
            xml_parts.append(pfrag_xml)
            step_meta["product_ids"].append(pfrag_id)
            all_frag_ids.append(pfrag_id)
            if prod_frag.ref.label:
                lbl_xml, _ = _build_label_element(
                    prod_frag.ref.label, prod_frag.cx, prod_frag.bbox[3] + 14.0, ids)
                xml_parts.append(lbl_xml)
            cursor_x = prod_frag.bbox[2] + inter_gap

        previous_bottom = branch_y + below_h
        lowest_y = max(lowest_y, previous_bottom)
        step_metadata.append(step_meta)

    if len(tail_ys) > 1:
        line_id, line_z = ids.next(), ids.next()
        xml_parts.append(
            f'<graphic id="{line_id}" Z="{line_z}" GraphicType="Line" '
            f'BoundingBox="{tail_x:.2f} {tail_ys[0]:.2f} {tail_x:.2f} {tail_ys[-1]:.2f}"/>')

    # Steps that do not start from the shared substrate are drawn as their own
    # rows underneath, rather than dropped.
    previous_bottom = _drawn_extent("\n".join(xml_parts))[1]
    for step in other_steps:
        row_scheme = SchemeDescriptor(structures=scheme.structures, steps=[step],
                                      layout="linear", gap_scale=scheme.gap_scale)

        def row_layout(at_y, row_ids, row_scheme=row_scheme):
            return _layout_linear(row_scheme, row_ids, source_data=source_data,
                                  reference_mol=reference_mol, arrow_y=at_y)

        row_y = _next_arrow_y(row_layout, arrow_y, previous_bottom, _DIVERGENT_BRANCH_GAP)
        row_xml, row_lowest = row_layout(row_y, ids)
        xml_parts.append(row_xml)
        previous_bottom = _drawn_extent(row_xml)[1]
        lowest_y = max(lowest_y, row_lowest, previous_bottom)

    # --- Build <scheme> ---
    scheme_id = ids.next()
    scheme_parts = [f'<scheme id="{scheme_id}">']
    for meta in step_metadata:
        step_id = ids.next()
        attrs = [f'id="{step_id}"']
        if meta["reactant_ids"]:
            attrs.append(f'ReactionStepReactants="{" ".join(str(x) for x in meta["reactant_ids"])}"')
        if meta["product_ids"]:
            attrs.append(f'ReactionStepProducts="{" ".join(str(x) for x in meta["product_ids"])}"')
        attrs.append(f'ReactionStepArrows="{meta["arrow_id"]}"')
        if meta["above_ids"]:
            attrs.append(f'ReactionStepObjectsAboveArrow="{" ".join(str(x) for x in meta["above_ids"])}"')
        if meta["below_ids"]:
            attrs.append(f'ReactionStepObjectsBelowArrow="{" ".join(str(x) for x in meta["below_ids"])}"')
        scheme_parts.append(f'<step {" ".join(attrs)}/>')
    scheme_parts.append('</scheme>')
    xml_parts.append("\n".join(scheme_parts))

    return "\n".join(xml_parts), lowest_y


_FAN_DIRECTIONS = ("left", "right", "down", "up")


def _layout_fan(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """One starting material in the middle, an arrow leaving it for each product: left, right,
    then down and up (BJOC 22-103 figure 1 panel c). Each branch keeps its conditions, and its
    co-substrates stand over its arrow (beside a vertical one) and count as reactants.

    Routes that are not one starting material with two to four products, that have steps not
    starting from it, or that would be wider than ``max_row_width``, are drawn by the divergent
    layout instead.
    """
    from collections import Counter

    steps = scheme.steps
    counts = Counter(sid for step in steps for sid in step.substrates)
    shared_id = counts.most_common(1)[0][0] if counts else None
    shared = [s for s in steps if shared_id in s.substrates]
    if shared_id is None or not 2 <= len(shared) <= len(_FAN_DIRECTIONS) or len(shared) != len(steps):
        scheme.drawn_layout = "divergent"
        return _layout_divergent(scheme, ids, source_data=source_data, reference_mol=reference_mol)

    def get_ref(sid: str) -> StructureRef:
        return scheme.structures.get(sid) or StructureRef(id=sid)

    def resolve(sid: str) -> _ResolvedFragment:
        return _resolve_structure(get_ref(sid), source_data=source_data, reference_mol=reference_mol)

    bond_len = ACS_BOND_LENGTH
    frag_gap = bond_len * LAYOUT_FRAG_GAP_BONDS * scheme.gap_scale
    inter_gap = max(bond_len * LAYOUT_INTER_GAP_BONDS, LAYOUT_PLUS_GAP) * scheme.gap_scale
    min_arrow_len = MIN_ARROW_BONDS * bond_len
    label_room = 14.0 + _DESCENT

    branches = []
    for step, direction in zip(shared, _FAN_DIRECTIONS):
        above = [normalize_condition_text(t) for t in step.above_arrow.text] if step.above_arrow else []
        below = [normalize_condition_text(t) for t in step.below_arrow.text] if step.below_arrow else []
        if step.yield_:
            below.append(step.yield_)
        over = [(resolve(sid), "reactant") for sid in step.substrates if sid != shared_id]
        if step.above_arrow:
            over += [(resolve(sid), "above") for sid in step.above_arrow.structures]
        products = [resolve(sid) for sid in step.products]
        over_w = sum(_bbox_width(f.bbox) for f, _ in over) + inter_gap * max(0, len(over) - 1)
        if direction in ("left", "right"):
            length = max(_estimate_text_width(above) + 10.0, _estimate_text_width(below) + 10.0,
                         over_w + 10.0, min_arrow_len)
        else:
            lines = len(above) + len(below)
            length = max(min_arrow_len, lines * _LINE_ADVANCE + 10.0)
        branches.append({"step": step, "dir": direction, "above": above, "below": below, "over": over,
                         "products": products, "len": length, "over_w": over_w})

    sm = resolve(shared_id)
    sm_w, sm_h = _bbox_width(sm.bbox), _bbox_height(sm.bbox)

    def reach(branch) -> float:
        """How far a branch reaches from the starting material's edge along its direction."""
        prods = branch["products"]
        if branch["dir"] in ("left", "right"):
            size = sum(_bbox_width(f.bbox) for f in prods) + inter_gap * max(0, len(prods) - 1)
        else:
            size = max((_bbox_height(f.bbox) + (label_room if f.ref.label else 0.0) for f in prods), default=0.0)
        return frag_gap * 1.3 + branch["len"] + size

    left = max((reach(b) for b in branches if b["dir"] == "left"), default=0.0)
    right = max((reach(b) for b in branches if b["dir"] == "right"), default=0.0)
    if scheme.max_row_width and left + sm_w + right > scheme.max_row_width:
        # Wider than the row budget: the stacked divergent layout fits where the fan does not.
        scheme.drawn_layout = "divergent"
        return _layout_divergent(scheme, ids, source_data=source_data, reference_mol=reference_mol)
    up = max((reach(b) for b in branches if b["dir"] == "up"), default=0.0)
    over_tall = max((_bbox_height(f.bbox) + 40.0 for b in branches if b["dir"] in ("left", "right")
                     for f, _ in b["over"]), default=0.0)
    cx = 80.0 + left + sm_w / 2.0
    cy = 120.0 + max(up, over_tall) + sm_h / 2.0
    _shift_atoms(sm.atoms, cx - sm.cx, cy - sm.cy)
    sm.bbox = _fragment_bbox(sm.atoms)
    sm.cx, sm.cy = _bbox_center(sm.bbox)

    xml_parts: List[str] = []
    sm_xml, _, sm_id = _build_fragment(sm.atoms, sm.bonds, ids)
    xml_parts.append(sm_xml)
    sm_bottom = sm.bbox[3]
    if sm.ref.label:
        lbl_xml, _ = _build_label_element(sm.ref.label, sm.cx, sm.bbox[3] + 14.0, ids)
        xml_parts.append(lbl_xml)
        sm_bottom = sm.bbox[3] + label_room
    lowest_y = sm_bottom

    def place(frag, x, y, meta_key, meta):
        _shift_atoms(frag.atoms, x - frag.cx, y - frag.cy)
        frag.bbox = _fragment_bbox(frag.atoms)
        frag.cx, frag.cy = _bbox_center(frag.bbox)
        fxml, _, fid = _build_fragment(frag.atoms, frag.bonds, ids)
        xml_parts.append(fxml)
        meta[meta_key].append(fid)
        bottom = frag.bbox[3]
        if frag.ref.label:
            lxml, _ = _build_label_element(frag.ref.label, frag.cx, frag.bbox[3] + 14.0, ids)
            xml_parts.append(lxml)
            bottom += label_room
        return bottom

    metadata: List[Dict[str, Any]] = []
    for b in branches:
        meta: Dict[str, Any] = {"reactant_ids": [sm_id], "product_ids": [], "arrow_id": 0,
                                "above_ids": [], "below_ids": []}
        step = b["step"]
        if b["dir"] in ("left", "right"):
            sign = -1.0 if b["dir"] == "left" else 1.0
            y = sm.cy
            tail = (sm.bbox[0] - frag_gap * 0.3) if sign < 0 else (sm.bbox[2] + frag_gap * 0.3)
            head = tail + sign * b["len"]
            tail_xy, head_xy = (tail, y), (head, y)
            mid = (tail + head) / 2.0
            text_h = 0.0
            if b["above"]:
                text_h = max(0, len(b["above"]) - 1) * _LINE_ADVANCE + _DESCENT + _CAP_HEIGHT
                txt, tid = _build_text_element(b["above"], mid, y - LAYOUT_BELOW_GAP - (text_h - _CAP_HEIGHT), ids,
                                               use_formatting=True)
                xml_parts.append(txt)
                meta["above_ids"].append(tid)
            if b["below"]:
                txt, tid = _build_text_element(b["below"], mid, y + LAYOUT_BELOW_GAP + _CAP_HEIGHT, ids)
                xml_parts.append(txt)
                meta["below_ids"].append(tid)
                lowest_y = max(lowest_y, y + LAYOUT_BELOW_GAP + _CAP_HEIGHT
                               + (len(b["below"]) - 1) * _LINE_ADVANCE + _DESCENT)
            base = (LAYOUT_BELOW_GAP + text_h + 4.0) if b["above"] else LAYOUT_ABOVE_GAP
            ox = mid - b["over_w"] / 2.0
            for frag, role in b["over"]:
                fw, fh = _bbox_width(frag.bbox), _bbox_height(frag.bbox)
                lift = _LABEL_ROOM_ABOVE if frag.ref.label else 0.0
                place(frag, ox + fw / 2.0, y - base - lift - fh / 2.0,
                      "reactant_ids" if role == "reactant" else "above_ids", meta)
                ox += fw + inter_gap
            cursor = head + sign * frag_gap
            for frag in b["products"]:
                fw = _bbox_width(frag.bbox)
                lowest_y = max(lowest_y, place(frag, cursor + sign * fw / 2.0, y, "product_ids", meta))
                cursor += sign * (fw + inter_gap)
        else:
            sign = 1.0 if b["dir"] == "down" else -1.0
            x = sm.cx
            tail = (sm_bottom + frag_gap * 0.3) if sign > 0 else (sm.bbox[1] - frag_gap * 0.3)
            head = tail + sign * b["len"]
            tail_xy, head_xy = (x, tail), (x, head)
            mid = (tail + head) / 2.0
            lines = b["above"] + b["below"]
            if lines:
                block = (len(lines) - 1) * _LINE_ADVANCE + _CAP_HEIGHT + _DESCENT
                txt, tid = _build_text_element(lines, x + 6.0, mid - block / 2.0 + _CAP_HEIGHT, ids,
                                               justification="Left", use_formatting=True)
                xml_parts.append(txt)
                meta["above_ids"].append(tid)
            ox = x - 8.0
            for frag, role in b["over"]:
                fw = _bbox_width(frag.bbox)
                place(frag, ox - fw / 2.0, mid, "reactant_ids" if role == "reactant" else "above_ids", meta)
                ox -= fw + inter_gap
            for frag in b["products"]:
                fh = _bbox_height(frag.bbox)
                if sign > 0:
                    y = head + frag_gap + fh / 2.0
                else:
                    y = head - frag_gap - (label_room if frag.ref.label else 0.0) - fh / 2.0
                lowest_y = max(lowest_y, place(frag, x, y, "product_ids", meta))
        arrow_xml, arrow_id = _build_arrow(tail_xy[0], tail_xy[1], head_xy[0], head_xy[1], ids,
                                           dashed=(step.arrow_style == "dashed"),
                                           nogo=(step.arrow_style == "failed"),
                                           kind=_arrow_kind(step.arrow_style))
        xml_parts.append(arrow_xml)
        meta["arrow_id"] = arrow_id
        lowest_y = max(lowest_y, head_xy[1])
        metadata.append(meta)

    scheme_parts = [f'<scheme id="{ids.next()}">']
    for meta in metadata:
        attrs = [f'id="{ids.next()}"',
                 f'ReactionStepReactants="{" ".join(str(x) for x in meta["reactant_ids"])}"',
                 f'ReactionStepProducts="{" ".join(str(x) for x in meta["product_ids"])}"',
                 f'ReactionStepArrows="{meta["arrow_id"]}"']
        if meta["above_ids"]:
            attrs.append(f'ReactionStepObjectsAboveArrow="{" ".join(str(x) for x in meta["above_ids"])}"')
        if meta["below_ids"]:
            attrs.append(f'ReactionStepObjectsBelowArrow="{" ".join(str(x) for x in meta["below_ids"])}"')
        scheme_parts.append(f'<step {" ".join(attrs)}/>')
    scheme_parts.append('</scheme>')
    xml_parts.append("\n".join(scheme_parts))
    return "\n".join(xml_parts), lowest_y


def _layout_divergent_stacked(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """Fallback: stack divergent steps vertically when no shared substrate."""
    xml_parts: List[str] = []
    arrow_y = 300.0
    lowest_y = arrow_y

    for step in scheme.steps:
        sub_scheme = SchemeDescriptor(
            structures=scheme.structures,
            steps=[step],
            layout="linear",
        )
        row_xml, row_lowest = _layout_linear(sub_scheme, ids, arrow_y=arrow_y,
                                             source_data=source_data,
                                             reference_mol=reference_mol)
        xml_parts.append(row_xml)
        lowest_y = row_lowest
        arrow_y = row_lowest + _DIVERGENT_BRANCH_GAP

    return "\n".join(xml_parts), lowest_y


# ---------------------------------------------------------------------------
# Layout engine: stacked-rows (multiple independent sub-schemes)
# ---------------------------------------------------------------------------

SECTION_GAP = 3.0 * ACS_BOND_LENGTH  # white space between independent sections
SECTION_LABEL_X = 40.0               # left margin for section labels


def _layout_stacked_rows(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
    section_gap: float = SECTION_GAP,
) -> Tuple[str, float]:
    """
    Layout multiple independent sub-schemes stacked vertically.

    Each section is rendered as its own row (linear or sequential).
    Section labels like "(i)", "(ii)" are placed at the left margin.

    Returns (inner_xml, lowest_y).
    """
    sections = scheme.sections
    if not sections:
        # Fall back to linear if no sections defined
        return _layout_linear(scheme, ids, source_data=source_data,
                              reference_mol=reference_mol)

    # --- Partition run_arrows by section ---
    # Global step numbers are 1-indexed across all sections.  Build a map
    # from global step number → (section_index, local_step_number).
    global_step = 1
    sec_run_arrows: Dict[int, List[StepRunArrows]] = {}
    for sec_idx, sec in enumerate(sections):
        for local_step in range(1, len(sec.steps) + 1):
            # Map this global step to this section
            for sra in scheme.run_arrows:
                if sra.step == global_step:
                    sec_run_arrows.setdefault(sec_idx, []).append(
                        StepRunArrows(step=local_step, runs=sra.runs))
            global_step += 1

    xml_parts: List[str] = []
    arrow_y = 300.0
    lowest_y = arrow_y
    previous_bottom: Optional[float] = None
    # Content starts to the right of section labels
    content_start_x = 100.0

    for sec_idx, sec in enumerate(sections):
        # --- Render this section's steps ---
        # Create a temporary sub-scheme sharing the structure definitions
        sub_scheme = SchemeDescriptor(
            structures=scheme.structures,
            steps=sec.steps,
            layout=sec.layout or "linear",
            gap_scale=scheme.gap_scale,
        )
        sec_ref_mol = _product_mol_for_steps(
            sec.steps, scheme.structures, source_data,
        )
        if sec_ref_mol is None:
            sec_ref_mol = reference_mol  # fallback to global
        section_ra = sec_run_arrows.get(sec_idx, [])

        # A section wider than the budget wraps like a sequential scheme; run
        # arrows are numbered within the section, so they keep it on one row.
        section_rows = [sec.steps]
        if scheme.max_row_width and len(sec.steps) > 1 and not section_ra:
            per_row = _auto_steps_per_row(
                sub_scheme, scheme.max_row_width, source_data, sec_ref_mol)
            section_rows = _split_into_rows(sec.steps, per_row)

        def layout(at_y: float, row_ids: _IDGen, row_steps=None, sub_scheme=sub_scheme,
                   section_ra=section_ra, sec_ref_mol=sec_ref_mol):
            return _layout_steps_row(
                sub_scheme, row_steps, row_ids,
                start_x=content_start_x, arrow_y=at_y,
                source_data=source_data,
                run_arrows=section_ra,
                reference_mol=sec_ref_mol,
            )

        arrow_y = _next_arrow_y(partial(layout, row_steps=section_rows[0]), arrow_y,
                                previous_bottom, section_gap)

        # --- Section label ---
        if sec.label:
            # Place label at the left margin, vertically centered on arrow_y
            # Baseline is below center, so use arrow_y + small offset
            label_xml, _ = _build_text_element(
                [sec.label], SECTION_LABEL_X, arrow_y + 4.0, ids,
                justification="Left", use_formatting=False,
            )
            xml_parts.append(label_xml)

        for row_index, row_steps in enumerate(section_rows):
            if row_index:
                arrow_y = _next_arrow_y(partial(layout, row_steps=row_steps), arrow_y,
                                        previous_bottom, LAYOUT_ROW_GAP)
            row_xml, row_lowest, _, _ = layout(arrow_y, ids, row_steps=row_steps)
            xml_parts.append(row_xml)
            previous_bottom = _drawn_extent(row_xml)[1]
            lowest_y = max(row_lowest, previous_bottom)

    return "\n".join(xml_parts), lowest_y


# ---------------------------------------------------------------------------
# Core row layout: positions steps in a single row (LTR or RTL)
# ---------------------------------------------------------------------------

STACKED_PLUS_GAP = 16.0   # vertical room for the "+" between stacked substrates
LABEL_ROOM = 20.0         # a compound number below a structure: baseline 14 pt + descent


def _place_stacked_substrates(
    subs: List[ResolvedFragment],
    cursor_x: float,
    arrow_y: float,
    inter_gap: float,
    is_rtl: bool,
    ids: _IDGen,
    xml_parts: List[str],
    step_meta: Dict[str, Any],
    all_frag_ids: List[int],
) -> Tuple[float, float, float]:
    """Draw a step's substrates in one column, a "+" between them, centred on the arrow.

    Returns the cursor after the column, the column's bottom and its centre x.
    """
    col_w = max(_bbox_width(s.bbox) for s in subs)
    heights = [_bbox_height(s.bbox) + (LABEL_ROOM if s.ref.label else 0.0) for s in subs]
    total = sum(heights) + STACKED_PLUS_GAP * (len(subs) - 1)
    col_cx = cursor_x - col_w / 2.0 if is_rtl else cursor_x + col_w / 2.0
    top = arrow_y - total / 2.0
    for k, sub in enumerate(subs):
        if k:
            plus_xml, _ = _build_text_element(['+'], col_cx, top + STACKED_PLUS_GAP / 2.0 + 3.5, ids)
            xml_parts.append(plus_xml)
            top += STACKED_PLUS_GAP
        _shift_atoms(sub.atoms, col_cx - sub.cx, top + _bbox_height(sub.bbox) / 2.0 - sub.cy)
        sub.bbox = _fragment_bbox(sub.atoms)
        sub.cx, sub.cy = _bbox_center(sub.bbox)
        frag_xml, _, frag_id = _build_fragment(sub.atoms, sub.bonds, ids)
        sub.xml, sub.frag_id = frag_xml, frag_id
        xml_parts.append(frag_xml)
        step_meta["reactant_ids"].append(frag_id)
        all_frag_ids.append(frag_id)
        if sub.ref.label:
            lbl_xml, _ = _build_label_element(sub.ref.label, sub.cx, sub.bbox[3] + 14.0, ids)
            xml_parts.append(lbl_xml)
        top += heights[k]
    next_x = col_cx - col_w / 2.0 - inter_gap if is_rtl else col_cx + col_w / 2.0 + inter_gap
    return next_x, top, col_cx


def _layout_steps_row(
    scheme: SchemeDescriptor,
    steps: List[StepDescriptor],
    ids: _IDGen,
    start_x: float = 100.0,
    arrow_y: float = 300.0,
    source_data: Optional[Dict[str, Dict]] = None,
    run_arrows: Optional[List[StepRunArrows]] = None,
    direction: str = "ltr",
    skip_first_substrate_id: Optional[int] = None,
    skip_first_substrate_cursor_x: Optional[float] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float, float, Dict]:
    """
    Layout multiple steps in a single row.

    Parameters
    ----------
    run_arrows : optional
        Run arrows for this row (step numbers are 1-indexed within this row).
        If None, uses scheme.run_arrows with original step numbers.
    direction : str
        "ltr" for left-to-right (default), "rtl" for right-to-left.
        For RTL, start_x is the RIGHT edge; cursor moves leftward.
        Arrows point left, substrates on right, products on left.
    skip_first_substrate_id : optional int
        If set, the first step's substrate is already placed (e.g. by a
        vertical arrow in serpentine mode). This frag_id is used in the
        scheme metadata but the structure is not drawn or positioned.
    skip_first_substrate_cursor_x : optional float
        When skip_first_substrate_id is set, the cursor position after the
        pre-placed substrate. For LTR this is the substrate's right_x + gap;
        for RTL this is the substrate's left_x - gap.

    Returns (xml_string, lowest_y, row_edge_x, row_info) where:
        lowest_y = bottom extent of this row including run arrows
        row_edge_x = rightmost x for LTR, leftmost x for RTL
        row_info = dict with extra metadata:
            "last_product_cx" : float — center x of the last product placed
            "last_product_bottom" : float — bottom y extent of last product + label
            "first_substrate_cx" : float — center x of the first substrate
    """
    is_rtl = (direction == "rtl")
    bond_len = ACS_BOND_LENGTH
    frag_gap = bond_len * LAYOUT_FRAG_GAP_BONDS * scheme.gap_scale    # gap between fragment and arrow
    inter_gap = max(bond_len * LAYOUT_INTER_GAP_BONDS, LAYOUT_PLUS_GAP) * scheme.gap_scale   # gap between adjacent fragments

    # --- Phase 1: Resolve all structures ---
    def _get_ref(sid: str) -> StructureRef:
        """Get StructureRef from declared structures or create a bare one for source lookup."""
        if sid in scheme.structures:
            return scheme.structures[sid]
        # Not declared — create a bare ref (will be resolved via source_data)
        return StructureRef(id=sid)

    resolved_steps: List[ResolvedStep] = []
    for step in steps:
        subs = [_resolve_structure(_get_ref(sid), source_data=source_data,
                                   reference_mol=reference_mol)
                for sid in step.substrates]
        prods = [_resolve_structure(_get_ref(pid), source_data=source_data,
                                    reference_mol=reference_mol)
                 for pid in step.products]

        above_structs = []
        above_text = []
        below_text = []
        below_structs = []

        if step.above_arrow:
            for sid in step.above_arrow.structures:
                above_structs.append(
                    _resolve_structure(_get_ref(sid), source_data=source_data,
                                      reference_mol=reference_mol))
            above_text = [normalize_condition_text(t) for t in step.above_arrow.text]
        if step.below_arrow:
            for sid in step.below_arrow.structures:
                below_structs.append(
                    _resolve_structure(_get_ref(sid), source_data=source_data,
                                      reference_mol=reference_mol))
            below_text = [normalize_condition_text(t) for t in step.below_arrow.text]

        # Add yield to below-arrow text if present
        if step.yield_:
            below_text.append(step.yield_)

        resolved_steps.append(ResolvedStep(
            descriptor=step,
            substrates=subs,
            products=prods,
            above_structures=above_structs,
            above_text=above_text,
            below_text=below_text,
            below_structures=below_structs,
        ))

    # --- Phase 2: Compute arrow lengths from content ---
    min_arrow_len = MIN_ARROW_BONDS * bond_len
    use_letter_conditions = bool(scheme.condition_key)

    for rs in resolved_steps:
        # Check if this step uses letter conditions
        is_letter = (use_letter_conditions
                     and _is_letter_condition(rs.below_text))
        rs._is_letter_cond = is_letter

        # Width of above-arrow content
        above_width = 0.0
        for af in rs.above_structures:
            above_width += _bbox_width(af.bbox) + inter_gap
        above_text_w = _estimate_text_width(rs.above_text)
        above_width = max(above_width, above_text_w)

        # Width of below-arrow content (skip if using letter label)
        below_width = 0.0
        if not is_letter:
            for bf in rs.below_structures:
                below_width += _bbox_width(bf.bbox) + inter_gap
            below_text_w = _estimate_text_width(rs.below_text)
            below_width = max(below_width, below_text_w)

        content_width = max(above_width, below_width)
        rs._arrow_len = max(content_width + 10.0, min_arrow_len)

    # --- Phase 3: Position everything ---
    # For LTR: cursor_x = left edge of next element, moves right.
    # For RTL: cursor_x = right edge of next element, moves left.
    cursor_x = start_x

    xml_parts: List[str] = []
    all_frag_ids: List[int] = []
    step_metadata: List[Dict] = []  # for <scheme><step> elements

    # Track extra info for serpentine/multi-row callers
    _last_product_cx: float = start_x
    _last_product_bottom: float = arrow_y
    _first_substrate_cx: float = start_x
    _stack_bottom: float = arrow_y
    _last_product_id: Optional[int] = None

    for step_idx, rs in enumerate(resolved_steps):
        step_meta: Dict[str, Any] = {
            "reactant_ids": [],
            "product_ids": [],
            "arrow_id": 0,
            "above_ids": [],
            "below_ids": [],
        }

        # -- Substrates --
        # Skip first substrate if: (a) pre-placed by vertical arrow (serpentine), or
        # (b) shared intermediate from the previous step's product.
        start_sub_idx = 0
        if step_idx == 0 and skip_first_substrate_id is not None:
            start_sub_idx = len(rs.substrates)
            step_meta["reactant_ids"].append(skip_first_substrate_id)
            if skip_first_substrate_cursor_x is not None:
                cursor_x = skip_first_substrate_cursor_x
        elif step_idx > 0:
            prev_products = resolved_steps[step_idx - 1].descriptor.products
            curr_substrates = rs.descriptor.substrates
            if (len(curr_substrates) >= 1 and len(prev_products) >= 1
                    and curr_substrates[0] == prev_products[-1]):
                start_sub_idx = 1
                # Use the previous step's product fragment ID
                prev_meta = step_metadata[step_idx - 1]
                if prev_meta["product_ids"]:
                    step_meta["reactant_ids"].append(prev_meta["product_ids"][-1])

        stacked = (getattr(scheme, "substrate_layout", "inline") == "stacked_plus"
                   and len(rs.substrates) - start_sub_idx >= 2)
        if stacked:
            cursor_x, column_bottom, column_cx = _place_stacked_substrates(
                rs.substrates[start_sub_idx:], cursor_x, arrow_y, inter_gap, is_rtl, ids,
                xml_parts, step_meta, all_frag_ids)
            _stack_bottom = max(_stack_bottom, column_bottom)
            if step_idx == 0 and start_sub_idx == 0:
                _first_substrate_cx = column_cx
        for i in range(len(rs.substrates) if stacked else start_sub_idx, len(rs.substrates)):
            sub = rs.substrates[i]
            if i:
                plus_x = cursor_x + inter_gap / 2 if is_rtl else cursor_x - inter_gap / 2
                plus_xml, _ = _build_text_element(['+'], plus_x, arrow_y + 3.5, ids)
                xml_parts.append(plus_xml)
            bbox = sub.bbox
            w = _bbox_width(bbox)
            # Shift to position
            if is_rtl:
                cx_target = cursor_x - w / 2.0
            else:
                cx_target = cursor_x + w / 2.0
            cy_target = arrow_y
            dx = cx_target - sub.cx
            dy = cy_target - sub.cy
            _shift_atoms(sub.atoms, dx, dy)
            sub.bbox = _fragment_bbox(sub.atoms)
            sub.cx, sub.cy = _bbox_center(sub.bbox)

            frag_xml, _, frag_id = _build_fragment(sub.atoms, sub.bonds, ids)
            sub.xml = frag_xml
            sub.frag_id = frag_id
            xml_parts.append(frag_xml)
            step_meta["reactant_ids"].append(frag_id)
            all_frag_ids.append(frag_id)

            # Track first substrate center
            if step_idx == 0 and i == 0:
                _first_substrate_cx = sub.cx

            # Compound label below structure
            if sub.ref.label:
                label_x = sub.cx
                label_y = sub.bbox[3] + 14.0  # below structure
                lbl_xml, lbl_id = _build_label_element(
                    sub.ref.label, label_x, label_y, ids,
                )
                xml_parts.append(lbl_xml)

            if is_rtl:
                cursor_x = sub.bbox[0] - inter_gap
            else:
                cursor_x = sub.bbox[2] + inter_gap

        # -- Arrow --
        if is_rtl:
            # RTL: tail on right (near substrate), head on left (near product)
            arrow_tail_x = cursor_x - frag_gap * 0.3
            arrow_head_x = arrow_tail_x - rs._arrow_len
        else:
            arrow_tail_x = cursor_x + frag_gap * 0.3
            arrow_head_x = arrow_tail_x + rs._arrow_len
        arrow_mid_x = (arrow_tail_x + arrow_head_x) / 2.0

        dashed = (rs.descriptor.arrow_style == "dashed")
        failed = (rs.descriptor.arrow_style == "failed")
        arrow_xml, arrow_id = _build_arrow(
            arrow_tail_x, arrow_y, arrow_head_x, arrow_y, ids,
            dashed=dashed, nogo=failed,
            kind=_arrow_kind(rs.descriptor.arrow_style),
        )
        xml_parts.append(arrow_xml)
        step_meta["arrow_id"] = arrow_id

        rs.arrow_tail_x = arrow_tail_x
        rs.arrow_tail_y = arrow_y
        rs.arrow_head_x = arrow_head_x
        rs.arrow_head_y = arrow_y

        # -- Above-arrow content --
        # Structures above arrow
        above_cursor_x = arrow_mid_x
        if rs.above_structures:
            total_above_w = sum(_bbox_width(af.bbox) for af in rs.above_structures)
            total_above_w += inter_gap * max(0, len(rs.above_structures) - 1)
            above_cursor_x = arrow_mid_x - total_above_w / 2.0

            # If there's text below above-arrow structures (e.g. "(1.2 eq)"),
            # push structures higher to make room for the text between
            # structure bottom and arrow.
            # Gap layout: struct -(4pt)- text -(2pt)- arrow
            above_text_height = 0.0
            if rs.above_text:
                n_abt = len(rs.above_text)
                # Visual text height + struct-to-text gap (4pt)
                above_text_height = (
                    _CAP_HEIGHT + max(0, n_abt - 1) * _LINE_ADVANCE + _DESCENT + 4.0
                )

            for af in rs.above_structures:
                af_w = _bbox_width(af.bbox)
                af_h = _bbox_height(af.bbox)
                label_height = 20.0 if af.ref.label else 0.0
                # Position above the arrow, with extra room for text below structure
                target_cx = above_cursor_x + af_w / 2.0
                if rs.above_text:
                    # text sits LAYOUT_BELOW_GAP above arrow; struct sits above text
                    target_cy = arrow_y - LAYOUT_BELOW_GAP - above_text_height - label_height - af_h / 2.0
                else:
                    # no text: struct sits LAYOUT_ABOVE_GAP above arrow
                    target_cy = arrow_y - LAYOUT_ABOVE_GAP - label_height - af_h / 2.0
                dx = target_cx - af.cx
                dy = target_cy - af.cy
                _shift_atoms(af.atoms, dx, dy)
                af.bbox = _fragment_bbox(af.atoms)
                af.cx, af.cy = _bbox_center(af.bbox)

                frag_xml, _, frag_id = _build_fragment(af.atoms, af.bonds, ids)
                af.xml = frag_xml
                af.frag_id = frag_id
                xml_parts.append(frag_xml)
                step_meta["above_ids"].append(frag_id)
                all_frag_ids.append(frag_id)

                # Label for above-arrow structure
                if af.ref.label:
                    lbl_xml, _ = _build_label_element(
                        af.ref.label, af.cx, af.bbox[3] + 14.0, ids,
                    )
                    xml_parts.append(lbl_xml)

                above_cursor_x += af_w + inter_gap

        # Above-arrow text (equiv text for above structures, or condition text)
        # p.y is the BASELINE of the first line.
        # Position so that the visual bottom (last baseline + descent) sits
        # LAYOUT_BELOW_GAP above the arrow (same gap as below-arrow content).
        if rs.above_text:
            n_abt = len(rs.above_text)
            # Distance from first baseline to visual bottom
            text_below_baseline = max(0, n_abt - 1) * _LINE_ADVANCE + _DESCENT
            text_y = arrow_y - LAYOUT_BELOW_GAP - text_below_baseline
            if rs.above_structures:
                # Text sits between above-arrow structures and the arrow.
                # Structures were already pushed higher to make room.
                pass  # text_y is already correct
            txt_xml, txt_id = _build_text_element(
                rs.above_text, arrow_mid_x, text_y, ids,
                use_formatting=True,
            )
            xml_parts.append(txt_xml)
            step_meta["above_ids"].append(txt_id)

        # -- Below-arrow content --
        # p.y is the BASELINE of the first line.
        # Visual top = p.y - _CAP_HEIGHT, so set p.y so that
        # (p.y - _CAP_HEIGHT) = arrow_y + LAYOUT_BELOW_GAP.
        below_y = arrow_y + LAYOUT_BELOW_GAP + _CAP_HEIGHT

        # Below-arrow structures first
        if rs.below_structures:
            below_struct_y = arrow_y + LAYOUT_BELOW_GAP
            total_below_w = sum(_bbox_width(bf.bbox) for bf in rs.below_structures)
            total_below_w += inter_gap * max(0, len(rs.below_structures) - 1)
            below_cursor_x = arrow_mid_x - total_below_w / 2.0

            for bf in rs.below_structures:
                bf_w = _bbox_width(bf.bbox)
                bf_h = _bbox_height(bf.bbox)
                target_cx = below_cursor_x + bf_w / 2.0
                target_cy = below_struct_y + bf_h / 2.0 + 4.0
                dx = target_cx - bf.cx
                dy = target_cy - bf.cy
                _shift_atoms(bf.atoms, dx, dy)
                bf.bbox = _fragment_bbox(bf.atoms)
                bf.cx, bf.cy = _bbox_center(bf.bbox)

                frag_xml, _, frag_id = _build_fragment(bf.atoms, bf.bonds, ids)
                bf.xml = frag_xml
                bf.frag_id = frag_id
                xml_parts.append(frag_xml)
                step_meta["below_ids"].append(frag_id)
                all_frag_ids.append(frag_id)

                if bf.ref.label:
                    lbl_xml, _ = _build_label_element(
                        bf.ref.label, bf.cx, bf.bbox[3] + 14.0, ids,
                    )
                    xml_parts.append(lbl_xml)

                below_cursor_x += bf_w + inter_gap

            # Adjust text Y to be below the structures
            below_y = max(bf.bbox[3] for bf in rs.below_structures) + 14.0

        # Below-arrow text (or letter label for condition key mode)
        if rs.below_text:
            if rs._is_letter_cond:
                # Letter condition: render small italic label above arrow
                letter_text = rs.below_text[0].strip()
                lbl_xml, lbl_id = _build_letter_label(
                    letter_text, arrow_mid_x, arrow_y, ids,
                )
                xml_parts.append(lbl_xml)
                step_meta["above_ids"].append(lbl_id)
            else:
                txt_xml, txt_id = _build_text_element(
                    rs.below_text, arrow_mid_x, below_y, ids,
                )
                xml_parts.append(txt_xml)
                step_meta["below_ids"].append(txt_id)

        if is_rtl:
            cursor_x = arrow_head_x - frag_gap * 1.0
        else:
            cursor_x = arrow_head_x + frag_gap * 1.0

        # -- Products --
        for i, prod in enumerate(rs.products):
            if i:
                plus_x = cursor_x + inter_gap / 2 if is_rtl else cursor_x - inter_gap / 2
                plus_xml, _ = _build_text_element(['+'], plus_x, arrow_y + 3.5, ids)
                xml_parts.append(plus_xml)
            bbox = prod.bbox
            w = _bbox_width(bbox)
            if is_rtl:
                cx_target = cursor_x - w / 2.0
            else:
                cx_target = cursor_x + w / 2.0
            cy_target = arrow_y
            dx = cx_target - prod.cx
            dy = cy_target - prod.cy
            _shift_atoms(prod.atoms, dx, dy)
            prod.bbox = _fragment_bbox(prod.atoms)
            prod.cx, prod.cy = _bbox_center(prod.bbox)

            frag_xml, _, frag_id = _build_fragment(prod.atoms, prod.bonds, ids)
            prod.xml = frag_xml
            prod.frag_id = frag_id
            xml_parts.append(frag_xml)
            step_meta["product_ids"].append(frag_id)
            all_frag_ids.append(frag_id)

            prod_bottom = prod.bbox[3]
            if prod.ref.label:
                lbl_xml, _ = _build_label_element(
                    prod.ref.label, prod.cx, prod.bbox[3] + 14.0, ids,
                )
                xml_parts.append(lbl_xml)
                prod_bottom = prod.bbox[3] + 14.0 + 6.0  # label baseline + descent

            # Track last product info for serpentine callers
            _last_product_cx = prod.cx
            _last_product_bottom = prod_bottom
            _last_product_id = getattr(prod, "frag_id", None)

            if is_rtl:
                cursor_x = prod.bbox[0] - inter_gap
            else:
                cursor_x = prod.bbox[2] + inter_gap

        step_metadata.append(step_meta)

    # --- Phase 4: Build <scheme> with <step> elements ---
    scheme_id = ids.next()
    scheme_parts = [f'<scheme id="{scheme_id}">']
    for meta in step_metadata:
        step_id = ids.next()
        attrs = [f'id="{step_id}"']
        if meta["reactant_ids"]:
            attrs.append(f'ReactionStepReactants="{" ".join(str(x) for x in meta["reactant_ids"])}"')
        if meta["product_ids"]:
            attrs.append(f'ReactionStepProducts="{" ".join(str(x) for x in meta["product_ids"])}"')
        attrs.append(f'ReactionStepArrows="{meta["arrow_id"]}"')
        if meta["above_ids"]:
            attrs.append(f'ReactionStepObjectsAboveArrow="{" ".join(str(x) for x in meta["above_ids"])}"')
        if meta["below_ids"]:
            attrs.append(f'ReactionStepObjectsBelowArrow="{" ".join(str(x) for x in meta["below_ids"])}"')
        scheme_parts.append(f'<step {" ".join(attrs)}/>')
    scheme_parts.append('</scheme>')
    xml_parts.append("\n".join(scheme_parts))

    # --- Phase 5: Run arrows ---
    # Compute lowest_y from row content (always needed for return value)
    lowest_y = arrow_y + 60.0  # baseline estimate
    for rs in resolved_steps:
        if rs.below_text:
            n_blt = len(rs.below_text)
            # below_y baseline + distance to visual bottom + margin
            below_vis_bottom = (
                arrow_y + LAYOUT_BELOW_GAP + _CAP_HEIGHT  # first baseline
                + max(0, n_blt - 1) * _LINE_ADVANCE + _DESCENT  # to visual bottom
                + 4.0  # margin
            )
            lowest_y = max(lowest_y, below_vis_bottom)
        for bf in rs.below_structures:
            lowest_y = max(lowest_y, bf.bbox[3] + 20.0)
    lowest_y = max(lowest_y, _stack_bottom + 4.0)

    # Determine which run arrows to render
    effective_run_arrows = run_arrows if run_arrows is not None else scheme.run_arrows

    if effective_run_arrows:
        base_y = lowest_y + 6.0  # tighter gap between conditions text and run arrows
        lowest_y = base_y
        # Each run takes the highest free line below the conditions; runs of
        # neighbouring steps share a line unless their labels would meet.
        taken: List[Tuple[float, float, float, float]] = []

        for sra in effective_run_arrows:
            step_idx = sra.step - 1  # 0-indexed
            if step_idx < 0 or step_idx >= len(resolved_steps):
                continue
            rs = resolved_steps[step_idx]

            # Run arrow matches the reaction arrow exactly (same tail/head X),
            # just translated vertically below the scheme content.
            run_tail_x = rs.arrow_tail_x
            run_head_x = rs.arrow_head_x
            run_y = base_y

            for run_entry in sra.runs:
                note_height = 10.0 if run_entry.note else 0.0
                left = run_tail_x - 5.0 - _estimate_text_width([run_entry.input_label])
                right = run_head_x + 5.0 + _estimate_text_width([run_entry.output_label or ""])
                if run_entry.note:
                    half = _estimate_text_width([run_entry.note]) / 2.0
                    mid = (run_tail_x + run_head_x) / 2.0
                    left, right = min(left, mid - half), max(right, mid + half)
                while any(left < x2 + RUN_LABEL_GAP and x1 < right + RUN_LABEL_GAP
                          and run_y < y2 and y1 < run_y + note_height + 18.0
                          for x1, y1, x2, y2 in taken):
                    run_y += 18.0
                taken.append((left, run_y, right, run_y + note_height + 18.0))
                if run_entry.note:
                    # Text centered above this specific run arrow
                    note_y = run_y - 1.0
                    note_xml, _ = _build_text_element(
                        [run_entry.note],
                        (run_tail_x + run_head_x) / 2.0,
                        note_y, ids,
                        justification="Center",
                        use_formatting=False,
                    )
                    xml_parts.append(note_xml)
                    run_y += note_height  # space for note text
                run_xml, _ = _build_run_arrow(
                    run_tail_x, run_head_x,
                    run_y,
                    run_entry.input_label,
                    run_entry.output_label,
                    ids,
                )
                xml_parts.append(run_xml)
                run_y += 18.0  # stack multiple runs
                lowest_y = max(lowest_y, run_y)

    # Compute row edge (rightmost for LTR, leftmost for RTL)
    row_edge_x = cursor_x

    row_info = {
        "last_product_cx": _last_product_cx,
        "last_product_bottom": _last_product_bottom,
        "first_substrate_cx": _first_substrate_cx,
        "last_product_id": _last_product_id,
    }
    return "\n".join(xml_parts), lowest_y, row_edge_x, row_info


# ---------------------------------------------------------------------------
# Document assembly
# ---------------------------------------------------------------------------

def _format_header(bbox: str) -> str:
    """Format the full CDXML header with ACS style."""
    return CDXML_HEADER.format(
        bbox=bbox,
        label_font=ACS_LABEL_FONT,
        label_size=ACS_LABEL_SIZE,
        label_face=ACS_LABEL_FACE,
        caption_size=ACS_CAPTION_SIZE,
        hash_spacing=ACS_HASH_SPACING,
        margin_width=ACS_MARGIN_WIDTH,
        line_width=ACS_LINE_WIDTH,
        bold_width=ACS_BOLD_WIDTH,
        bond_length=ACS_BOND_LENGTH_STR,
        bond_spacing=ACS_BOND_SPACING,
        chain_angle=ACS_CHAIN_ANGLE_STR,
    )


_PAGE_OPEN = (
    '<page id="{page_id}" BoundingBox="0 0 1620 2160" '
    'HeaderPosition="36" FooterPosition="36" '
    'PrintTrimMarks="yes" HeightPages="3" WidthPages="3">'
)
_PAGE_CLOSE = "</page>"


def _product_mol_for_steps(
    steps: List[StepDescriptor],
    structures: Dict[str, StructureRef],
    source_data: Optional[Dict[str, Dict]] = None,
    pick: str = "last",
) -> Optional[Any]:
    """
    Build an RDKit Mol with 2D coords for the product of given steps.

    Parameters
    ----------
    steps : list of StepDescriptor
    structures : dict mapping structure IDs to StructureRef
    source_data : optional reaction_parser JSON species lookup
    pick : "last" (default) = last step's last product;
           "first" = first step's first product (for divergent).

    Returns an RDKit Mol with 2D conformer, or None.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None

    pid: Optional[str] = None
    if pick == "first":
        if steps and steps[0].products:
            pid = steps[0].products[0]
    else:
        if steps and steps[-1].products:
            pid = steps[-1].products[-1]

    if pid is None:
        return None

    product_smiles: Optional[str] = None
    ref = structures.get(pid)
    if ref and ref.smiles:
        product_smiles = ref.smiles
    elif source_data:
        sp = source_data.get(pid) or source_data.get(pid.lower())
        if sp:
            product_smiles = sp.get("smiles")

    if not product_smiles:
        return None

    mol = Chem.MolFromSmiles(product_smiles)
    if mol is None:
        return None

    AllChem.Compute2DCoords(mol)
    return mol


def _identify_product_mol(
    scheme: SchemeDescriptor,
    source_data: Optional[Dict[str, Dict]] = None,
) -> Optional[Any]:
    """
    Identify the main product and create an RDKit Mol with 2D coords.

    This mol serves as the reference for MCS alignment of all other
    structures in the scheme — shared scaffolds are oriented to match
    the product's layout, producing visually consistent schemes.

    Reference product selection:
      - linear / sequential: product of the LAST step
      - divergent: product of the FIRST (horizontal) step
      - stacked-rows: product of the first section's last step

    Returns an RDKit Mol with 2D conformer, or None if RDKit is
    unavailable or no product SMILES can be resolved.
    """
    if scheme.layout == "divergent":
        return _product_mol_for_steps(
            scheme.steps, scheme.structures, source_data, pick="first",
        )
    elif scheme.layout == "stacked-rows":
        if scheme.sections:
            return _product_mol_for_steps(
                scheme.sections[0].steps, scheme.structures, source_data,
            )
        return None
    else:
        # linear / sequential — last step's last product
        return _product_mol_for_steps(
            scheme.steps, scheme.structures, source_data,
        )


def _layout_convergent(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """Layout a convergent multi-branch route as stacked arm rows merging into the final row."""
    rows = SchemeDescriptor(
        structures=scheme.structures,
        layout="stacked-rows",
        wrap=scheme.wrap,
        max_row_width=scheme.max_row_width,
        gap_scale=scheme.gap_scale,
        run_arrows=scheme.run_arrows,
        sections=[SectionDescriptor(steps=steps) for steps in convergent_sections(scheme.steps)],
    )
    return _layout_stacked_rows(
        rows,
        ids,
        source_data=source_data,
        reference_mol=reference_mol,
        section_gap=LAYOUT_ROW_GAP,
    )


def _boxes_overlap(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> bool:
    return min(a[2], b[2]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[1], b[1])


def _ray_box_exit_t(
    px: float,
    py: float,
    dx: float,
    dy: float,
    box: Tuple[float, float, float, float],
    pad: float,
) -> float:
    """Return smallest positive t in (0, 1) where P(t) = (px + t*dx, py + t*dy) exits padded *box*."""
    bx0, by0, bx1, by1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    tx = ((bx1 - px) / dx if dx > 0 else (bx0 - px) / dx) if abs(dx) > 1e-6 else float("inf")
    ty = ((by1 - py) / dy if dy > 0 else (by0 - py) / dy) if abs(dy) > 1e-6 else float("inf")
    t = min(tx, ty)
    if not math.isfinite(t):
        return 0.25
    return max(0.08, min(0.45, t))


def _ray_box_exit_unclamped(
    px: float,
    py: float,
    dx: float,
    dy: float,
    box: Tuple[float, float, float, float],
    pad: float,
) -> float:
    """Return the true ``t`` in (0, 1] where ``P(t)`` leaves the padded *box*.

    Unlike :func:`_ray_box_exit_t` this does not clamp to ``[0.08, 0.45]``, so the
    result really sits on the box edge — needed to start a connector arrow at a
    structure's border instead of inside it.
    """
    bx0, by0, bx1, by1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    tx = ((bx1 - px) / dx if dx > 0 else (bx0 - px) / dx) if abs(dx) > 1e-6 else float("inf")
    ty = ((by1 - py) / dy if dy > 0 else (by0 - py) / dy) if abs(dy) > 1e-6 else float("inf")
    t = min(tx, ty)
    if not math.isfinite(t) or t <= 0:
        return 0.0
    return min(1.0, t)


def _text_box_at(lines: List[str], x: float, y: float) -> Tuple[float, float, float, float]:
    max_chars = max((len(ln) for ln in lines), default=1)
    w = max_chars * _CHAR_WIDTH
    top = y - _CAP_HEIGHT
    bottom = y + max(0, len(lines) - 1) * _LINE_ADVANCE + _DESCENT
    return (x - w / 2.0, top, x + w / 2.0, bottom)


def _overlay_box(atoms: List[Dict], overlays: Any) -> Optional[Tuple[float, float, float, float]]:
    """Box of the dashed outlines a brief will draw over this structure, or None.

    ``overlays`` are [{atoms, shape}] with 0-based SMILES atom indices; atom dicts carry the
    1-based ``index`` of the same order. The box holds the outline as it will be drawn.
    """
    by_index = {a.get("index"): a for a in atoms}
    box: Optional[Tuple[float, float, float, float]] = None
    for item in overlays or ():
        if item.get("shape") == "ring":
            # A generic ring: an arc through two bonded atoms, away from the rest of the molecule.
            ends = [(float(by_index[i + 1]["x"]), float(by_index[i + 1]["y"]))
                    for i in item.get("atoms") or () if i + 1 in by_index]
            rest = [(float(by_index[i + 1]["x"]), float(by_index[i + 1]["y"]))
                    for i in item.get("away") or () if i + 1 in by_index]
            if len(ends) == 2:
                away = (sum(p[0] for p in rest) / len(rest), sum(p[1] for p in rest) / len(rest)) if rest else None
                part = generic_ring_box(ends[0], ends[1], away, item.get("members", 6))
                box = part if box is None else _union_box(box, part)
            continue
        points = [(float(by_index[i + 1]["x"]), float(by_index[i + 1]["y"]))
                  for i in item.get("atoms") or () if i + 1 in by_index]
        if not points:
            continue
        part = overlay_box(list(dict.fromkeys(points)), shape=item.get("shape", "circle"))
        box = part if box is None else _union_box(box, part)
    return box


def _union_box(a: Tuple[float, float, float, float], b: Optional[Tuple[float, float, float, float]]
               ) -> Tuple[float, float, float, float]:
    if b is None:
        return a
    return (min(a[0], b[0]), min(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3]))


def _cycle_arc_span(cx: float, cy: float, radius: float, a0: float, a1: float,
                    box0: Tuple[float, float, float, float], box1: Tuple[float, float, float, float],
                    pad: float, min_length: float) -> Optional[Tuple[float, float]]:
    """Angles of the arc from the species at a0 to the next at a1 (a1 > a0, both on the circle),
    clipped where the circle leaves each species' box and ``pad`` beyond; None if too short."""
    def inside(box, t):
        x, y = cx + radius * math.cos(t), cy + radius * math.sin(t)
        return box[0] - pad < x < box[2] + pad and box[1] - pad < y < box[3] + pad

    step = math.radians(0.25)
    start, end = a0, a1
    while start < a1 and inside(box0, start):
        start += step
    while end > start and inside(box1, end):
        end -= step
    if (end - start) * radius < min_length:
        return None
    return start, end


def _build_cycle_arc(cx: float, cy: float, radius: float, start: float, end: float, ids: _IDGen,
                     dashed: bool = False) -> Tuple[str, int, List[Tuple[float, float]]]:
    """A step of a cycle as ChemDraw's own arc arrow, running clockwise on the page from the
    angle ``start`` to ``end``. Returns (xml, id, points from tail to head).

    ChemDraw runs an arc anticlockwise from Tail3D through AngularSize degrees, so a clockwise
    arrow is written from its far end with the head on the tail end, as figure_mechanism writes
    curved arrows (toolkit/tests/data/native-arc-arrows.cdxml).
    """
    pts = [(cx + radius * math.cos(start + (end - start) * k / 32), cy + radius * math.sin(start + (end - start) * k / 32))
           for k in range(33)]
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    aid, z = ids.next(), ids.next()
    first, last = pts[0], pts[-1]
    attrs = [
        f'id="{aid}"', f'BoundingBox="{min(xs):.2f} {min(ys):.2f} {max(xs):.2f} {max(ys):.2f}"', f'Z="{z}"',
        'FillType="None"', f'LineWidth="{ACS_LINE_WIDTH}"', 'ArrowheadHead="None"', 'ArrowheadTail="Full"',
        'ArrowheadType="Solid"', 'HeadSize="1000"', 'ArrowheadCenterSize="875"', 'ArrowheadWidth="250"',
        f'AngularSize="{math.degrees(end - start):.2f}"',
        f'Head3D="{first[0]:.2f} {first[1]:.2f} 0"', f'Tail3D="{last[0]:.2f} {last[1]:.2f} 0"',
        f'Center3D="{cx:.2f} {cy:.2f} 0"', f'MajorAxisEnd3D="{cx + radius:.2f} {cy:.2f} 0"',
        f'MinorAxisEnd3D="{cx:.2f} {cy - radius:.2f} 0"',
    ]
    if dashed:
        attrs.append('LineType="Dashed"')
    return f'<arrow {" ".join(attrs)}/>', aid, pts


def _along(path: List[Tuple[float, float]], frac: float) -> Tuple[float, float]:
    """The point a fraction of the way along a polyline, by length."""
    lengths = [math.dist(a, b) for a, b in zip(path, path[1:])]
    goal = sum(lengths) * frac
    for (a, b), length in zip(zip(path, path[1:]), lengths):
        if goal <= length or length == lengths[-1] and (a, b) == (path[-2], path[-1]):
            t = 0.0 if length == 0 else min(1.0, goal / length)
            return a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t
        goal -= length
    return path[-1]


def _place_radial_text(
    lines: List[str],
    mx: float,
    my: float,
    dir_x: float,
    dir_y: float,
    placed_boxes: List[Tuple[float, float, float, float]],
    max_dist: Optional[float] = None,
) -> Optional[Tuple[float, float, Tuple[float, float, float, float]]]:
    """Find the closest radial offset along (dir_x, dir_y) from (mx, my) with zero box overlap."""
    limit = max_dist if max_dist is not None else 220.0
    d = 12.0
    while d <= limit:
        cx_t = mx + dir_x * d
        cy_t = my + dir_y * d
        # Shift baseline y so visual vertical center of text block sits at cy_t
        block_h = max(0, len(lines) - 1) * _LINE_ADVANCE + _CAP_HEIGHT + _DESCENT
        baseline_y = cy_t - block_h / 2.0 + _CAP_HEIGHT
        tbox = _text_box_at(lines, cx_t, baseline_y)
        if not any(_boxes_overlap(tbox, ob) for ob in placed_boxes):
            return cx_t, baseline_y, tbox
        d += 2.0
    return None


def _layout_catalytic_cycle(
    scheme: SchemeDescriptor,
    ids: _IDGen,
    source_data: Optional[Dict[str, Dict]] = None,
    reference_mol: Optional[Any] = None,
) -> Tuple[str, float]:
    """
    Layout a closed catalytic cycle around a circle centered at (280.0, 240.0)
    with radius R = max(90.0, n_steps * 32.0), chord reaction arrows between
    consecutive intermediates, and radial condition labels outside/inside each arrow.
    """
    from .route_graph import cycle_order

    if not scheme.steps:
        return "", 240.0

    try:
        order = cycle_order(scheme.steps)
        steps = [scheme.steps[i] for i in order]
    except Exception:
        steps = list(scheme.steps)

    n_steps = len(steps)
    cx, cy = 280.0, 240.0
    gap_scale = max(float(getattr(scheme, "gap_scale", 1.0) or 1.0), 0.5)
    radius = max(90.0, n_steps * 32.0) * gap_scale

    # Prepare one intermediate fragment per step around the ring
    inter_frags: List[_ResolvedFragment] = []
    for k, step in enumerate(steps):
        sub_id = step.substrates[0] if step.substrates else None
        if not sub_id or sub_id not in scheme.structures:
            raise ValueError(f"Catalytic cycle step {k + 1} has no valid substrate structure.")
        frag = _resolve_structure(
            scheme.structures[sub_id],
            center_x=cx,
            center_y=cy,
            source_data=source_data,
            reference_mol=reference_mol,
        )
        if frag is None:
            raise ValueError(f"Could not prepare structure '{sub_id}' in catalytic cycle.")
        inter_frags.append(frag)
    # A dashed overlay drawn over an intermediate later (its fragment marked in every one) is
    # part of what the arrows and captions must clear.
    inter_overlays = [scheme.structures[s.substrates[0]].overlays if getattr(scheme.structures[s.substrates[0]],
                      "overlays", ()) else () for s in steps]
    inter_rel = [_overlay_box(frag.atoms, ov) for frag, ov in zip(inter_frags, inter_overlays)]

    pad = 8.0 * gap_scale
    min_arrow_len = 28.0

    # Ensure radius is large enough that structures do not overlap and chord arrows fit cleanly
    for _ in range(12):
        angles = [-math.pi / 2.0 + 2.0 * math.pi * k / n_steps for k in range(n_steps)]
        centers = [(cx + radius * math.cos(a), cy + radius * math.sin(a)) for a in angles]
        boxes: List[Tuple[float, float, float, float]] = []
        for (ix, iy), frag, rel in zip(centers, inter_frags, inter_rel):
            hw = _bbox_width(frag.bbox) / 2.0
            hh = _bbox_height(frag.bbox) / 2.0
            lbl_extra = (14.0 + _DESCENT) if frag.ref.label else 0.0
            if rel is not None:
                dx_, dy_ = ix - frag.cx, iy - frag.cy
                rel_box = (rel[0] + dx_, rel[1] + dy_, rel[2] + dx_, rel[3] + dy_)
                lbl_extra += max(0.0, rel_box[3] - (iy + hh))
                boxes.append(_union_box((ix - hw, iy - hh, ix + hw, iy + hh + lbl_extra), rel_box))
            else:
                boxes.append((ix - hw, iy - hh, ix + hw, iy + hh + lbl_extra))

        any_struct_overlap = any(
            _boxes_overlap(
                (boxes[i][0] - pad, boxes[i][1] - pad, boxes[i][2] + pad, boxes[i][3] + pad),
                (boxes[j][0] - pad, boxes[j][1] - pad, boxes[j][2] + pad, boxes[j][3] + pad),
            )
            for i in range(n_steps)
            for j in range(i + 1, n_steps)
        )
        arrows_ok = True
        for k in range(n_steps):
            nk = (k + 1) % n_steps
            x0, y0 = centers[k]
            x1, y1 = centers[nk]
            dx, dy = x1 - x0, y1 - y0
            chord = math.hypot(dx, dy)
            t0 = _ray_box_exit_t(x0, y0, dx, dy, boxes[k], pad)
            t1 = 1.0 - _ray_box_exit_t(x1, y1, -dx, -dy, boxes[nk], pad)
            if chord * (t1 - t0) < min_arrow_len:
                arrows_ok = False
                break
        if not any_struct_overlap and arrows_ok:
            break
        radius += 8.0

    xml_parts: List[str] = []
    placed_boxes: List[Tuple[float, float, float, float]] = []
    frag_ids: List[int] = []
    final_boxes: List[Tuple[float, float, float, float]] = []
    angles = [-math.pi / 2.0 + 2.0 * math.pi * k / n_steps for k in range(n_steps)]
    centers = [(cx + radius * math.cos(a), cy + radius * math.sin(a)) for a in angles]
    lowest_y = cy + radius

    # 1. Place intermediate structures around the ring
    for k, ((ix, iy), frag) in enumerate(zip(centers, inter_frags)):
        _shift_atoms(frag.atoms, ix - frag.cx, iy - frag.cy)
        frag.bbox = _fragment_bbox(frag.atoms)
        frag.cx, frag.cy = _bbox_center(frag.bbox)
        frag_xml, _, fid = _build_fragment(frag.atoms, frag.bonds, ids)
        xml_parts.append(frag_xml)
        frag_ids.append(fid)

        ring = _overlay_box(frag.atoms, inter_overlays[k])
        box_y2 = frag.bbox[3]
        if frag.ref.label:
            # The number goes under the overlay when the outline reaches below the structure.
            lbl_y = max(frag.bbox[3], ring[3] if ring else frag.bbox[3]) + 14.0
            lbl_xml, _ = _build_label_element(frag.ref.label, frag.cx, lbl_y, ids)
            xml_parts.append(lbl_xml)
            box_y2 = lbl_y + _DESCENT

        fbox = _union_box((frag.bbox[0] - 2.0, frag.bbox[1] - 2.0, frag.bbox[2] + 2.0, box_y2 + 2.0), ring)
        final_boxes.append(fbox)
        placed_boxes.append(fbox)
        lowest_y = max(lowest_y, box_y2)

    # 2. Place the step arrows between consecutive intermediates: arcs of the circle the species
    # sit on (a cycle reads as a ring, BJOC 22-125 scheme 4), or a straight chord for an
    # equilibrium or failed step, and where the arc left between two species is too short.
    arrow_records: List[Tuple[float, float, float, float, int]] = []
    arrow_paths: List[List[Tuple[float, float]]] = []
    for k, step in enumerate(steps):
        nk = (k + 1) % n_steps
        span = None
        if n_steps >= 2 and step.arrow_style in (None, "solid", "dashed"):
            span = _cycle_arc_span(cx, cy, radius, angles[k], angles[k] + 2.0 * math.pi / n_steps,
                                   final_boxes[k], final_boxes[nk], pad, min_arrow_len)
        if span is not None:
            arc_xml, arc_id, path = _build_cycle_arc(cx, cy, radius, span[0], span[1], ids,
                                                     dashed=(step.arrow_style == "dashed"))
            xml_parts.append(arc_xml)
            arrow_records.append((path[0][0], path[0][1], path[-1][0], path[-1][1], arc_id))
            arrow_paths.append(path)
            for a, b in zip(path[::4], path[4::4] + path[-1:]):
                placed_boxes.append((min(a[0], b[0]) - 3.0, min(a[1], b[1]) - 3.0,
                                     max(a[0], b[0]) + 3.0, max(a[1], b[1]) + 3.0))
            lowest_y = max(lowest_y, max(p[1] for p in path) + 3.0)
            continue
        x0, y0 = centers[k]
        x1, y1 = centers[nk]
        dx, dy = x1 - x0, y1 - y0
        t0 = _ray_box_exit_t(x0, y0, dx, dy, final_boxes[k], pad)
        t1 = 1.0 - _ray_box_exit_t(x1, y1, -dx, -dy, final_boxes[nk], pad)
        if t1 <= t0 + 0.14:
            mid_t = (t0 + t1) / 2.0
            t0 = max(0.12, mid_t - 0.08)
            t1 = min(0.88, mid_t + 0.08)
        tx = x0 + dx * t0
        ty = y0 + dy * t0
        hx = x0 + dx * t1
        hy = y0 + dy * t1
        arrow_xml, arrow_id = _build_arrow(
            tx,
            ty,
            hx,
            hy,
            ids,
            dashed=(step.arrow_style == "dashed"),
            nogo=(step.arrow_style == "failed"),
            kind=_arrow_kind(step.arrow_style),
        )
        xml_parts.append(arrow_xml)
        arrow_records.append((tx, ty, hx, hy, arrow_id))
        arrow_paths.append([(tx, ty), (hx, hy)])
        abox = (min(tx, hx) - 1.5, min(ty, hy) - 3.0, max(tx, hx) + 1.5, max(ty, hy) + 3.0)
        placed_boxes.append(abox)
        lowest_y = max(lowest_y, abox[3])

    # 2b. Place entering reactants and departing products on radial connectors.
    # The ring intermediate drawn after step k is steps[k+1].substrates[0]; any other
    # species listed for the step enters or leaves the cycle there. Without this the
    # cycle only shows its intermediates and silently drops ArX / ArB(OH)2 / product.
    extra_in: Dict[int, List[int]] = {}
    extra_out: Dict[int, List[int]] = {}
    for k, step in enumerate(steps):
        nk = (k + 1) % n_steps
        ring_here = step.substrates[0] if step.substrates else None
        ring_next = steps[nk].substrates[0] if steps[nk].substrates else None
        entering = [s for s in step.substrates[1:] if s != ring_here]
        leaving = [p for p in step.products if p != ring_next]
        if not entering and not leaving:
            continue

        tx, ty, hx, hy, _aid = arrow_records[k]
        mx_, my_ = _along(arrow_paths[k], 0.5)
        rx_, ry_ = mx_ - cx, my_ - cy
        r_norm = math.hypot(rx_, ry_) or 1.0
        rx_, ry_ = rx_ / r_norm, ry_ / r_norm

        for kind, sids in (("in", entering), ("out", leaving)):
            for j, sid in enumerate(sids):
                struct = scheme.structures.get(sid)
                if struct is None:
                    continue
                # stagger along the chord so an entering and a leaving species
                # belonging to the same step do not land on top of each other
                frac = 0.16 if kind == "in" else 0.84
                ax, ay = _along(arrow_paths[k], frac)
                dist = radius * 0.48 + 34.0 + 30.0 * j
                px_ = ax + rx_ * dist
                py_ = ay + ry_ * dist

                frag = _resolve_structure(
                    struct,
                    center_x=px_,
                    center_y=py_,
                    source_data=source_data,
                    reference_mol=reference_mol,
                )
                if frag is None:
                    continue
                _shift_atoms(frag.atoms, px_ - frag.cx, py_ - frag.cy)
                frag.bbox = _fragment_bbox(frag.atoms)
                frag.cx, frag.cy = _bbox_center(frag.bbox)
                frag_xml, _, fid = _build_fragment(frag.atoms, frag.bonds, ids)
                xml_parts.append(frag_xml)

                ring = _overlay_box(frag.atoms, getattr(struct, "overlays", ()))
                box_y2 = frag.bbox[3]
                if frag.ref.label:
                    lbl_y = max(frag.bbox[3], ring[3] if ring else frag.bbox[3]) + 14.0
                    lbl_xml, _ = _build_label_element(frag.ref.label, frag.cx, lbl_y, ids)
                    xml_parts.append(lbl_xml)
                    box_y2 = lbl_y + _DESCENT

                fbox = _union_box((frag.bbox[0] - 2.0, frag.bbox[1] - 2.0, frag.bbox[2] + 2.0, box_y2 + 2.0), ring)
                placed_boxes.append(fbox)
                lowest_y = max(lowest_y, box_y2)

                # short connector arrow between the species and the chord arrow
                ddx, ddy = ax - frag.cx, ay - frag.cy
                edge_t = _ray_box_exit_unclamped(frag.cx, frag.cy, ddx, ddy, fbox, 3.0)
                ex = frag.cx + ddx * edge_t
                ey = frag.cy + ddy * edge_t
                if kind == "in":
                    conn_xml, conn_id = _build_arrow(ex, ey, ax, ay, ids)
                else:
                    conn_xml, conn_id = _build_arrow(ax, ay, ex, ey, ids)
                xml_parts.append(conn_xml)
                placed_boxes.append(
                    (min(ex, ax) - 1.5, min(ey, ay) - 3.0, max(ex, ax) + 1.5, max(ey, ay) + 3.0)
                )

                if kind == "in":
                    extra_in.setdefault(k, []).append(fid)
                else:
                    extra_out.setdefault(k, []).append(fid)

    # 3. Place above / below condition labels radially outside / inside each arrow
    step_metadata: List[Dict[str, Any]] = []
    for k, step in enumerate(steps):
        nk = (k + 1) % n_steps
        tx, ty, hx, hy, arrow_id = arrow_records[k]
        mx, my = _along(arrow_paths[k], 0.5)
        rx, ry = mx - cx, my - cy
        r_norm = math.hypot(rx, ry)
        if r_norm > 1e-6:
            rx, ry = rx / r_norm, ry / r_norm
        else:
            dx, dy = hx - tx, hy - ty
            c_len = math.hypot(dx, dy) or 1.0
            rx, ry = dy / c_len, -dx / c_len

        above_lines: List[str] = []
        if step.above_arrow and step.above_arrow.text:
            for raw in step.above_arrow.text:
                above_lines.append(normalize_condition_text(raw))

        below_lines: List[str] = []
        if step.below_arrow and step.below_arrow.text:
            for raw in step.below_arrow.text:
                below_lines.append(normalize_condition_text(raw))
        if step.yield_:
            y_str = f"({step.yield_})" if not step.yield_.startswith("(") else step.yield_
            below_lines.append(y_str)

        above_ids: List[int] = []
        below_ids: List[int] = []

        if above_lines:
            placed = _place_radial_text(above_lines, mx, my, rx, ry, placed_boxes)
            if placed is not None:
                px_t, py_t, tbox = placed
                placed_boxes.append((tbox[0] - 2.0, tbox[1] - 2.0, tbox[2] + 2.0, tbox[3] + 2.0))
                txt_xml, txt_id = _build_text_element(
                    above_lines, px_t, py_t, ids, use_formatting=True
                )
                xml_parts.append(txt_xml)
                above_ids.append(txt_id)
                lowest_y = max(lowest_y, tbox[3])

        if below_lines:
            # Try radially inside first; fall back to radially outside if interior is tight
            placed = _place_radial_text(
                below_lines, mx, my, -rx, -ry, placed_boxes, max_dist=radius * 0.58
            )
            if placed is None:
                placed = _place_radial_text(below_lines, mx, my, rx, ry, placed_boxes)
            if placed is not None:
                px_t, py_t, tbox = placed
                placed_boxes.append((tbox[0] - 2.0, tbox[1] - 2.0, tbox[2] + 2.0, tbox[3] + 2.0))
                txt_xml, txt_id = _build_text_element(
                    below_lines, px_t, py_t, ids, use_formatting=True
                )
                xml_parts.append(txt_xml)
                below_ids.append(txt_id)
                lowest_y = max(lowest_y, tbox[3])

        step_metadata.append({
            "reactant_ids": [frag_ids[k]] + extra_in.get(k, []),
            "product_ids": [frag_ids[nk]] + extra_out.get(k, []),
            "arrow_id": arrow_id,
            "above_ids": above_ids,
            "below_ids": below_ids,
        })

    # 4. Build <scheme> with <step> elements
    scheme_id = ids.next()
    scheme_parts = [f'<scheme id="{scheme_id}">']
    for meta in step_metadata:
        step_id = ids.next()
        attrs = [
            f'id="{step_id}"',
            f'ReactionStepReactants="{" ".join(str(x) for x in meta["reactant_ids"])}"',
            f'ReactionStepProducts="{" ".join(str(x) for x in meta["product_ids"])}"',
            f'ReactionStepArrows="{meta["arrow_id"]}"',
        ]
        if meta["above_ids"]:
            attrs.append(
                f'ReactionStepObjectsAboveArrow="{" ".join(str(x) for x in meta["above_ids"])}"'
            )
        if meta["below_ids"]:
            attrs.append(
                f'ReactionStepObjectsBelowArrow="{" ".join(str(x) for x in meta["below_ids"])}"'
            )
        scheme_parts.append(f'<step {" ".join(attrs)}/>')
    scheme_parts.append("</scheme>")
    xml_parts.append("\n".join(scheme_parts))

    return "\n".join(xml_parts), lowest_y


def _render_impl(scheme: SchemeDescriptor, yaml_dir: Optional[str] = None) -> str:
    """
    Render a SchemeDescriptor to a CDXML document string.

    Parameters
    ----------
    scheme : SchemeDescriptor
        Parsed scheme from YAML.
    yaml_dir : str, optional
        Directory of the YAML file (for resolving relative source paths).

    Returns
    -------
    str
        Complete CDXML document.
    """
    ids = _IDGen(1000)

    # Load source JSON if specified
    source_data = None
    if scheme.source:
        source_path = scheme.source
        if not os.path.isabs(source_path) and yaml_dir:
            source_path = os.path.join(yaml_dir, source_path)
        source_data = _load_source_json(source_path)

    # --- Identify product reference Mol for alignment ---
    # The product of the last step is the reference for MCS alignment.
    # All other structures are aligned to match its scaffold orientation.
    ref_mol = _identify_product_mol(scheme, source_data)

    # Choose layout
    if scheme.layout == "linear":
        inner_xml, lowest_y = _layout_linear(scheme, ids, source_data=source_data,
                                              reference_mol=ref_mol)
    elif scheme.layout == "sequential":
        inner_xml, lowest_y = _layout_sequential(scheme, ids, source_data=source_data,
                                                  reference_mol=ref_mol)
    elif scheme.layout == "divergent":
        inner_xml, lowest_y = _layout_divergent(scheme, ids, source_data=source_data,
                                                 reference_mol=ref_mol)
    elif scheme.layout == "fan":
        inner_xml, lowest_y = _layout_fan(scheme, ids, source_data=source_data, reference_mol=ref_mol)
    elif scheme.layout == "stacked-rows":
        inner_xml, lowest_y = _layout_stacked_rows(scheme, ids, source_data=source_data,
                                                    reference_mol=ref_mol)
    elif scheme.layout == "convergent":
        inner_xml, lowest_y = _layout_convergent(scheme, ids, source_data=source_data,
                                                  reference_mol=ref_mol)
    elif scheme.layout in ("catalytic-cycle", "cycle"):
        inner_xml, lowest_y = _layout_catalytic_cycle(scheme, ids, source_data=source_data,
                                                       reference_mol=ref_mol)
    else:
        raise NotImplementedError(
            f"Layout '{scheme.layout}' is not yet implemented. "
            f"Supported: linear, sequential, divergent, stacked-rows, convergent, catalytic-cycle"
        )

    # Condition key block below the scheme
    if scheme.condition_key:
        key_y = lowest_y + 20.0
        key_xml, _ = _build_condition_key(
            {k: normalize_condition_text(v) for k, v in scheme.condition_key.items()},
            80.0, key_y, ids,
        )
        inner_xml += "\n" + key_xml

    return _wrap_document(inner_xml, ids)


def _wrap_document(inner_xml: str, ids: _IDGen) -> str:
    """Put page content into a complete CDXML document."""
    page_id = ids.next()

    # Use a generous bounding box
    bbox = "0 0 1620 2160"

    doc_parts = [
        _format_header(bbox),
        _PAGE_OPEN.format(page_id=page_id),
        inner_xml,
        _PAGE_CLOSE,
        CDXML_FOOTER,
    ]
    return "\n".join(doc_parts)


class _ValidatedCDXML(str):
    pass


def render_with(build) -> str:
    """Build page content with ``build(ids) -> inner_xml`` and wrap it like a rendered scheme.

    For figures that place their own objects (mechanisms) but still want the
    ACS document settings and the per-fragment validation receipts of ``render``.
    """
    def document() -> str:
        ids = _IDGen(1000)
        return _wrap_document(build(ids), ids)
    return _validated(document)


def render(scheme: SchemeDescriptor, yaml_dir: Optional[str] = None) -> str:
    """Render a scheme and retain per-fragment output-derived validation receipts."""
    return _validated(lambda: _render_impl(scheme, yaml_dir))


def _validated(produce) -> str:
    receipts = []
    token = _fragment_receipts.set(receipts)
    try:
        result = _ValidatedCDXML(produce())
        from xml.etree import ElementTree as ET
        root = ET.fromstring(result)
        fragment_ids = {node.get('id') for node in root.iter('fragment')}
        validated_ids = {str(receipt['fragment_id']) for receipt in receipts}
        result.chemistry_validation = {
            'status': 'preserved' if fragment_ids and fragment_ids == validated_ids else 'unverified',
            'method': 'final_coordinate_fragment_roundtrips', 'fragments': receipts}
        return result
    finally:
        _fragment_receipts.reset(token)


def render_to_file(
    scheme: SchemeDescriptor,
    output_path: str,
    yaml_dir: Optional[str] = None,
) -> None:
    """Render and write to a file."""
    cdxml = render(scheme, yaml_dir=yaml_dir)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(cdxml)
