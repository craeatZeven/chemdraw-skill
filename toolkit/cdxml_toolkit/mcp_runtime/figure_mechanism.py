"""Reaction mechanisms for draw_figure: structures with the curved arrows the brief lists.

The brief names every species, the steps in the order they happen and, for
each step, the electron moves as (from atom or bond, to atom or bond). This
module places the structures, turns each move into a curved arrow from the
drawn 2D coordinates, checks that no arrow runs through another atom or bond (one
automatic re-layout if one does), and reads the saved file back. The arrows are
ChemDraw's own arc arrows (an <arrow> with AngularSize), each in the colour its move
asks for. It never
decides which electrons move: a brief without moves is refused, and a move
that looks wrong is drawn as given and reported.
"""
from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rdkit import Chem

from ..cdxml_check import overlays_collide
from ..constants import (
    ACS_BOND_LENGTH,
    ACS_LINE_WIDTH,
    LAYOUT_BELOW_GAP,
    LAYOUT_FRAG_GAP_BONDS,
    LAYOUT_INTER_GAP_BONDS,
    LAYOUT_PLUS_GAP,
)
from ..render import renderer as R
from ..render.schema import StructureRef
from ..text_formatting import normalize_condition_text
from . import artifact_safety, figure_draw, figure_house, figure_lint

MECH_KEYS = {"kind", "species", "steps", "max_width", "title", "style", "show_atom_numbers"}
SPECIES_KEYS = {"id", "smiles", "name", "label", "newman", "generic_rings"}
NEWMAN_KEYS = {"front", "rear", "conformation", "dihedral_offset"}
NEWMAN_SPOKE = 1.4 * ACS_BOND_LENGTH   # a projection's spokes, sized beside 14.4 pt bonds
# Closest a rear spoke may come to a front one; ChemDraw-style eclipsed drawings use 14.
NEWMAN_MIN_APART = 10.0
STEP_KEYS = {"substrates", "products", "above", "below", "yield", "arrow", "moves"}
MOVE_KEYS = {"from", "to", "electrons", "color"}
END_KEYS = {"species", "atom", "bond"}
MECH_ARROWS = ("solid", "reversible", "resonance")

X0 = 100.0                 # left edge of every row
FIRST_ROW_TOP = 120.0      # page y of the top of the first row
ROW_GAP = 18.0
TITLE_ROOM = 18.0
BULGE_MIN, BULGE_MAX = 9.0, 28.0
RELAYOUT_BULGE = 1.35
SAMPLES = 24
CARBON_REACH, LABEL_REACH = 3.5, 6.5   # how close a curve may pass an atom
CARBON_STANDOFF, LABEL_STANDOFF = 3.0, 7.0   # where an arrow starts or stops from an atom
HEIGHT_STEPS = (1.0, 1.4, 1.9, 2.5)   # arc heights tried, relative to the default
LABEL_OBSTACLE = 8.0
BOND_CLEARANCE = 2.0       # how close an arrow may pass a bond it does not start or end at
ARROW_CLEARANCE = 3.0      # how close two arrows may run, away from their tails
HEAD_CLEARANCE = 6.0       # how close two arrowheads may come
TAIL_SHARE = 0.2           # the part of an arrow, from its tail, that may meet another's tail


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _end_errors(where: str, end: Any, ids: set[str]) -> list[str]:
    if not isinstance(end, dict):
        return [f"{where}: give the species and exactly one of atom or bond."]
    errors: list[str] = []
    extra = sorted(set(end) - END_KEYS)
    if extra:
        errors.append(f"{where}: unknown keys {extra}.")
    species = end.get("species")
    if species is None:
        errors.append(f"{where}: species is missing.")
    elif species not in ids:
        errors.append(f"{where}: species {species!r} is not in species.")
    given = [key for key in ("atom", "bond") if end.get(key) is not None]
    if len(given) != 1:
        errors.append(f"{where}: give exactly one of atom or bond.")
    elif given[0] == "atom" and not (_is_int(end["atom"]) and end["atom"] >= 0):
        errors.append(f"{where}: atom must be a whole number from 0.")
    elif given[0] == "bond":
        bond = end["bond"]
        if not (isinstance(bond, (list, tuple)) and len(bond) == 2 and all(_is_int(i) and i >= 0 for i in bond)
                and bond[0] != bond[1]):
            errors.append(f"{where}: bond must be two different atom numbers, such as [1, 2].")
    return errors


def _validate(brief: dict) -> list[str]:
    errors: list[str] = []
    unknown = sorted(set(brief) - MECH_KEYS)
    if unknown:
        errors.append(f"Unknown brief keys: {unknown}")
    style_error = figure_draw.style_preset_error(brief.get("style"))
    if style_error:
        errors.append(style_error)
    width = brief.get("max_width")
    if width is not None and not (isinstance(width, (int, float)) and not isinstance(width, bool) and width > 0):
        errors.append("max_width must be a positive number of points.")
    if brief.get("show_atom_numbers") is not None and not isinstance(brief["show_atom_numbers"], bool):
        errors.append("show_atom_numbers must be true or false.")
    if brief.get("title") is not None and not isinstance(brief["title"], str):
        errors.append("title must be text.")

    species = brief.get("species")
    ids: set[str] = set()
    newman_ids: set[str] = set()
    if not isinstance(species, list) or not species:
        errors.append("species: list every structure that is drawn, each with an id and a smiles or a name.")
        species = []
    for entry in species:
        if not isinstance(entry, dict) or not str(entry.get("id") or "").strip():
            errors.append("Every species needs an id.")
            continue
        ident = str(entry["id"])
        if ident in ids:
            errors.append(f"Species {ident}: duplicate id.")
        ids.add(ident)
        extra = sorted(set(entry) - SPECIES_KEYS)
        if extra:
            errors.append(f"Species {ident}: unknown keys {extra}.")
        if "newman" in entry:
            newman_ids.add(ident)
            errors += _newman_errors(ident, entry["newman"])
        elif not entry.get("smiles") and not entry.get("name"):
            errors.append(f"Species {ident} needs a smiles or a name (or a newman projection).")
        errors += figure_house.generic_ring_errors(entry)

    steps = brief.get("steps")
    if not isinstance(steps, list) or not steps:
        errors.append("steps: give at least one step with substrates, products and moves.")
        steps = []
    produced: dict[str, int] = {}
    for number, step in enumerate(steps, 1):
        if not isinstance(step, dict):
            errors.append(f"Step {number} must be an object.")
            continue
        extra = sorted(set(step) - STEP_KEYS)
        if extra:
            errors.append(f"Step {number}: unknown keys {extra}.")
        members: dict[str, list] = {}
        for side in ("substrates", "products"):
            members[side] = _as_ids(step.get(side))
            if not members[side]:
                errors.append(f"Step {number} needs {side}.")
            for member in members[side]:
                if member not in ids:
                    errors.append(f"Step {number} names {member!r} in {side}, which is not in species.")
        for product in members["products"]:
            if product in produced:
                errors.append(f"Species {product} is produced by step {produced[product]} and step {number}; "
                              "draw a second copy under another id if both are meant.")
            produced.setdefault(product, number)
        if step.get("arrow", "solid") not in MECH_ARROWS:
            errors.append(f"Step {number}: arrow must be one of {list(MECH_ARROWS)}.")
        moves = step.get("moves")
        conformers = bool(newman_ids & set(members["substrates"] + members["products"]))
        if conformers and moves in (None, []):
            continue  # a change of conformation, or a step drawn on projections, has no arrows to list
        if not isinstance(moves, list) or not moves:
            errors.append(f"Step {number} needs moves: list each electron move as from and to "
                          "(an atom number or a bond). A mechanism is never inferred.")
            continue
        for index, move in enumerate(moves, 1):
            errors += _move_errors(number, index, move, members["substrates"], ids)
            for side in ("from", "to"):
                end = move.get(side) if isinstance(move, dict) else None
                if isinstance(end, dict) and str(end.get("species")) in newman_ids:
                    errors.append(f"Step {number} move {index} {side}: a Newman projection has no numbered "
                                  "atoms for an arrow; point it at a structure given by smiles or name.")
    return errors


def _newman_errors(ident: str, spec: Any) -> list[str]:
    if not isinstance(spec, dict):
        return [f"Species {ident}: newman must be an object with front and rear substituents."]
    errors = []
    extra = sorted(set(spec) - NEWMAN_KEYS)
    if extra:
        errors.append(f"Species {ident}: unknown newman keys {extra}.")
    for side in ("front", "rear"):
        group = spec.get(side)
        if not (isinstance(group, list) and len(group) == 3 and all(isinstance(g, str) and g for g in group)):
            errors.append(f"Species {ident}: newman {side} must list three substituents, e.g. [\"Me\", \"H\", \"H\"].")
    if spec.get("conformation", "staggered") not in ("staggered", "eclipsed"):
        errors.append(f"Species {ident}: newman conformation must be staggered or eclipsed.")
    offset = spec.get("dihedral_offset")
    if offset is not None and (isinstance(offset, bool) or not isinstance(offset, (int, float))):
        errors.append(f"Species {ident}: newman dihedral_offset must be a number of degrees.")
    elif offset is not None:
        # The three front spokes are 120 degrees apart, so a rear set turned by a multiple of
        # 120 lies on them, and each rear label is printed over a front one.
        apart = abs((float(offset) + 60.0) % 120.0 - 60.0)
        if apart < NEWMAN_MIN_APART:
            errors.append(f"Species {ident}: newman dihedral_offset {offset:g} puts the rear spokes on top of the "
                          f"front ones ({apart:g} degrees apart); it is the rear set's angle from the front set, "
                          "60 for staggered, 14 for a visible eclipsed one. For another conformer, reorder rear.")
    return errors


def _newman_frag(ident: str, spec: dict, label: str | None) -> Any:
    """A Newman projection as a structure the row layout can place: atoms, bonds, box, label."""
    from ..templates.conformers import newman_projection
    drawn = newman_projection(conformation=spec.get("conformation", "staggered"),
                              dihedral_offset=spec.get("dihedral_offset"), front=tuple(spec["front"]),
                              rear=tuple(spec["rear"]), cx=0.0, cy=0.0, bond_length=NEWMAN_SPOKE)
    atoms = [dict(a) for a in drawn["atoms"]]
    frag = SimpleNamespace(atoms=atoms, bonds=[dict(b) for b in drawn["bonds"]], newman=dict(spec),
                           ref=StructureRef(id=ident, label=label))
    frag.bbox = R._fragment_bbox(atoms)
    frag.cx, frag.cy = R._bbox_center(frag.bbox)
    return frag


def _newman_xml(frag: Any, next_id) -> list[str]:
    """The projection drawn where the layout put it, with ids from the page's own sequence.

    The template orders its parts by Z (rear spokes, then the filled disc, then the front
    spokes); new Z values keep that order, so the disc still hides the rear spokes' roots.
    """
    from ..templates.conformers import newman_projection
    centre = frag.atoms[0]  # the front carbon, at the centre of the projection
    spec = frag.newman
    drawn = newman_projection(conformation=spec.get("conformation", "staggered"),
                              dihedral_offset=spec.get("dihedral_offset"), front=tuple(spec["front"]),
                              rear=tuple(spec["rear"]), cx=centre["x"], cy=centre["y"], bond_length=NEWMAN_SPOKE)
    elements = [ET.fromstring(drawn["fragment_xml"])] + [ET.fromstring(g) for g in drawn["graphics_xml"]]
    every = [node for element in elements for node in element.iter()]
    renumber = {node.get("id"): str(next_id()) for node in every if node.get("id")}
    for node in every:
        if node.get("id"):
            node.set("id", renumber[node.get("id")])
        for key in ("B", "E"):
            if node.get(key) in renumber:
                node.set(key, renumber[node.get(key)])
    for node in sorted((n for n in every if n.get("Z")), key=lambda n: int(n.get("Z"))):
        node.set("Z", str(next_id()))
    return [ET.tostring(element, encoding="unicode") for element in elements]


def _newman_inventory(spec: dict) -> list[str]:
    """What the readback finds for this projection: the template read the way pages are read."""
    import tempfile
    from ..chemistry_semantics import document_inventory
    from ..templates.conformers import newman_projection
    drawn = newman_projection(conformation=spec.get("conformation", "staggered"),
                              dihedral_offset=spec.get("dihedral_offset"), front=tuple(spec["front"]),
                              rear=tuple(spec["rear"]), bond_length=NEWMAN_SPOKE)
    with tempfile.TemporaryDirectory(prefix="newman-") as scratch:
        path = Path(scratch) / "newman.cdxml"
        path.write_text(drawn["cdxml"], encoding="utf-8")
        try:
            return list(document_inventory(path).elements())
        except ValueError:
            return []


def _as_ids(value: Any) -> list[str]:
    return [str(i) for i in figure_draw._as_list(value)]


def _move_errors(number: int, index: int, move: Any, substrates: list[str], ids: set[str]) -> list[str]:
    where = f"Step {number} move {index}"
    if not isinstance(move, dict):
        return [f"{where} must be an object with from and to."]
    errors: list[str] = []
    extra = sorted(set(move) - MOVE_KEYS)
    if extra:
        errors.append(f"{where}: unknown keys {extra}.")
    if move.get("electrons", 2) not in (1, 2) or isinstance(move.get("electrons"), bool):
        errors.append(f"{where}: electrons must be 1 or 2.")
    if "color" in move and not (isinstance(move["color"], str) and re.fullmatch(r"#[0-9A-Fa-f]{6}", move["color"])):
        errors.append(f"{where}: color must be a hex colour like #2E7D32.")
    for side in ("from", "to"):
        errors += _end_errors(f"{where} {side}", move.get(side), ids)
        end = move.get(side)
        if isinstance(end, dict) and end.get("species") in ids and end["species"] not in substrates:
            errors.append(f"{where} {side}: species {end['species']!r} is not one of the step's substrates; "
                          "electrons move within the species that react.")
    if not errors and _end_key(move["from"]) == _end_key(move["to"]):
        errors.append(f"{where}: the arrow starts and ends on the same atom or bond.")
    return errors


def _end_key(end: dict) -> tuple:
    if end.get("atom") is not None:
        return (end["species"], "atom", end["atom"])
    return (end["species"], "bond", tuple(sorted(end["bond"])))


def _index_errors(brief: dict, mols: dict[str, Any]) -> list[str]:
    """Atom and bond numbers must exist in the grounded structures."""
    errors: list[str] = []
    for number, step in enumerate(brief["steps"], 1):
        for index, move in enumerate(step.get("moves") or [], 1):
            for side in ("from", "to"):
                end = move[side]
                mol = mols[end["species"]]
                where = f"Step {number} move {index} {side}"
                if end.get("atom") is not None:
                    if end["atom"] >= mol.GetNumAtoms():
                        errors.append(f"{where}: atom {end['atom']} does not exist in {end['species']} "
                                      f"({mol.GetNumAtoms()} atoms, numbered from 0).")
                else:
                    i, j = end["bond"]
                    if max(i, j) >= mol.GetNumAtoms() or mol.GetBondBetweenAtoms(i, j) is None:
                        errors.append(f"{where}: bond {[i, j]} does not exist in {end['species']}.")
    return errors


# ---------------------------------------------------------------------------
# Move checks (reported, never blocking)
# ---------------------------------------------------------------------------

def _describe(end: dict, mol) -> str:
    if end.get("atom") is not None:
        return f"{end['species']} atom {end['atom']} ({mol.GetAtomWithIdx(end['atom']).GetSymbol()})"
    i, j = end["bond"]
    return f"{end['species']} bond {i}-{j}"


def _lone_pairs(atom) -> int:
    outer = Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())
    nonbonding = outer - atom.GetFormalCharge() - atom.GetTotalValence() - atom.GetNumRadicalElectrons()
    return max(nonbonding // 2, 0)


def _move_checks(move: dict, mol) -> list[str]:
    checks: list[str] = []
    source = move["from"]
    if source.get("atom") is None:
        return checks
    atom = mol.GetAtomWithIdx(source["atom"])
    name = _describe(source, mol)
    if move.get("electrons", 2) == 2:
        if atom.GetFormalCharge() >= 0 and _lone_pairs(atom) == 0:
            checks.append(f"{name} has no lone pair or negative charge to donate; "
                          "start the arrow on a bond or on an atom that has one.")
    elif atom.GetNumRadicalElectrons() == 0:
        checks.append(f"{name} has no unpaired electron for a one-electron arrow.")
    return checks


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

@dataclass
class _Inst:
    uid: int
    species: str
    frag: Any
    row: int = 0


@dataclass
class _Drawing:
    cdxml: Any
    lint_cdxml: Any
    rows: list
    curves: list
    atoms: list
    bonds: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def _reach(atom: dict) -> float:
    plain = atom.get("symbol", "C") == "C" and not atom.get("charge") and not atom.get("isotope")
    return CARBON_REACH if plain else LABEL_REACH


def _place(inst: _Inst, left: float, row: int) -> float:
    """Put the structure's left edge at ``left`` and its centre on the arrow line; return its right edge."""
    frag = inst.frag
    R._shift_atoms(frag.atoms, left + R._bbox_width(frag.bbox) / 2.0 - frag.cx, -frag.cy)
    frag.bbox = R._fragment_bbox(frag.atoms)
    frag.cx, frag.cy = R._bbox_center(frag.bbox)
    inst.row = row
    return frag.bbox[2]


def _text_lines(step: dict) -> tuple[list[str], list[str]]:
    above = [normalize_condition_text(str(t)) for t in figure_draw._as_list(step.get("above"))]
    below = [normalize_condition_text(str(t)) for t in figure_draw._as_list(step.get("below"))]
    if step.get("yield") not in (None, ""):
        below.append(str(step["yield"]))
    return above, below


def _arrange(brief: dict, grounded: list[dict], width: float, gap_scale: float):
    """Place every structure and step arrow in rows; returns (rows, per-step instances, notes)."""
    bond = ACS_BOND_LENGTH
    frag_gap = bond * LAYOUT_FRAG_GAP_BONDS * gap_scale
    inter_gap = max(bond * LAYOUT_INTER_GAP_BONDS, LAYOUT_PLUS_GAP) * gap_scale
    limit = X0 + width
    by_id = {g["id"]: g for g in grounded}
    insts: list[_Inst] = []
    notes: list[str] = []

    def structure(species: str, label: str | None = None) -> Any:
        item = by_id[species]
        if item.get("newman"):
            return _newman_frag(species, item["newman"], label)
        return R._resolve_structure(StructureRef(id=species, smiles=item["smiles"], label=label))

    def new(species: str) -> _Inst:
        item = by_id[species]
        label = str(item["label"]) if item["label"] else None
        inst = _Inst(len(insts), species, structure(species, label))
        insts.append(inst)
        return inst

    def place_group(row_items: list, ids: list[str], cursor: float, row: int) -> tuple[list[_Inst], float]:
        placed: list[_Inst] = []
        for position, species in enumerate(ids):
            if position:
                row_items.append({"t": "plus", "x": cursor - inter_gap / 2.0})
            inst = new(species)
            cursor = _place(inst, cursor, row) + inter_gap
            row_items.append({"t": "frag", "inst": inst})
            placed.append(inst)
        return placed, cursor

    rows: list[list[dict]] = []
    step_insts: list[dict] = []
    drawn: set[str] = set()
    previous_products: list[str] = []
    last_product: _Inst | None = None

    for number, step in enumerate(brief["steps"], 1):
        substrates, products = _as_ids(step["substrates"]), _as_ids(step["products"])
        above, below = _text_lines(step)
        arrow_len = max(max(R._estimate_text_width(above), R._estimate_text_width(below)) + 10.0, 5.0 * bond)
        product_widths = []
        for species in products:
            probe = structure(species)
            product_widths.append(R._bbox_width(probe.bbox))
        after_arrow = 0.3 * frag_gap + arrow_len + frag_gap + sum(product_widths) + inter_gap * (len(products) - 1)

        chained = number > 1 and substrates == previous_products and len(previous_products) == 1 and last_product
        sub_insts: list[_Inst]
        if chained:
            row_items = rows[-1]
            cursor = last_product.frag.bbox[2] + inter_gap
            if cursor + after_arrow > limit and any(i["t"] == "arrow" for i in row_items):
                rows.append([])
                row_items = rows[-1]
                sub_insts, cursor = place_group(row_items, substrates, X0, len(rows) - 1)
                notes.append(f"species {substrates[0]} is drawn again at the start of row {len(rows)} "
                             "because the route does not fit on one row.")
            else:
                sub_insts = [last_product]
        else:
            rows.append([])
            row_items = rows[-1]
            for species in substrates:
                if species in drawn:
                    notes.append(f"species {species} is drawn again in step {number} (row {len(rows)}) because "
                                 "that step has more than one substrate or does not directly follow the previous one.")
            sub_insts, cursor = place_group(row_items, substrates, X0, len(rows) - 1)
        drawn.update(substrates)

        tail = cursor + 0.3 * frag_gap
        head = tail + arrow_len
        row_items.append({"t": "arrow", "tail": tail, "head": head, "above": above, "below": below,
                          "style": step.get("arrow", "solid")})
        prod_insts, cursor = place_group(row_items, products, head + frag_gap, len(rows) - 1)
        drawn.update(products)
        step_insts.append({"subs": dict(zip(substrates, sub_insts)), "prods": dict(zip(products, prod_insts))})
        previous_products = products
        last_product = prod_insts[0] if len(prod_insts) == 1 else None
    return rows, step_insts, notes


def _row_extent(row: list[dict]) -> tuple[float, float]:
    top = bottom = 0.0
    for item in row:
        if item["t"] == "frag":
            frag = item["inst"].frag
            top = min(top, frag.bbox[1])
            bottom = max(bottom, frag.bbox[3] + (20.0 if frag.ref.label else 0.0))
        elif item["t"] == "arrow":
            if item["above"]:
                lines = len(item["above"])
                top = min(top, -LAYOUT_BELOW_GAP - ((lines - 1) * R._LINE_ADVANCE + R._DESCENT) - R._CAP_HEIGHT)
            if item["below"]:
                lines = len(item["below"])
                bottom = max(bottom, LAYOUT_BELOW_GAP + R._CAP_HEIGHT + (lines - 1) * R._LINE_ADVANCE + R._DESCENT)
    return top, bottom


# ---------------------------------------------------------------------------
# Curved arrows
# ---------------------------------------------------------------------------

def _anchor(inst: _Inst, end: dict) -> tuple[tuple[float, float], float, set]:
    """Page position, how far an arrow stays off the atom, and the atoms the endpoint belongs to."""
    by_index = dict(enumerate(inst.frag.atoms))  # the renderer numbers atoms from 1; moves count from 0
    if end.get("atom") is not None:
        atom = by_index[end["atom"]]
        return (atom["x"], atom["y"]), (LABEL_STANDOFF if _reach(atom) == LABEL_REACH else CARBON_STANDOFF),             {(inst.uid, end["atom"])}
    i, j = end["bond"]
    a, b = by_index[i], by_index[j]
    return ((a["x"] + b["x"]) / 2.0, (a["y"] + b["y"]) / 2.0), 0.0, {(inst.uid, i), (inst.uid, j)}


def _unit(dx: float, dy: float) -> tuple[float, float]:
    length = math.hypot(dx, dy) or 1.0
    return dx / length, dy / length


def _up_angle(point, centre) -> float:
    """Angle of ``point`` about ``centre`` with y up: anticlockwise as seen on the page."""
    return math.atan2(-(point[1] - centre[1]), point[0] - centre[0])


def _arc(src, src_r: float, dst, dst_r: float, side: float, scale: float) -> dict:
    """A circular arc from ``src`` to ``dst`` that bulges to ``side`` of the chord.

    The bulge (the sagitta) is 0.45 of the chord, kept between BULGE_MIN and BULGE_MAX and
    multiplied by ``scale``; ``side`` +1 or -1 picks the normal (-uy, ux) or its opposite. The
    arc then stops ``src_r`` and ``dst_r`` short of the atoms, along the circle itself, so it
    lands on an atom's edge from the side it arrives on. ``pts`` samples the arc from source
    to target; ``sweep`` is signed, positive when the arc turns anticlockwise on the page.
    """
    chord = max(math.dist(src, dst), 1.0)
    ux, uy = _unit(dst[0] - src[0], dst[1] - src[1])
    nx, ny = -uy * side, ux * side
    height = min(max(0.45 * chord, BULGE_MIN), BULGE_MAX) * scale
    radius = chord * chord / (8.0 * height) + height / 2.0
    mx, my = (src[0] + dst[0]) / 2.0, (src[1] + dst[1]) / 2.0
    centre = (mx - nx * (radius - height), my - ny * (radius - height))
    total = 4.0 * math.atan(2.0 * height / chord)          # the angle the whole arc subtends
    begin = _up_angle(src, centre)
    apex = (mx + nx * height, my + ny * height)

    def at(angle):
        return centre[0] + radius * math.cos(angle), centre[1] - radius * math.sin(angle)
    # Turn the way that passes the apex half-way: anticlockwise (+1) or clockwise (-1) on the page.
    turn = 1.0 if math.dist(at(begin + total / 2.0), apex) <= math.dist(at(begin - total / 2.0), apex) else -1.0
    trim_src = 2.0 * math.asin(min(1.0, src_r / (2.0 * radius)))
    trim_dst = 2.0 * math.asin(min(1.0, dst_r / (2.0 * radius)))
    sweep = max(total - trim_src - trim_dst, 0.15 * total)
    first = begin + turn * trim_src
    pts = [at(first + turn * sweep * k / SAMPLES) for k in range(SAMPLES + 1)]
    return {"centre": centre, "radius": radius, "sweep": turn * sweep, "pts": pts}


def _shift_arc(arc: dict, dy: float) -> dict:
    return {**arc, "centre": (arc["centre"][0], arc["centre"][1] + dy),
            "pts": [(x, y + dy) for x, y in arc["pts"]]}


def _hits(points, atoms: list[dict], skip: set) -> list[dict]:
    found: dict[Any, dict] = {}
    samples = points[1:-1]
    for atom in atoms:
        if atom["key"] in skip:
            continue
        if any(math.hypot(x - atom["x"], y - atom["y"]) < atom["reach"] for x, y in samples):
            found[atom["key"]] = atom
    return list(found.values())


def _segment_distance(p, q, a, b) -> float:
    """Shortest distance between the segments p-q and a-b."""
    def point(c, s, t):
        dx, dy = t[0] - s[0], t[1] - s[1]
        length = dx * dx + dy * dy
        u = 0.0 if not length else max(0.0, min(1.0, ((c[0] - s[0]) * dx + (c[1] - s[1]) * dy) / length))
        return math.hypot(c[0] - s[0] - u * dx, c[1] - s[1] - u * dy)

    def side(o, s, t):
        return (t[0] - s[0]) * (o[1] - s[1]) - (t[1] - s[1]) * (o[0] - s[0])
    if side(a, p, q) * side(b, p, q) < 0 and side(p, a, b) * side(q, a, b) < 0:
        return 0.0
    return min(point(p, a, b), point(q, a, b), point(a, p, q), point(b, p, q))


def _bond_hits(points, bonds: list[dict], skip: set) -> list[dict]:
    """Bonds an arrow cuts or grazes, other than those at the atoms it starts or ends on."""
    found = []
    for bond in bonds:
        if bond["atoms"] & skip:
            continue
        if any(_segment_distance(p, q, bond["a"], bond["b"]) < BOND_CLEARANCE for p, q in zip(points, points[1:])):
            found.append(bond)
    return found


def _arrows_clash(first, second) -> bool:
    """Whether two arrows (sampled source to target) run into each other.

    Their heads may not come within HEAD_CLEARANCE, nor their bodies within ARROW_CLEARANCE;
    the first TAIL_SHARE of each is left out, as two arrows may leave the same atom.
    """
    if math.dist(first[-1], second[-1]) < HEAD_CLEARANCE:
        return True
    a, b = first[int(len(first) * TAIL_SHARE):], second[int(len(second) * TAIL_SHARE):]
    return any(_segment_distance(p, q, r, t) < ARROW_CLEARANCE for p, q in zip(a, a[1:]) for r, t in zip(b, b[1:]))


def _choose_arc(src, src_r: float, dst, dst_r: float, *, atoms: list[dict], bonds: list[dict], skip: set,
                placed: list, bulge: float, centre) -> dict:
    """The arc that meets the fewest atoms, bonds and arrows already drawn; then the lowest; then
    the preferred side (outside the molecule for a move within one, else higher on the page)."""
    best = None
    for height, scale in enumerate(HEIGHT_STEPS):
        for side in (-1.0, 1.0):
            arc = _arc(src, src_r, dst, dst_r, side, scale * bulge)
            points = arc["pts"]
            hits = (len(_hits(points, atoms, skip)) + len(_bond_hits(points, bonds, skip))
                    + sum(_arrows_clash(points, other) for other in placed))
            mid = points[SAMPLES // 2]
            if centre is not None:
                preference = -math.hypot(mid[0] - centre[0], mid[1] - centre[1])
            else:
                preference = mid[1]
            key = (hits, height, preference)
            if best is None or key < best[0]:
                best = (key, arc)
    return best[1]


def _row_bonds(row: list[dict]) -> list[dict]:
    """The drawn bonds of a row, as segments between atom positions."""
    found = []
    for item in row:
        if item["t"] != "frag":
            continue
        atoms = item["inst"].frag.atoms
        number = {a["index"]: k for k, a in enumerate(atoms)}
        for bond in item["inst"].frag.bonds:
            i, j = number.get(bond.get("atom1")), number.get(bond.get("atom2"))
            if i is None or j is None:
                continue
            found.append({"key": ("bond", item["inst"].uid, min(i, j), max(i, j)),
                          "label": f"the {item['inst'].species} bond {min(i, j)}-{max(i, j)}",
                          "a": (atoms[i]["x"], atoms[i]["y"]), "b": (atoms[j]["x"], atoms[j]["y"]),
                          "atoms": {(item["inst"].uid, i), (item["inst"].uid, j)}})
    return found


def _label_obstacles(row: list[dict]) -> list[dict]:
    """The compound numbers under the structures, so an arrow does not run through them."""
    found = []
    for item in row:
        if item["t"] == "frag" and item["inst"].frag.ref.label:
            frag = item["inst"].frag
            found.append({"key": ("label", item["inst"].uid), "label": f"the label {frag.ref.label}",
                          "x": (frag.bbox[0] + frag.bbox[2]) / 2.0, "y": frag.bbox[3] + 11.0,
                          "reach": LABEL_OBSTACLE})
    return found


def _row_curves(brief: dict, step_insts: list[dict], rows: list[list[dict]], bulge: float) -> dict[int, list[dict]]:
    """Curved arrows in row-local coordinates, keyed by row."""
    atoms_by_row: dict[int, list[dict]] = {}
    for row_index, row in enumerate(rows):
        atoms_by_row[row_index] = [
            {"key": (item["inst"].uid, number), "label": f"{item['inst'].species} atom {number}",
             "x": a["x"], "y": a["y"], "reach": _reach(a)}
            for item in row if item["t"] == "frag" for number, a in enumerate(item["inst"].frag.atoms)
        ] + _label_obstacles(row)
    bonds_by_row = {row_index: _row_bonds(row) for row_index, row in enumerate(rows)}
    curves: dict[int, list[dict]] = {i: [] for i in range(len(rows))}
    for number, step in enumerate(brief["steps"], 1):
        subs = step_insts[number - 1]["subs"]
        for index, move in enumerate(step.get("moves") or [], 1):
            src_inst, dst_inst = subs[move["from"]["species"]], subs[move["to"]["species"]]
            (sx, sy), src_r, src_keys = _anchor(src_inst, move["from"])
            (tx, ty), dst_r, dst_keys = _anchor(dst_inst, move["to"])
            row = src_inst.row
            skip = src_keys | dst_keys
            centre = None
            if src_inst is dst_inst:
                pts = src_inst.frag.atoms
                centre = (sum(a["x"] for a in pts) / len(pts), sum(a["y"] for a in pts) / len(pts))
            # Arrows are routed in the order given, each around those already drawn in its row.
            arc = _choose_arc((sx, sy), src_r, (tx, ty), dst_r, atoms=atoms_by_row[row], bonds=bonds_by_row[row],
                              skip=skip, placed=[c["pts"] for c in curves[row]], bulge=bulge, centre=centre)
            curves[row].append({"step": number, "move": index, "electrons": move.get("electrons", 2),
                                "color": move.get("color"), "arc": arc, "pts": arc["pts"], "skip": skip})
    return curves


def _crossings(curves: list[dict], atoms: list[dict], bonds: list[dict] = ()) -> list[dict]:
    """Atoms, bonds and other arrows a curved arrow runs into (other than at its own endpoints)."""
    found = []
    for number, curve in enumerate(curves):
        for atom in _hits(curve["pts"], atoms, curve["skip"]):
            found.append({"step": curve["step"], "move": curve["move"], "with": atom["label"]})
        for bond in _bond_hits(curve["pts"], list(bonds), curve["skip"]):
            found.append({"step": curve["step"], "move": curve["move"], "with": bond["label"]})
        for other in curves[:number]:
            if _arrows_clash(curve["pts"], other["pts"]):
                found.append({"step": curve["step"], "move": curve["move"],
                              "with": f"the arrow of step {other['step']} move {other['move']}"})
    return found


def _arc_xml(curve: dict, ids) -> str:
    """The arc as ChemDraw writes a curved arrow (see toolkit/tests/data/native-arc-arrows.cdxml).

    ChemDraw runs an arc from Tail3D anticlockwise on the page through AngularSize degrees to
    Head3D. An arrow that turns that way has its tail at the source and the head on Head3D; one
    that turns clockwise is the same arc from the target, with the head on the tail end.
    """
    arc = curve["arc"]
    (cx, cy), r, pts = arc["centre"], arc["radius"], arc["pts"]
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    head = "HalfLeft" if curve["electrons"] == 1 else "Full"
    # A half head is the full head's 1 : 0.875 : 0.25 shape at 0.6 of its size, as ChemDraw writes
    # its own: the thin barb that reads as one electron.
    shape = ('HeadSize="600" ArrowheadCenterSize="525" ArrowheadWidth="150"' if curve["electrons"] == 1
             else 'HeadSize="1000" ArrowheadCenterSize="875" ArrowheadWidth="250"')
    if arc["sweep"] >= 0:
        tail, tip, ends = pts[0], pts[-1], f'ArrowheadHead="{head}"'
    else:
        tail, tip, ends = pts[-1], pts[0], f'ArrowheadTail="{head}" ArrowheadHead="None"'
    return (f'<arrow id="{ids()}" BoundingBox="{min(xs):.2f} {min(ys):.2f} {max(xs):.2f} {max(ys):.2f}" '
            f'Z="{ids()}" FillType="None" LineWidth="{ACS_LINE_WIDTH}" {ends} ArrowheadType="Solid" '
            f'{shape} '
            f'AngularSize="{abs(math.degrees(arc["sweep"])):.2f}" '
            f'Head3D="{tip[0]:.2f} {tip[1]:.2f} 0" Tail3D="{tail[0]:.2f} {tail[1]:.2f} 0" '
            f'Center3D="{cx:.2f} {cy:.2f} 0" MajorAxisEnd3D="{cx + r:.2f} {cy:.2f} 0" '
            f'MinorAxisEnd3D="{cx:.2f} {cy - r:.2f} 0"/>')


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _draw(brief: dict, grounded: list[dict], width: float, gap_scale: float, bulge: float,
          numbers: bool) -> _Drawing:
    rows, step_insts, notes = _arrange(brief, grounded, width, gap_scale)
    local_curves = _row_curves(brief, step_insts, rows, bulge)

    # Row-local layout -> page coordinates: stack the rows, leaving room for the loops.
    title = brief.get("title")
    cursor = FIRST_ROW_TOP + (TITLE_ROOM if title else 0.0)
    offsets: list[float] = []
    for row_index, row in enumerate(rows):
        top, bottom = _row_extent(row)
        for curve in local_curves[row_index]:
            ys = [y for _, y in curve["pts"]]
            top, bottom = min(top, min(ys)), max(bottom, max(ys))
        offset = cursor - top
        offsets.append(offset)
        cursor = offset + bottom + ROW_GAP
        for item in row:
            if item["t"] == "frag":
                R._shift_atoms(item["inst"].frag.atoms, 0.0, offset)
                item["inst"].frag.bbox = R._fragment_bbox(item["inst"].frag.atoms)
        for curve in local_curves[row_index]:
            curve["arc"] = _shift_arc(curve["arc"], offset)
            curve["pts"] = curve["arc"]["pts"]

    curves = [c for row_index in range(len(rows)) for c in local_curves[row_index]]
    atoms = [
        {"key": (item["inst"].uid, number), "label": f"{item['inst'].species} atom {number}",
         "x": a["x"], "y": a["y"], "reach": _reach(a)}
        for row in rows for item in row if item["t"] == "frag" for number, a in enumerate(item["inst"].frag.atoms)
    ] + [o for row in rows for o in _label_obstacles(row)]
    bonds = [b for row in rows for b in _row_bonds(row)]

    def build(with_numbers: bool):
        def make(ids) -> str:
            parts: list[str] = []
            if title:
                xml, _ = R._build_text_element([str(title)], X0, FIRST_ROW_TOP - 4.0, ids, justification="Left",
                                               use_formatting=False)
                parts.append(xml)
            for row_index, row in enumerate(rows):
                y = offsets[row_index]
                for item in row:
                    if item["t"] == "frag":
                        frag = item["inst"].frag
                        if getattr(frag, "newman", None):
                            parts += _newman_xml(frag, ids.next)
                        else:
                            xml, _, _ = R._build_fragment(frag.atoms, frag.bonds, ids)
                            parts.append(xml)
                        if frag.ref.label:
                            label, _ = R._build_label_element(frag.ref.label, frag.cx, frag.bbox[3] + 14.0, ids)
                            parts.append(label)
                    elif item["t"] == "plus":
                        xml, _ = R._build_text_element(["+"], item["x"], y + 3.5, ids)
                        parts.append(xml)
                    else:
                        parts += _arrow_xml(item, y, ids)
            if with_numbers:
                for row in rows:
                    for item in row:
                        if item["t"] == "frag":
                            for number, a in enumerate(item["inst"].frag.atoms):
                                xml, _ = R._build_text_element([str(number)], a["x"] + 4.0, a["y"] - 4.0, ids,
                                                               justification="Left", use_formatting=False,
                                                               font_size=7)
                                parts.append(xml)
            parts += [_arc_xml(c, ids.next) for c in curves]
            return "\n".join(parts)
        return make

    # Shifted atoms are final; recompute fragment centres used for labels.
    for row in rows:
        for item in row:
            if item["t"] == "frag":
                frag = item["inst"].frag
                frag.cx, frag.cy = R._bbox_center(frag.bbox)

    cdxml = R.render_with(build(numbers))
    clean = R.render_with(build(False)) if numbers else cdxml
    return _Drawing(cdxml=cdxml, lint_cdxml=clean, rows=rows, curves=curves, atoms=atoms, bonds=bonds, notes=notes)


def _arrow_xml(item: dict, y: float, ids) -> list[str]:
    parts: list[str] = []
    style = item["style"]
    xml, _ = R._build_arrow(item["tail"], y, item["head"], y, ids, kind=R._arrow_kind(style))
    parts.append(xml)
    middle = (item["tail"] + item["head"]) / 2.0
    if item["above"]:
        below_baseline = (len(item["above"]) - 1) * R._LINE_ADVANCE + R._DESCENT
        xml, _ = R._build_text_element(item["above"], middle, y - LAYOUT_BELOW_GAP - below_baseline, ids)
        parts.append(xml)
    if item["below"]:
        xml, _ = R._build_text_element(item["below"], middle, y + LAYOUT_BELOW_GAP + R._CAP_HEIGHT, ids)
        parts.append(xml)
    return parts


def _color_moves(cdxml: Any, curves: list[dict]) -> Any:
    """Give each curved arrow the colour its move names; it wins over a house colour."""
    if not any(c.get("color") for c in curves):
        return cdxml
    root = ET.fromstring(str(cdxml).split("?>", 1)[-1])
    colors = figure_house._Colors(root)
    arcs = [a for a in root.iter("arrow") if a.get("AngularSize")]
    for arrow, curve in zip(arcs, curves):  # written in the same order
        if curve.get("color"):
            arrow.set("color", colors.index(curve["color"]))
    wrapped = R._ValidatedCDXML('<?xml version="1.0" encoding="UTF-8" ?>\n' + ET.tostring(root, encoding="unicode"))
    wrapped.chemistry_validation = getattr(cdxml, "chemistry_validation", None)
    return wrapped


def _move_report(brief: dict, mols: dict[str, Any]) -> tuple[list[dict], list[str]]:
    report: list[dict] = []
    warnings: list[str] = []
    for number, step in enumerate(brief["steps"], 1):
        for index, move in enumerate(step.get("moves") or [], 1):
            checks = _move_checks(move, mols[move["from"]["species"]])
            report.append({
                "step": number, "move": index,
                "from": _describe(move["from"], mols[move["from"]["species"]]),
                "to": _describe(move["to"], mols[move["to"]["species"]]),
                "electrons": move.get("electrons", 2),
                "checks": checks,
            })
            warnings += [f"Step {number} move {index}: {text}" for text in checks]
    return report, warnings


def draw_mechanism(brief: dict, output_path: str | None, render: bool) -> dict:
    """Draw a mechanism; called by draw_figure for ``kind: mechanism``."""
    errors = _validate(brief)
    if errors:
        return figure_draw._failure(errors, kind="mechanism")
    # A Newman projection is drawn from its substituents, not from a structure to ground.
    projections = [e for e in brief["species"] if "newman" in e]
    grounded, errors = figure_draw._ground([e for e in brief["species"] if "newman" not in e])
    if errors:
        return figure_draw._failure(errors, kind="mechanism")
    grounded += [{"id": str(e["id"]), "label": e.get("label"), "name": None, "smiles": None, "source": "newman",
                  "newman": e["newman"], "inventory": _newman_inventory(e["newman"])} for e in projections]
    mols = {g["id"]: Chem.MolFromSmiles(g["smiles"]) for g in grounded if g.get("smiles")}
    errors = _index_errors(brief, mols)
    if errors:
        return figure_draw._failure(errors, kind="mechanism")
    smiles = {g["id"]: g["smiles"] for g in grounded if g.get("smiles")}

    balance, warnings = figure_draw._check_balance(brief, smiles)
    moves, move_warnings = _move_report(brief, mols)
    warnings += move_warnings
    notes: list[str] = []
    used = {i for step in brief["steps"] for side in ("substrates", "products") for i in _as_ids(step[side])}
    for item in grounded:
        if item["id"] not in used:
            notes.append(f"species {item['id']} is not used in any step, so it is not drawn.")
    numbers = bool(brief.get("show_atom_numbers"))
    if numbers:
        warnings.append("The atom numbers are drawn so the moves can be checked; "
                        "remove show_atom_numbers for the final figure.")

    requested_width = float(brief.get("max_width") or figure_draw.DEFAULT_MAX_WIDTH)
    # Laid out at ACS geometry; the style preset rescales the saved file afterwards.
    slug, scale = figure_draw.brief_style(brief)
    budget = requested_width / scale
    width, gap_scale, bulge = budget, 1.0, 1.0
    try:
        drawing = _draw(brief, grounded, width, gap_scale, bulge, numbers)
        issues = figure_lint.lint_cdxml(drawing.lint_cdxml, max_width=budget)
        crossings = _crossings(drawing.curves, drawing.atoms, drawing.bonds)
        relayout = None
        if issues or crossings:
            first = {"issues": issues, "crossings": crossings}
            new_width, new_gap, changed, _ = figure_draw._relayout(issues, width, gap_scale)
            if crossings:
                if new_gap == gap_scale:
                    new_gap = round(gap_scale * figure_draw.RELAYOUT_GAP_SCALE, 2)
                    changed.append(f"widened the gaps between structures by x{figure_draw.RELAYOUT_GAP_SCALE:g}")
                bulge = RELAYOUT_BULGE
                changed.append(f"raised the arrow loops by x{RELAYOUT_BULGE:g}")
            if changed:
                drawing = _draw(brief, grounded, new_width, new_gap, bulge, numbers)
                issues = figure_lint.lint_cdxml(drawing.lint_cdxml, max_width=budget)
                crossings = _crossings(drawing.curves, drawing.atoms, drawing.bonds)
                relayout = {"changed": changed, "first_issues": first["issues"],
                            "first_crossings": first["crossings"]}
    except Exception as exc:
        return figure_draw._failure([f"The layout failed: {exc}"], kind="mechanism", warnings=warnings)

    try:
        ringed, ring_applied, _ = figure_house.apply_generic_rings(drawing.cdxml, brief, grounded)
    except ValueError as exc:
        return figure_draw._failure([f"The layout failed: {exc}"], kind="mechanism", warnings=warnings)
    if ring_applied:
        warnings += [f"Generic ring: {issue}." for issue in overlays_collide(str(ringed))]
    house = figure_draw.brief_house(brief)
    final = figure_draw.styled(ringed, slug)
    house_applied: list[str] = []
    house_skipped = ["the reaction-scheme conventions (labels, yields, abbreviations) apply to schemes only"]
    if house is not None:
        final, house_applied, _ = figure_house.apply_measured_geometry(final, house, slug)
        final, applied, more = figure_house.apply_figure_type_conventions(final, house, "mechanism", slug)
        house_applied += applied
        house_skipped += more
    final = _color_moves(final, drawing.curves)
    house_applied += ring_applied
    try:
        target = artifact_safety.resolve_destination(
            source=None, output_path=output_path, tag="figure", suffix=".cdxml")
        with artifact_safety.staging_file(target) as stage:
            stage.write_text(str(final), encoding="utf-8", newline="")
            artifact_safety.publish_file(stage, target)
    except ValueError as exc:  # includes OutputExistsError
        return figure_draw._failure([str(exc)], kind="mechanism", warnings=warnings)

    notes += drawing.notes
    for issue in issues:
        warnings.append(f"Layout check: {issue['message']}")
    for crossing in crossings:
        warnings.append(f"Step {crossing['step']} move {crossing['move']}: the arrow passes through "
                        f"{crossing['with']}; move the species apart or reorder them.")
    readback = figure_draw._readback(target, grounded, used, getattr(drawing.cdxml, "chemistry_validation", None))
    if readback["status"] != "matched":
        warnings.append(f"The saved file does not match the brief: missing {readback['missing']}, "
                        f"{readback['unexpected']} unexpected structure(s).")
    if readback["fragment_roundtrip"] != "preserved":
        warnings.append(f"Per-structure round-trip is {readback['fragment_roundtrip']}: "
                        "a structure may have changed while drawing.")
    notes.append("The curved arrows are drawn exactly as listed in moves; the tool checks lone pairs and atom "
                 "balance but does not decide whether the mechanism is chemically right.")

    native: dict[str, Any]
    produced = [target]
    if not render:
        native = {"status": "not_run", "reason": "native rendering was not requested"}
        notes.append("The drawing has not been rendered natively; look at a render before treating it as final.")
    else:
        png = figure_draw._png_destination(target)
        try:
            figure_draw._render_png(str(target), str(png))
            native = {"status": "rendered", "png_path": str(png)}
            produced.append(png)
        except Exception as exc:
            native = {"status": "not_run", "reason": str(exc)}
            notes.append("Native rendering did not run, so the drawing is an unrendered draft.")

    species = []
    for g in grounded:
        entry = {k: g[k] for k in ("id", "label", "name", "smiles", "source")}
        if g["id"] in mols:
            entry["atoms"] = ", ".join(f"{a.GetIdx()} {a.GetSymbol()}" for a in mols[g["id"]].GetAtoms())
        else:
            entry["newman"] = g["newman"]  # drawn from its substituents; no atoms to number
        species.append(entry)
    report = {
        "ok": True,
        "operation": "draw_figure",
        "kind": "mechanism",
        "output_path": str(target),
        "species": species,
        "steps": balance,
        "moves": moves,
        "arrows": {"count": len(drawing.curves), "crossings": crossings, "relayout": relayout},
        "style": figure_draw.style_report(slug, house, house_applied, house_skipped),
        "lint": {"issues": issues, "relayout": relayout, "max_width": requested_width},
        "readback": readback,
        "native": native,
        "checks_passed": (not issues and not crossings and readback["status"] == "matched"
                          and readback["fragment_roundtrip"] == "preserved"),
        "notes": notes,
        "warnings": warnings,
    }
    return artifact_safety.with_artifacts(report, produced)


__all__ = ["MECH_KEYS", "draw_mechanism"]
