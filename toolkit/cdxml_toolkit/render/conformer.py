"""Conformer templates as resolved fragments in scheme and mechanism layouts."""
from __future__ import annotations

import copy
import math
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Mapping

from ..constants import ACS_BOND_LENGTH
from ..templates.conformers import _charge_errors, _substituent, build_conformer
from .schema import StructureRef

TEMPLATES = ("chair", "newman", "haworth", "fischer", "norbornane", "bicyclo_222", "adamantane")
NEWMAN_MIN_APART = 10.0
NEWMAN_SPOKE = 1.4 * ACS_BOND_LENGTH
CONFORMER_FILL_COLOR = '<color r="0.8392" g="0.9176" b="0.9725"/>'
CONFORMER_FILL_INDEX = "9"

_TEMPLATE_KEYS = {
    "chair": {"template", "flip", "substituents", "bond_length"},
    "newman": {"template", "front", "rear", "conformation", "dihedral_offset", "attach", "charges", "bond_length"},
    "haworth": {"template", "anomer", "substituents", "bond_length"},
    "fischer": {"template", "centers", "top", "bottom", "bond_length"},
    "norbornane": {"template", "substituents", "bond_length"},
    "bicyclo_222": {"template", "substituents", "bond_length"},
    "adamantane": {"template", "substituents", "bond_length"},
}


def _label_errors(where: str, value: Any) -> list[str]:
    try:
        _substituent(value)
    except (TypeError, ValueError) as exc:
        return [f"{where}: {exc}."]
    return []


def _position_errors(where: str, value: Any, positions: range, directions: set[str]) -> list[str]:
    if not isinstance(value, Mapping):
        return [f"{where} must map ring positions to substituents."]
    errors: list[str] = []
    for position, entries in value.items():
        try:
            parsed = int(position)
        except (TypeError, ValueError):
            parsed = -1
        if parsed not in positions:
            errors.append(f"{where} positions must be in {positions.start}..{positions.stop - 1}, not {position!r}.")
            continue
        if not isinstance(entries, Mapping) or set(entries) - directions:
            errors.append(f"{where}[{position}] allows only {sorted(directions)}.")
            continue
        for direction, entry in entries.items():
            errors += _label_errors(f"{where}[{position}].{direction}", entry)
    return errors


def _newman_errors(where: str, spec: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    for side in ("front", "rear"):
        group = spec.get(side)
        if not isinstance(group, list) or len(group) != 3:
            errors.append(f"{where} {side} must list three substituents.")
        else:
            for index, entry in enumerate(group):
                errors += _label_errors(f"{where} {side}[{index}]", entry)
    if spec.get("conformation", "staggered") not in ("staggered", "eclipsed"):
        errors.append(f"{where} conformation must be staggered or eclipsed.")
    offset = spec.get("dihedral_offset")
    if offset is not None and (isinstance(offset, bool) or not isinstance(offset, (int, float))):
        errors.append(f"{where} dihedral_offset must be a number of degrees.")
    elif offset is not None:
        apart = abs((float(offset) + 60.0) % 120.0 - 60.0)
        if apart < NEWMAN_MIN_APART:
            errors.append(f"{where} dihedral_offset {offset:g} puts the rear spokes on top of the front ones.")
    front = spec.get("front") if isinstance(spec.get("front"), list) else []
    rear = spec.get("rear") if isinstance(spec.get("rear"), list) else []
    if len(front) == len(rear) == 3 and spec.get("charges") is not None:
        try:
            _charge_errors(spec["charges"], front, rear)
        except (TypeError, ValueError) as exc:
            errors.append(f"{where} charges: {exc}.")
    attach = spec.get("attach")
    labels = []
    for entry in front:
        try:
            labels.append(_substituent(entry)[0])
        except (TypeError, ValueError):
            pass
    if attach is None:
        if "*" in labels:
            errors.append(f'{where} attach is required when a front substituent is "*".')
    else:
        if not isinstance(attach, Mapping) or set(attach) != {"smiles"} or not isinstance(attach.get("smiles"), str):
            errors.append(f'{where} attach must be {{"smiles": "..."}} with one [*].')
        elif attach["smiles"].count("*") != 1:
            errors.append(f"{where} attach SMILES must contain exactly one [*].")
        if labels.count("*") != 1:
            errors.append(f'{where} attach replaces exactly one front substituent written "*".')
        try:
            if "*" in [_substituent(entry)[0] for entry in rear]:
                errors.append(f"{where} attach meets the front carbon only.")
        except (TypeError, ValueError):
            pass
    return errors


def conformer_errors(ident: str, spec: Any) -> list[str]:
    where = f"Species {ident}: conformer"
    if not isinstance(spec, Mapping):
        return [f"{where} must be an object with a template."]
    template = str(spec.get("template") or "").strip().lower()
    if template not in TEMPLATES:
        return [f"{where} template must be one of {list(TEMPLATES)}, not {spec.get('template')!r}."]
    errors: list[str] = []
    extra = sorted(set(spec) - _TEMPLATE_KEYS[template])
    if extra:
        errors.append(f"{where} has unknown keys {extra}.")
    length = spec.get("bond_length")
    if length is not None and (isinstance(length, bool) or not isinstance(length, (int, float)) or length <= 0):
        errors.append(f"{where} bond_length must be a positive number of points.")
    if template == "chair":
        if not isinstance(spec.get("flip", False), bool):
            errors.append(f"{where} flip must be true or false.")
        if spec.get("substituents") is not None:
            errors += _position_errors(f"{where} substituents", spec["substituents"], range(1, 7),
                                       {"ax", "axial", "eq", "equatorial"})
    elif template == "newman":
        errors += _newman_errors(where, spec)
    elif template == "haworth":
        if spec.get("anomer", "beta") not in ("alpha", "beta"):
            errors.append(f"{where} anomer must be alpha or beta.")
        if spec.get("substituents") is not None:
            errors += _position_errors(f"{where} substituents", spec["substituents"], range(1, 6), {"up", "down"})
    elif template == "fischer":
        centers = spec.get("centers", [{"left": "H", "right": "OH"}])
        if not isinstance(centers, list) or not centers:
            errors.append(f"{where} centers must be a non-empty list.")
        else:
            for index, center in enumerate(centers):
                if not isinstance(center, Mapping) or set(center) - {"left", "right"}:
                    errors.append(f"{where} centers[{index}] allows only left and right.")
                    continue
                for direction in ("left", "right"):
                    if direction in center:
                        errors += _label_errors(f"{where} centers[{index}].{direction}", center[direction])
        for direction in ("top", "bottom"):
            if direction in spec:
                errors += _label_errors(f"{where} {direction}", spec[direction])
    elif spec.get("substituents") is not None:
        size = {"norbornane": 7, "bicyclo_222": 8, "adamantane": 10}[template]
        substituents = spec["substituents"]
        if not isinstance(substituents, Mapping):
            errors.append(f"{where} substituents must map ring positions to labels.")
        else:
            for position, entry in substituents.items():
                try:
                    parsed = int(position)
                except (TypeError, ValueError):
                    parsed = -1
                if parsed not in range(1, size + 1):
                    errors.append(f"{where} substituent positions must be in 1..{size}, not {position!r}.")
                else:
                    errors += _label_errors(f"{where} substituents[{position}]", entry)
    return errors


_ATTACHED: dict[str, dict] = {}


def attach_group(smiles: str) -> dict:
    if smiles in _ATTACHED:
        return _ATTACHED[smiles]
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"the attached group {smiles!r} is not valid SMILES")
    stars = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(stars) != 1 or stars[0].GetDegree() != 1:
        raise ValueError(f"the attached group {smiles!r} needs one [*] bonded to one atom")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError(f"the attached group {smiles!r} carries a radical, which a projection cannot draw yet")
    if any(abs(a.GetFormalCharge()) > 1 for a in mol.GetAtoms()):
        raise ValueError(f"the attached group {smiles!r} has a charge beyond plus or minus one")
    Chem.Kekulize(mol, clearAromaticFlags=True)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = {a.GetIdx(): (conf.GetAtomPosition(a.GetIdx()).x, -conf.GetAtomPosition(a.GetIdx()).y)
          for a in mol.GetAtoms()}
    lengths = [math.dist(xy[b.GetBeginAtomIdx()], xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    scale = ACS_BOND_LENGTH / (sum(lengths) / len(lengths))
    star = stars[0].GetIdx()
    root = stars[0].GetNeighbors()[0].GetIdx()
    angle = math.atan2(xy[root][1] - xy[star][1], xy[root][0] - xy[star][0])
    cos, sin = math.cos(-angle), math.sin(-angle)
    from ..templates.conformers import _HYDRIDE_LABELS, _SIMPLE_ELEMENTS

    atoms = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() == star:
            continue
        symbol, hydrogens = atom.GetSymbol(), atom.GetTotalNumHs()
        label = symbol if symbol == "C" or not hydrogens else symbol + "H" + (str(hydrogens) if hydrogens > 1 else "")
        if symbol not in _SIMPLE_ELEMENTS or (label != symbol and label not in _HYDRIDE_LABELS):
            raise ValueError(f"the attached group {smiles!r} has an atom a projection cannot label ({label})")
        dx, dy = ((xy[atom.GetIdx()][axis] - xy[root][axis]) * scale for axis in (0, 1))
        atoms.append({"index": atom.GetIdx() + 1, "symbol": label, "charge": atom.GetFormalCharge(),
                      "num_hydrogens": hydrogens, "show_label": symbol != "C",
                      "x": dx * cos - dy * sin, "y": dx * sin + dy * cos})
    bonds = [{"atom1": b.GetBeginAtomIdx() + 1, "atom2": b.GetEndAtomIdx() + 1,
              "order": int(b.GetBondTypeAsDouble())}
             for b in mol.GetBonds() if star not in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())]
    _ATTACHED[smiles] = {"root": root + 1, "atoms": atoms, "bonds": bonds}
    return _ATTACHED[smiles]


def drawing(spec: Mapping[str, Any], cx: float = 0.0, cy: float = 0.0) -> dict[str, Any]:
    normalized = copy.deepcopy(dict(spec))
    normalized.setdefault("template", "newman")
    normalized["cx"], normalized["cy"] = float(cx), float(cy)
    normalized.setdefault("bond_length", NEWMAN_SPOKE if normalized["template"] == "newman" else ACS_BOND_LENGTH)
    attach = normalized.get("attach")
    if isinstance(attach, Mapping) and "smiles" in attach:
        normalized["attach"] = attach_group(attach["smiles"])
    return build_conformer(normalized)


def fragment(ident: str, spec: Mapping[str, Any], label: str | None) -> Any:
    from . import renderer as R

    from ..label_geometry import label_box

    normalized = copy.deepcopy(dict(spec))
    normalized.setdefault("template", "newman")
    drawn = drawing(normalized)
    atoms = [dict(a) for a in drawn["atoms"]]
    for atom in atoms:
        reach = list(atom.get("reach") or ())
        if atom.get("charge_box"):
            x1, y1, x2, y2 = atom.pop("charge_box")
            reach.extend([(x1 - atom["x"], y1 - atom["y"]), (x2 - atom["x"], y2 - atom["y"])])
        sym = str(atom.get("symbol", "C"))
        if atom.get("show_label", sym != "C"):
            bx1, by1, bx2, by2 = label_box(sym, atom["x"], atom["y"], str(atom.get("justification", "Left")))
            reach.extend([(bx1 - atom["x"], by1 - atom["y"]), (bx2 - atom["x"], by2 - atom["y"])])
        if reach:
            atom["reach"] = reach
    ref = StructureRef(id=ident, label=label, conformer=normalized)
    frag = R.ResolvedFragment(ref=ref, atoms=atoms, bonds=[dict(b) for b in drawn["bonds"]])
    frag.bbox = R._fragment_bbox(atoms)
    frag.cx, frag.cy = R._bbox_center(frag.bbox)
    frag.conformer = normalized
    frag.conformer_offset = (drawn["atoms"][0]["x"] - drawn["cx"], drawn["atoms"][0]["y"] - drawn["cy"])
    if normalized["template"] == "newman":
        frag.newman = dict(normalized)
    return frag


def inventory(spec: Mapping[str, Any]) -> list[str]:
    from ..chemistry_semantics import document_inventory

    drawn = drawing(spec)
    with tempfile.TemporaryDirectory(prefix="conformer-") as scratch:
        path = Path(scratch) / "conformer.cdxml"
        path.write_text(drawn["cdxml"], encoding="utf-8")
        try:
            return list(document_inventory(path).elements())
        except ValueError:
            return []


def _receipt(frag: Any, fragment_element: ET.Element, renumber: Mapping[str, str]) -> None:
    from . import renderer as R

    collector = R._fragment_receipts.get()
    if collector is None:
        return
    fid = renumber.get(fragment_element.get("id"), fragment_element.get("id"))
    generic = [n.get("GenericNickname") for n in fragment_element.iter("n") if n.get("NodeType") == "GenericNickname"]
    if generic:
        collector.append({"fragment_id": fid, "status": "not_checkable",
                          "reason": f"generic label {', '.join(sorted(set(generic)))} cannot be read back"})
        return
    from ..chemistry_semantics import semantic_key
    from ..publication.layout import fragment_mol

    try:
        drawn = semantic_key(fragment_mol(fragment_element))
    except Exception:
        collector.append({"fragment_id": fid, "status": "unreadable"})
        return
    expected = inventory(frag.conformer)
    collector.append({"fragment_id": fid, "status": "preserved" if drawn in expected else "changed",
                      "method": "conformer_template_readback"})


def xml(frag: Any, next_id) -> list[str]:
    dx, dy = frag.conformer_offset
    centre = frag.atoms[0]
    drawn = drawing(frag.conformer, centre["x"] - dx, centre["y"] - dy)
    elements = [ET.fromstring(drawn["fragment_xml"])] + [ET.fromstring(g) for g in drawn["graphics_xml"]]
    every = [node for element in elements for node in element.iter()]
    renumber = {node.get("id"): str(next_id()) for node in every if node.get("id")}
    _receipt(frag, elements[0], renumber)
    for node in every:
        if node.get("id"):
            node.set("id", renumber[node.get("id")])
        for key in ("B", "E", "object"):
            if node.get(key) in renumber:
                node.set(key, renumber[node.get(key)])
        if node.get("CrossingBonds"):
            node.set("CrossingBonds", " ".join(renumber[value] for value in node.get("CrossingBonds").split()))
    for node in sorted((n for n in every if n.get("Z")), key=lambda n: int(n.get("Z"))):
        node.set("Z", str(next_id()))
    for element in elements[1:]:
        if element.get("OvalType") == "Circle Filled":
            element.set("color", CONFORMER_FILL_INDEX)
        elif element.get("OvalType") == "Circle":
            element.set("color", "0")
    return [ET.tostring(element, encoding="unicode") for element in elements]


def palette_xml(inner_xml: str) -> str:
    return CONFORMER_FILL_COLOR if f'OvalType="Circle Filled" color="{CONFORMER_FILL_INDEX}"' in inner_xml else ""
