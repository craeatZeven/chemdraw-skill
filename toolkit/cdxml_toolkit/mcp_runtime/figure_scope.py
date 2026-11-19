"""Substrate-scope panels for draw_figure: a list of compounds with labels and values.

The brief carries the compounds and the numbers the user gave. This module
grounds each structure, writes the table the publication layout reads and lets
``publication_figure`` pack, label and check the pages. A value the brief does
not contain is drawn as a dash and named in the notes; nothing is filled in.
"""
from __future__ import annotations

import math

import copy
import csv
import json
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem

from ..chemistry_semantics import document_inventory, semantic_key
from . import artifact_safety
from . import figure_house
from .figure_draw import (
    _as_list, _failure, _ground, _render_png, brief_house, brief_style, style_preset_error, style_report,
)
from ..chemdraw.style_presets import DEFAULT_STYLE_PRESET, style_preset_as_manifest_style

SCOPE_KEYS = {"kind", "species", "fields", "conditions", "footnotes", "title", "style", "width_mm", "height_mm",
              "highlight_varying", "header", "columns", "depiction_rotation", "rmark_in_ring"}
HEADER_KEYS = {"substrates", "products", "above", "below"}
HIGHLIGHT = "#FF6600"  # the house orange, already in the default check palette
SPECIES_STYLE = ("color", "atom_colors", "bond_colors", "abbreviate")
SPECIES_KEYS = {"id", "smiles", "name", "label", "group", "values", "note", "definition", "secondary_of",
                "drawn_with", *SPECIES_STYLE}
NOTE_COLUMN = "note"
# A line naming what a structure stands for ("Alk = C7H15"), drawn under its number.
DEFINITION_COLUMN = "definition"
SECONDARY_COLUMN = "secondary_of"
# A compound listed under another's structure (BJOC 22-127: 3q under 3p's NH-Alk).
DRAWN_WITH_COLUMN = "drawn_with"
MAX_SPECIES = 1000


def _fields(brief: dict) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    given = brief.get("fields")
    if given is None:
        keys = list(dict.fromkeys(k for e in brief["species"] if isinstance(e, dict) for k in (e.get("values") or {})))
        return [{"key": k, "label": k, "unit": None, "show_label": False, "hide_label": False, "optional": False}
                for k in keys], errors
    if not isinstance(given, list):
        return [], ["fields must be a list of keys or {key, label, unit} objects."]
    fields: list[dict] = []
    for item in given:
        if isinstance(item, str):
            item = {"key": item}
        if not isinstance(item, dict) or not str(item.get("key") or "").strip():
            errors.append("Every field needs a key.")
            continue
        key = str(item["key"])
        if any(f["key"] == key for f in fields):
            errors.append(f"Field {key}: duplicate key.")
            continue
        label = item.get("label", key)
        for flag in ("show_label", "optional"):
            if flag in item and not isinstance(item[flag], bool):
                errors.append(f"Field {key}: {flag} must be true or false.")
        # show_label false keeps a field's line unnamed even beside other fields; optional
        # leaves the line out where the value is empty instead of drawing a dash.
        fields.append({"key": key, "label": "" if label is None else str(label), "unit": item.get("unit") or None,
                       "show_label": item.get("show_label") is True, "hide_label": item.get("show_label") is False,
                       "optional": item.get("optional") is True,
                        "format": item.get("format") or None})
    return fields, errors


def _header_errors(header) -> list[str]:
    """Problems with a scope's header scheme (the general reaction over the grid)."""
    if not isinstance(header, dict):
        return ["header must be an object: {substrates, products, above?, below?}."]
    errors = []
    unknown = sorted(set(header) - HEADER_KEYS)
    if unknown:
        errors.append(f"header has unknown keys {unknown}; use {sorted(HEADER_KEYS)}.")
    for side in ("substrates", "products"):
        items = header.get(side)
        if not isinstance(items, list) or not items:
            errors.append(f"header {side}: list at least one structure, each {{smiles or name, label?}}.")
            continue
        for number, item in enumerate(items, 1):
            if not isinstance(item, dict) or not (item.get("smiles") or item.get("name")):
                if isinstance(item, dict) and item.get('text'):
                    # A reagent written as its abbreviation ("S8") rather than drawn.
                    unknown = sorted(set(item) - {'text', 'label'})
                    if unknown:
                        errors.append(f"header {side}[{number}] text item has unknown keys {unknown}.")
                    continue
                errors.append(f"header {side}[{number}] needs a smiles, a name or a text; write R groups as [*:1], [*:2].")
            elif set(item) - {"smiles", "name", "label", "rotation"}:
                errors.append(f"header {side}[{number}] has unknown keys {sorted(set(item) - {'smiles', 'name', 'label', 'rotation'})}.")
    for side in ("above", "below"):
        if side in header and not (isinstance(header[side], list) and all(isinstance(t, str) for t in header[side])):
            errors.append(f"header {side} must be a list of condition lines.")
    return errors


def _validate(brief: dict) -> list[str]:
    errors: list[str] = []
    unknown = sorted(set(brief) - SCOPE_KEYS - figure_house.RECORD_KEYS)
    if unknown:
        errors.append(f"Unknown brief keys: {unknown}")
    errors += figure_house.record_errors(brief)
    if brief.get("header") is not None:
        errors += _header_errors(brief["header"])
    style_error = style_preset_error(brief.get("style"))
    if style_error:
        errors.append(style_error)
    for entry in brief.get("species") or []:
        if isinstance(entry, dict):
            errors += figure_house.species_style_errors(entry)
    for key, low, high in (("width_mm", 60, 400), ("height_mm", 60, 500), ("depiction_rotation", -360, 360)):
        value = brief.get(key)
        if value is not None and not (isinstance(value, (int, float)) and low <= value <= high):
            errors.append(f"{key} must be a number from {low} to {high} millimetres.")
    columns = brief.get("columns")
    if columns is not None and not (isinstance(columns, int) and not isinstance(columns, bool) and 1 <= columns <= 8):
        errors.append("columns must be a whole number from 1 to 8: the compounds to a row.")
    footnotes = brief.get("footnotes")
    if footnotes is not None and not (isinstance(footnotes, list) and all(isinstance(f, str) for f in footnotes)):
        errors.append("footnotes must be a list of text lines.")
    for key in ("title", "conditions"):
        if brief.get(key) is not None and not isinstance(brief[key], str):
            errors.append(f"{key} must be text.")
    highlight = brief.get("highlight_varying")
    if highlight is not None and not (highlight is True or highlight is False or (
            isinstance(highlight, str) and re.fullmatch(r"#[0-9A-Fa-f]{6}", highlight))):
        errors.append("highlight_varying must be true, false or a hex colour like #FF6600.")

    species = brief.get("species")
    if not isinstance(species, list) or not species:
        return errors + ["species: list every compound, each with an id and a smiles or a name."]
    if len(species) > MAX_SPECIES:
        errors.append(f"A scope panel holds at most {MAX_SPECIES} compounds.")
    ids: set[str] = set()
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
        if not entry.get("smiles") and not entry.get("name"):
            errors.append(f"Species {ident} needs a smiles or a name.")
        if not isinstance(entry.get("values") or {}, dict):
            errors.append(f"Species {ident}: values must be an object of field -> value.")
        if entry.get("secondary_of") is not None:
            every_id = {str(e.get("id")) for e in brief.get("species") or [] if isinstance(e, dict)}
            main = str(entry["secondary_of"])
            if main == str(ident) or main not in every_id:
                errors.append(f"Species {ident}: secondary_of must name another species of this brief.")
        if entry.get("drawn_with") is not None:
            by_id = {str(e.get("id")): e for e in brief.get("species") or [] if isinstance(e, dict)}
            main = by_id.get(str(entry["drawn_with"]))
            if main is None or main is entry:
                errors.append(f"Species {ident}: drawn_with must name another species of this brief.")
            elif main.get("drawn_with") is not None or main.get("secondary_of") is not None:
                errors.append(f"Species {ident}: drawn_with must name a species drawn with its own structure.")
            elif entry.get("secondary_of") is not None:
                errors.append(f"Species {ident}: a species is a secondary product or drawn with another, not both.")
            elif _text(main.get("group")) != _text(entry.get("group")):
                errors.append(f"Species {ident}: drawn_with must name a species of the same group.")
        if entry.get("definition") is not None and not isinstance(entry["definition"], str):
            errors.append(f'Species {ident}: definition must be text (a line such as "Alk = C7H15").')
        if entry.get("note") is not None and not isinstance(entry["note"], str):
            errors.append(f"Species {ident}: note must be text (a line under the compound, e.g. 'not formed').")

    fields, field_errors = _fields(brief)
    errors += field_errors
    if any(f["key"] == NOTE_COLUMN for f in fields):
        errors.append(f"{NOTE_COLUMN!r} is reserved for the per-compound note; name the field otherwise.")
    if any(f["key"] == DEFINITION_COLUMN for f in fields):
        errors.append(f"{DEFINITION_COLUMN!r} is reserved for the per-compound definition line; name the field "
                      "otherwise.")
    known = {f["key"] for f in fields}
    for entry in species:
        if isinstance(entry, dict) and isinstance(entry.get("values") or {}, dict):
            for key in entry.get("values") or {}:
                if key not in known:
                    errors.append(f"Species {entry.get('id')}: a value for {key!r}, which is not in fields.")
    return errors


def _free_directory() -> Path:
    base = Path(tempfile.gettempdir()) / "scope-figure"
    candidate, index = base, 2
    while candidate.exists():
        candidate = base.with_name(f"{base.name}_{index}")
        index += 1
    return candidate


def _text(value: Any) -> str:
    return "" if value is None else str(value)


def _varying_atoms(grounded: list[dict]) -> tuple[dict[str, list[int]], str | None]:
    """For each compound, the atoms outside the core all of them share (their SMILES indices).

    The core is the maximum common substructure of every compound, matched by element and
    bond order, whole rings only, so a ring is never split into a highlighted half.
    """
    from rdkit.Chem import rdFMCS
    mols = {g["id"]: Chem.MolFromSmiles(g["smiles"]) for g in grounded}
    if len(mols) < 2:
        return {}, "highlight_varying needs at least two compounds to compare."
    found = rdFMCS.FindMCS(list(mols.values()), timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True,
                           atomCompare=rdFMCS.AtomCompare.CompareElements,
                           bondCompare=rdFMCS.BondCompare.CompareOrder)
    if found.canceled or found.numAtoms < 3:
        return {}, ("highlight_varying found no core of three or more atoms that every compound shares, "
                    "so nothing was highlighted.")
    core = Chem.MolFromSmarts(found.smartsString)
    matched = {ident: set(mol.GetSubstructMatch(core)) for ident, mol in mols.items()}
    # One awkard compound must not shrink the core for everyone: an aryl ring shared by a clear
    # majority (two thirds) of the compounds belongs to the scaffold, so a table whose last rows
    # are alkyl amides still highlights only the substituent of its aryl rows.
    matched = _extend_core_by_shared_rings(mols, matched)
    varying = {}
    for ident, mol in mols.items():
        keep = matched[ident]
        varying[ident] = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in keep]
    return varying, None


def _ring_systems(mol) -> list[set[int]]:
    """Atom index sets of a molecule's fused rings, one set per ring system."""
    systems: list[set[int]] = []
    for ring in mol.GetRingInfo().AtomRings():
        merged = set(ring)
        for other in [s for s in systems if s & merged]:
            merged |= other
            systems.remove(other)
        systems.append(merged)
    return systems


def _ring_query(mol, atoms: set[int]):
    """(SMARTS, query) for one ring system, so the same ring can be counted in the others."""
    bonds = [b.GetIdx() for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms]
    if not bonds:
        return None, None
    smart = Chem.MolToSmarts(Chem.PathToSubmol(mol, bonds))
    query = Chem.MolFromSmarts(smart)
    return (smart, query) if query is not None else (None, None)


def _extend_core_by_shared_rings(mols: dict, matched: dict, share: float = 2 / 3) -> dict:
    """Add every ring system most compounds share to each compound's core atoms."""
    per_mol = {}
    for ident, mol in mols.items():
        systems = []
        for atoms in _ring_systems(mol):
            if len(atoms) < 6:
                continue
            smart, query = _ring_query(mol, atoms)
            if query is not None:
                systems.append((atoms, smart, query))
        per_mol[ident] = systems
    candidates = {}
    for systems in per_mol.values():
        for _, smart, query in systems:
            candidates[smart] = query
    wanted = max(2, math.ceil(len(mols) * share))
    shared = {smart for smart, query in candidates.items()
              if sum(1 for mol in mols.values() if mol.HasSubstructMatch(query)) >= wanted}
    if not shared:
        return matched
    out = {}
    for ident, mol in mols.items():
        keep = set(matched[ident])
        for atoms, smart, _ in per_mol[ident]:
            if smart in shared:
                keep |= atoms
        out[ident] = keep
    return out


def _highlight_brief(brief: dict, grounded: list[dict]) -> tuple[dict, dict[str, int], list[str]]:
    """The brief with the varying group of each compound coloured, as atom and bond colours.

    A compound that sets its own color, atom_colors or bond_colors keeps them. The bond that
    joins a group to the core is coloured with it, so a methyl drawn as a bare line shows.
    """
    wanted = brief.get("highlight_varying")
    if not wanted:
        return brief, {}, []
    colour = wanted if isinstance(wanted, str) else HIGHLIGHT
    varying, problem = _varying_atoms(grounded)
    if problem:
        return brief, {}, [problem]
    smiles = {g["id"]: g["smiles"] for g in grounded}
    styled = copy.deepcopy(brief)
    counts: dict[str, int] = {}
    for entry in styled["species"]:
        ident = str(entry["id"])
        atoms = varying.get(ident) or []
        if not atoms or any(k in entry for k in ("color", "atom_colors", "bond_colors")):
            continue
        mol = Chem.MolFromSmiles(smiles[ident])
        chosen = set(atoms)
        entry["atom_colors"] = [{"atoms": atoms, "color": colour}]
        entry["bond_colors"] = [{"bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()
                                           if {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} & chosen],
                                 "color": colour}]
        counts[ident] = len(atoms)
    return styled, counts, []


def _cell_fields(fields: list[dict], notes: bool, definitions: bool = False) -> list[dict]:
    """The lines under each structure, as the publication layout reads them.

    One field reads as the corpus scope tables do, on the number's line ("3a, 92%", the
    number bold); say what the number is in a footnote or the conditions, or set
    ``show_label``. With several fields each has its own centred line, naming its field. A note, when any
    compound has one, is an extra unlabelled line that is left out where it is empty. A
    definition ("Alk = C7H15") joins the number's line: "3p, Alk = C7H15, 63%".
    A field with ``show_label`` false is unnamed even beside others; an ``optional`` one is
    left out where its value is empty rather than drawn as a dash.
    """
    columns = {f["key"]: f"v{i}" for i, f in enumerate(fields, 1)}
    compact = len(fields) == 1 and not fields[0].get("show_label")
    cells = [{"column": columns[f["key"]], "label": "" if compact or f.get("hide_label") else f["label"],
              **({"unit": f["unit"]} if f["unit"] else {}), **({"optional": True} if f.get("optional") else {}),
              **({"format": f["format"]} if f.get("format") else {}),
              **({"inline": True} if compact else {})}
             for f in fields]
    if definitions:
        cells.insert(0, {"column": DEFINITION_COLUMN, "label": "", "optional": True, "inline": True})
    if notes:
        cells.append({"column": NOTE_COLUMN, "label": "", "optional": True})
    return cells


def _ground_header(header: dict | None) -> tuple[dict | None, list[str]]:
    """The header with every structure as a checked SMILES (names resolved like species)."""
    if not header:
        return None, []
    out = {k: list(header.get(k) or []) for k in ("above", "below")}
    errors: list[str] = []
    for side in ("substrates", "products"):
        items = [{"id": f"header {side} {n}", **item} for n, item in enumerate(header[side], 1)]
        # A reagent written as its abbreviation ("S8") is not grounded: it has no structure to check.
        structural = [i for i in items if not i.get("text")]
        grounded, more = _ground([{k: v for k, v in i.items() if k != "label"} for i in structural])
        errors += more
        rebuilt = []
        queue = list(grounded)
        for item in items:
            if item.get("text"):
                rebuilt.append({"text": item["text"], **({"label": str(item["label"])} if item.get("label") else {})})
            elif queue:
                g = queue.pop(0)
                entry = {"smiles": g["smiles"], **({"label": str(item["label"])} if item.get("label") else {})}
                if item.get("rotation") is not None:
                    entry["rotation"] = item["rotation"]
                rebuilt.append(entry)
        out[side] = rebuilt
    return out, errors


def _write_inputs(work: Path, brief: dict, grounded: list[dict], fields: list[dict], render: bool,
                  header: dict | None = None, header_mark_color: str | None = None) -> Path:
    columns = {f["key"]: f"v{i}" for i, f in enumerate(fields, 1)}
    notes = any(isinstance(e, dict) and e.get("note") for e in brief["species"])
    definitions = any(isinstance(e, dict) and e.get("definition") for e in brief["species"])
    secondary = any(isinstance(e, dict) and e.get("secondary_of") for e in brief["species"])
    family = any(isinstance(e, dict) and e.get("drawn_with") for e in brief["species"])
    with (work / "table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["compound_id", "smiles", "label", "group", *columns.values(), *([NOTE_COLUMN] if notes else []),
                         *([DEFINITION_COLUMN] if definitions else []),
                         *([SECONDARY_COLUMN] if secondary else []),
                         *([DRAWN_WITH_COLUMN] if family else [])])
        by_id = {str(e["id"]): e for e in brief["species"]}
        for item in grounded:
            entry = by_id[item["id"]]
            values = entry.get("values") or {}
            writer.writerow([item["id"], item["smiles"], _text(entry.get("label")), _text(entry.get("group")),
                             *(_text(values.get(f["key"])) for f in fields),
                             *([_text(entry.get("note"))] if notes else []),
                             *([_text(entry.get("definition"))] if definitions else []),
                             *([_text(entry.get("secondary_of"))] if secondary else []),
                             *([_text(entry.get("drawn_with"))] if family else [])])
    style = {k: brief[k] for k in ("width_mm", "height_mm") if brief.get(k) is not None}
    slug = brief_style(brief)[0]
    if slug != DEFAULT_STYLE_PRESET:
        # The panel keeps its own text size; the preset sets the structures' geometry and font.
        preset = style_preset_as_manifest_style(slug)
        preset.pop("font_size")
        style.update(preset)
    spec = {
        "mode": "scope", "table": "table.csv", "label_column": "label", "group_column": "group",
        **({"secondary_column": SECONDARY_COLUMN} if secondary else {}),
        **({"drawn_with_column": DRAWN_WITH_COLUMN} if family else {}),
        # A definition such as "Alk = C7H15" also draws that group as its label.
        **({"label_groups_from": DEFINITION_COLUMN} if definitions else {}),
        "fields": _cell_fields(fields, notes, definitions),
        # An explicit title wins, including "" which asks for no title at all (a journal figure
        # carries none; the caption does that job).
        "title": brief["title"] if "title" in brief else "Substrate scope",
        "conditions": brief.get("conditions") or "",
        "footnotes": brief.get("footnotes") or [],
        "style": style,
        "native_render": render,
        **({"header": header} if header else {}),
        **({"columns": brief["columns"]} if brief.get("columns") else {}),
        # Turns the reference depiction of the whole table, so a benzothiazole can sit on the left
        # where the source figure has it rather than on the right where RDKit puts it.
        **({"depiction_rotation": brief["depiction_rotation"]} if brief.get("depiction_rotation") else {}),
        **({"rmark_in_ring": True} if brief.get("rmark_in_ring") else {}),
        **({"header_mark_color": header_mark_color} if header_mark_color else {}),
    }
    path = work / "spec.json"
    path.write_text(json.dumps(spec, ensure_ascii=False), encoding="utf-8")
    return path


def _readback(pages: list[Path], grounded: list[dict], header: dict | None = None,
              drawn_with: set[str] | None = None) -> dict:
    wanted: Counter = Counter()
    for item in grounded:
        if item["id"] in (drawn_with or set()):
            continue  # shown by the structure it is drawn with
        wanted.update(semantic_key(m) for m in Chem.GetMolFrags(Chem.MolFromSmiles(item["smiles"]), asMols=True))
    if header:
        from ..publication.layout import header_identities
        wanted.update(header_identities({"header": header}))
    found: Counter = Counter()
    for page in pages:
        found.update(document_inventory(page, like=wanted))
    missing = sorted((wanted - found).elements())
    extra = sorted((found - wanted).elements())
    return {
        "status": "matched" if not missing and not extra else "mismatch",
        "structures_drawn": sum(found.values()),
        "missing": len(missing),
        "unexpected": len(extra),
    }


def draw_scope(brief: dict, output_path: str | None, render: bool) -> dict:
    """Draw a scope panel; called by draw_figure for ``kind: scope``."""
    from . import publication_figure as publication

    errors = _validate(brief)
    if errors:
        return _failure(errors, kind="scope")
    grounded, errors = _ground(brief["species"])
    if errors:
        return _failure(errors, kind="scope")
    header, errors = _ground_header(brief.get("header"))
    if errors:
        return _failure(errors, kind="scope")
    fields, _ = _fields(brief)
    # Species listed under another's structure have no fragment of their own on the page.
    members = {str(e["id"]) for e in brief["species"] if e.get("drawn_with") is not None}
    drawn = [g for g in grounded if g["id"] not in members]

    target = Path(output_path).expanduser().resolve() if output_path else _free_directory()
    if target.exists():
        return _failure([f"Refusing to overwrite an existing folder: {target}"], kind="scope")

    house = brief_house(brief)
    try:
        styled, highlighted, highlight_warnings = _highlight_brief(brief, drawn)
    except Exception as exc:  # a failed comparison leaves the panel uncoloured, and says so
        styled, highlighted, highlight_warnings = brief, {}, [f"highlight_varying could not compare the compounds: {exc}"]
    restyle = house is not None or any(isinstance(e, dict) and any(k in e for k in SPECIES_STYLE)
                                       for e in styled["species"])
    wanted = brief.get("highlight_varying")
    # The source figures draw the general reaction's R marks in the same black as the rest of the
    # scheme; painting them in the highlight colour adds a second visual focus. An explicit
    # `header_mark_color` on the brief still wins, and so does an explicit highlight colour.
    mark_color = (wanted if isinstance(wanted, str) else HIGHLIGHT) if header and highlighted else None
    try:
        with tempfile.TemporaryDirectory(prefix="scope-inputs-") as scratch:
            # A page that is edited afterwards is rendered after the edit, not by the builder.
            spec = _write_inputs(Path(scratch), brief, grounded, fields, render and not restyle, header, mark_color)
            result = publication.publication_figure("create", str(spec), str(target))
    except Exception as exc:
        return _failure([f"The scope panel could not be drawn: {exc}"], kind="scope")
    if not result.get("ok"):
        return _failure([f"The panel has {result.get('conflicts', 'some')} data-binding conflict(s); "
                         f"see {result.get('output_dir')}."], kind="scope")

    project = json.loads((target / "project.json").read_text(encoding="utf-8"))
    per_page = Counter(c["page"] for c in project["compounds"].values())
    pages = [{"path": str(target / p["file"]), "file": p["file"], "compounds": per_page[i],
              **({"grid": p["grid"]} if p.get("grid") else {})}
             for i, p in enumerate(project["pages"])]
    checks = result["checks"]
    applied: list[str] = []
    skipped: list[str] = []
    style_warnings: list[str] = []
    if restyle:
        slug = brief_style(brief)[0]
        for page in pages:
            path = Path(page["path"])
            text, page_applied, skipped, warned = figure_house.apply_scope_conventions(
                path.read_text(encoding="utf-8"), styled, drawn, house, slug)
            path.write_text(str(text), encoding="utf-8")
            applied += [a for a in page_applied if a not in applied]
            style_warnings += [w for w in warned if w not in style_warnings]
    if mark_color:
        # The general reaction's R marks take the colour of what varies below (BJOC 22-127).
        for page in pages:
            path = Path(page["path"])
            text, painted = figure_house.paint_header_marks(path.read_text(encoding="utf-8"), mark_color)
            if painted:
                path.write_text(str(text), encoding="utf-8")
                applied.append("header_marks_highlighted")
    # One font declaration per face: two drawing paths each declare Arial, which leaves the page
    # with a duplicate entry.
    for page in pages:
        path = Path(page["path"])
        merged, dropped = figure_house.merge_duplicate_fonts(path.read_text(encoding="utf-8"))
        if dropped:
            path.write_text(str(merged), encoding="utf-8")
    if render:
        # Rendering is not a styling step: a plain panel needs its PNG too. This used
        # to sit inside ``if restyle`` so an unstyled scope panel with render=true
        # reported native_preview "unavailable" and wrote no png/ folder.
        png_dir = target / "png"
        png_dir.mkdir(exist_ok=True)
        try:
            for page in pages:
                _render_png(page["path"], str(png_dir / (Path(page["file"]).stem + ".png")))
            checks.setdefault("native_preview", {})["status"] = "rendered_pending_visual_review"
        except Exception as exc:
            checks.setdefault("native_preview", {})["error"] = str(exc)
    readback = _readback([Path(p["path"]) for p in pages], grounded, header, members)

    notes: list[str] = []
    warnings: list[str] = []
    for entry in brief["species"]:
        values = entry.get("values") or {}
        blank = [f["key"] for f in fields if values.get(f["key"]) in (None, "")]
        if blank:
            notes.append(f"species {entry['id']} has no {', '.join(blank)}; shown as a dash, not as zero.")
    warnings += style_warnings + highlight_warnings
    for item in checks["layout"]["issues"]:
        warnings.append(f"Layout check: {item}")
    for item in checks["data_binding"]["errors"]:
        warnings.append(f"Data binding: {item}")
    for item in project.get("warnings", []):
        warnings.append(f"{item.get('compound_id')}: {item.get('message')}")
    if readback["status"] != "matched":
        warnings.append(f"The saved pages do not match the brief: {readback['missing']} missing, "
                        f"{readback['unexpected']} unexpected structure(s).")

    preview = checks.get("native_preview", {})
    produced: list[Path] = [Path(p["path"]) for p in pages]
    if render and preview.get("status") == "rendered_pending_visual_review":
        pngs = sorted((target / "png").glob("*.png"))
        native: dict[str, Any] = {"status": "rendered", "png_paths": [str(p) for p in pngs]}
        produced += pngs
    else:
        reason = preview.get("error") or ("native rendering was not requested" if not render else "native rendering did not run")
        native = {"status": "not_run", "reason": reason}
        notes.append("The panel has not been rendered natively; look at a render before treating it as final."
                     if not render else "Native rendering did not run, so the panel is an unrendered draft.")

    report = {
        "ok": True,
        "operation": "draw_figure",
        "kind": "scope",
        "output_dir": str(target),
        "pages": pages,
        "species": [{k: g[k] for k in ("id", "label", "name", "smiles", "source")} for g in grounded],
        "fields": fields,
        "highlighted": highlighted,
        "style": style_report(brief_style(brief)[0], house, applied, skipped),
        "checks": {"data_binding": checks["data_binding"], "chemistry": checks["chemistry"]},
        "lint": {"issues": checks["layout"]["issues"], "relayout": None},
        "readback": readback,
        "native": native,
        "checks_passed": (checks["data_binding"]["status"] == "passed" and checks["chemistry"]["status"] == "passed"
                          and not checks["layout"]["issues"] and readback["status"] == "matched"),
        "notes": notes,
        "warnings": warnings,
    }
    return artifact_safety.with_artifacts(report, produced)


__all__ = ["SCOPE_KEYS", "draw_scope"]
