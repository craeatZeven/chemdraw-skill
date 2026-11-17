"""Several schemes on one page, each under a boxed title, with a legend (``kind: "panels"``).

Journal figures often set two or three schemes one above the other, each under a boxed one-line
title ("(a) examples of ...", "(b) previous work", "(c) this work"), with a legend of variable
groups under the last (BJOC 22-103 figure 1). Each panel is one of:

- ``brief``: an ordinary scheme or mechanism brief, drawn and checked by ``draw_figure``;
- ``gallery``: structures with their names and no arrows (panel (a), the drugs that share a
  motif), packed into rows and each row spread to the widest;
- ``sections``: several briefs, each opened by a boxed line on a dashed rule across the figure
  (panel (b), one section per earlier report).

The parts are put on one page with ``cdxml_compose.merge_cdxml``, which keeps every part's
colours, fonts and chemistry.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from collections import Counter
from typing import Any

from . import artifact_safety, figure_house

PANEL_BRIEF_KEYS = {"kind", "panels", "legend", "style", "max_width"}
PANEL_KEYS = {"title", "brief", "gallery", "sections"}
PANEL_CONTENT = ("brief", "gallery", "sections")
SECTION_KEYS = {"title", "brief"}
GALLERY_KEYS = {"id", "smiles", "name", "caption", *figure_house.SPECIES_STYLE_KEYS}
MAX_PANELS = 8
MAX_SECTIONS = 8
MAX_GALLERY = 40
TITLE_HEIGHT = 16.0   # a 10 pt line in its box
TITLE_GAP = 8.0       # between a title box and its panel
PANEL_GAP = 16.0      # between one panel and the next title
TITLE_INSET = 4.0     # the title's text from the box's left edge
SECTION_LEAD = 8.0    # the dashed piece before a section's box
SECTION_GAP = 4.0     # between a dashed piece and the box
SECTION_SHARE = 0.6   # a section box spans at least this share of the figure (the witness's ~60%)
SECTION_SPACING = 12.0  # between one section's drawing and the next section's rule
GALLERY_GAP = 24.0    # between gallery cells, across and down
GALLERY_WIDTH = 480.0  # the default row width of a gallery, in points


def _validate(brief: dict) -> list[str]:
    errors: list[str] = []
    unknown = sorted(set(brief) - PANEL_BRIEF_KEYS)
    if unknown:
        errors.append(f"Unknown brief keys: {unknown}")
    panels = brief.get("panels")
    if not isinstance(panels, list) or not 1 <= len(panels) <= MAX_PANELS:
        errors.append(f"panels: give 1 to {MAX_PANELS} panels, each {{title, brief}}.")
        panels = []
    for number, panel in enumerate(panels, 1):
        where = f"panels[{number}]"
        if not isinstance(panel, dict):
            errors.append(f"{where} must be an object {{title, brief}}.")
            continue
        extra = sorted(set(panel) - PANEL_KEYS)
        if extra:
            errors.append(f"{where} has unknown keys {extra}.")
        if not isinstance(panel.get("title"), str) or not panel["title"].strip():
            errors.append(f"{where} needs a title, the one line in its box (\"(a) previous work\").")
        given = [k for k in PANEL_CONTENT if k in panel]
        if len(given) != 1:
            errors.append(f"{where} needs exactly one of brief, gallery or sections.")
        elif given[0] == "brief":
            errors += _inner_errors(panel["brief"], where)
        elif given[0] == "gallery":
            errors += _gallery_errors(panel["gallery"], where)
        else:
            errors += _sections_errors(panel["sections"], where)
    legend = brief.get("legend")
    if legend is not None:
        if not isinstance(legend, list):
            errors.append("legend must be a list of lines (text, a run {text, color, bold, italic}, or a list of runs).")
        else:
            for number, line in enumerate(legend, 1):
                errors += figure_house.condition_line_errors(line, f"legend line {number}")
    return errors


def _inner_errors(inner: Any, where: str) -> list[str]:
    if not isinstance(inner, dict) or inner.get("kind", "scheme") not in ("scheme", "mechanism"):
        return [f"{where}: brief must be a scheme or a mechanism brief (a scope panel is a page of its own)."]
    return []


def _gallery_errors(gallery: Any, where: str) -> list[str]:
    if not isinstance(gallery, list) or not 1 <= len(gallery) <= MAX_GALLERY:
        return [f"{where}: gallery must list 1 to {MAX_GALLERY} structures, each {{id, smiles or name, caption?}}."]
    errors: list[str] = []
    seen: set[str] = set()
    for number, entry in enumerate(gallery, 1):
        at = f"{where} gallery[{number}]"
        if not isinstance(entry, dict):
            errors.append(f"{at} must be an object {{id, smiles or name, caption?}}.")
            continue
        extra = sorted(set(entry) - GALLERY_KEYS)
        if extra:
            errors.append(f"{at} has unknown keys {extra}.")
        ident = str(entry.get("id") or "")
        if not ident or ident in seen:
            errors.append(f"{at} needs an id of its own.")
        seen.add(ident)
        if not entry.get("smiles") and not entry.get("name"):
            errors.append(f"{at} needs a smiles or a name.")
        caption = entry.get("caption")
        if caption is not None and not isinstance(caption, str) and not (
                isinstance(caption, list) and all(isinstance(line, str) for line in caption)):
            errors.append(f"{at}: caption must be text (lines split by \\n) or a list of lines.")
        errors += [f"{at}: {e}" for e in figure_house.species_style_errors({**entry, "id": ident})]
    return errors


def _sections_errors(sections: Any, where: str) -> list[str]:
    if not isinstance(sections, list) or not 1 <= len(sections) <= MAX_SECTIONS:
        return [f"{where}: sections must list 1 to {MAX_SECTIONS} sections, each {{title, brief}}."]
    errors: list[str] = []
    for number, section in enumerate(sections, 1):
        at = f"{where} sections[{number}]"
        if not isinstance(section, dict):
            errors.append(f"{at} must be an object {{title, brief}}.")
            continue
        extra = sorted(set(section) - SECTION_KEYS)
        if extra:
            errors.append(f"{at} has unknown keys {extra}.")
        if not isinstance(section.get("title"), str) or not section["title"].strip():
            errors.append(f"{at} needs a title, the line in its box (\"L. Grimaud 2009 [35]\").")
        errors += _inner_errors(section.get("brief"), at)
    return errors


def _gallery_rows(sizes: list[tuple[float, float]], max_width: float,
                  gap: float = GALLERY_GAP) -> tuple[list[tuple[float, float]], float, float]:
    """Top-left corners of gallery cells: rows filled left to right up to max_width, then each
    row of two or more spread to the widest row's width. Returns (corners, width, height)."""
    rows: list[list[int]] = [[]]
    used = 0.0
    for index, (w, _) in enumerate(sizes):
        if rows[-1] and used + gap + w > max_width:
            rows.append([])
            used = 0.0
        used = w if not rows[-1] else used + gap + w
        rows[-1].append(index)
    natural = [sum(sizes[i][0] for i in row) + gap * (len(row) - 1) for row in rows]
    width = max(natural)
    corners: list[tuple[float, float]] = [(0.0, 0.0)] * len(sizes)
    y = 0.0
    for row, used in zip(rows, natural):
        spread = gap + ((width - used) / (len(row) - 1) if len(row) > 1 else 0.0)
        x = 0.0
        for i in row:
            corners[i] = (x, y)
            x += sizes[i][0] + spread
        y += max(sizes[i][1] for i in row) + gap
    return corners, width, y - gap


def _section_rule(title: str, width: float, text_width: float) -> list[dict]:
    """A section's opening line: a short dashed piece, the boxed title, a dashed rule to the
    figure's right edge, all through the middle of the box (BJOC 22-103 figure 1, panel (b))."""
    mid = TITLE_HEIGHT / 2
    start = SECTION_LEAD + SECTION_GAP
    room = width - start - SECTION_GAP - SECTION_LEAD
    box = min(room, max(SECTION_SHARE * width, text_width + 2 * TITLE_INSET))
    end = start + box
    return [{"type": "line", "start": [0.0, mid], "end": [SECTION_LEAD, mid], "style": "dashed"},
            {"type": "rectangle", "start": [start, 0.0], "end": [end, TITLE_HEIGHT]},
            {"type": "text", "text": title, "position": [start + TITLE_INSET, TITLE_HEIGHT - 4.5]},
            {"type": "line", "start": [end + SECTION_GAP, mid], "end": [width, mid], "style": "dashed"}]


def _caption_lines(entry: dict) -> list[str]:
    caption = entry.get("caption")
    if caption is None:
        return []
    lines = caption.splitlines() if isinstance(caption, str) else list(caption)
    return [line for line in lines if line.strip()]


def _draw_gallery(entries: list[dict], path: Path, max_width: float | None) -> dict:
    """Structures with their caption lines and no arrows, in rows; checked like a scheme."""
    from rdkit import Chem

    from ..cdxml_check import check_cdxml
    from ..chemistry_semantics import document_inventory, semantic_key
    from ..publication.layout import CELL, caption_runs, molecule_coordinates, style_settings, text_width
    from .figure_draw import _ground

    grounded, errors = _ground(entries)
    if errors:
        return {"ok": False, "errors": errors}
    style = style_settings({})
    size = style["font_size"]
    pad = size * CELL["pad"]
    advance = size * CELL["line_advance"] + CELL["line_gap"]
    cells = []
    for entry, item in zip(entries, grounded):
        _, xy, _ = molecule_coordinates(Chem.MolFromSmiles(item["smiles"]), style["bond_length"])
        xs, ys = [p[0] for p in xy], [p[1] for p in xy]
        lines = [caption_runs(line) for line in _caption_lines(entry)]
        molecule_w, molecule_h = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
        caption_w = max([text_width(plain, style) for plain, _ in lines], default=0.0)
        caption_h = (size * CELL["caption_baseline"] - pad + advance * (len(lines) - 1) + size * 0.3) if lines else 0.0
        cells.append({"item": item, "xy": xy, "x0": min(xs), "y0": min(ys), "lines": lines,
                      "mw": molecule_w, "mh": molecule_h, "w": max(molecule_w, caption_w), "h": molecule_h + caption_h})
    corners, _, _ = _gallery_rows([(c["w"], c["h"]) for c in cells], max_width or GALLERY_WIDTH)
    objects: list[dict] = []
    for number, (cell, (x, y)) in enumerate(zip(cells, corners)):
        cx = x + cell["w"] / 2
        dx = cx - cell["mw"] / 2 + pad - cell["x0"]
        dy = y + pad - cell["y0"]
        objects.append({"type": "molecule", "id": f"g{number}", "smiles": cell["item"]["smiles"],
                        "coordinates": [[px + dx, py + dy] for px, py in cell["xy"]]})
        baseline = y + cell["mh"] - pad + size * CELL["caption_baseline"]
        for plain, runs in cell["lines"]:
            objects.append({"type": "text", "text": plain, "align": "Center", "size": size,
                            "position": [cx, baseline], **({"runs": runs} if runs else {})})
            baseline += advance
    _compose(objects, path)
    text, applied, warnings = figure_house.apply_species_styling(
        path.read_text(encoding="utf-8"), {"species": entries}, grounded, None)
    if applied:
        path.write_text(str(text), encoding="utf-8")
    wanted: Counter = Counter()
    for item in grounded:
        wanted.update(semantic_key(m) for m in Chem.GetMolFrags(Chem.MolFromSmiles(item["smiles"]), asMols=True))
    found = document_inventory(str(path), like=wanted)
    check = check_cdxml(str(path), palette=None)
    return {"ok": True, "output_path": str(path), "layout": "gallery",
            "readback": {"status": "matched" if found == wanted else "mismatch"},
            "checks_passed": not check["errors"], "warnings": list(warnings)}


def _stack(sizes: list[tuple[float, float]], title_height: float = TITLE_HEIGHT, gap: float = TITLE_GAP,
           panel_gap: float = PANEL_GAP) -> tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]:
    """Top-left corners of each title box and each panel, one under the other, every title box
    as wide as the widest panel. Returns (titles, panels, width, bottom)."""
    width = max(w for w, _ in sizes)
    titles, panels, y = [], [], 0.0
    for index, (_, height) in enumerate(sizes):
        if index:
            y += panel_gap
        titles.append((0.0, y))
        y += title_height + gap
        panels.append((0.0, y))
        y += height
    return titles, panels, width, y


def _runs(line: Any) -> list[dict]:
    if isinstance(line, str):
        return [{"text": line}]
    if isinstance(line, dict):
        return [dict(line)]
    return [dict(run) for run in line]


def _compose(objects: list[dict], path: Path) -> Path:
    from .figure_tools import compose_chemical_figure

    manifest = path.with_suffix(".json")
    manifest.write_text(json.dumps({"objects": objects}), encoding="utf-8")
    compose_chemical_figure(str(manifest), str(path))
    return path


def draw_panels(brief: dict, output_path: str | None, render: bool) -> dict:
    """Draw each panel with draw_figure, then stack them under boxed titles on one page."""
    from ..cdxml_check import check_cdxml
    from ..cdxml_compose import content_box, merge_cdxml, parse_cdxml
    from .figure_draw import _failure, _png_destination, _render_png, draw_figure

    errors = _validate(brief)
    if errors:
        return _failure(errors, kind="panels")
    reports: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="panels-") as scratch:
        work = Path(scratch)

        def size_of(path: Path) -> tuple[float, float]:
            x1, y1, x2, y2 = content_box(list(parse_cdxml(str(path)).find("page")))
            return x2 - x1, y2 - y1

        def drawn(inner_brief: dict, name: str, where: str) -> tuple[dict | None, list[str]]:
            inner = dict(inner_brief)
            for key in ("style", "max_width"):
                if brief.get(key) is not None:
                    inner.setdefault(key, brief[key])
            report = draw_figure(inner, output_path=str(work / name), render=False)
            if not report.get("ok"):
                return None, [f"{where}: {e}" for e in report.get("errors", [])]
            return report, []

        # Each panel is a list of blocks drawn one under another: (rule title or None, file, size).
        blocks: list[list[tuple[str | None, Path, tuple[float, float]]]] = []
        for number, panel in enumerate(brief["panels"], 1):
            where = f"Panel {number} ({panel['title']})"
            if "gallery" in panel:
                report = _draw_gallery(panel["gallery"], work / f"panel-{number}.cdxml", brief.get("max_width"))
                if not report.get("ok"):
                    return _failure([f"{where}: {e}" for e in report["errors"]], kind="panels")
                path = Path(report["output_path"])
                blocks.append([(None, path, size_of(path))])
            elif "sections" in panel:
                parts_of: list[dict] = []
                block = []
                for k, section in enumerate(panel["sections"], 1):
                    report, failed = drawn(section["brief"], f"panel-{number}-{k}.cdxml",
                                           f"{where}, section {k} ({section['title']})")
                    if failed:
                        return _failure(failed, kind="panels")
                    parts_of.append(report)
                    path = Path(report["output_path"])
                    block.append((section["title"], path, size_of(path)))
                blocks.append(block)
                report = {"layout": "sections",
                          "readback": {"status": "matched" if all(r["readback"]["status"] == "matched"
                                                                  for r in parts_of) else "mismatch"},
                          "checks_passed": all(r.get("checks_passed") for r in parts_of),
                          "warnings": [w for r in parts_of for w in r.get("warnings", [])]}
            else:
                report, failed = drawn(panel["brief"], f"panel-{number}.cdxml", where)
                if failed:
                    return _failure(failed, kind="panels")
                path = Path(report["output_path"])
                blocks.append([(None, path, size_of(path))])
            reports.append(report)
        rule_height = TITLE_HEIGHT + TITLE_GAP
        sizes = [(max(w for _, _, (w, _) in block),
                  sum(h + (rule_height if title else 0) for title, _, (_, h) in block)
                  + SECTION_SPACING * (len(block) - 1))
                 for block in blocks]
        titles, panels, width, bottom = _stack(sizes)
        from ..publication.layout import style_settings, text_width
        text_style = style_settings({})
        parts = []
        for number, (panel, block) in enumerate(zip(brief["panels"], blocks), 1):
            title = _compose([{"type": "rectangle", "start": [0, 0], "end": [width, TITLE_HEIGHT]},
                              {"type": "text", "text": panel["title"], "position": [TITLE_INSET, TITLE_HEIGHT - 4.5]}],
                             work / f"title-{number}.cdxml")
            parts.append({"file": str(title), "at": list(titles[number - 1]), "anchor": "top-left"})
            y = panels[number - 1][1]
            for k, (rule, path, (_, height)) in enumerate(block, 1):
                if rule:
                    drawn_rule = _compose(_section_rule(rule, width, text_width(rule, text_style)),
                                          work / f"rule-{number}-{k}.cdxml")
                    parts.append({"file": str(drawn_rule), "at": [0.0, y], "anchor": "top-left"})
                    y += rule_height
                parts.append({"file": str(path), "at": [0.0, y], "anchor": "top-left"})
                y += height + SECTION_SPACING
        if brief.get("legend"):
            lines = [{"type": "text", "text": "".join(r["text"] for r in _runs(line)), "runs": _runs(line),
                      "position": [0, 10 + 12 * n]} for n, line in enumerate(brief["legend"])]
            legend = _compose(lines, work / "legend.cdxml")
            parts.append({"file": str(legend), "at": [0.0, bottom + PANEL_GAP], "anchor": "top-left"})
        merged = merge_cdxml(parts)
        target = artifact_safety.resolve_destination(source=None, output_path=output_path, tag="figure",
                                                     suffix=".cdxml")
        with artifact_safety.staging_file(target) as stage:
            stage.write_text(merged["cdxml"], encoding="utf-8", newline="")
            # staging_file only lends a temporary path; without publishing, closing the block
            # deletes it and the target never exists — which is what check_cdxml then tripped on.
            artifact_safety.publish_file(stage, target)
    produced = [target]
    native = {"status": "not_run", "reason": "render was not requested"}
    if render:
        try:
            png = _png_destination(target)
            _render_png(str(target), str(png))
            native = {"status": "rendered", "png_path": str(png)}
            produced.append(png)
        except Exception as exc:
            native = {"status": "not_run", "reason": str(exc)}
    check = check_cdxml(str(target), palette=None)
    summary = [{"title": p["title"], "layout": r.get("layout"), "readback": r["readback"]["status"],
                "checks_passed": r.get("checks_passed"), "warnings": r.get("warnings", [])}
               for p, r in zip(brief["panels"], reports)]
    report = {
        "ok": True,
        "operation": "draw_figure",
        "kind": "panels",
        "output_path": str(target),
        "panels": summary,
        "check": {"errors": check["errors"], "warnings": check["warnings"]},
        "native": native,
        "checks_passed": all(p["checks_passed"] for p in summary) and not check["errors"],
        "warnings": merged["warnings"],
    }
    return artifact_safety.with_artifacts(report, produced)
