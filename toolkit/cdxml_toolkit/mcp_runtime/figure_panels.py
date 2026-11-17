"""Several schemes on one page, each under a boxed title, with a legend (``kind: "panels"``).

Journal figures often set two or three schemes one above the other, each under a boxed one-line
title ("(a) examples of ...", "(b) previous work", "(c) this work"), with a legend of variable
groups under the last (BJOC 22-103 figure 1). Each panel is an ordinary scheme or mechanism
brief, drawn and checked by ``draw_figure``; the panels are then put on one page with
``cdxml_compose.merge_cdxml``, which keeps every part's colours, fonts and chemistry.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from . import artifact_safety, figure_house

PANEL_BRIEF_KEYS = {"kind", "panels", "legend", "style", "max_width"}
PANEL_KEYS = {"title", "brief"}
MAX_PANELS = 8
TITLE_HEIGHT = 16.0   # a 10 pt line in its box
TITLE_GAP = 8.0       # between a title box and its panel
PANEL_GAP = 16.0      # between one panel and the next title
TITLE_INSET = 4.0     # the title's text from the box's left edge


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
        inner = panel.get("brief")
        if not isinstance(inner, dict) or inner.get("kind", "scheme") not in ("scheme", "mechanism"):
            errors.append(f"{where}: brief must be a scheme or a mechanism brief (a scope panel is a page of its own).")
    legend = brief.get("legend")
    if legend is not None:
        if not isinstance(legend, list):
            errors.append("legend must be a list of lines (text, a run {text, color, bold, italic}, or a list of runs).")
        else:
            for number, line in enumerate(legend, 1):
                errors += figure_house.condition_line_errors(line, f"legend line {number}")
    return errors


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
        files, sizes = [], []
        for number, panel in enumerate(brief["panels"], 1):
            inner = dict(panel["brief"])
            for key in ("style", "max_width"):
                if brief.get(key) is not None:
                    inner.setdefault(key, brief[key])
            report = draw_figure(inner, output_path=str(work / f"panel-{number}.cdxml"), render=False)
            if not report.get("ok"):
                return _failure([f"Panel {number} ({panel['title']}): {e}" for e in report.get("errors", [])],
                                kind="panels")
            reports.append(report)
            path = Path(report["output_path"])
            x1, y1, x2, y2 = content_box(list(parse_cdxml(str(path)).find("page")))
            files.append(path)
            sizes.append((x2 - x1, y2 - y1))
        titles, panels, width, bottom = _stack(sizes)
        parts = []
        for number, (panel, path) in enumerate(zip(brief["panels"], files), 1):
            title = _compose([{"type": "rectangle", "start": [0, 0], "end": [width, TITLE_HEIGHT]},
                              {"type": "text", "text": panel["title"], "position": [TITLE_INSET, TITLE_HEIGHT - 4.5]}],
                             work / f"title-{number}.cdxml")
            parts.append({"file": str(path), "at": list(panels[number - 1]), "anchor": "top-left"})
            parts.append({"file": str(title), "at": list(titles[number - 1]), "anchor": "top-left"})
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
