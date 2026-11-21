"""Optimization and screening tables for ``draw_figure`` (``kind: "table"``).

Draws publication-grade booktabs (three-line) tables with an optional general
reaction scheme above, Arial-metric column sizing so rules never strike
through text, automatic chemical formula subscripting and italic-superscript
footnote markers, optimal-row highlighting, and wrapped footnotes.
"""
from __future__ import annotations

import json
import math
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

from ..chemistry_semantics import document_inventory
from ..chemdraw.style_presets import DEFAULT_STYLE_PRESET, style_preset, style_preset_as_manifest_style
from ..publication.house_styles import HouseStyle
from ..publication.layout import (
    caption_runs,
    header_identities,
    header_objects,
    text_width,
    wrap,
)
from ..publication.palettes import resolve_palette
from . import artifact_safety
from .figure_draw import (
    DEFAULT_MAX_WIDTH,
    _failure,
    _png_destination,
    _render_png,
    brief_house,
    brief_style,
    style_preset_error,
    style_report,
)
from .figure_scope import _ground_header, _header_errors
from .figure_tools import compose_chemical_figure

TABLE_KEYS = {
    "kind",
    "title",
    "header",
    "columns",
    "rows",
    "footnotes",
    "style",
    "table_style",
    "max_width",
    "highlight_color",
}
COLUMN_KEYS = {"key", "label", "align"}
ALIGNS = ("left", "center", "right")
TABLE_STYLES = ("horizontal_rules", "grid", "none")
HIGHLIGHT_MODES = ("bold", "color", "shaded", "boxed")
_HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")
MAX_ROWS = 200
MAX_COLS = 20


def _normalise_columns(raw: Any) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    if not isinstance(raw, list) or not raw:
        return [], ["columns: list at least one column (name or {key, label?, align?})."]
    if len(raw) > MAX_COLS:
        errors.append(f"A table holds at most {MAX_COLS} columns.")
    cols: list[dict[str, str]] = []
    seen: set[str] = set()
    for idx, item in enumerate(raw, 1):
        if isinstance(item, str):
            item = {"key": item, "label": item}
        if not isinstance(item, dict) or not str(item.get("key") or "").strip():
            errors.append(f"columns[{idx}] must be a column name or {{key, label?, align?}}.")
            continue
        extra = sorted(set(item) - COLUMN_KEYS)
        if extra:
            errors.append(f"columns[{idx}] has unknown keys {extra}.")
        key = str(item["key"])
        if key in seen:
            errors.append(f"columns[{idx}]: duplicate key {key!r}.")
            continue
        seen.add(key)
        label = str(item.get("label") if item.get("label") is not None else key)
        align = str(item.get("align") or "center").lower()
        if align not in ALIGNS:
            errors.append(f"columns[{idx}]: align must be one of {list(ALIGNS)}.")
            align = "center"
        cols.append({"key": key, "label": label, "align": align})
    return cols, errors


def _validate(brief: dict) -> list[str]:
    errors: list[str] = []
    unknown = sorted(set(brief) - TABLE_KEYS)
    if unknown:
        errors.append(f"Unknown table brief keys: {unknown}")
    style_err = style_preset_error(brief.get("style"))
    if style_err:
        errors.append(style_err)
    if brief.get("header") is not None:
        errors += _header_errors(brief["header"])
    if brief.get("title") is not None and not isinstance(brief["title"], str):
        errors.append("title must be text.")
    ts = brief.get("table_style")
    if ts is not None and ts not in TABLE_STYLES:
        errors.append(f"table_style must be one of {list(TABLE_STYLES)}.")
    hc = brief.get("highlight_color")
    if hc is not None and not (isinstance(hc, str) and _HEX_RE.match(hc)):
        errors.append("highlight_color must be a hex colour like #800000.")
    width = brief.get("max_width")
    if width is not None and not (isinstance(width, (int, float)) and width >= 120):
        errors.append("max_width must be a number of points >= 120.")
    footnotes = brief.get("footnotes")
    if footnotes is not None and not (isinstance(footnotes, list) and all(isinstance(f, str) for f in footnotes)):
        errors.append("footnotes must be a list of text lines.")

    cols, col_errors = _normalise_columns(brief.get("columns"))
    errors += col_errors
    rows = brief.get("rows")
    if not isinstance(rows, list) or not rows:
        errors.append("rows: list at least one data row.")
        return errors
    if len(rows) > MAX_ROWS:
        errors.append(f"A table holds at most {MAX_ROWS} rows.")
    for idx, row in enumerate(rows, 1):
        if isinstance(row, list):
            if cols and len(row) != len(cols):
                errors.append(f"rows[{idx}] has {len(row)} values, expected {len(cols)}.")
        elif isinstance(row, dict):
            hl = row.get("highlight")
            if hl is not None and not (
                isinstance(hl, bool)
                or (isinstance(hl, str) and (hl in HIGHLIGHT_MODES or _HEX_RE.match(hl)))
            ):
                errors.append(
                    f"rows[{idx}].highlight must be true, false, one of {list(HIGHLIGHT_MODES)}, or #RRGGBB."
                )
        else:
            errors.append(f"rows[{idx}] must be an object or a list of cell values.")
    return errors


def _format_cell(value: Any, *, bold: bool = False, color: str | None = None) -> tuple[str, list[dict] | None]:
    text = "-" if value is None or str(value) == "" else str(value)
    plain, runs = caption_runs(text, bold=bold)
    if color:
        if runs is None:
            runs = [{"text": plain, **({"bold": True} if bold else {}), "color": color}]
        else:
            runs = [{**r, "color": color} for r in runs]
    return plain, runs


def _measure_runs(plain: str, runs: list[dict] | None, style: dict, bold: bool = False) -> float:
    size = float(style["font_size"])
    if not runs:
        return text_width(plain, style, size, bold=bold)
    total = 0.0
    for r in runs:
        r_size = size * 0.7 if (r.get("subscript") or r.get("superscript")) else size
        total += text_width(r["text"], style, r_size, bold=bool(r.get("bold") or bold))
    return total


def _align_name(align: str) -> str:
    return {"left": "Left", "center": "Center", "right": "Right"}[align]


def _resolve_highlight_style(
    row_hl: Any,
    house: HouseStyle | None,
    palette: dict[str, Any],
    override_color: str | None,
) -> tuple[bool, str | None, str | None]:
    """Return (bold, text colour, rectangle fill type) for a row's highlight setting."""
    if row_hl is True:
        entry_hl = house.figure_convention("optimization_table", "entry_highlight") if house else None
        row_hl = {"bold_row": "bold", "shaded_row": "shaded", "boxed_row": "boxed", "none": None}.get(entry_hl, "color")
    if not row_hl:
        return False, None, None
    accent = override_color or palette["primary"]
    if isinstance(row_hl, str) and _HEX_RE.match(row_hl):
        return True, row_hl.upper(), None
    if row_hl == "bold":
        return True, None, None
    if row_hl == "shaded":
        return True, accent, "solid"
    if row_hl == "boxed":
        return True, accent, "none"
    return True, accent, None


def draw_table(brief: dict, output_path: str | None, render: bool) -> dict:
    """Draw a publication-grade optimization or screening table (``kind: "table"``)."""
    errors = _validate(brief)
    if errors:
        return _failure(errors, kind="table")

    header, header_errors = _ground_header(brief.get("header"))
    if header_errors:
        return _failure(header_errors, kind="table")

    cols, _ = _normalise_columns(brief["columns"])
    slug, _ = brief_style(brief)
    house = brief_house(brief)
    palette = resolve_palette(house if house is not None else slug)
    preset = style_preset(slug)
    attrs = preset["attributes"]
    manifest_style = style_preset_as_manifest_style(slug)

    # Preserve measured ratios when a house style is requested
    applied: list[str] = []
    skipped: list[str] = []
    if house is not None:
        measured = house.measured or {}
        bond = float(manifest_style["bond_length"])
        lw = measured.get("line_width_over_bond_length")
        if isinstance(lw, (int, float)) and lw > 0:
            manifest_style["line_width"] = float(f"{min(max(float(lw) * bond, 0.3), 3.0):.3g}")
            applied.append(f"measured.line_width={manifest_style['line_width']:g} pt")
        bw = measured.get("bold_width_over_bond_length")
        if isinstance(bw, (int, float)) and bw > 0:
            manifest_style["bold_width"] = float(f"{min(max(float(bw) * bond, 0.8), 6.0):.3g}")
            applied.append(f"measured.bold_width={manifest_style['bold_width']:g} pt")

    size = float(manifest_style["font_size"])
    layout_style = {
        "font": manifest_style.get("font", "Arial"),
        "font_size": size,
        "bond_length": float(manifest_style["bond_length"]),
    }

    table_style = brief.get("table_style")
    if table_style is None and house is not None:
        hs_ts = house.figure_convention("optimization_table", "table_style")
        if hs_ts in TABLE_STYLES:
            table_style = hs_ts
            applied.append(f"optimization_table.table_style={hs_ts}")
    table_style = table_style or "horizontal_rules"

    if house is not None:
        hs_hl = house.figure_convention("optimization_table", "entry_highlight")
        if hs_hl:
            applied.append(f"optimization_table.entry_highlight={hs_hl}")

    # Normalize row records
    norm_rows: list[tuple[dict[str, Any], Any]] = []
    for raw_row in brief["rows"]:
        if isinstance(raw_row, list):
            norm_rows.append(({c["key"]: v for c, v in zip(cols, raw_row)}, False))
        else:
            norm_rows.append((raw_row, raw_row.get("highlight", False)))

    # Format headers & cells and measure column widths from Arial metrics
    col_pad = max(16.0, size * 1.6)
    margin = 16.0
    header_cells = [_format_cell(c["label"], bold=True) for c in cols]
    col_widths = [max(24.0, _measure_runs(p, r, layout_style, bold=True)) for p, r in header_cells]

    formatted_rows: list[tuple[list[tuple[str, list[dict] | None]], str | None, str | None]] = []
    override_hl_color = brief.get("highlight_color")
    for row_map, row_hl in norm_rows:
        bold, color, fill_type = _resolve_highlight_style(row_hl, house, palette, override_hl_color)
        cells: list[tuple[str, list[dict] | None]] = []
        for idx, col in enumerate(cols):
            plain, runs = _format_cell(row_map.get(col["key"]), bold=bold, color=color)
            cells.append((plain, runs))
            col_widths[idx] = max(col_widths[idx], _measure_runs(plain, runs, layout_style, bold=bold))
        formatted_rows.append((cells, color, fill_type))

    # Measure optional header reaction scheme first so table width aligns cleanly with it
    objects: list[dict[str, Any]] = []
    y = margin
    right_edge = margin

    if brief.get("title"):
        title_plain, title_runs = _format_cell(brief["title"], bold=True)
        objects.append({
            "type": "text",
            "id": "table_title",
            "text": title_plain,
            "position": [margin, y + size],
            "size": size,
            "align": "Left",
            **({"runs": title_runs} if title_runs else {}),
        })
        right_edge = max(right_edge, margin + _measure_runs(title_plain, title_runs, layout_style, bold=True))
        y += size * 1.8

    header_width = 0.0
    if header:
        drawn_header, header_bottom, header_right = header_objects({"header": header}, layout_style, margin, y)
        objects.extend(drawn_header)
        header_width = header_right - margin
        right_edge = max(right_edge, header_right)
        y = header_bottom + size * 0.9

    # Compute final column geometry
    inner_pad = 6.0
    raw_table_w = sum(col_widths) + col_pad * (len(cols) - 1) + 2 * inner_pad
    max_w = float(brief.get("max_width") or DEFAULT_MAX_WIDTH)
    table_w = max(raw_table_w, header_width)
    if table_w > raw_table_w and len(cols) > 1:
        extra_per_gap = (table_w - raw_table_w) / (len(cols) - 1)
        col_gap = col_pad + extra_per_gap
    else:
        col_gap = col_pad

    col_lefts: list[float] = []
    cx = margin + inner_pad
    for w in col_widths:
        col_lefts.append(cx)
        cx += w + col_gap
    table_left = margin
    table_right = margin + table_w
    right_edge = max(right_edge, table_right)

    # Booktabs spacing: guarantees generous clearance between rules and text cap/descender
    top_rule_y = y
    if table_style in ("horizontal_rules", "grid"):
        objects.append({
            "type": "line",
            "start": [table_left, top_rule_y],
            "end": [table_right, top_rule_y],
            "style": "line",
        })

    header_baseline = top_rule_y + size * 1.35
    for idx, (col, (plain, runs)) in enumerate(zip(cols, header_cells)):
        align = col["align"]
        if align == "left":
            tx = col_lefts[idx]
        elif align == "right":
            tx = col_lefts[idx] + col_widths[idx]
        else:
            tx = col_lefts[idx] + col_widths[idx] / 2.0
        objects.append({
            "type": "text",
            "id": f"col_header_{idx}",
            "text": plain,
            "position": [tx, header_baseline],
            "size": size,
            "align": _align_name(align),
            **({"runs": runs} if runs else {}),
        })

    mid_rule_y = header_baseline + size * 0.60
    if table_style in ("horizontal_rules", "grid"):
        objects.append({
            "type": "line",
            "start": [table_left, mid_rule_y],
            "end": [table_right, mid_rule_y],
            "style": "line",
        })

    row_pitch = size * 1.65
    row_baseline = mid_rule_y + size * 1.35
    for r_idx, (cells, color, fill_type) in enumerate(formatted_rows):
        if fill_type:
            # Frame strictly clear of cap-top and descenders so text_over_ink stays zero
            box_top = row_baseline - size * 1.15
            box_bottom = row_baseline + size * 0.48
            objects.append({
                "type": "rectangle",
                "start": [table_left + 1.0, box_top],
                "end": [table_right - 1.0, box_bottom],
                "color": palette["highlight_fill"] if fill_type == "solid" else color,
                "fill_type": fill_type,
                "background": fill_type == "solid",
            })
        for c_idx, (col, (plain, runs)) in enumerate(zip(cols, cells)):
            align = col["align"]
            if align == "left":
                tx = col_lefts[c_idx]
            elif align == "right":
                tx = col_lefts[c_idx] + col_widths[c_idx]
            else:
                tx = col_lefts[c_idx] + col_widths[c_idx] / 2.0
            objects.append({
                "type": "text",
                "id": f"cell_{r_idx}_{c_idx}",
                "text": plain,
                "position": [tx, row_baseline],
                "size": size,
                "align": _align_name(align),
                **({"runs": runs} if runs else {}),
            })
        if table_style == "grid" and r_idx < len(formatted_rows) - 1:
            sep_y = row_baseline + size * 0.55
            objects.append({
                "type": "line",
                "start": [table_left, sep_y],
                "end": [table_right, sep_y],
                "style": "line",
            })
        row_baseline += row_pitch

    last_baseline = row_baseline - row_pitch
    bottom_rule_y = last_baseline + size * 0.65
    if table_style in ("horizontal_rules", "grid"):
        objects.append({
            "type": "line",
            "start": [table_left, bottom_rule_y],
            "end": [table_right, bottom_rule_y],
            "style": "line",
        })

    cursor_y = bottom_rule_y
    footnotes = brief.get("footnotes") or []
    if footnotes:
        foot_size = min(size, 8.0)
        foot_style = {**layout_style, "font_size": foot_size}
        cursor_y += foot_size * 1.5
        wrap_w = max(table_w, min(max_w, 360.0))
        for f_idx, raw_line in enumerate(footnotes):
            wrapped_lines = wrap(raw_line, wrap_w, foot_style).splitlines()
            for w_idx, wline in enumerate(wrapped_lines):
                plain, runs = _format_cell(wline)
                objects.append({
                    "type": "text",
                    "id": f"footnote_{f_idx}_{w_idx}",
                    "text": plain,
                    "position": [table_left, cursor_y],
                    "size": foot_size,
                    "align": "Left",
                    **({"runs": runs} if runs else {}),
                })
                right_edge = max(right_edge, table_left + _measure_runs(plain, runs, foot_style))
                cursor_y += foot_size * 1.45

    page_w = math.ceil(right_edge + margin)
    page_h = math.ceil(cursor_y + margin)

    try:
        target = artifact_safety.resolve_destination(
            source=None, output_path=output_path, tag="table", suffix=".cdxml"
        )
        with tempfile.TemporaryDirectory(prefix="table-figure-") as scratch:
            manifest_path = Path(scratch) / "table.manifest.json"
            staged_cdxml = Path(scratch) / "table.cdxml"
            manifest = {
                "version": 1,
                "style": {
                    k: manifest_style[k]
                    for k in ("font", "font_size", "bond_length", "line_width", "bold_width",
                              "hash_spacing", "margin_width")
                    if k in manifest_style
                },
                "objects": objects,
            }
            manifest_path.write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")
            compose_chemical_figure(str(manifest_path), str(staged_cdxml))

            # Crop page bounding box tightly to table content
            tree = ET.parse(staged_cdxml)
            root = tree.getroot()
            # Keep measured LabelSize if house style declared one
            if house is not None:
                cap = (house.measured or {}).get("bond_length_over_cap_height")
                if isinstance(cap, (int, float)) and cap > 0:
                    bond = float(manifest_style["bond_length"])
                    lsize = f"{min(max(bond / float(cap) / 0.716, 6.0), 14.0):.3g}"
                    root.set("LabelSize", lsize)
            page = root.find("page")
            bbox = f"0 0 {page_w:g} {page_h:g}"
            root.set("BoundingBox", bbox)
            if page is not None:
                page.set("BoundingBox", bbox)
                page.set("Width", f"{page_w:g}")
                page.set("Height", f"{page_h:g}")
                page.set("WidthPages", "1")
                page.set("HeightPages", "1")
            xml_bytes = ET.tostring(root, encoding="utf-8", xml_declaration=True)
            # Ensure uppercase UTF-8 declaration for native ChemDraw compatibility
            xml_text = xml_bytes.decode("utf-8").replace("encoding='utf-8'", 'encoding="UTF-8"', 1)
            with artifact_safety.staging_file(target) as stage:
                stage.write_text(xml_text, encoding="utf-8", newline="")
                artifact_safety.publish_file(stage, target)
    except ValueError as exc:
        return _failure([str(exc)], kind="table")
    except Exception as exc:
        return _failure([f"The table figure could not be drawn: {exc}"], kind="table")

    # Chemistry readback of optional header scheme
    if header:
        wanted = Counter(header_identities({"header": header}))
        found = document_inventory(target, like=wanted)
        missing = len(list((wanted - found).elements()))
        extra = len(list((found - wanted).elements()))
        readback = {
            "status": "matched" if missing == 0 and extra == 0 else "mismatch",
            "structures_drawn": sum(found.values()),
            "missing": missing,
            "unexpected": extra,
        }
    else:
        readback = {"status": "matched", "structures_drawn": 0, "missing": 0, "unexpected": 0}

    # Layout check via text_over_ink
    from ..cdxml_check import text_over_ink
    crowded = text_over_ink(ET.parse(target).getroot())
    lint_issues = [{"code": "text_over_ink", "message": msg} for msg in crowded]

    warnings: list[str] = []
    notes: list[str] = []
    if table_w > max_w:
        warnings.append(f"Table width {table_w:.1f} pt exceeds max_width {max_w:.1f} pt.")
    for item in lint_issues:
        warnings.append(f"Layout check: {item['message']}")
    if readback["status"] != "matched":
        warnings.append(
            f"The saved header structures do not match the brief: "
            f"{readback['missing']} missing, {readback['unexpected']} unexpected."
        )

    produced = [target]
    if not render:
        native: dict[str, Any] = {"status": "not_run", "reason": "native rendering was not requested"}
        notes.append("The drawing has not been rendered natively; look at a render before treating it as final.")
    else:
        png = _png_destination(target)
        try:
            _render_png(str(target), str(png))
            native = {"status": "rendered", "png_path": str(png)}
            produced.append(png)
        except Exception as exc:
            native = {"status": "not_run", "reason": str(exc)}
            notes.append("Native rendering did not run, so the drawing is an unrendered draft.")

    report = {
        "ok": True,
        "operation": "draw_figure",
        "kind": "table",
        "output_path": str(target),
        "table_style": table_style,
        "columns": cols,
        "row_count": len(formatted_rows),
        "style": style_report(slug, house, applied, skipped),
        "lint": {"issues": lint_issues, "relayout": None, "max_width": max_w, "width_pt": page_w},
        "readback": readback,
        "native": native,
        "checks_passed": not lint_issues and readback["status"] == "matched",
        "notes": notes,
        "warnings": warnings,
    }
    return artifact_safety.with_artifacts(report, produced)


__all__ = ["TABLE_KEYS", "draw_table"]
