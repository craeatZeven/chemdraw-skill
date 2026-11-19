"""Lay several scientific drawings out as one page: rows, labels, arrows, curved arrows.

Each part is a CDXML produced by another scientific generator (``lewis``, ``mo_diagram``,
``energy_profile`` ...), so this is the bridge between a single structure and a textbook figure.
Textbook conventions are the point: every row has one label baseline, an arrow sits on the bond
line rather than the middle of the box (lone pairs hang below it), a resonance sign is a
double-headed arrow and an equilibrium sign is two half-headed arrows.
"""
from __future__ import annotations

import copy
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Mapping

_IDS_BASE = 1000


def _bbox(page):
    xs, ys = [], []
    for element in page.iter():
        box = element.get("BoundingBox")
        if box and len(box.split()) == 4:
            v = [float(t) for t in box.split()]
            xs += [v[0], v[2]]
            ys += [v[1], v[3]]
        point = element.get("p")
        if point:
            xs.append(float(point.split()[0]))
            ys.append(float(point.split()[1]))
    if not xs:
        return (0.0, 0.0, 0.0, 0.0)
    return (min(xs), min(ys), max(xs), max(ys))


def _spine_y(page):
    """y of the longest horizontal bond: where an arrow between parts belongs."""
    best = None
    for graphic in page.iter("graphic"):
        box = graphic.get("BoundingBox")
        if not box:
            continue
        x0, y0, x1, y1 = [float(v) for v in box.split()]
        if abs(y1 - y0) < 1.5 and (best is None or abs(x1 - x0) > best[0]):
            best = (abs(x1 - x0), (y0 + y1) / 2.0)
    if best:
        return best[1]
    x0, y0, x1, y1 = _bbox(page)
    return (y0 + y1) / 2.0


class _Ids:
    def __init__(self):
        self.n = _IDS_BASE

    def next(self):
        self.n += 1
        return self.n


def _shift(element, dx, dy, ids):
    for node in element.iter():
        for key in ("p", "Tail3D", "Head3D", "Center3D", "MajorAxisEnd3D", "MinorAxisEnd3D"):
            value = node.get(key)
            if value and len(value.split()) >= 2:
                parts = value.split()
                parts[0] = str(float(parts[0]) + dx)
                parts[1] = str(float(parts[1]) + dy)
                node.set(key, " ".join(parts))
        box = node.get("BoundingBox")
        if box and len(box.split()) == 4:
            v = [float(t) for t in box.split()]
            node.set("BoundingBox", "%s %s %s %s" % (v[0] + dx, v[1] + dy, v[2] + dx, v[3] + dy))
        for key in ("id", "Z"):
            value = node.get(key)
            if value and value.isdigit():
                node.set(key, str(ids.next()))


def runs(text):
    """[(piece, face)]: 32 on an atom count after a letter (CH3), 64 on a charge written ^- or ^2-."""
    pieces, index = [], 0
    for match in re.finditer(r"\^(\d*[-+])", text):
        if match.start() > index:
            pieces.append((text[index:match.start()], "sub"))
        pieces.append((match.group(1), "sup"))
        index = match.end()
    if index < len(text):
        pieces.append((text[index:], "sub"))
    out = []
    for piece, mode in pieces:
        if mode == "sup":
            out.append((piece, "64"))
            continue
        start, runs_out = 0, []
        for match in re.finditer(r"[A-Za-z)](\d+)", piece):
            if match.start() + 1 > start:
                runs_out.append((piece[start:match.start() + 1], "0"))
            runs_out.append((match.group(1), "32"))
            start = match.end()
        if start < len(piece):
            runs_out.append((piece[start:], "0"))
        out += [item for item in runs_out if item[0]]
    return out or [(text, "0")]


def _text(page, ids, x, y, text, size, anchor="Center"):
    node = ET.SubElement(page, "t", {"id": str(ids.next()), "Z": str(ids.next()),
                                     "p": "%.2f %.2f" % (x, y),
                                     "InterpretChemically": "no", "Justification": anchor,
                                     "CaptionJustification": anchor})
    for piece, face in runs(text):
        ET.SubElement(node, "s", {"font": "3", "size": str(size), "face": face}).text = piece
    return node


def compose_page(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Lay out ``rows`` of CDXML ``parts`` as one page.

    Spec: ``{"rows": [{"parts": [path, ...], "gap": 46, "margin_label": "(a)",
    "labels": ["(a)", ""], "arrows": [[0, 1], ["equilibrium", 0, 1]],
    "curves": [{"centre": [x, y], "radius": r, "from": deg, "to": deg}]}], "row_gap": 22}``.
    """
    rows = spec.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError('compose page needs "rows": [{"parts": [file, ...]}]')
    size = float(spec.get("font_size", 10.0))
    row_gap = float(spec.get("row_gap", 22.0))
    label_offset = float(spec.get("label_offset", 12.0))
    ids = _Ids()
    first = ET.parse(rows[0]["parts"][0]).getroot()
    root = ET.Element(first.tag, dict(first.attrib))
    for child in first:
        if child.tag != "page":
            root.append(copy.deepcopy(child))
    page = ET.Element("page", dict(first.find("page").attrib))
    root.append(page)
    y_cursor, report, warnings = 0.0, [], []
    for row in rows:
        parts = row.get("parts") or []
        if not parts:
            raise ValueError("a row needs at least one part")
        gap = float(row.get("gap", 46.0))
        boxes, spines, x = [], [], 56.0 if row.get("margin_label") else 14.0
        for path in parts:
            source = ET.parse(path).getroot().find("page")
            box = _bbox(source)
            dx, dy = x - box[0], y_cursor - box[1]
            for child in list(source):
                clone = copy.deepcopy(child)
                _shift(clone, dx, dy, ids)
                page.append(clone)
            boxes.append((box[0] + dx, box[1] + dy, box[2] + dx, box[3] + dy))
            spines.append(_spine_y(source) + dy)
            x = box[2] + dx + gap
        baseline = max(b[3] for b in boxes) + label_offset
        if row.get("margin_label"):
            _text(page, ids, 26.0, (spines[0] + spines[-1]) / 2.0, row["margin_label"], size, "Left")
        for box, text in zip(boxes, row.get("labels") or []):
            if text:
                _text(page, ids, (box[0] + box[2]) / 2.0, baseline, text, size)
        margin = max(10.0, 1.2 * size)
        for item in row.get("arrows") or []:
            kind, a, b = ("resonance", int(item[0]), int(item[1])) if len(item) == 2 else (item[0], int(item[1]), int(item[2]))
            left, right = boxes[a], boxes[b]
            if not (0 <= a < len(boxes) and 0 <= b < len(boxes)):
                raise ValueError("an arrow names a part the row does not have")
            y = (spines[a] + spines[b]) / 2.0
            x0, x1 = left[2] + margin, right[0] - margin
            if x1 <= x0:
                warnings.append("the gap between two parts is too small for an arrow; widen gap")
                continue

            def arrow(at, head, tail):
                ET.SubElement(page, "arrow", {
                    "id": str(ids.next()), "Z": str(ids.next()),
                    "Tail3D": "%.2f %.2f 0" % (x0, at), "Head3D": "%.2f %.2f 0" % (x1, at),
                    "ArrowheadHead": head, "ArrowheadTail": tail, "ArrowheadType": "Solid",
                    "HeadSize": "700", "ArrowheadWidth": "180", "FillType": "None",
                    "LineWidth": "0.6"})
            if kind == "equilibrium":
                arrow(y - 4.0, "HalfLeft", "None")
                arrow(y + 4.0, "None", "HalfLeft")
            elif kind == "resonance":
                arrow(y, "Full", "Full")
            else:
                raise ValueError('an arrow kind is "resonance" or "equilibrium", not %r' % (kind,))
        for curve in row.get("curves") or []:
            cx, cy = curve["centre"]
            radius = float(curve["radius"])
            start, end = math.radians(float(curve["from"])), math.radians(float(curve["to"]))
            electrons = int(curve.get("electrons", 2))
            ET.SubElement(page, "arrow", {
                "id": str(ids.next()), "Z": str(ids.next()),
                "Tail3D": "%.2f %.2f 0" % (cx + radius * math.cos(start), cy + radius * math.sin(start)),
                "Head3D": "%.2f %.2f 0" % (cx + radius * math.cos(end), cy + radius * math.sin(end)),
                "Center3D": "%.2f %.2f 0" % (cx, cy),
                "AngularSize": "%.2f" % (float(curve["to"]) - float(curve["from"])),
                "ArrowheadHead": "HalfLeft" if electrons == 1 else "Full", "ArrowheadTail": "None",
                "ArrowheadType": "Solid",
                "HeadSize": "600" if electrons == 1 else "1000",
                "ArrowheadCenterSize": "525" if electrons == 1 else "875",
                "ArrowheadWidth": "150" if electrons == 1 else "250",
                "FillType": "None", "LineWidth": "0.6"})
        report.append({"panels": [list(map(float, b)) for b in boxes], "baseline": baseline})
        y_cursor = baseline + row_gap
    from .energy import _save_text
    text = ET.tostring(root, encoding="unicode")
    out = {"diagram": "compose_page", "rows": report, "warnings": warnings, "cdxml": text}
    if output is not None:
        out["output_path"] = _save_text(text, output)
    return out
