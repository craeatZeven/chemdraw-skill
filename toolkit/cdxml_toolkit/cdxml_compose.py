"""Combine CDXML drawings and add aromatic ring circles, without ChemDraw or RDKit.

``merge_cdxml`` places several documents on one page (an energy curve plus the structures
of its stationary points, panels of a composite figure). Object ids, id references, colour
indices and font ids are renumbered, so every part keeps its own colours, fonts and
chemistry. ``add_ring_circles`` turns Kekule aromatic rings into the circle depiction.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from xml.etree import ElementTree as ET

# Attributes whose value is one id or a list of ids of other objects.
ID_REFERENCES = (
    'B', 'E', 'Attachments', 'BracketedObjectIDs', 'GraphicID', 'CrossingBonds', 'BondOrdering', 'object',
    'ReactionStepReactants', 'ReactionStepProducts', 'ReactionStepArrows', 'ReactionStepPlusses',
    'ReactionStepObjectsAboveArrow', 'ReactionStepObjectsBelowArrow', 'ReactionStepAtomMap', 'SupersededBy',
    'BasisObjects', 'ConstraintType_objects', 'ChemicalPropertyBasisObjects',
)
COLOR_ATTRS = ('color', 'bgcolor', 'fgcolor', 'HighlightColor')
FONT_ATTRS = ('font', 'LabelFont', 'CaptionFont')
PAIR_ATTRS = ('p', 'BoundingBox', 'CurvePoints')
TRIPLE_ATTRS = ('Head3D', 'Tail3D', 'Center3D', 'MajorAxisEnd3D', 'MinorAxisEnd3D')
STYLE_ATTRS = ('BondLength', 'LineWidth', 'BoldWidth', 'LabelSize', 'CaptionSize', 'HashSpacing', 'MarginWidth')
DOCTYPE = '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >'
ANCHORS = {'top-left': (0.0, 0.0), 'top': (0.5, 0.0), 'top-right': (1.0, 0.0), 'left': (0.0, 0.5),
           'center': (0.5, 0.5), 'right': (1.0, 0.5), 'bottom-left': (0.0, 1.0), 'bottom': (0.5, 1.0),
           'bottom-right': (1.0, 1.0)}


def parse_cdxml(source) -> ET.Element:
    """Root element of a CDXML file path or CDXML text."""
    if isinstance(source, ET.Element):
        return source
    text = str(source)
    if not text.lstrip().startswith('<'):
        text = Path(text).read_text(encoding='utf-8-sig')
    text = re.sub(r'<!DOCTYPE[^>]*>', '', text, count=1)
    return ET.fromstring(text.encode('utf-8'))


def serialize(root: ET.Element) -> str:
    return '<?xml version="1.0" encoding="UTF-8" ?>\n' + DOCTYPE + '\n' + ET.tostring(root, encoding='unicode')


def _numbers(value):
    return [float(v) for v in value.split()]


def content_box(elements: Iterable[ET.Element]):
    """Bounding box ``(x1, y1, x2, y2)`` of everything positioned inside ``elements``."""
    xs, ys = [], []
    for element in elements:
        for node in element.iter():
            for attr in PAIR_ATTRS:
                if node.get(attr):
                    values = _numbers(node.get(attr))
                    xs += values[0::2]; ys += values[1::2]
            for attr in TRIPLE_ATTRS:
                if node.get(attr):
                    values = _numbers(node.get(attr))
                    xs += values[0::3]; ys += values[1::3]
    if not xs:
        raise ValueError('The drawing has no positioned objects')
    return min(xs), min(ys), max(xs), max(ys)


def translate(elements: Iterable[ET.Element], dx: float, dy: float) -> None:
    """Move everything inside ``elements`` by ``(dx, dy)``."""
    for element in elements:
        for node in element.iter():
            for attr, step in [(a, 2) for a in PAIR_ATTRS] + [(a, 3) for a in TRIPLE_ATTRS]:
                value = node.get(attr)
                if not value:
                    continue
                values = _numbers(value)
                for i in range(0, len(values) - 1, step):
                    values[i] += dx
                    values[i + 1] += dy
                node.set(attr, ' '.join(f'{v:.2f}'.rstrip('0').rstrip('.') for v in values))


def _color_table(root):
    """Colour index -> ``(r, g, b)``: 0 is black, 1 white, the table starts at 2."""
    table = {0: (0.0, 0.0, 0.0), 1: (1.0, 1.0, 1.0)}
    colortable = root.find('colortable')
    if colortable is not None:
        for i, color in enumerate(colortable.findall('color'), start=2):
            table[i] = tuple(round(float(color.get(k, '0')), 4) for k in 'rgb')
    return table


def _font_table(root):
    fonttable = root.find('fonttable')
    return {} if fonttable is None else {font.get('id'): (font.get('name'), font.get('charset'))
                                         for font in fonttable.findall('font')}


def _fmt(value):
    return f'{value:.4f}'.rstrip('0').rstrip('.') if value not in (0.0, 1.0) else str(int(value))


def merge_cdxml(parts: Sequence[Mapping[str, Any]], output=None, *, gap: float = 20.0, margin: float = 12.0) -> dict[str, Any]:
    """Put several CDXML drawings on one page.

    Each part is ``{"file": path or "cdxml": text, "at": [x, y], "anchor": "top-left"}``:
    the part's ``anchor`` point (top-left, top, center, bottom, ...) lands on ``at``. A part
    without ``at`` is placed to the right of the previous one. Document settings come from
    the first part; parts drawn with other settings are reported, not rescaled.
    """
    if not parts:
        raise ValueError('Give at least one part to merge')
    merged = None
    page = None
    colors: dict[tuple, int] = {}
    fonts: dict[tuple, str] = {}
    next_id = 2          # the page is object 1
    z_base = 0
    cursor_x = 0.0
    warnings: list[str] = []
    placed = []
    for index, part in enumerate(parts, start=1):
        source = part.get('cdxml') if part.get('cdxml') is not None else part.get('file')
        if source is None:
            raise ValueError(f'part {index}: give "file" or "cdxml"')
        root = parse_cdxml(source)
        if root.tag != 'CDXML':
            raise ValueError(f'part {index}: not a CDXML document')
        children = [child for pg in root.findall('page') for child in list(pg)]
        if not children:
            raise ValueError(f'part {index}: the document has no objects')
        if merged is None:
            merged = ET.Element('CDXML', {k: v for k, v in root.attrib.items() if k != 'BoundingBox'})
            ET.SubElement(merged, 'colortable')
            ET.SubElement(merged, 'fonttable')
            page = ET.SubElement(merged, 'page', {'id': '1'})
            colors.update({(1.0, 1.0, 1.0): 2, (0.0, 0.0, 0.0): 3})
            for rgb in ((1.0, 1.0, 1.0), (0.0, 0.0, 0.0)):
                ET.SubElement(merged.find('colortable'), 'color', dict(zip('rgb', map(_fmt, rgb))))
        else:
            differ = [f'{k} {root.get(k)} (first part: {merged.get(k)})' for k in STYLE_ATTRS
                      if root.get(k) and merged.get(k) and abs(float(root.get(k)) - float(merged.get(k))) > 1e-6]
            if differ:
                warnings.append(f'part {index} was drawn with other settings: ' + ', '.join(differ))

        # ---- colours and fonts keep their meaning under new indices
        table = _color_table(root)
        font_names = _font_table(root)
        for node in (n for child in children for n in child.iter()):
            for attr in COLOR_ATTRS:
                value = node.get(attr)
                if value is None or not value.isdigit():
                    continue
                old = int(value)
                if old < 2 or old not in table:
                    continue
                rgb = table[old]
                if rgb not in colors:
                    colors[rgb] = len(colors) + 2
                    ET.SubElement(merged.find('colortable'), 'color', dict(zip('rgb', map(_fmt, rgb))))
                node.set(attr, str(colors[rgb]))
            for attr in FONT_ATTRS:
                value = node.get(attr)
                if value is None:
                    continue
                key = font_names.get(value, (None, None))
                if key[0] is None:
                    continue
                if key not in fonts:
                    wanted = value if value not in fonts.values() else str(max(int(v) for v in fonts.values()) + 1)
                    fonts[key] = wanted
                    ET.SubElement(merged.find('fonttable'), 'font', {'id': wanted, 'charset': key[1] or 'iso-8859-1',
                                                                    'name': key[0]})
                node.set(attr, fonts[key])
        if index == 1:
            for attr in ('LabelFont', 'CaptionFont'):
                key = font_names.get(root.get(attr, ''), (None, None))
                if key[0] is not None:
                    if key not in fonts:
                        fonts[key] = root.get(attr)
                        ET.SubElement(merged.find('fonttable'), 'font', {'id': fonts[key], 'name': key[0],
                                                                        'charset': key[1] or 'iso-8859-1'})
                    merged.set(attr, fonts[key])

        # ---- ids
        mapping = {}
        nodes = [n for child in children for n in child.iter()]
        for node in nodes:
            if node.get('id') is not None:
                mapping[node.get('id')] = str(next_id)
                next_id += 1
        for node in nodes:
            if node.get('id') is not None:
                node.set('id', mapping[node.get('id')])
            for attr in ID_REFERENCES:
                value = node.get(attr)
                if value is None:
                    continue
                tokens = value.split()
                missing = [t for t in tokens if t != '0' and t not in mapping]
                if missing:
                    warnings.append(f'part {index}: {node.tag} {attr} refers to missing object(s) {" ".join(missing)}')
                node.set(attr, ' '.join(mapping.get(t, t) for t in tokens))
            if node.get('Z') is not None:
                try:
                    node.set('Z', str(z_base + int(float(node.get('Z')))))
                except ValueError:
                    pass
        z_values = [int(float(n.get('Z'))) for n in nodes if n.get('Z') and re.fullmatch(r'-?\d+(\.\d+)?', n.get('Z'))]
        z_base = max(z_values + [z_base]) + 10

        # ---- position
        x1, y1, x2, y2 = content_box(children)
        anchor = str(part.get('anchor', 'top-left')).lower()
        if anchor not in ANCHORS:
            raise ValueError(f'part {index}: anchor must be one of {", ".join(ANCHORS)}')
        fx, fy = ANCHORS[anchor]
        if part.get('at') is not None:
            tx, ty = float(part['at'][0]), float(part['at'][1])
        else:
            tx, ty = cursor_x + fx * (x2 - x1), fy * (y2 - y1)
        dx = tx - (x1 + fx * (x2 - x1))
        dy = ty - (y1 + fy * (y2 - y1))
        translate(children, dx, dy)
        box = (x1 + dx, y1 + dy, x2 + dx, y2 + dy)
        cursor_x = max(cursor_x, box[2] + gap)
        placed.append({'part': index, 'box': [round(v, 2) for v in box], 'objects': len(children)})
        for child in children:
            page.append(child)

    # ---- frame the page
    children = list(page)
    x1, y1, x2, y2 = content_box(children)
    translate(children, margin - x1, margin - y1)
    width, height = x2 - x1 + 2 * margin, y2 - y1 + 2 * margin
    for item in placed:
        item['box'] = [round(item['box'][0] + margin - x1, 2), round(item['box'][1] + margin - y1, 2),
                       round(item['box'][2] + margin - x1, 2), round(item['box'][3] + margin - y1, 2)]
    box = f'0 0 {width:.2f} {height:.2f}'
    merged.set('BoundingBox', box)
    page.set('BoundingBox', box)
    page.set('Width', f'{width:.2f}')
    page.set('Height', f'{height:.2f}')
    for a, b in ((a, b) for i, a in enumerate(placed) for b in placed[i + 1:]):
        if (a['box'][0] < b['box'][2] and b['box'][0] < a['box'][2]
                and a['box'][1] < b['box'][3] and b['box'][1] < a['box'][3]):
            warnings.append(f'parts {a["part"]} and {b["part"]} overlap; check the render or move one with "at"')
    text = serialize(merged)
    saved = None
    if output is not None:
        target = Path(output)
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite an existing file: {target}')
        target.write_text(text, encoding='utf-8')
        saved = str(target)
    return {'output_path': saved, 'cdxml': text, 'parts': placed, 'size': [round(width, 2), round(height, 2)],
            'warnings': warnings}


# ---------------------------------------------------------------------------- aromatic circles
def _rings(adjacency, sizes=(5, 6)):
    """Simple cycles of the given sizes, each once, as tuples of node ids."""
    found = {}
    for start in adjacency:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            for nxt in adjacency[node]:
                if nxt == start and len(path) in sizes:
                    found.setdefault(frozenset(path), tuple(path))
                elif nxt not in path and len(path) < max(sizes) and nxt > start:
                    stack.append((nxt, path + [nxt]))
    return list(found.values())


def add_ring_circles(source, output=None, *, fragment_ids: Iterable[str] | None = None, ratio: float = 0.62) -> dict[str, Any]:
    """Redraw Kekule aromatic rings with a circle.

    A six-membered ring qualifies when every ring atom carries a double bond inside the fused
    ring system; a five-membered ring when it has two ring double bonds and one N, O or S
    that has none. The ring bonds become single and a circle is added at the ring centre.
    The result is a depiction: read the chemistry from the source file, because a reader
    that ignores the circle sees a saturated ring.
    """
    root = parse_cdxml(source)
    wanted = None if fragment_ids is None else {str(v) for v in fragment_ids}
    next_id = max([int(n.get('id')) for n in root.iter() if (n.get('id') or '').isdigit()] + [0]) + 1
    circles = 0
    changed_bonds = 0
    report = []
    for page in root.iter('page'):
        for fragment in list(page.iter('fragment')):
            if wanted is not None and fragment.get('id') not in wanted:
                continue
            atoms = {n.get('id'): n for n in fragment.findall('n') if n.get('p')}
            bonds = [b for b in fragment.findall('b') if b.get('B') in atoms and b.get('E') in atoms]
            adjacency = {aid: set() for aid in atoms}
            order = {}
            for bond in bonds:
                a, b = bond.get('B'), bond.get('E')
                adjacency[a].add(b); adjacency[b].add(a)
                order[frozenset((a, b))] = bond
            adjacency = {k: sorted(v) for k, v in adjacency.items()}

            def is_double(a, b):
                return order[frozenset((a, b))].get('Order') in ('2', '1.5')

            def ring_bonds(ring):
                return [frozenset((ring[i], ring[(i + 1) % len(ring)])) for i in range(len(ring))]

            rings = _rings(adjacency)
            ring_atoms = {a for ring in rings for a in ring}
            aromatic = []
            for ring in rings:
                doubles = [pair for pair in ring_bonds(ring) if is_double(*pair)]
                if len(ring) == 6:
                    ok = all(any(is_double(a, n) and n in ring_atoms for n in adjacency[a]) for a in ring) and len(doubles) >= 2
                else:
                    bare = [a for a in ring if not any(is_double(a, n) for n in adjacency[a])]
                    ok = (len(doubles) == 2 and len(bare) == 1
                          and atoms[bare[0]].get('Element') in ('7', '8', '16'))
                if ok:
                    aromatic.append(ring)
            if not aromatic:
                continue
            for ring in aromatic:
                for pair in ring_bonds(ring):
                    bond = order[pair]
                    if bond.get('Order') in ('2', '1.5'):
                        del bond.attrib['Order']
                        bond.attrib.pop('DoublePosition', None)
                        changed_bonds += 1
                points = [tuple(_numbers(atoms[a].get('p'))[:2]) for a in ring]
                cx = sum(p[0] for p in points) / len(points)
                cy = sum(p[1] for p in points) / len(points)
                inradius = min(math.dist((cx, cy), ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2))
                               for p, q in zip(points, points[1:] + points[:1]))
                r = inradius * ratio
                graphic = ET.Element('graphic', {
                    'id': str(next_id), 'GraphicType': 'Oval', 'OvalType': 'Circle',
                    'BoundingBox': f'{cx - r:.2f} {cy - r:.2f} {cx + r:.2f} {cy + r:.2f}',
                    'Center3D': f'{cx:.2f} {cy:.2f} 0', 'MajorAxisEnd3D': f'{cx + r:.2f} {cy:.2f} 0',
                    'MinorAxisEnd3D': f'{cx:.2f} {cy + r:.2f} 0'})
                next_id += 1
                parent = next(p for p in root.iter() if fragment in list(p))
                parent.insert(list(parent).index(fragment) + 1, graphic)
                circles += 1
            report.append({'fragment': fragment.get('id'), 'rings': len(aromatic)})
    text = serialize(root)
    saved = None
    if output is not None:
        target = Path(output)
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite an existing file: {target}')
        target.write_text(text, encoding='utf-8')
        saved = str(target)
    return {'output_path': saved, 'cdxml': text, 'circles': circles, 'bonds_made_single': changed_bonds,
            'fragments': report,
            'note': 'depiction only: the circled rings read back as saturated; keep the source for chemistry'}
