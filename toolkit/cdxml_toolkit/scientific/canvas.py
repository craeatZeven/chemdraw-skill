"""Drawing surface shared by the diagram generators.

Extends :class:`Drawing` with what a labelled diagram needs: rich text with
``{sub:}``/``{sup:}``/``{i:}``/``{b:}`` markup, CJK runs in SimSun, straight arrows,
circles, rectangles, sampled smooth paths and wavy lines. Every primitive is written the
way an existing toolkit workflow already writes it (lines, degenerate-handle curves,
``<arrow>``, ``Oval`` and ``Rectangle`` graphics), so nothing new has to be proven in
native ChemDraw beyond the layout itself.

Coordinates are points, x right, y down. ``finish()`` moves the content to the page
margin and sets the document and page bounding boxes.
"""
from __future__ import annotations

import math
import re
from typing import Iterable, Sequence
from xml.etree import ElementTree as ET

from .common import Drawing

BLACK = '#000000'
BLUE = '#336699'
ORANGE = '#FF6600'
WHITE = '#FFFFFF'
# Named colours a spec may use; anything else must be #RRGGBB.
PALETTE = {'black': BLACK, 'blue': BLUE, 'orange': ORANGE, 'white': WHITE}
# Default colour of the first, second and third series (pathway, species, phase).
SERIES_COLORS = (BLACK, BLUE, ORANGE)
FONT_SIZE = 10.0
CJK_FONT = '5'
LATIN_FONT = '3'

_TOKEN = re.compile(r'\{(i|b|bi|sub|sup):([^{}]*)\}')
_FACE = {'b': 1, 'i': 2, 'bi': 3, 'sub': 32, 'sup': 64}
_CJK = re.compile(r'[\u3000-\u303f\u3400-\u9fff\uff00-\uffef]')
_SPLIT = re.compile(r'[\u3000-\u303f\u3400-\u9fff\uff00-\uffef]+|[^\u3000-\u303f\u3400-\u9fff\uff00-\uffef]+')
_NARROW = set("ijlI|.,:;'!() []")
_WIDE = set('mwMW@%→←⇌')
# Attributes that hold absolute coordinates and move with the content.
_POINT_ATTRS = ('p', 'BoundingBox', 'Head3D', 'Tail3D', 'Center3D', 'MajorAxisEnd3D', 'MinorAxisEnd3D', 'CurvePoints')
_TRIPLE_ATTRS = ('Head3D', 'Tail3D', 'Center3D', 'MajorAxisEnd3D', 'MinorAxisEnd3D')


def resolve_color(value, default=BLACK):
    """A palette name or ``#RRGGBB`` (or ``#RGB``) as upper-case ``#RRGGBB``."""
    if value is None:
        return default
    text = str(value).strip()
    if text.lower() in PALETTE:
        return PALETTE[text.lower()]
    if re.fullmatch(r'#[0-9a-fA-F]{3}', text):
        text = '#' + ''.join(ch * 2 for ch in text[1:])
    if not re.fullmatch(r'#[0-9a-fA-F]{6}', text):
        raise ValueError(f'Unknown colour {value!r}: use {", ".join(PALETTE)} or #RRGGBB')
    return text.upper()


def parse_markup(markup):
    """``'E{sub:a}'`` -> ``[('E', 0), ('a', 32)]`` (text, CDXML face bits)."""
    out, pos = [], 0
    text = str(markup)
    for match in _TOKEN.finditer(text):
        if match.start() > pos:
            out.append((text[pos:match.start()], 0))
        out.append((match.group(2), _FACE[match.group(1)]))
        pos = match.end()
    if pos < len(text):
        out.append((text[pos:], 0))
    return [(t, f) for t, f in out if t]


def plain(markup):
    """The text of a markup string without its formatting tokens."""
    return ''.join(t for t, _ in parse_markup(markup))


def text_width(markup, size=FONT_SIZE):
    """Approximate width of one line of Arial/SimSun text, in points."""
    widest = 0.0
    for line in str(markup).split('\n'):
        width = 0.0
        for text, face in parse_markup(line):
            scale = 0.75 if face & 96 else 1.0
            for ch in text:
                if _CJK.match(ch):
                    unit = 1.0
                elif ch in _NARROW:
                    unit = 0.28
                elif ch in _WIDE:
                    unit = 0.85
                elif ch.isupper() or ch.isdigit():
                    unit = 0.62 if ch.isupper() else 0.56
                else:
                    unit = 0.52
                width += unit * size * scale * (1.06 if face & 1 else 1.0)
        widest = max(widest, width)
    return widest


def text_height(markup, size=FONT_SIZE):
    return (str(markup).count('\n') + 1) * size * 1.2


def boxes_overlap(a, b, gap=0.0):
    return a[0] < b[2] + gap and b[0] < a[2] + gap and a[1] < b[3] + gap and b[1] < a[3] + gap


class Canvas(Drawing):
    """A :class:`Drawing` that tracks its extent and writes rich captions."""

    def __init__(self, font_size=FONT_SIZE):
        super().__init__(10, 10)
        self.font_size = float(font_size)
        self.root.set('ChainAngle', '120')
        self.root.set('HashSpacing', '2.5')
        self.root.set('MarginWidth', '1.6')
        self.root.set('BondSpacing', '18')
        self.root.set('CaptionJustification', 'Left')
        ET.SubElement(self.root.find('fonttable'), 'font', {'id': CJK_FONT, 'charset': 'Unicode', 'name': 'SimSun'})
        self._extent = [math.inf, math.inf, -math.inf, -math.inf]
        self.text_boxes: list[tuple[float, float, float, float, str]] = []
        self.obstacles: list[tuple[float, float, float, float]] = []

    # ------------------------------------------------------------------ bookkeeping
    def see(self, x, y, pad=0.0):
        e = self._extent
        e[0] = min(e[0], x - pad); e[1] = min(e[1], y - pad)
        e[2] = max(e[2], x + pad); e[3] = max(e[3], y + pad)

    def see_box(self, box):
        self.see(box[0], box[1]); self.see(box[2], box[3])

    @property
    def extent(self):
        if self._extent[0] > self._extent[2]:
            return (0.0, 0.0, 0.0, 0.0)
        return tuple(self._extent)

    def _stroke(self, attrs, color, dashed):
        index = self._resolve_color(resolve_color(color)) if color is not None else None
        if index is not None:
            attrs['color'] = index
        if dashed:
            attrs['LineType'] = 'Dashed'
        return attrs

    # ------------------------------------------------------------------ primitives
    def segment(self, a, b, color=None, dashed=False, parent=None):
        """A straight line from ``a`` to ``b``."""
        a = (round(float(a[0]), 2), round(float(a[1]), 2)); b = (round(float(b[0]), 2), round(float(b[1]), 2))
        self.see(*a); self.see(*b)
        return self.line(a, b, parent=parent, color=resolve_color(color) if color is not None else None, dashed=dashed)

    def path(self, points: Iterable[Sequence[float]], color=None, dashed=False, closed=False, filled=False, parent=None):
        """A polyline through ``points`` (an editable curve with degenerate handles)."""
        pts = [(round(float(x), 2), round(float(y), 2)) for x, y in points]
        if len(pts) < 2:
            raise ValueError('A path needs at least two points')
        for p in pts:
            self.see(*p)
        stroke = resolve_color(color) if color is not None else None
        if closed or len(pts) <= 250:
            return [self.curve(pts, parent=parent, closed=closed, filled=filled, color=stroke, dashed=dashed)]
        # Long open paths are split: ChemDraw limits the control points of one curve.
        return [self.curve(pts[i:i + 251], parent=parent, color=stroke, dashed=dashed)
                for i in range(0, len(pts) - 1, 250)]

    def arrow(self, tail, head, color=None, dashed=False, head_style='Full', tail_style=None, head_size=1000,
              parent=None):
        """A straight arrow. ``head_style``/``tail_style``: ``Full``, ``HalfLeft``, ``HalfRight`` or ``None``."""
        x1, y1 = round(float(tail[0]), 2), round(float(tail[1]), 2)
        x2, y2 = round(float(head[0]), 2), round(float(head[1]), 2)
        if math.hypot(x2 - x1, y2 - y1) < 1e-6:
            raise ValueError('An arrow needs two different points')
        self.see(x1, y1, 2); self.see(x2, y2, 2)
        size = int(head_size)
        attrs = {'FillType': 'None', 'ArrowheadHead': head_style or 'None', 'ArrowheadType': 'Solid',
                 'HeadSize': size, 'ArrowheadCenterSize': int(size * 0.875), 'ArrowheadWidth': int(size * 0.25),
                 'Head3D': f'{x2} {y2} 0', 'Tail3D': f'{x1} {y1} 0',
                 'BoundingBox': f'{min(x1, x2)} {min(y1, y2)} {max(x1, x2)} {max(y1, y2)}'}
        if tail_style:
            attrs['ArrowheadTail'] = tail_style
        return self.add('arrow', self._stroke(attrs, color, dashed), parent)

    def circle(self, center, radius, fill=None, outline=BLACK, dashed=False, parent=None):
        """A circle; ``fill`` paints the disc first so it hides what lies behind it."""
        cx, cy, r = round(float(center[0]), 2), round(float(center[1]), 2), round(float(radius), 2)
        self.see(cx - r, cy - r); self.see(cx + r, cy + r)
        geometry = {'BoundingBox': f'{cx - r} {cy - r} {cx + r} {cy + r}', 'Center3D': f'{cx} {cy} 0',
                    'MajorAxisEnd3D': f'{cx + r} {cy} 0', 'MinorAxisEnd3D': f'{cx} {cy + r} 0'}
        nodes = []
        if fill is not None:
            attrs = {'GraphicType': 'Oval', 'OvalType': 'Circle Filled', **geometry}
            nodes.append(self.add('graphic', self._stroke(attrs, fill, False), parent))
        if outline is not None:
            attrs = {'GraphicType': 'Oval', 'OvalType': 'Circle', **geometry}
            nodes.append(self.add('graphic', self._stroke(attrs, outline, dashed), parent))
        return nodes

    def rect(self, x1, y1, x2, y2, color=None, dashed=False, rounded=False, fill=None, parent=None):
        """A rectangle outline, optionally over a solid fill."""
        x1, x2 = sorted((round(float(x1), 2), round(float(x2), 2)))
        y1, y2 = sorted((round(float(y1), 2), round(float(y2), 2)))
        self.see(x1, y1); self.see(x2, y2)
        nodes = []
        if fill is not None:
            attrs = {'GraphicType': 'Rectangle', 'BoundingBox': f'{x1} {y1} {x2} {y2}', 'FillType': 'Solid'}
            if rounded:
                attrs['RectangleType'] = 'RoundEdge'
            nodes.append(self.add('graphic', self._stroke(attrs, fill, False), parent))
        attrs = {'GraphicType': 'Rectangle', 'BoundingBox': f'{x1} {y1} {x2} {y2}'}
        if rounded:
            attrs['RectangleType'] = 'RoundEdge'
        nodes.append(self.add('graphic', self._stroke(attrs, color, dashed), parent))
        return nodes

    def wavy(self, a, b, color=None, amplitude=2.2, wavelength=7.0, arrow=True, parent=None):
        """A wavy line from ``a`` to ``b`` (a non-radiative transition), ending in an arrowhead."""
        x1, y1, x2, y2 = float(a[0]), float(a[1]), float(b[0]), float(b[1])
        length = math.hypot(x2 - x1, y2 - y1)
        if length < 1e-6:
            raise ValueError('A wavy line needs two different points')
        ux, uy = (x2 - x1) / length, (y2 - y1) / length
        tip = min(6.0, length * 0.3) if arrow else 0.0
        body = length - tip
        waves = max(1, round(body / wavelength))
        steps = waves * 12
        points = []
        for i in range(steps + 1):
            t = body * i / steps
            off = amplitude * math.sin(2 * math.pi * waves * i / steps)
            points.append((x1 + ux * t - uy * off, y1 + uy * t + ux * off))
        nodes = self.path(points, color=color, parent=parent)
        if arrow:
            nodes.append(self.arrow((x1 + ux * (body - 0.5), y1 + uy * (body - 0.5)), (x2, y2), color=color,
                                    head_size=800, parent=parent))
        return nodes

    def label(self, xy, markup, size=None, align='Left', valign='baseline', color=None, bold=False, parent=None):
        """Rich caption at ``xy``. Returns its estimated box ``(x1, y1, x2, y2)``.

        ``align`` is ``Left``, ``Center`` or ``Right``; ``valign`` places ``y`` on the first
        ``baseline``, the ``middle`` of the block or its ``top``. CJK runs use SimSun.
        """
        size = float(size or self.font_size)
        markup = str(markup)
        x, y = float(xy[0]), float(xy[1])
        lines = markup.count('\n') + 1
        if valign == 'middle':
            y += size * 0.36 - (lines - 1) * size * 0.6
        elif valign == 'top':
            y += size * 0.8
        elif valign != 'baseline':
            raise ValueError("valign must be 'baseline', 'middle' or 'top'")
        if align not in ('Left', 'Center', 'Right'):
            raise ValueError("align must be 'Left', 'Center' or 'Right'")
        x, y = round(x, 2), round(y, 2)
        attrs = {'p': f'{x} {y}', 'InterpretChemically': 'no', 'Justification': align, 'CaptionJustification': align}
        index = self._resolve_color(resolve_color(color)) if color is not None else None
        if index is not None:
            attrs['color'] = index
        node = self.add('t', attrs, parent)
        base = 1 if bold else 0
        for text, face in parse_markup(markup):
            for piece in _SPLIT.findall(text):
                run = {'font': CJK_FONT if _CJK.match(piece[0]) else LATIN_FONT, 'size': f'{size:g}',
                       'face': str(face | base)}
                if index is not None:
                    run['color'] = index
                ET.SubElement(node, 's', run).text = piece
        width = text_width(markup, size)
        left = x if align == 'Left' else x - width / 2 if align == 'Center' else x - width
        box = (left, y - size * 0.8, left + width, y + size * 0.25 + (lines - 1) * size * 1.2)
        self.see_box(box)
        self.text_boxes.append((*box, plain(markup)))
        return box

    # ------------------------------------------------------------------ label placement
    def measure(self, xy, markup, size=None, align='Left', valign='baseline'):
        """The box :meth:`label` would report, without drawing anything."""
        size = float(size or self.font_size)
        x, y = float(xy[0]), float(xy[1])
        lines = str(markup).count('\n') + 1
        if valign == 'middle':
            y += size * 0.36 - (lines - 1) * size * 0.6
        elif valign == 'top':
            y += size * 0.8
        width = text_width(markup, size)
        left = x if align == 'Left' else x - width / 2 if align == 'Center' else x - width
        return (left, y - size * 0.8, left + width, y + size * 0.25 + (lines - 1) * size * 1.2)

    def block(self, box, pad=1.0):
        """Reserve ``box`` so that :meth:`place` keeps captions out of it."""
        x1, x2 = sorted((float(box[0]), float(box[2])))
        y1, y2 = sorted((float(box[1]), float(box[3])))
        self.obstacles.append((x1 - pad, y1 - pad, x2 + pad, y2 + pad))

    def block_path(self, points, pad=1.5):
        """Reserve the ink of a polyline (one small box per segment)."""
        pts = list(points)
        for a, b in zip(pts, pts[1:]):
            self.block((a[0], a[1], b[0], b[1]), pad)

    def place(self, candidates, markup, size=None, color=None, bold=False, parent=None):
        """Draw ``markup`` at the first candidate ``(xy, align, valign)`` that is free.

        Returns ``(box, clear)``; ``clear`` is False when every candidate met reserved ink
        or an earlier caption and the first one was used anyway.
        """
        chosen, clear = None, True
        for xy, align, valign in candidates:
            box = self.measure(xy, markup, size, align, valign)
            if not any(boxes_overlap(box, other) for other in self.obstacles):
                chosen = (xy, align, valign)
                break
        if chosen is None:
            chosen, clear = candidates[0], False
        box = self.label(chosen[0], markup, size=size, align=chosen[1], valign=chosen[2], color=color, bold=bold,
                         parent=parent)
        self.block(box, 0.5)
        return box, clear

    # ------------------------------------------------------------------ output
    def overlapping_labels(self, gap=0.5):
        """Pairs of captions whose estimated boxes overlap: a layout warning, not a proof."""
        hits = []
        for i, a in enumerate(self.text_boxes):
            for b in self.text_boxes[i + 1:]:
                if boxes_overlap(a, b, -gap):
                    hits.append((a[4], b[4]))
        return hits

    def finish(self, margin=12.0):
        """Move the content to ``margin`` from the page origin and size the page around it."""
        x1, y1, x2, y2 = self.extent
        dx, dy = margin - x1, margin - y1
        for node in self.page.iter():
            for attr in _POINT_ATTRS:
                value = node.get(attr)
                if value is None:
                    continue
                numbers = [float(v) for v in value.split()]
                step = 3 if attr in _TRIPLE_ATTRS else 2
                for i in range(0, len(numbers), step):
                    numbers[i] += dx
                    numbers[i + 1] += dy
                node.set(attr, ' '.join(f'{v:.2f}'.rstrip('0').rstrip('.') for v in numbers))
        width, height = x2 - x1 + 2 * margin, y2 - y1 + 2 * margin
        box = f'0 0 {width:.2f} {height:.2f}'
        self.root.set('BoundingBox', box)
        self.page.set('BoundingBox', box)
        self.page.set('Width', f'{width:.2f}')
        self.page.set('Height', f'{height:.2f}')
        self.text_boxes = [(a + dx, b + dy, c + dx, d + dy, t) for a, b, c, d, t in self.text_boxes]
        self._extent = [margin, margin, width - margin, height - margin]
        return width, height

    def to_string(self):
        return ('<?xml version="1.0" encoding="UTF-8" ?>\n'
                '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n'
                + ET.tostring(self.root, encoding='unicode'))

    def write(self, output=None, margin=12.0):
        """Finish the page and save it; returns ``(saved path or None, CDXML text)``."""
        self.finish(margin)
        return (self.save(output) if output is not None else None), self.to_string()
