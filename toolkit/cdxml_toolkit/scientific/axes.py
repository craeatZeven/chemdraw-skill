"""A framed x-y plot area on a :class:`Canvas`: data-to-page mapping, ticks and axis titles."""
from __future__ import annotations

import math

from .canvas import Canvas
from .energy import MINUS


def nice_ticks(low, high, count=5):
    """Round tick values covering ``[low, high]`` (1, 2, 2.5 or 5 times a power of ten)."""
    if high <= low:
        raise ValueError('An axis needs a range with high > low')
    raw = (high - low) / max(count - 1, 1)
    power = 10 ** math.floor(math.log10(raw))
    step = min((m * power for m in (1, 2, 2.5, 5, 10)), key=lambda v: abs(v - raw))
    first = math.ceil(low / step - 1e-9) * step
    ticks = []
    value = first
    while value <= high + step * 1e-9:
        ticks.append(round(value, 10))
        value += step
    return ticks, step


def tick_text(value, step):
    decimals = max(0, -math.floor(math.log10(step) + 1e-9)) if step < 1 else 0
    if step == 2.5 * 10 ** math.floor(math.log10(step)) and step < 10:
        decimals += 1
    text = f'{value:.{decimals}f}'
    if text in ('-0', '-0.0', '-0.00'):
        text = text[1:]
    return text.replace('-', MINUS)


class Plot:
    """Map ``(x, y)`` data into a ``width`` by ``height`` frame whose top-left corner is the origin."""

    def __init__(self, canvas: Canvas, x_range, y_range, width=220.0, height=160.0, x_ticks=None, y_ticks=None):
        self.c = canvas
        self.x0, self.x1 = (float(v) for v in x_range)
        self.y0, self.y1 = (float(v) for v in y_range)
        if self.x1 <= self.x0 or self.y1 <= self.y0:
            raise ValueError('An axis range must increase')
        self.w, self.h = float(width), float(height)
        self.x_ticks = x_ticks
        self.y_ticks = y_ticks

    def at(self, x, y):
        return ((x - self.x0) / (self.x1 - self.x0) * self.w, self.h - (y - self.y0) / (self.y1 - self.y0) * self.h)

    def inside(self, x, y):
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1

    def frame(self, x_title='', y_title='', ticks=True, size=8):
        c = self.c
        c.rect(0.0, 0.0, self.w, self.h)
        for edge in ((0, 0, self.w, 0), (0, self.h, self.w, self.h), (0, 0, 0, self.h), (self.w, 0, self.w, self.h)):
            c.block(edge, 1.5)
        if ticks:
            xt, xs = nice_ticks(self.x0, self.x1) if self.x_ticks is None else (self.x_ticks, _step(self.x_ticks))
            yt, ys = nice_ticks(self.y0, self.y1) if self.y_ticks is None else (self.y_ticks, _step(self.y_ticks))
            for value in xt:
                px, _ = self.at(value, self.y0)
                c.segment((px, self.h), (px, self.h - 3.0))
                c.label((px, self.h + 10.0), tick_text(value, xs), align='Center', size=size)
            for value in yt:
                _, py = self.at(self.x0, value)
                c.segment((0.0, py), (3.0, py))
                c.label((-4.0, py), tick_text(value, ys), align='Right', valign='middle', size=size)
        if x_title:
            c.label((self.w / 2, self.h + 24.0), x_title, align='Center')
        if y_title:
            left = min([b[0] for b in c.text_boxes] + [0.0]) - 6.0
            c.label((left, -10.0), y_title, align='Left')

    def polyline(self, points, dashed=False, color=None, clip=True):
        """Draw data points as one or more lines, cut where they leave the frame."""
        runs, run = [], []
        for x, y in points:
            if clip and not self.inside(x, y):
                if len(run) > 1:
                    runs.append(run)
                run = []
                continue
            run.append(self.at(x, y))
        if len(run) > 1:
            runs.append(run)
        for r in runs:
            self.c.path(r, dashed=dashed, color=color)
            self.c.block_path(r, 1.5)
        return runs

    def segment(self, a, b, dashed=False, color=None):
        """A data-space straight line clipped to the frame (Liang-Barsky)."""
        clipped = clip_segment(a, b, (self.x0, self.y0, self.x1, self.y1))
        if clipped is None:
            return None
        pa, pb = self.at(*clipped[0]), self.at(*clipped[1])
        self.c.segment(pa, pb, dashed=dashed, color=color)
        self.c.block_path([pa, pb], 1.5)
        return clipped


def _step(ticks):
    return abs(ticks[1] - ticks[0]) if len(ticks) > 1 else 1.0


def clip_segment(a, b, box):
    (x0, y0), (x1, y1) = a, b
    xmin, ymin, xmax, ymax = box
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - xmin), (dx, xmax - x0), (-dy, y0 - ymin), (dy, ymax - y0)):
        if abs(p) < 1e-15:
            if q < 0:
                return None
            continue
        r = q / p
        if p < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return None
    return (x0 + t0 * dx, y0 + t0 * dy), (x0 + t1 * dx, y0 + t1 * dy)


def clip_polygon(polygon, a, b, c):
    """Keep the part of a convex polygon where ``a x + b y <= c`` (Sutherland-Hodgman, one edge)."""
    out = []
    n = len(polygon)
    for i in range(n):
        p, q = polygon[i], polygon[(i + 1) % n]
        fp, fq = a * p[0] + b * p[1] - c, a * q[0] + b * q[1] - c
        if fp <= 1e-12:
            out.append(p)
        if (fp < -1e-12 < fq) or (fq < -1e-12 < fp):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out
