"""Potential-energy-surface contour plots with their stationary points located.

The surface is a grid of values (from a calculation) or the Mueller-Brown model potential.
Contours come from marching squares; minima and first-order saddle points (transition states)
are found from the grid and refined by Newton steps on a bicubic-free local quadratic fit, then
reported with their energies, so the drawing and the numbers agree.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

from .axes import Plot
from .canvas import BLUE, ORANGE, Canvas
from .common import number
from .energy import _result

_MB = ((-200.0, -1.0, 0.0, -10.0, 1.0, 0.0), (-100.0, -1.0, 0.0, -10.0, 0.0, 0.5),
       (-170.0, -6.5, 11.0, -6.5, -0.5, 1.5), (15.0, 0.7, 0.6, 0.7, -1.0, 1.0))


def mueller_brown(x, y):
    """The Mueller-Brown potential (Theor. Chim. Acta 1979, 53, 75): three minima, two saddles."""
    return sum(A * math.exp(a * (x - x0) ** 2 + b * (x - x0) * (y - y0) + c * (y - y0) ** 2)
               for A, a, b, c, x0, y0 in _MB)


def _interp(p, q, vp, vq, level):
    t = (level - vp) / (vq - vp) if vq != vp else 0.5
    return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))


def contour_segments(xs, ys, z, level):
    """Line segments of one contour level (marching squares, saddle cells split by the centre value)."""
    segments = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            corners = [(xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]), (xs[i], ys[j + 1])]
            values = [z[j][i], z[j][i + 1], z[j + 1][i + 1], z[j + 1][i]]
            inside = [v < level for v in values]
            crossings = []
            for k in range(4):
                a, b = k, (k + 1) % 4
                if inside[a] != inside[b]:
                    crossings.append(_interp(corners[a], corners[b], values[a], values[b], level))
            if len(crossings) == 2:
                segments.append(tuple(crossings))
            elif len(crossings) == 4:
                centre = sum(values) / 4
                if (centre < level) == inside[0]:
                    segments += [(crossings[0], crossings[1]), (crossings[2], crossings[3])]
                else:
                    segments += [(crossings[3], crossings[0]), (crossings[1], crossings[2])]
    return segments


def _chain(segments, tol=1e-9):
    """Join contour segments into polylines."""
    remaining = list(segments)
    lines = []
    while remaining:
        line = list(remaining.pop())
        grown = True
        while grown:
            grown = False
            for k, (a, b) in enumerate(remaining):
                if math.dist(a, line[-1]) < tol:
                    line.append(b)
                elif math.dist(b, line[-1]) < tol:
                    line.append(a)
                elif math.dist(b, line[0]) < tol:
                    line.insert(0, a)
                elif math.dist(a, line[0]) < tol:
                    line.insert(0, b)
                else:
                    continue
                remaining.pop(k)
                grown = True
                break
        lines.append(line)
    return lines


def _stationary(xs, ys, z):
    """Grid minima and saddle points (interior points only), classified by the discrete Hessian."""
    found = []
    for i in range(1, len(xs) - 1):
        for j in range(1, len(ys) - 1):
            v = z[j][i]
            ring = [z[j - 1][i - 1], z[j - 1][i], z[j - 1][i + 1], z[j][i + 1], z[j + 1][i + 1], z[j + 1][i],
                    z[j + 1][i - 1], z[j][i - 1]]
            if all(v < w for w in ring):
                found.append(('minimum', i, j))
                continue
            # A saddle: going round the neighbours, the sign of (neighbour - v) changes four times.
            signs = [w > v for w in ring]
            changes = sum(1 for k in range(8) if signs[k] != signs[(k + 1) % 8])
            if changes >= 4:
                hx, hy = xs[i + 1] - xs[i], ys[j + 1] - ys[j]
                fxx = (z[j][i + 1] - 2 * v + z[j][i - 1]) / hx ** 2
                fyy = (z[j + 1][i] - 2 * v + z[j - 1][i]) / hy ** 2
                fxy = (z[j + 1][i + 1] - z[j + 1][i - 1] - z[j - 1][i + 1] + z[j - 1][i - 1]) / (4 * hx * hy)
                gx = (z[j][i + 1] - z[j][i - 1]) / (2 * hx)
                gy = (z[j + 1][i] - z[j - 1][i]) / (2 * hy)
                det = fxx * fyy - fxy ** 2
                if det < 0:
                    # Newton step to the stationary point of the local quadratic: stay within the cell.
                    dx = -(fyy * gx - fxy * gy) / det
                    dy = -(fxx * gy - fxy * gx) / det
                    if abs(dx) <= hx * 0.75 and abs(dy) <= hy * 0.75:
                        found.append(('saddle', i, j))
    return found


def _refine(fn, x, y, kind, h=1e-4):
    """Newton iterations on an analytic surface."""
    for _ in range(40):
        gx = (fn(x + h, y) - fn(x - h, y)) / (2 * h)
        gy = (fn(x, y + h) - fn(x, y - h)) / (2 * h)
        fxx = (fn(x + h, y) - 2 * fn(x, y) + fn(x - h, y)) / h ** 2
        fyy = (fn(x, y + h) - 2 * fn(x, y) + fn(x, y - h)) / h ** 2
        fxy = (fn(x + h, y + h) - fn(x + h, y - h) - fn(x - h, y + h) + fn(x - h, y - h)) / (4 * h * h)
        det = fxx * fyy - fxy ** 2
        if abs(det) < 1e-14:
            break
        dx = -(fyy * gx - fxy * gy) / det
        dy = -(fxx * gy - fxy * gx) / det
        x, y = x + dx, y + dy
        if math.hypot(dx, dy) < 1e-10:
            break
    fxx = (fn(x + h, y) - 2 * fn(x, y) + fn(x - h, y)) / h ** 2
    fyy = (fn(x, y + h) - 2 * fn(x, y) + fn(x, y - h)) / h ** 2
    fxy = (fn(x + h, y + h) - fn(x + h, y - h) - fn(x - h, y + h) + fn(x - h, y - h)) / (4 * h * h)
    det = fxx * fyy - fxy ** 2
    actual = 'minimum' if det > 0 and fxx > 0 else 'saddle' if det < 0 else 'other'
    return x, y, actual


def potential_surface(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Contour plot of a 2D potential-energy surface.

    ``{"preset": "mueller_brown"}`` or ``{"grid": {"x": [...], "y": [...], "z": [[z for each x]
    for each y]}}`` (e.g. a relaxed scan). ``levels``: a list of values or a count (default
    14); ``cap`` drops levels above a value; ``x_label``, ``y_label``; ``mark``: label the
    stationary points (default true); ``path``: ``[[x, y], ...]`` drawn on top (e.g. an IRC).
    """
    analytic = None
    if spec.get('preset') is not None:
        if str(spec['preset']).lower() not in ('mueller_brown', 'muller_brown'):
            raise ValueError('The surface preset is "mueller_brown"')
        analytic = mueller_brown
        n = int(spec.get('resolution', 90))
        x_range, y_range = spec.get('x', [-1.6, 1.2]), spec.get('y', [-0.4, 2.0])
        xs = [x_range[0] + (x_range[1] - x_range[0]) * k / (n - 1) for k in range(n)]
        ys = [y_range[0] + (y_range[1] - y_range[0]) * k / (n - 1) for k in range(n)]
        z = [[analytic(x, y) for x in xs] for y in ys]
        cap = spec.get('cap', 50.0)
    else:
        grid = spec.get('grid') or {}
        xs = [number(v, 'grid x', -1e9, 1e9) for v in grid.get('x', [])]
        ys = [number(v, 'grid y', -1e9, 1e9) for v in grid.get('y', [])]
        z = [[number(v, 'grid z', -1e12, 1e12) for v in row] for row in grid.get('z', [])]
        if len(xs) < 3 or len(ys) < 3 or len(z) != len(ys) or any(len(row) != len(xs) for row in z):
            raise ValueError('"grid" needs x (n values), y (m values) and z (m rows of n values)')
        if any(b <= a for a, b in zip(xs, xs[1:])) or any(b <= a for a, b in zip(ys, ys[1:])):
            raise ValueError('grid x and y must increase')
        cap = spec.get('cap')
    flat = [v for row in z for v in row]
    low = min(flat)
    high = min(max(flat), float(cap)) if cap is not None else max(flat)
    levels = spec.get('levels', 14)
    if isinstance(levels, int):
        levels = [low + (high - low) * (k + 0.5) / levels for k in range(levels)]
    levels = sorted(float(v) for v in levels)

    c = Canvas()
    plot = Plot(c, (xs[0], xs[-1]), (ys[0], ys[-1]), float(spec.get('width', 220.0)), float(spec.get('height', 190.0)))
    plot.frame(spec.get('x_label', 'q{sub:1}'), spec.get('y_label', 'q{sub:2}'))
    drawn = 0
    for level in levels:
        for line in _chain(contour_segments(xs, ys, z, level)):
            if len(line) > 1:
                pts = [plot.at(x, y) for x, y in line]
                c.path(pts)
                drawn += 1
    if spec.get('path'):
        plot.polyline([(float(p[0]), float(p[1])) for p in spec['path']], color=ORANGE)
    points = []
    for kind, i, j in _stationary(xs, ys, z):
        x, y, value = xs[i], ys[j], z[j][i]
        if analytic is not None:
            x, y, kind = _refine(analytic, x, y, kind)
            value = analytic(x, y)
        if kind == 'other' or not plot.inside(x, y):
            continue
        if any(math.dist((x, y), (p['x'], p['y'])) < 1e-3 for p in points):
            continue
        points.append({'kind': 'transition state' if kind == 'saddle' else 'minimum', 'x': round(x, 4),
                       'y': round(y, 4), 'energy': round(value, 3)})
    points.sort(key=lambda p: p['energy'])
    minima = [p for p in points if p['kind'] == 'minimum']
    saddles = [p for p in points if p['kind'] != 'minimum']
    for k, p in enumerate(minima, start=1):
        p['name'] = f'M{k}'
    for k, p in enumerate(saddles, start=1):
        p['name'] = f'TS{k}'
    if spec.get('mark', True):
        for p in points:
            px, py = plot.at(p['x'], p['y'])
            if p['kind'] == 'minimum':
                c.circle((px, py), 1.8, fill=BLUE, outline=None)
            else:
                c.segment((px - 2.5, py - 2.5), (px + 2.5, py + 2.5), color=ORANGE)
                c.segment((px - 2.5, py + 2.5), (px + 2.5, py - 2.5), color=ORANGE)
            name = p['name'] + ('{sup:‡}' if p['kind'] != 'minimum' else '')
            c.place([((px + 4.0, py - 4.0), 'Left', 'baseline'), ((px - 4.0, py - 4.0), 'Right', 'baseline'),
                     ((px + 4.0, py + 10.0), 'Left', 'baseline'), ((px - 4.0, py + 10.0), 'Right', 'baseline')],
                    name, size=8)
    report = {'diagram': 'potential_surface', 'levels': len(levels), 'contour_lines': drawn,
              'stationary_points': points, 'warnings': [] if drawn else ['no contour line falls inside the plot']}
    if minima and saddles:
        report['barriers_from_global_minimum'] = {s['name']: round(s['energy'] - minima[0]['energy'], 3) for s in saddles}
    return _result(c, output, report)
