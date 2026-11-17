"""Dashed outlines that mark one fragment of a drawn structure (shared by the scheme renderer,
which keeps arrows and captions off them, and the brief tools, which draw them)."""

from __future__ import annotations

import math

from .constants import ACS_BOND_LENGTH

# Clearance between a marked atom and its overlay: clear of a one-letter atom label, inside
# the gap to a neighbouring structure.
OVERLAY_PAD = 0.5 * ACS_BOND_LENGTH
# Handle length of a cubic Bezier quarter of a 45-degree arc, as ChemDraw writes a circle.
_ARC_KAPPA = 4.0 / 3.0 * math.tan(math.pi / 16.0)
# Half the dashed line and a little air: what a box around an outline adds beyond it.
OVERLAY_LINE = 1.0


def _closed_curve(anchors: list[tuple[float, float]], tangents: list[tuple[float, float]]) -> str:
    """(in-handle, anchor, out-handle) triples of a closed curve, the way ChemDraw writes one:
    no Closed attribute, the last out-handle leading back into the first anchor."""
    pts = []
    for (x, y), (tx, ty) in zip(anchors, tangents):
        pts += [(x - tx, y - ty), (x, y), (x + tx, y + ty)]
    return " ".join(f"{x:.2f} {y:.2f}" for x, y in pts)


def _circle(cx: float, cy: float, r: float) -> str:
    angles = [math.pi / 4.0 * k for k in range(8)]
    anchors = [(cx + r * math.cos(a), cy + r * math.sin(a)) for a in angles]
    tangents = [(-_ARC_KAPPA * r * math.sin(a), _ARC_KAPPA * r * math.cos(a)) for a in angles]
    return _closed_curve(anchors, tangents)


def _hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list = []
    upper: list = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 1e-9:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 1e-9:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def overlay_curve_points(points: list[tuple[float, float]], shape: str = "circle",
                         radius: float | None = None) -> str:
    """CurvePoints of a closed outline around ``points``.

    circle: centred on the points' centroid, OVERLAY_PAD outside the farthest one (or of the
    given radius). outline: the convex hull pushed OVERLAY_PAD outwards, its corners rounded
    into arcs; fewer than three distinct corners make a circle.
    """
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    hull = _hull(points) if shape == "outline" else []
    if len(hull) < 3:
        if radius is None:
            radius = max(math.hypot(x - cx, y - cy) for x, y in points) + OVERLAY_PAD
        return _circle(cx, cy, radius)
    anchors: list[tuple[float, float]] = []
    tangents: list[tuple[float, float]] = []
    n = len(hull)
    edges = []
    for i in range(n):
        (x1, y1), (x2, y2) = hull[i], hull[(i + 1) % n]
        length = math.hypot(x2 - x1, y2 - y1)
        edges.append(((x2 - x1) / length, (y2 - y1) / length))
    # Walk the hull: at each corner an arc of radius OVERLAY_PAD from the incoming edge's
    # outward normal to the outgoing one, then straight on along the next edge.
    sign = 1.0 if sum((hull[i][0] * hull[(i + 1) % n][1] - hull[(i + 1) % n][0] * hull[i][1])
                      for i in range(n)) > 0 else -1.0
    for i in range(n):
        x, y = hull[i]
        a_in = math.atan2(edges[i - 1][1], edges[i - 1][0]) - sign * math.pi / 2.0
        a_out = math.atan2(edges[i][1], edges[i][0]) - sign * math.pi / 2.0
        sweep = (a_out - a_in) * sign % (2.0 * math.pi) * sign
        steps = max(1, math.ceil(abs(sweep) / (math.pi / 4.0) - 1e-9))
        k = 4.0 / 3.0 * math.tan(abs(sweep) / steps / 4.0) * OVERLAY_PAD
        for j in range(steps + 1):
            a = a_in + sweep * j / steps
            anchors.append((x + OVERLAY_PAD * math.cos(a), y + OVERLAY_PAD * math.sin(a)))
            dx, dy = -math.sin(a) * sign, math.cos(a) * sign
            tangents.append((dx * k, dy * k))
    # Between corners the curve runs straight on: an arc's last handle and the next arc's first
    # both lie along the shared edge, so that segment's four control points are collinear.
    return _closed_curve(anchors, tangents)


def overlay_box(points: list[tuple[float, float]], shape: str = "circle") -> tuple[float, float, float, float]:
    """Box that holds the outline drawn around ``points``, its line included.

    The control points of a Bezier chain enclose the curve, so their box encloses it too.
    """
    v = [float(x) for x in overlay_curve_points(points, shape=shape).split()]
    xs, ys = v[0::2], v[1::2]
    return (min(xs) - OVERLAY_LINE, min(ys) - OVERLAY_LINE, max(xs) + OVERLAY_LINE, max(ys) + OVERLAY_LINE)
