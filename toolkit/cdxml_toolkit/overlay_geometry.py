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
# ChemDraw's CurveType flags: 1 closed, 2 dashed (128 filled). A dashed outline is both.
CLOSED_DASHED_CURVE = {"CurveType": "3", "Closed": "yes", "LineType": "Dashed"}
# Half the dashed line and a little air: what a box around an outline adds beyond it.
OVERLAY_LINE = 1.0


def _closed_curve(anchors: list[tuple[float, float]], tangents: list[tuple[float, float]]) -> str:
    """(in-handle, anchor, out-handle) triples of a closed curve, the first one repeated last.

    ChemDraw draws one cubic segment between each pair of consecutive anchors and nothing from
    the last back to the first: its own dashed circle (CurveType 2, no closed bit) renders with
    a gap. So the first triple is written again at the end, and the curve is also marked closed
    (``CLOSED_DASHED_CURVE``) for readers that honour the flag.
    """
    pts = []
    for (x, y), (tx, ty) in zip(anchors, tangents):
        pts += [(x - tx, y - ty), (x, y), (x + tx, y + ty)]
    pts += pts[:3]
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


def _ring_frame(a, b, away, members):
    """Centre and radius of the ring of ``members`` atoms through the bond a-b, on the side of
    the bond away from ``away`` (the rest of the molecule)."""
    chord = math.dist(a, b)
    if chord < 1e-6:
        raise ValueError("a generic ring needs two atoms at different places")
    r = chord / (2.0 * math.sin(math.pi / members))
    h = math.sqrt(max(r * r - chord * chord / 4.0, 0.0))
    mx, my = (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0
    nx, ny = -(b[1] - a[1]) / chord, (b[0] - a[0]) / chord
    if away is not None and nx * (mx - away[0]) + ny * (my - away[1]) < 0:
        nx, ny = -nx, -ny
    return (mx + nx * h, my + ny * h), r


def generic_ring_points(a, b, away=None, members: int = 6, clear_a: float = 0.0, clear_b: float = 0.0) -> str:
    """CurvePoints of the dashed arc that stands for an unspecified ring fused at the bond a-b.

    The arc lies on the circle of a regular ring of ``members`` atoms through a and b, centred
    on the side away from ``away``, and runs the long way round from a to b: the short arc,
    where the solid bond is, is left out (BJOC 22-125 scheme 4). ``clear_a`` and ``clear_b``
    stop it that many points short of a labelled atom. It is an open curve.
    """
    (cx, cy), r = _ring_frame(a, b, away, members)
    ta = math.atan2(a[1] - cy, a[0] - cx)
    tb = math.atan2(b[1] - cy, b[0] - cx)
    delta = (tb - ta) % (2.0 * math.pi)
    if delta < math.pi:          # that is the short way: go round the other side
        delta -= 2.0 * math.pi
    sign = 1.0 if delta > 0 else -1.0
    start = ta + sign * clear_a / r
    sweep = delta - sign * (clear_a + clear_b) / r
    steps = max(1, math.ceil(abs(sweep) / (math.pi / 4.0) - 1e-9))
    k = 4.0 / 3.0 * math.tan(sweep / steps / 4.0) * r
    pts = []
    for i in range(steps + 1):
        t = start + sweep * i / steps
        x, y = cx + r * math.cos(t), cy + r * math.sin(t)
        tx, ty = -math.sin(t) * k, math.cos(t) * k
        into = (x, y) if i == 0 else (x - tx, y - ty)        # an open curve: no handle past its ends
        out = (x, y) if i == steps else (x + tx, y + ty)
        pts += [into, (x, y), out]
    return " ".join(f"{x:.2f} {y:.2f}" for x, y in pts)


def generic_ring_box(a, b, away=None, members: int = 6) -> tuple[float, float, float, float]:
    """Box that holds a generic ring's arc, its line included."""
    v = [float(x) for x in generic_ring_points(a, b, away, members).split()]
    xs, ys = v[0::2], v[1::2]
    return (min(xs) - OVERLAY_LINE, min(ys) - OVERLAY_LINE, max(xs) + OVERLAY_LINE, max(ys) + OVERLAY_LINE)


def sample_curve(points: str, per_segment: int = 8) -> list[tuple[float, float]]:
    """Points along an open (in, anchor, out) curve, as ChemDraw draws it."""
    v = [float(x) for x in points.split()]
    pts = list(zip(v[0::2], v[1::2]))
    triples = [pts[i:i + 3] for i in range(0, len(pts) - 2, 3)]
    out = [triples[0][1]] if triples else []
    for (_, a, h1), (h2, b, _) in zip(triples, triples[1:]):
        for n in range(1, per_segment + 1):
            u = n / per_segment
            out.append(tuple((1 - u) ** 3 * a[i] + 3 * (1 - u) ** 2 * u * h1[i] + 3 * (1 - u) * u * u * h2[i]
                             + u ** 3 * b[i] for i in (0, 1)))
    return out


def _segment_distance(p, a, b) -> float:
    dx, dy = b[0] - a[0], b[1] - a[1]
    length = dx * dx + dy * dy
    t = 0.0 if not length else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


# The least room an arc keeps from the molecule's bonds, away from its two ends.
GENERIC_RING_CLEARANCE = 3.0
_END_ZONE = 6.0   # samples this close to an end atom are not counted: the arc starts there


def arc_clearance(points: str, bonds, ends) -> float:
    """How close the arc comes to any bond, leaving out the stretch next to its two end atoms."""
    samples = [p for p in sample_curve(points) if all(math.dist(p, e) > _END_ZONE for e in ends)]
    if not samples or not bonds:
        return float("inf")
    return min(_segment_distance(p, a, b) for p in samples for a, b in bonds)


def choose_generic_ring(a, b, away=None, members: int = 6, bonds=(), clear_a: float = 0.0,
                        clear_b: float = 0.0) -> tuple[str, float, int]:
    """The generic-ring arc through a and b that keeps clear of the molecule's other bonds.

    Both sides of the bond are tried, the one away from ``away`` first; the side that keeps
    further from every bond (``arc_clearance``) wins. If neither keeps GENERIC_RING_CLEARANCE,
    smaller rings are tried, down to four atoms. Returns (CurvePoints, clearance, members);
    a clearance under GENERIC_RING_CLEARANCE means the arc still meets a bond.
    """
    mx, my = (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0
    if away is None:
        sides = [None]
    else:
        # The same point reflected through the bond's midpoint puts the arc on the other side.
        sides = [away, (2 * mx - away[0], 2 * my - away[1])]
    best = None
    for size in range(int(members), 3, -1):
        for side in sides:
            points = generic_ring_points(a, b, side, size, clear_a, clear_b)
            room = arc_clearance(points, bonds, (a, b))
            if best is None or room > best[1] + 1e-6:
                best = (points, room, size)
        if best[1] >= GENERIC_RING_CLEARANCE:
            return best
    return best
