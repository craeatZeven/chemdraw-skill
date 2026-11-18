"""Steric clash marks: a short double zigzag across the gap between two crowding groups.

Journals mark a steric clash in a conformer or transition state with a small red double
zigzag between the two groups (BJOC 22-125 scheme 4: between Ar and the ring in V-1, between
the ring and [Rh] in V-2 and V-3). The mark sits in the middle of the clear gap between the
groups, runs across the line joining them and is drawn as two open ChemDraw curves whose
control points lie on their anchors, so each is a polyline of straight segments.
"""

from __future__ import annotations

import math

LENGTH = 9.0        # how far each zigzag runs across the line between the groups
TEETH = 4           # segments per zigzag
AMPLITUDE = 1.0     # how far a zigzag swings either way along that line
SEPARATION = 2.4    # between the two zigzags' centre lines
CLEARANCE = 1.0     # left between the mark and each group's reach
MIN_GAP = SEPARATION + 2 * AMPLITUDE + 2 * CLEARANCE
COLOR = "#FF0000"   # the red of the witness


def clash_mark(a: tuple[float, float], b: tuple[float, float], reach_a: float, reach_b: float
               ) -> tuple[list[list[tuple[float, float]]], float]:
    """The mark between atoms at ``a`` and ``b``, whose drawings reach ``reach_a`` and
    ``reach_b`` from their centres. Returns (two polylines, the clear gap between the reaches);
    a gap under MIN_GAP leaves no room for the mark without touching a group."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    distance = math.hypot(dx, dy) or 1.0
    ux, uy = dx / distance, dy / distance
    nx, ny = -uy, ux
    gap = distance - reach_a - reach_b
    middle = reach_a + gap / 2.0
    mx, my = a[0] + ux * middle, a[1] + uy * middle
    lines = []
    for side in (-1.0, 1.0):
        cx, cy = mx + ux * side * SEPARATION / 2.0, my + uy * side * SEPARATION / 2.0
        line = []
        for j in range(TEETH + 1):
            across = -LENGTH / 2.0 + LENGTH * j / TEETH
            along = AMPLITUDE if j % 2 else -AMPLITUDE
            line.append((cx + nx * across + ux * along, cy + ny * across + uy * along))
        lines.append(line)
    return lines, gap


def curve_points(line: list[tuple[float, float]]) -> str:
    """CurvePoints for a polyline: each point as a (control in, anchor, control out) triple
    whose controls sit on the anchor, so ChemDraw draws straight segments between them."""
    return " ".join(f"{x:.2f} {y:.2f} {x:.2f} {y:.2f} {x:.2f} {y:.2f}" for x, y in line)
