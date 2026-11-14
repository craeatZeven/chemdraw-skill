"""Binary (temperature-composition) phase diagrams: eutectic and isomorphous systems.

The diagram is built from the melting points, the eutectic point and the solid solubility
limits. Marked points are located in their phase field and, in a two-phase field, the tie line
is drawn and the lever rule gives the composition and amount of each phase. Curves are
schematic (straight between the given points unless ``curvature`` is set), so the lever-rule
numbers are read from the drawn diagram, not from measured data.
"""
from __future__ import annotations

from typing import Any, Mapping

from .canvas import Canvas, boxes_overlap as _overlap
from .common import number
from .energy import _result

PLOT_W, PLOT_H = 220.0, 160.0


def _curve(x0, t0, x1, t1, bend):
    """Liquidus or solidus between two points, bowed by ``bend`` (0 is straight)."""
    def temperature(x):
        f = (x - x0) / (x1 - x0)
        return t0 + (t1 - t0) * f + bend * (t0 - t1) * f * (1 - f) * (1 if t0 > t1 else -1)
    return temperature


def _solve(fn, lo, hi, target):
    """x in [lo, hi] with fn(x) = target (fn monotonic there)."""
    flo = fn(lo) - target
    for _ in range(80):
        mid = (lo + hi) / 2
        fmid = fn(mid) - target
        if (fmid > 0) == (flo > 0):
            lo, flo = mid, fmid
        else:
            hi = mid
    return (lo + hi) / 2


def binary_phase_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Temperature-composition diagram of components A and B (x is the mole or mass fraction of B).

    ``{"type": "eutectic", "components": ["Pb", "Sn"], "melting": [327, 232],
    "eutectic": [0.62, 183], "solubility": [0.18, 0.97]}`` (``solubility``: the largest
    fraction of B in alpha and of B where beta ends, at the eutectic temperature; omit it for
    immiscible solids). ``{"type": "isomorphous", "components": ["Cu", "Ni"],
    "melting": [1085, 1455]}`` draws a lens with one solid solution. ``points``:
    ``[{"x": 0.4, "T": 250, "label": "P"}]`` are located and analysed with the lever rule.
    ``unit`` (default "°C"), ``basis`` ("mole" or "mass"), ``curvature`` (0 to 0.5).
    """
    kind = str(spec.get('type', 'eutectic')).lower()
    if kind not in ('eutectic', 'isomorphous'):
        raise ValueError('type must be "eutectic" or "isomorphous"')
    names = [str(v) for v in spec.get('components', ['A', 'B'])]
    if len(names) != 2:
        raise ValueError('"components" names A and B')
    melting = spec.get('melting')
    if not isinstance(melting, (list, tuple)) or len(melting) != 2:
        raise ValueError('"melting" is [T_m(A), T_m(B)]')
    ta, tb = (number(v, 'melting point', -1e5, 1e5) for v in melting)
    bend = number(spec.get('curvature', 0.18), 'curvature', 0.0, 0.5)
    unit = str(spec.get('unit', '°C'))
    basis = str(spec.get('basis', 'mole')).lower()
    if basis not in ('mole', 'mass'):
        raise ValueError('basis must be "mole" or "mass"')

    point_temps = [number(item.get('T'), f'point {i + 1} T', -1e5, 1e5) for i, item in enumerate(spec.get('points') or [])]

    def widen(low, high):
        """Temperature range of the plot, stretched so that every marked point is inside."""
        pad = (high - low) * 0.06
        return min([low] + [t - pad for t in point_temps]), max([high] + [t + pad for t in point_temps])
    if kind == 'eutectic':
        eut = spec.get('eutectic')
        if not isinstance(eut, (list, tuple)) or len(eut) != 2:
            raise ValueError('An eutectic diagram needs "eutectic": [x_e, T_e]')
        xe, te = number(eut[0], 'eutectic composition', 0.0, 1.0), number(eut[1], 'eutectic temperature', -1e5, 1e5)
        if not 0 < xe < 1:
            raise ValueError('The eutectic composition lies strictly between pure A and pure B')
        if not te < min(ta, tb):
            raise ValueError('The eutectic temperature lies below both melting points')
        sol = spec.get('solubility')
        xa, xb = (0.0, 1.0) if sol is None else (number(sol[0], 'solubility', 0.0, 1.0), number(sol[1], 'solubility', 0.0, 1.0))
        if not (0 <= xa < xe < xb <= 1):
            raise ValueError('The solubility limits must lie on either side of the eutectic: 0 <= alpha < x_e < beta <= 1')
        liq_a, liq_b = _curve(0.0, ta, xe, te, bend), _curve(1.0, tb, xe, te, bend)
        sol_a = _curve(0.0, ta, xa, te, 0.0) if xa > 0 else None
        sol_b = _curve(1.0, tb, xb, te, 0.0) if xb < 1 else None
        t_bottom = te - (min(ta, tb) - te) * 0.8
        t_bottom = widen(t_bottom, te)[0]
        # Solvus: solubility falls on cooling, to 40 % of its eutectic value at the bottom of the plot.
        solv_a = (lambda t: xa * (0.4 + 0.6 * (t - t_bottom) / (te - t_bottom))) if xa > 0 else None
        solv_b = (lambda t: 1 - (1 - xb) * (0.4 + 0.6 * (t - t_bottom) / (te - t_bottom))) if xb < 1 else None
        alpha = '\u03b1' if xa > 0 else f'{names[0]}(s)'
        beta = '\u03b2' if xb < 1 else f'{names[1]}(s)'

        def edge_a(t):
            """Largest B fraction alpha holds at t (0 for an immiscible solid A)."""
            if xa == 0:
                return 0.0
            return solv_a(t) if t <= te else _solve(sol_a, 0.0, xa, t) if t < ta else 0.0

        def edge_b(t):
            if xb == 1:
                return 1.0
            return solv_b(t) if t <= te else _solve(sol_b, xb, 1.0, t) if t < tb else 1.0

        def liquidus_x(t):
            return (_solve(liq_a, 0.0, xe, t) if t < ta else 0.0), (_solve(liq_b, xe, 1.0, t) if t < tb else 1.0)

        def locate(x, t):
            la, lb = liquidus_x(t) if t > te else (xe, xe)
            if t > te and la <= x <= lb:
                return 'L', [('L', x)]
            if xa > 0 and x <= edge_a(t):
                return alpha, [(alpha, x)]
            if xb < 1 and x >= edge_b(t):
                return beta, [(beta, x)]
            if t > te:
                if x < la:
                    return f'L + {alpha}', [(alpha, edge_a(t)), ('L', la)]
                return f'L + {beta}', [('L', lb), (beta, edge_b(t))]
            return f'{alpha} + {beta}', [(alpha, edge_a(t)), (beta, edge_b(t))]
        t_top = widen(te, max(ta, tb) + (max(ta, tb) - te) * 0.25)[1]
        report_extra = {'eutectic': {'x': xe, 'T': te,
                                     'reaction': f'L → {alpha} + {beta} at {te:g} {unit}'}}
    else:
        if abs(ta - tb) < 1e-9:
            raise ValueError('An isomorphous lens needs two different melting points')
        lens = number(spec.get('lens', 0.5), 'lens', 0.05, 0.9)
        span = abs(ta - tb)
        line = _curve(0.0, ta, 1.0, tb, 0.0)

        def measured(key):
            """Piecewise-linear curve through measured [x, T] points, pinned to the melting points."""
            raw = spec.get(key)
            if raw is None:
                return None
            pts = sorted({(0.0, ta), (1.0, tb), *((number(x, key, 0.0, 1.0), number(t, key, -1e5, 1e5))
                                                  for x, t in raw)})
            temps = [t for _x, t in pts]
            if temps != sorted(temps) and temps != sorted(temps, reverse=True):
                raise ValueError(f'"{key}" must rise or fall steadily from one melting point to the other')

            def fn(x):
                for (x0, t0), (x1, t1) in zip(pts, pts[1:]):
                    if x <= x1:
                        return t0 + (t1 - t0) * (x - x0) / ((x1 - x0) or 1.0)
                return pts[-1][1]
            return fn
        liquidus = measured('liquidus') or (lambda x: line(x) + lens * span * x * (1 - x))
        solidus = measured('solidus') or (lambda x: line(x) - lens * span * x * (1 - x))
        for k in range(1, 50):
            if liquidus(k / 50) <= solidus(k / 50):
                raise ValueError('The liquidus must lie above the solidus between the pure components')
        alpha = str(spec.get('solid_label', '\u03b1'))
        increasing = tb > ta

        def locate(x, t):
            if t >= liquidus(x):
                return 'L', [('L', x)]
            if t <= solidus(x):
                return alpha, [(alpha, x)]
            xl = _solve(liquidus, 0.0, 1.0, t)
            xs = _solve(solidus, 0.0, 1.0, t)
            order = [('L', xl), (alpha, xs)] if increasing else [(alpha, xs), ('L', xl)]
            return f'L + {alpha}', order
        t_bottom, t_top = widen(min(ta, tb) - span * 0.35, max(ta, tb) + span * 0.2)
        report_extra = {}

    def at(x, t):
        return x * PLOT_W, PLOT_H - (t - t_bottom) / (t_top - t_bottom) * PLOT_H

    c = Canvas()
    c.rect(0.0, 0.0, PLOT_W, PLOT_H)

    def draw(fn_of_x, x0, x1, steps=40, dashed=False):
        pts = [at(x0 + (x1 - x0) * k / steps, fn_of_x(x0 + (x1 - x0) * k / steps)) for k in range(steps + 1)]
        c.path(pts, dashed=dashed)
        c.block_path(pts, 2.0)

    def draw_t(fn_of_t, t0, t1, steps=30):
        pts = [at(fn_of_t(t0 + (t1 - t0) * k / steps), t0 + (t1 - t0) * k / steps) for k in range(steps + 1)]
        c.path(pts)
        c.block_path(pts, 2.0)

    if kind == 'eutectic':
        draw(liq_a, 0.0, xe)
        draw(liq_b, xe, 1.0)
        if sol_a:
            draw(sol_a, 0.0, xa)
        if sol_b:
            draw(sol_b, xb, 1.0)
        c.segment(at(xa, te), at(xb, te))
        c.block_path([at(xa, te), at(xb, te)], 2.0)
        if solv_a:
            draw_t(solv_a, t_bottom, te)
        if solv_b:
            draw_t(solv_b, t_bottom, te)
        c.circle(at(xe, te), 1.5, fill='#000000', outline=None)
    else:
        draw(liquidus, 0.0, 1.0)
        draw(solidus, 0.0, 1.0)
    for x in (0.0, 1.0):
        c.block_path([at(x, t_bottom), at(x, t_top)], 1.0)
    c.block_path([at(0.0, t_bottom), at(1.0, t_bottom)], 1.0)
    c.block_path([at(0.0, t_top), at(1.0, t_top)], 1.0)

    # ---- axes
    for k in range(6):
        x = k / 5
        px, py = at(x, t_bottom)
        c.segment((px, py), (px, py - 3.0))
        c.label((px, py + 11.0), f'{x:.1f}'.rstrip('0').rstrip('.') if k not in (0, 5) else ('0' if k == 0 else '1'),
                align='Center', size=8)
    c.label((0.0, PLOT_H + 24.0), names[0], align='Center')
    c.label((PLOT_W, PLOT_H + 24.0), names[1], align='Center')
    c.label((PLOT_W / 2, PLOT_H + 24.0), spec.get('x_label', f'{"x" if basis == "mole" else "w"}{{sub:{names[1]}}}'),
            align='Center')
    c.label((-6.0, -8.0), spec.get('y_label', f'T / {unit}'), align='Left')
    for t, name in ((ta, names[0]), (tb, names[1])):
        px, py = at(0.0 if name == names[0] else 1.0, t)
        side = -5.0 if name == names[0] else 5.0
        c.label((px + side, py), f'{t:g}', align='Right' if side < 0 else 'Left', valign='middle', size=8)
    if kind == 'eutectic':
        px, py = at(1.0, te)
        c.label((px + 5.0, py), f'{te:g}', valign='middle', size=8)
        c.segment((px, py), (px + 3.0, py))

    # ---- field labels: the centre of each field, found by sampling the plot
    warnings: list[str] = []
    samples: dict[str, list] = {}
    for i in range(1, 60):
        for j in range(1, 60):
            x, t = i / 60, t_bottom + (t_top - t_bottom) * j / 60
            samples.setdefault(locate(x, t)[0], []).append((x, t))
    for name, pts in samples.items():
        xs = sorted(p[0] for p in pts)
        ts = sorted(p[1] for p in pts)
        cx, ct = xs[len(xs) // 2], ts[len(ts) // 2]
        candidates = [(at(*p), 'Center', 'middle') for p in sorted(pts, key=lambda p: (p[0] - cx) ** 2 * 4 +
                                                                  ((p[1] - ct) / (t_top - t_bottom)) ** 2)[:90:3]]
        inside = next((cand for cand in candidates
                       if not any(_overlap(c.measure(cand[0], name, None, cand[1], cand[2]), o) for o in c.obstacles)),
                      None)
        if inside is not None:
            c.place([inside], name)
            continue
        # Too narrow for its name: label it beside the frame with a leader into the field.
        fx, fy = at(cx, ct)
        right = cx >= 0.5
        lx = PLOT_W + 14.0 if right else -14.0
        _box, clear = c.place([((lx, fy + dy), 'Left' if right else 'Right', 'middle') for dy in (0.0, -12.0, 12.0, -24.0, 24.0)],
                              name)
        c.segment((fx, fy), (lx - (2.0 if right else -2.0), fy))
        if not clear:
            warnings.append(f'the field label {name} overlaps other ink')

    # ---- marked points and the lever rule
    analysed = []
    for i, item in enumerate(spec.get('points') or []):
        x = number(item.get('x'), f'point {i + 1} x', 0.0, 1.0)
        t = number(item.get('T'), f'point {i + 1} T', t_bottom, t_top)
        field, phases = locate(x, t)
        entry: dict[str, Any] = {'label': item.get('label'), 'x': x, 'T': t, 'field': field}
        px, py = at(x, t)
        c.circle((px, py), 1.4, fill='#000000', outline=None)
        if len(phases) == 2:
            (p1, x1), (p2, x2) = phases
            c.segment(at(x1, t), at(x2, t), dashed=True)
            f2 = (x - x1) / (x2 - x1) if abs(x2 - x1) > 1e-12 else 0.0
            entry['phases'] = [{'phase': p1, 'composition': round(x1, 3), 'fraction': round(1 - f2, 3)},
                               {'phase': p2, 'composition': round(x2, 3), 'fraction': round(f2, 3)}]
        else:
            entry['phases'] = [{'phase': phases[0][0], 'composition': round(x, 3), 'fraction': 1.0}]
        if item.get('label'):
            c.place([((px + 4.0, py - 4.0), 'Left', 'baseline'), ((px - 4.0, py - 4.0), 'Right', 'baseline'),
                     ((px + 4.0, py + 10.0), 'Left', 'baseline')], str(item['label']))
        analysed.append(entry)
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report = {'diagram': 'binary_phase_diagram', 'type': kind, 'components': names, 'melting': [ta, tb],
              'fields': sorted(samples), 'points': analysed, **report_extra,
              'scale': 'schematic curves: lever-rule values are read from the drawn diagram', 'warnings': warnings}
    return _result(c, output, report)
