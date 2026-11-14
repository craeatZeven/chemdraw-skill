"""Schematic one-component phase diagram (pressure against temperature)."""
from __future__ import annotations

from typing import Any, Mapping

from .canvas import Canvas
from .energy import _result


def phase_diagram(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Pressure-temperature phase diagram of one substance, schematic.

    ``triple`` and ``critical`` are ``[T, P]`` in any units (only their order matters:
    the critical point lies above the triple point in both). ``melting_slope`` is
    ``"positive"`` (most substances) or ``"negative"`` (water, bismuth). Optional:
    ``labels`` for the three regions, ``marks``: ``[{"T", "P", "label"}]`` (e.g. 1 atm and
    the normal boiling point; a mark must lie inside the drawn range), ``x_label``,
    ``y_label``, ``supercritical`` (label the region beyond the critical point).
    """
    triple, critical = spec.get('triple'), spec.get('critical')
    if not (isinstance(triple, (list, tuple)) and isinstance(critical, (list, tuple))
            and len(triple) == 2 and len(critical) == 2):
        raise ValueError('Give "triple" and "critical" as [T, P]')
    (tt, tp), (ct, cp) = [float(v) for v in triple], [float(v) for v in critical]
    if not (ct > tt and cp > tp):
        raise ValueError('The critical point lies at a higher temperature and pressure than the triple point')
    slope = str(spec.get('melting_slope', 'positive')).lower()
    if slope not in ('positive', 'negative'):
        raise ValueError('melting_slope must be "positive" or "negative"')
    labels = {'solid': 'solid', 'liquid': 'liquid', 'gas': 'gas', **(spec.get('labels') or {})}

    width, height = 200.0, 150.0
    # Schematic, not to scale: the triple and critical points sit at fixed places and other values
    # are placed piecewise-linearly between them, so marks keep their order against both points.
    t_lo, t_hi = tt - (ct - tt) * 0.9, ct + (ct - tt) * 0.6
    p_lo, p_hi = 0.0, cp * 1.6
    t_knots = ((t_lo, 0.0), (tt, 0.36 * width), (ct, 0.76 * width), (t_hi, width))
    p_knots = ((p_lo, height), (tp, 0.70 * height), (cp, 0.24 * height), (p_hi, 0.0))

    def interpolate(value, knots):
        for (a, fa), (b, fb) in zip(knots, knots[1:]):
            if value <= b or (b, fb) == knots[-1]:
                return fa + (fb - fa) * (value - a) / ((b - a) or 1.0)
        return knots[-1][1]

    def at(t, p):
        return interpolate(t, t_knots), interpolate(p, p_knots)

    c = Canvas()
    c.arrow((0.0, height), (0.0, -10.0))
    c.arrow((0.0, height), (width + 10.0, height))
    c.label((-4.0, -16.0), spec.get('y_label', 'Pressure'))
    c.label((width + 10.0, height + 13.0), spec.get('x_label', 'Temperature'), align='Right')
    tx, ty = at(tt, tp)
    cx, cy = at(ct, cp)
    # Sublimation curve: from near the origin, rising ever more steeply to the triple point.
    sub = [(tx * k / 20, height - (height - ty) * (k / 20) ** 2.4) for k in range(21)]
    # Vaporisation curve: triple point to critical point, concave upward.
    vap = [(tx + (cx - tx) * k / 20, ty + (cy - ty) * (k / 20) ** 1.6) for k in range(21)]
    # Melting curve: nearly vertical from the triple point to the top of the frame.
    lean = (0.12 if slope == 'positive' else -0.08) * (ty - 0.0)
    melt = [(tx + lean * k / 20, ty - ty * k / 20) for k in range(21)]
    for curve in (sub, vap, melt):
        c.path(curve)
        c.block_path(curve, 2.0)
    c.block((0.0, -10.0, 0.0, height), 2.0)
    c.block((0.0, height, width + 10.0, height), 2.0)
    c.circle((tx, ty), 1.6, fill='#000000', outline=None)
    c.circle((cx, cy), 1.6, fill='#000000', outline=None)
    warnings = []

    def put(name, text, candidates):
        _box, clear = c.place(candidates, text)
        if not clear:
            warnings.append(f'the {name} label overlaps a curve or another label')

    put('triple point', spec.get('triple_label', 'triple point'),
        [((tx + 6.0, ty + 12.0), 'Left', 'baseline'), ((tx + 6.0, ty + 22.0), 'Left', 'baseline'),
         ((tx - 6.0, ty + 14.0), 'Right', 'baseline')])
    put('critical point', spec.get('critical_label', 'critical point'),
        [((cx + 5.0, cy + 10.0), 'Left', 'baseline'), ((cx + 5.0, cy - 5.0), 'Left', 'baseline'),
         ((cx - 5.0, cy - 6.0), 'Right', 'baseline')])
    put('solid', labels['solid'], [((tx * f, height * g), 'Center', 'middle') for g in (0.35, 0.25, 0.5)
                                    for f in (0.45, 0.35, 0.55)])
    put('liquid', labels['liquid'], [((tx + (cx - tx) * f, ty + (cy - ty) * g), 'Center', 'middle')
                                      for g in (0.55, 0.75, 0.35) for f in (0.3, 0.45, 0.2)])
    put('gas', labels['gas'], [((tx + (width - tx) * f, height - (height - ty) * g), 'Center', 'middle')
                                for g in (0.3, 0.2, 0.45) for f in (0.6, 0.45, 0.75)])
    if spec.get('supercritical'):
        put('supercritical', str(spec['supercritical']),
            [((cx + (width - cx) * f, cy - g), 'Center', 'middle') for g in (24.0, 34.0, 14.0) for f in (0.5, 0.3)])
    marks = []
    for i, mark in enumerate(spec.get('marks') or []):
        t, p = float(mark.get('T', tt)), float(mark.get('P', tp))
        if not (t_lo <= t <= t_hi and p_lo <= p <= p_hi):
            raise ValueError(f'mark {i + 1} lies outside the drawn range')
        x, y = at(t, p)
        if mark.get('P') is not None and mark.get('T') is None:
            c.segment((0.0, y), (width, y), dashed=True)
        elif mark.get('T') is not None and mark.get('P') is None:
            c.segment((x, height), (x, 0.0), dashed=True)
        else:
            c.circle((x, y), 1.3, fill='#000000', outline=None)
        if mark.get('label'):
            put(f'mark {i + 1}', str(mark['label']), [((x + 4.0, y - 4.0), 'Left', 'baseline'),
                                                       ((width - 2.0, y - 4.0), 'Right', 'baseline'),
                                                       ((x + 4.0, y + 11.0), 'Left', 'baseline')])
        marks.append({'T': t, 'P': p})
    title = spec.get('title')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')
    report = {'diagram': 'phase_diagram', 'triple': [tt, tp], 'critical': [ct, cp], 'melting_slope': slope,
              'scale': 'schematic: curve shapes are qualitative', 'marks': marks,
              'warnings': warnings}
    return _result(c, output, report)
