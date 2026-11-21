"""Editable ellipse lobes, hybrid-orbital panels and sigma/pi overlap drawings.

Coordinates are points, x right, y down. A lobe's ``size`` and ``width`` are its full
major and minor axis lengths; ``coefficient`` scales both. Phase is encoded by colour,
not by chemical connectivity: these panels contain no molecular fragments.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

from .canvas import Canvas, resolve_color
from .common import number, point
from .energy import _result

POSITIVE_PHASE = '#C0392B'
NEGATIVE_PHASE = '#0072B2'
_TETRAHEDRON = (
    (1.0, -1.0, 1.0),
    (-1.0, -1.0, -1.0),
    (-1.0, 1.0, 1.0),
    (1.0, 1.0, -1.0),
)


def _unequal_directions():
    h_h = math.radians(110.3)
    x = -math.sqrt((2 * math.cos(h_h) + 1) / 3)
    radius = math.sqrt(1 - x * x)
    return ((1.0, 0.0, 0.0),) + tuple(
        (x, radius * math.cos(math.radians(phi)), radius * math.sin(math.radians(phi)))
        for phi in (90, 210, 330)
    )


def _vector_angle(first, second, dimensions=3):
    a, b = first[:dimensions], second[:dimensions]
    norm_a = math.sqrt(sum(value * value for value in a))
    norm_b = math.sqrt(sum(value * value for value in b))
    cosine = sum(x * y for x, y in zip(a, b)) / (norm_a * norm_b)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _phase_color(phase, color):
    if phase not in (-1, 1):
        raise ValueError('phase must be 1 or -1')
    return resolve_color(color, POSITIVE_PHASE if phase == 1 else NEGATIVE_PHASE)


def _ellipse(canvas, origin, direction, size, width, color, oval_type, parent=None):
    norm = math.sqrt(sum(value * value for value in direction))
    if norm < 1e-9:
        raise ValueError('A lobe needs a nonzero direction')
    unit = tuple(value / norm for value in direction)
    center = tuple(origin[i] + unit[i] * size / 2 for i in range(3))
    major = tuple(center[i] + unit[i] * size / 2 for i in range(3))
    planar = math.hypot(unit[0], unit[1])
    perpendicular = ((-unit[1] / planar, unit[0] / planar, 0.0) if planar > 1e-9 else (1.0, 0.0, 0.0))
    minor = tuple(center[i] + perpendicular[i] * width / 2 for i in range(3))
    a, b = size / 2, width / 2
    ex = math.sqrt((a * unit[0]) ** 2 + (b * perpendicular[0]) ** 2)
    ey = math.sqrt((a * unit[1]) ** 2 + (b * perpendicular[1]) ** 2)
    box = (center[0] - ex, center[1] - ey, center[0] + ex, center[1] + ey)
    attrs = {
        'GraphicType': 'Oval',
        'BoundingBox': ' '.join(f'{value:.2f}' for value in box),
        'Center3D': ' '.join(f'{value:.2f}' for value in center),
        'MajorAxisEnd3D': ' '.join(f'{value:.2f}' for value in major),
        'MinorAxisEnd3D': ' '.join(f'{value:.2f}' for value in minor),
    }
    if oval_type is not None:
        attrs['OvalType'] = oval_type
    index = canvas._resolve_color(color)
    if index is not None:
        attrs['color'] = index
    node = canvas.add('graphic', attrs, parent)
    canvas.see(box[0], box[1]); canvas.see(box[2], box[3])
    return node


def _lobe_vector(canvas, at, direction, size, width, phase, fill, color, parent=None):
    origin3 = (at[0], at[1], 0.0)
    oval_type = {'filled': 'Shaded', 'hollow': None, 'solid': 'Filled'}.get(fill, False)
    if oval_type is False:
        raise ValueError('fill must be "filled" (phase shading), "hollow" (outline), or "solid"')
    return _ellipse(canvas, origin3, direction, size, width, _phase_color(phase, color), oval_type, parent)


def _dimensions(size, width, coefficient):
    scale = abs(number(coefficient, 'coefficient', 0.0, 1000.0))
    return (number(size, 'size', 0.01, 1000.0) * scale,
            number(width, 'width', 0.01, 1000.0) * scale)


def lobe(canvas, at, angle_deg=-90.0, size=36.0, width=16.0, coefficient=1.0,
         phase=1, fill='filled', color=None, parent=None):
    """Place one ellipse lobe, with its narrow end at ``at`` and axis at ``angle_deg``."""
    origin = point(at)
    angle = number(angle_deg, 'angle_deg', -36000, 36000)
    major, minor = _dimensions(size, width, coefficient)
    radians = math.radians(angle)
    return _lobe_vector(canvas, origin, (math.cos(radians), math.sin(radians), 0.0),
                        major, minor, phase, fill, color, parent)


def p_orbital(canvas, at, angle_deg=-90.0, size=36.0, width=16.0, coefficient=1.0,
              phase=1, fill='filled', color=None, parent=None, negative_color=None, axis='xy'):
    """Two opposite lobes with opposite phases; ``axis='z'`` is normal to the page."""
    if axis not in ('xy', 'z'):
        raise ValueError('p orbital axis must be "xy" or "z"')
    origin = point(at)
    major, minor = _dimensions(size, width, coefficient)
    if axis == 'z':
        return [
            _lobe_vector(canvas, origin, (0.0, 0.0, 1.0), major, minor, phase, fill, color, parent),
            _lobe_vector(canvas, origin, (0.0, 0.0, -1.0), major, minor, -phase, fill, negative_color, parent),
        ]
    angle = number(angle_deg, 'angle_deg', -36000, 36000)
    radians = math.radians(angle)
    return [
        _lobe_vector(canvas, origin, (math.cos(radians), math.sin(radians), 0.0),
                     major, minor, phase, fill, color, parent),
        _lobe_vector(canvas, origin, (-math.cos(radians), -math.sin(radians), 0.0),
                     major, minor, -phase, fill, negative_color, parent),
    ]


def _angle_arc(canvas, center, vectors, radius, label, label_level=0):
    first, second = vectors
    a0 = math.atan2(first[1], first[0])
    a1 = math.atan2(second[1], second[0])
    start = (center[0] + radius * math.cos(a0), center[1] + radius * math.sin(a0))
    end = (center[0] + radius * math.cos(a1), center[1] + radius * math.sin(a1))
    sweep = math.degrees((a0 - a1 + math.pi) % (2 * math.pi) - math.pi)
    canvas.add('arrow', {
        'Tail3D': f'{start[0]:.2f} {start[1]:.2f} 0',
        'Head3D': f'{end[0]:.2f} {end[1]:.2f} 0',
        'Center3D': f'{center[0]:.2f} {center[1]:.2f} 0',
        'AngularSize': f'{sweep:.2f}',
        'ArrowheadHead': 'Full', 'ArrowheadTail': 'None', 'ArrowheadType': 'Solid',
        'HeadSize': '400', 'ArrowheadCenterSize': '350', 'ArrowheadWidth': '100', 'FillType': 'None',
    })
    canvas.see(*start); canvas.see(*end)
    canvas.label((center[0], center[1] - radius * 2.5 - 20 - 14 * label_level),
                 label, align='Center', valign='middle')


def _shape_spec(spec):
    kind = str(spec.get('hybridisation', 'sp3')).lower().replace('³', '3').replace('²', '2')
    if kind not in ('sp3', 'sp2', 'sp', 'unequal'):
        raise ValueError('hybridisation must be "sp3", "sp2", "sp", or "unequal"')
    origin = point(spec.get('at', (80, 80)))
    size = number(spec.get('size', 36), 'size', 0.01, 1000)
    width = number(spec.get('width', 16), 'width', 0.01, 1000)
    fill = str(spec.get('fill', 'filled')).lower()
    positive = spec.get('positive_color')
    negative = spec.get('negative_color')
    return kind, origin, size, width, fill, positive, negative


def draw_orbital_shapes(canvas, spec: Mapping[str, Any]) -> dict[str, Any]:
    """Draw one hybridisation panel on an existing canvas."""
    kind, origin, size, width, fill, positive, negative = _shape_spec(spec)
    if kind in ('sp3', 'unequal'):
        directions = _TETRAHEDRON if kind == 'sp3' else _unequal_directions()
        angle = math.degrees(math.acos(-1 / 3)) if kind == 'sp3' else _vector_angle(directions[0], directions[1])
        default_title = 'unequal hybridisation (CH{sub:3}Cl)' if kind == 'unequal' else 'sp{sup:3}'
    elif kind == 'sp2':
        directions = tuple((math.cos(math.radians(a)), math.sin(math.radians(a)), 0.0) for a in (0, 120, 240))
        angle, default_title = 120.0, 'sp{sup:2}'
    else:
        directions = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
        angle, default_title = 180.0, 'sp'

    if kind == 'unequal':
        sizes = spec.get('lobe_sizes', [size * 1.3, size, size, size])
    else:
        sizes = spec.get('lobe_sizes', [size] * len(directions))
    if not isinstance(sizes, (list, tuple)) or len(sizes) != len(directions):
        raise ValueError(f'{kind} needs {len(directions)} lobe_sizes')
    sizes = [number(value, 'lobe size', 0.01, 1000) for value in sizes]
    show_axes = bool(spec.get('show_axes', True))
    for direction, lobe_size in zip(directions, sizes):
        if show_axes:
            canvas.segment(origin, (origin[0] + direction[0] * lobe_size,
                                    origin[1] + direction[1] * lobe_size))
        _lobe_vector(canvas, origin, direction, lobe_size, width, 1, fill, positive)
    p_count = 0
    if kind == 'sp2':
        p_orbital(canvas, origin, size=size * 0.8, width=width * 0.8, phase=1, fill=fill,
                  color=positive, negative_color=negative, axis='z')
        p_count = 2
    elif kind == 'sp':
        p_orbital(canvas, origin, -90, size * 0.75, width * 0.75, phase=1, fill=fill,
                  color=positive, negative_color=negative)
        p_orbital(canvas, origin, size=size * 0.65, width=width * 0.65, phase=1, fill=fill,
                  color=positive, negative_color=negative, axis='z')
        p_count = 4
    if spec.get('show_angle', True):
        if kind == 'unequal':
            arc_specs = (((directions[0], directions[1]), angle, 'H-C-Cl 3D'),
                         ((directions[1], directions[2]), 110.3, 'H-C-H 3D'))
        else:
            arc_specs = ((directions[:2], angle, '3D angle' if kind == 'sp3' else ''),)
        for label_level, (vectors, spatial, prefix) in enumerate(arc_specs):
            projected = _vector_angle(*vectors, dimensions=2)
            label = (f'{prefix + " " if prefix else ""}{spatial:.1f}°' if kind != 'sp3'
                     else f'3D angle {spatial:.1f}° (projection {projected:.1f}°)')
            if kind == 'unequal':
                label += f' (projection {projected:.1f}°)'
            _angle_arc(canvas, origin, vectors, min(sizes) * 0.42, label, label_level)
    title = spec.get('title', default_title)
    if title:
        canvas.label((origin[0], origin[1] + size + 30), str(title), align='Center')
    report = {'diagram': 'orbital_shapes', 'hybridisation': kind,
              'hybrid_lobes': len(directions), 'p_lobes': p_count,
              'lobe_count': len(directions) + p_count, 'angle_deg': round(angle, 4),
              'lobe_sizes': sizes, 'warnings': []}
    if kind == 'unequal':
        report['angles_deg'] = {'H-C-Cl': round(angle, 4), 'H-C-H': 110.3}
    return report


def orbital_shapes(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Draw sp3, sp2, sp or unequal hybrid lobes with axes and angle annotations."""
    canvas = Canvas()
    report = draw_orbital_shapes(canvas, spec)
    return _result(canvas, output, report)


def _overlap_options(spec):
    size = number(spec.get('size', 64), 'size', 0.01, 1000)
    width = number(spec.get('width', 22), 'width', 0.01, 1000)
    separation = number(spec.get('separation', 16), 'separation', 0.01, 1000)
    fill = str(spec.get('fill', 'filled')).lower()
    return size, width, separation, fill


def _inner_radii(x_full, y_full, center_distance):
    ratio = center_distance / x_full
    x_radius = (x_full - center_distance) / 2
    y_radius = y_full / 2 * math.sqrt(1 - ratio * ratio)
    scale = 0.98 * (-ratio + math.sqrt(ratio * ratio + 4 * ratio)) / 2
    return scale * x_radius, scale * y_radius


def _draw_sigma(canvas, origin, spec):
    size, width, separation, fill = _overlap_options(spec)
    if separation >= size:
        raise ValueError('sigma overlap needs separation smaller than size so the lobe boxes intersect')
    left = (origin[0] - separation / 2, origin[1])
    right = (origin[0] + separation / 2, origin[1])
    positive = spec.get('positive_color')
    canvas.segment((left[0] - 12, origin[1]), (right[0] + 12, origin[1]), dashed=True)
    overlap_major, overlap_minor = _inner_radii(size, width, size - separation)
    lobe(canvas, (origin[0] - overlap_major, origin[1]), 0, 2 * overlap_major, 2 * overlap_minor,
         phase=1, fill='solid', color=positive)
    lobe(canvas, left, 0, size, width, phase=1, fill=fill, color=positive)
    lobe(canvas, right, 180, size, width, phase=1, fill=fill, color=positive)
    canvas.label((origin[0], origin[1] + 48), 'σ overlap', align='Center')
    canvas.label((origin[0], origin[1] - 42), 'internuclear axis', align='Center')
    return {'kind': 'sigma', 'lobe_count': 2, 'overlap_regions': 1,
            'nuclei': [left, right]}


def _draw_pi(canvas, origin, spec):
    size, width, separation, fill = _overlap_options(spec)
    if separation >= width:
        raise ValueError('pi overlap needs separation smaller than width so same-phase lobe boxes intersect')
    left = (origin[0] - separation / 2, origin[1])
    right = (origin[0] + separation / 2, origin[1])
    positive, negative = spec.get('positive_color'), spec.get('negative_color')
    canvas.segment((left[0] - 12, origin[1]), (right[0] + 12, origin[1]), dashed=True)
    overlap_width, overlap_height = _inner_radii(width, size, separation)
    lobe(canvas, (origin[0], origin[1] - size / 2 + overlap_height), -90, 2 * overlap_height, 2 * overlap_width,
         phase=1, fill='solid', color=positive)
    lobe(canvas, (origin[0], origin[1] + size / 2 - overlap_height), 90, 2 * overlap_height, 2 * overlap_width,
         phase=-1, fill='solid', color=negative)
    for nucleus in (left, right):
        p_orbital(canvas, nucleus, -90, size, width, phase=1, fill=fill,
                  color=positive, negative_color=negative)
    canvas.label((origin[0], origin[1] + size + 30), 'π overlap', align='Center')
    canvas.label((origin[0], origin[1] - size - 12), 'internuclear axis', align='Center')
    return {'kind': 'pi', 'lobe_count': 4, 'overlap_regions': 2,
            'nuclei': [left, right], 'upper_phase': 1, 'lower_phase': -1}


def draw_orbital_overlap(canvas, spec: Mapping[str, Any]) -> dict[str, Any]:
    """Draw sigma, pi, or both overlap panels on an existing canvas."""
    kind = str(spec.get('kind', 'both')).lower()
    if kind not in ('sigma', 'pi', 'both'):
        raise ValueError('orbital overlap kind must be "sigma", "pi", or "both"')
    origin = point(spec.get('at', (90, 80)))
    panels = []
    if kind in ('sigma', 'both'):
        panels.append(_draw_sigma(canvas, origin, spec))
    if kind in ('pi', 'both'):
        pi_origin = (origin[0] + (180 if kind == 'both' else 0), origin[1])
        panels.append(_draw_pi(canvas, pi_origin, spec))
    return {'diagram': 'orbital_overlap', 'kind': kind, 'panels': panels, 'warnings': []}


def orbital_overlap(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Draw head-on sigma and side-by-side pi overlap with phase-coloured lobes."""
    canvas = Canvas()
    report = draw_orbital_overlap(canvas, spec)
    return _result(canvas, output, report)
