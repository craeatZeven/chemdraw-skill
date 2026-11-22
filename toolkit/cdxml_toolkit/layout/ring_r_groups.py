"""Resolve and place native ChemDraw variable-attachment R groups on a drawn ring."""
from __future__ import annotations

import math
import re
from typing import NamedTuple
from xml.etree import ElementTree as ET

_HEX_COLOR_RE = re.compile(r'^#[0-9A-Fa-f]{6}$')
_LATIN1_SUPERSCRIPTS = {'1': '\u00b9', '2': '\u00b2', '3': '\u00b3'}


class RingRGroupPlacement(NamedTuple):
    """Resolved geometry for one variable-attachment R group crossing a ring edge."""

    number: int
    sites: list[int]
    anchor: tuple[float, float]
    label: tuple[float, float]
    crossed: int
    alignment: str
    gap_margin: float | None
    color: str | None = None


def ring_r_group_sites(mol, groups):
    """Resolve each numbered R group to candidate atom indices, one perceived ring, and optional color."""
    if not isinstance(groups, list):
        raise ValueError('ring_r_groups must be a list of {number, atoms} objects')
    used_numbers = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0}
    rings = mol.GetRingInfo().AtomRings()
    sites_by_ring = []
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError('ring_r_groups entries must contain number and atoms')
        keys = set(group)
        if not {'number', 'atoms'} <= keys or keys - {'number', 'atoms', 'color'}:
            raise ValueError('ring_r_groups entries must contain number and atoms (and optional color)')
        number, atoms = group['number'], group['atoms']
        color = group.get('color')
        if color is not None and not (isinstance(color, str) and _HEX_COLOR_RE.fullmatch(color)):
            raise ValueError('ring_r_groups color must be a hex colour like #C81E1E')
        if type(number) is not int or number < 1 or number in used_numbers:
            raise ValueError('ring_r_groups numbers must be unique positive integers, distinct from fixed R groups')
        if (not isinstance(atoms, list) or len(atoms) < 2
                or any(type(i) is not int or not 0 <= i < mol.GetNumAtoms() for i in atoms)
                or len(set(atoms)) != len(atoms)):
            raise ValueError('ring_r_groups atoms must be at least two distinct zero-based atom indices')
        matching = [ring for ring in rings if set(atoms).issubset(ring)]
        if len(matching) != 1 or any(mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in atoms):
            raise ValueError('ring_r_groups atoms must identify candidate atoms on one unambiguous ring')
        used_numbers.add(number)
        sites_by_ring.append((number, list(atoms), matching[0], color.upper() if color else None))
    return sites_by_ring


def _side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def _intersects(a, b, c, d):
    return (all(max(min(a[k], b[k]), min(c[k], d[k])) <= min(max(a[k], b[k]), max(c[k], d[k]))
                for k in (0, 1))
            and _side(a, b, c) * _side(a, b, d) <= 0
            and _side(c, d, a) * _side(c, d, b) <= 0)


def ring_r_group_layout(mol, points, sites_by_ring, bond_length, font_size,
                        *, bond_spacing=0.18, margin_width=1.6, line_width=0.6):
    """Draw a short connector normal to one ring edge, with room for native label clipping."""
    if not sites_by_ring:
        return []
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    placements: list[RingRGroupPlacement] = []
    used_edges: set[int] = set()

    for number, sites, ring, color in sites_by_ring:
        center = tuple(sum(points[i][axis] for i in ring) / len(ring) for axis in (0, 1))
        candidates = []
        for i, j in zip(ring, ring[1:] + ring[:1]):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond is None:
                continue
            crossed = bond.GetIdx()
            if crossed in used_edges:
                continue
            a, b = points[i], points[j]
            mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
            edge_length = math.dist(a, b)
            side = _side(a, b, center)
            if edge_length < 0.01 or abs(side) < 0.01:
                continue

            # The centroid picks the outward side. The edge itself defines the normal.
            direction = math.copysign(1 / edge_length, side)
            dx, dy = (b[1] - a[1]) * direction, (a[0] - b[0]) * direction
            is_double = bond.GetBondTypeAsDouble() == 2.0
            inset = (bond_spacing + 0.145 if is_double else 0.20) * edge_length
            inset = min(inset, 0.6 * abs(side) / edge_length)
            alignment = 'Right' if dx < -0.2 else ('Left' if dx > 0.2 else 'Center')

            # ChemDraw trims the terminal bond at the padded label box. Estimate that
            # clipping along the connector, including the different R and superscript bearings.
            horizontal = margin_width + font_size * (0.06 if alignment == 'Right' else 0.28)
            vertical = margin_width + font_size * 0.35
            clip = min(horizontal / max(abs(dx), 1e-6), vertical / max(abs(dy), 1e-6))
            outset = clip + 0.15 * bond_length
            anchor = (mx - dx * inset, my - dy * inset)
            label = (mx + dx * outset, my + dy * outset)

            if any(_intersects(anchor, label, points[u], points[v])
                   for index, (u, v) in enumerate(bonds) if index != crossed):
                continue
            if any(math.dist(label, prev.label) < 0.85 * bond_length
                   or _intersects(anchor, label, prev.anchor, prev.label)
                   for prev in placements):
                continue

            site_hits = int(i in sites) + int(j in sites)
            probe = (mx + dx * 0.55 * bond_length, my + dy * 0.55 * bond_length)
            clearance = min(math.dist(probe, pt) for pt in points)
            prefer_left_first = 1 if (not placements and dx < -0.2) else 0
            score = (
                is_double,
                site_hits,
                round(clearance / max(bond_length, 1e-6), 2),
                round(abs(dx), 2),
                prefer_left_first,
                -crossed,
            )
            # Tighten the oblique gap on the crossed bond only. This keeps the two
            # strokes of a double bond visually paired around one interruption.
            gap_margin = min(margin_width, max(line_width, 0.69 * margin_width)) if abs(dy) > 0.35 else None
            candidates.append((score, anchor, label, crossed, alignment, gap_margin))

        if not candidates:
            raise ValueError(f'ring_r_groups R{number}: no clear edge for variable attachment')
        _, anchor, label, crossed, alignment, gap_margin = max(candidates)
        used_edges.add(crossed)
        placements.append(RingRGroupPlacement(
            number=number,
            sites=sites,
            anchor=anchor,
            label=label,
            crossed=crossed,
            alignment=alignment,
            gap_margin=gap_margin,
            color=color,
        ))
    return placements


def build_variable_r_label_node(
    node_id: str,
    z_id: str,
    text_id: str,
    label_xy: tuple[float, float],
    number: int,
    alignment: str,
    font_id: str,
    font_size: float,
    color_idx: str,
) -> ET.Element:
    """Build a ChemDraw-safe variable R-group label node that preserves ``R^n`` order.

    ChemDraw flips multi-run ``GenericNickname`` labels on left-pointing bonds into
    ``1R``/``2R`` or anchors ``R`` instead of the superscript over the ring edge.
    For ``R¹``–``R³`` (ISO-8859-1 superscripts), emitting a single non-chemical text
    run on ``NodeType="Unspecified"`` with matching ``LabelDisplay`` and
    ``LabelAlignment`` anchors the outer edge cleanly at ``label_xy`` and never
    reverses the label across native ChemDraw round-trips.
    """
    lx, ly = label_xy
    num_str = str(number)
    use_latin1_sup = all(ch in _LATIN1_SUPERSCRIPTS for ch in num_str)

    if use_latin1_sup:
        sup_text = ''.join(_LATIN1_SUPERSCRIPTS[ch] for ch in num_str)
        node = ET.Element('n', {
            'id': node_id,
            'p': f'{lx:.2f} {ly:.2f}',
            'Z': z_id,
            'NodeType': 'Unspecified',
            'LabelDisplay': alignment,
            'color': color_idx,
        })
        t_el = ET.SubElement(node, 't', {
            'id': text_id,
            'p': f'{lx:.2f} {ly + font_size * 0.35:.2f}',
            'LabelAlignment': alignment,
            'LabelJustification': alignment,
            'InterpretChemically': 'no',
        })
        run = ET.SubElement(t_el, 's', {
            'font': font_id,
            'size': f'{font_size:g}',
            'color': color_idx,
        })
        run.text = f'R{sup_text}'
        return node

    node = ET.Element('n', {
        'id': node_id,
        'p': f'{lx:.2f} {ly:.2f}',
        'Z': z_id,
        'NodeType': 'GenericNickname',
        'GenericNickname': f'R{num_str}',
        'LabelDisplay': 'Left',
        'color': color_idx,
    })
    t_el = ET.SubElement(node, 't', {
        'id': text_id,
        'p': f'{lx:.2f} {ly + font_size * 0.35:.2f}',
        'LabelJustification': 'Left',
        'InterpretChemically': 'no',
    })
    r_run = ET.SubElement(t_el, 's', {
        'font': font_id,
        'size': f'{font_size:g}',
        'color': color_idx,
        'face': '96',
    })
    r_run.text = 'R'
    n_run = ET.SubElement(t_el, 's', {
        'font': font_id,
        'size': f'{font_size:g}',
        'color': color_idx,
        'face': '64',
    })
    n_run.text = num_str
    return node
