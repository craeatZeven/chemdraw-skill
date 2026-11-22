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


def ring_r_group_layout(mol, points, sites_by_ring, bond_length, font_size):
    """Place each variable-attachment R group as a short crossing bond near a ring double bond.

    Enforces publication geometry:
    - Anchor sits just inside the crossed ring double bond, never at the ring center
      (``anchor_dist / ring_radius >= 0.45``).
    - Bond length is short relative to ring radius (``0.30 <= bond_len / ring_radius <= 0.55``).
    - Prefers crossing a ring double bond (``Order == 2``) between candidate attachment sites.
    """
    if not sites_by_ring:
        return []
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    placements: list[RingRGroupPlacement] = []
    used_edges: set[int] = set()

    for entry in sites_by_ring:
        if len(entry) == 3:
            number, sites, ring = entry
            color = None
        else:
            number, sites, ring, color = entry

        center = tuple(sum(points[i][axis] for i in ring) / len(ring) for axis in (0, 1))
        ring_radius = sum(math.dist(points[i], center) for i in ring) / len(ring)
        if ring_radius < 0.01:
            raise ValueError(f'ring_r_groups R{number}: degenerate ring coordinates')

        candidates = []
        for i, j in zip(ring, ring[1:] + ring[:1]):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond is None:
                continue
            crossed = bond.GetIdx()
            if crossed in used_edges:
                continue
            a, b = points[i], points[j]
            midpoint = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
            apothem = math.dist(midpoint, center)
            if apothem < 0.01 or abs(_side(a, b, center)) < 0.01:
                continue

            dx = (midpoint[0] - center[0]) / apothem
            dy = (midpoint[1] - center[1]) / apothem
            edge_len = math.dist(a, b)
            apothem_ratio = apothem / ring_radius
            # Inner line of a ChemDraw double bond sits at ~18% BondSpacing inside the edge.
            inner_line_ratio = (apothem - 0.18 * edge_len) / ring_radius

            if abs(dy) <= 0.35:
                if dx > 0.2:
                    # Right-facing labels lead with wide 'R', so push t_out slightly further out
                    # while keeping anchor inside the inner double-bond line and bond_ratio <= 0.54.
                    target_in = min(0.62, inner_line_ratio - 0.06)
                    target_bond = 0.54
                else:
                    target_in = min(0.60, inner_line_ratio - 0.08)
                    target_bond = 0.50
            elif dy < -0.35:
                # Upward-angled labels extend downward toward the ring edge; give extra outward room.
                target_in = min(0.64, inner_line_ratio - 0.045)
                target_bond = 0.54
            else:
                target_in = min(0.63, inner_line_ratio - 0.055)
                target_bond = 0.54

            # Strictly guarantee anchor_dist / ring_radius >= 0.45 and inside the edge,
            # and 0.30 <= bond_len / ring_radius <= 0.55 while crossing cleanly outside the edge.
            t_in = max(0.46, min(target_in, apothem_ratio - 0.12))
            min_cross_bond = (apothem_ratio - t_in) + 0.16
            bond_ratio = max(0.32, min(0.54, max(target_bond, min_cross_bond)))
            if t_in + bond_ratio <= apothem_ratio + 0.08:
                t_in = max(0.46, apothem_ratio - 0.34)
                bond_ratio = min(0.54, max(0.32, (apothem_ratio - t_in) + 0.18))
            t_out = t_in + bond_ratio

            anchor = (center[0] + dx * t_in * ring_radius, center[1] + dy * t_in * ring_radius)
            label = (center[0] + dx * t_out * ring_radius, center[1] + dy * t_out * ring_radius)

            if any(_intersects(anchor, label, points[u], points[v])
                   for index, (u, v) in enumerate(bonds) if index != crossed):
                continue
            if any(math.dist(label, prev.label) < 0.85 * bond_length
                   or _intersects(anchor, label, prev.anchor, prev.label)
                   for prev in placements):
                continue

            is_double = 1 if bond.GetBondTypeAsDouble() == 2.0 else 0
            site_hits = int(i in sites) + int(j in sites)
            probe = (center[0] + dx * (apothem + 0.55 * bond_length),
                     center[1] + dy * (apothem + 0.55 * bond_length))
            clearance = min(math.dist(probe, pt) for pt in points)
            prefer_left_first = 1 if (not placements and dx < -0.2) else 0
            alignment = 'Right' if dx < -0.2 else ('Left' if dx > 0.2 else 'Center')
            score = (
                is_double,
                site_hits,
                round(clearance / max(bond_length, 1e-6), 2),
                round(abs(dx), 2),
                prefer_left_first,
                -crossed,
            )
            candidates.append((score, anchor, label, crossed, alignment))

        if not candidates:
            raise ValueError(f'ring_r_groups R{number}: no clear edge for variable attachment')
        _, anchor, label, crossed, alignment = max(candidates)
        used_edges.add(crossed)
        placements.append(RingRGroupPlacement(
            number=number,
            sites=sites,
            anchor=anchor,
            label=label,
            crossed=crossed,
            alignment=alignment,
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
