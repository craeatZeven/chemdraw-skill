"""Resolve and place variable-position R groups on a drawn ring."""
from __future__ import annotations

import math


def ring_r_group_sites(mol, groups):
    """Resolve each numbered R group to candidate atom indices and one perceived ring."""
    if not isinstance(groups, list):
        raise ValueError('ring_r_groups must be a list of {number, atoms} objects')
    used_numbers = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0}
    rings = mol.GetRingInfo().AtomRings()
    sites_by_ring = []
    for group in groups:
        if not isinstance(group, dict) or set(group) != {'number', 'atoms'}:
            raise ValueError('ring_r_groups entries must contain number and atoms')
        number, atoms = group['number'], group['atoms']
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
        sites_by_ring.append((number, atoms, matching[0]))
    return sites_by_ring


def _side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def _intersects(a, b, c, d):
    return (all(max(min(a[k], b[k]), min(c[k], d[k])) <= min(max(a[k], b[k]), max(c[k], d[k]))
                for k in (0, 1))
            and _side(a, b, c) * _side(a, b, d) <= 0
            and _side(c, d, a) * _side(c, d, b) <= 0)


def ring_r_group_layout(mol, points, sites_by_ring, bond_length, font_size):
    """Return the number, candidate sites, anchor, label and crossed RDKit bond index."""
    if not sites_by_ring:
        return []
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    placements = []
    for number, sites, ring in sites_by_ring:
        center = tuple(sum(points[i][axis] for i in ring) / len(ring) for axis in (0, 1))
        candidates = []
        for i, j in zip(ring, ring[1:] + ring[:1]):
            a, b = points[i], points[j]
            midpoint = tuple((x + y) / 2 for x, y in zip(a, b))
            distance = math.dist(midpoint, center)
            if distance < 0.01 or abs(_side(a, b, center)) < 0.01:
                continue
            dx, dy = ((p - c) / distance for c, p in zip(center, midpoint))
            # The superscript is raised to the right of the R glyph. Allow room when
            # that part of the label faces back toward the crossed edge.
            label_clearance = font_size * (max(0, -dx) * (0.75 + 0.4 * len(str(number)))
                                           + max(0, dy) * 0.8)
            extension = max(0.6 * bond_length, label_clearance)
            label = midpoint[0] + dx * extension, midpoint[1] + dy * extension
            crossed = mol.GetBondBetweenAtoms(i, j).GetIdx()
            if any(_intersects(center, label, points[u], points[v])
                   for index, (u, v) in enumerate(bonds) if index != crossed):
                continue
            if any(math.dist(label, other_label) < bond_length
                   or _intersects(midpoint, label, other_center, other_label)
                   for _, _, other_center, other_label, _ in placements):
                continue
            nominal = midpoint[0] + dx * 0.6 * bond_length, midpoint[1] + dy * 0.6 * bond_length
            clearance = min(math.dist(nominal, point) for point in points)
            horizontal = abs(dx)
            score = (round(clearance / bond_length, 2), int(i in sites) + int(j in sites),
                     round(horizontal, 2), -crossed)
            candidates.append((score, label, crossed))
        if not candidates:
            raise ValueError(f'ring_r_groups R{number}: no clear bond from the ring center to its exterior')
        _, label, crossed = max(candidates)
        placements.append((number, sites, center, label, crossed))
    return placements
