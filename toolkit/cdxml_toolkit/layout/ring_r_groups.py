"""Place variable-position substituents through the edge of a drawn ring."""
from __future__ import annotations

import math


def ring_r_group_sites(mol, groups):
    """Resolve each numbered R group to its allowed atoms and a single perceived ring."""
    if not isinstance(groups, list):
        raise ValueError("ring_r_groups must be a list of {number, atoms} objects")
    used_numbers = {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomicNum() == 0}
    rings = mol.GetRingInfo().AtomRings()
    for group in groups:
        if not isinstance(group, dict) or set(group) != {"number", "atoms"}:
            raise ValueError("ring_r_groups entries must contain number and atoms")
        number, atoms = group["number"], group["atoms"]
        if type(number) is not int or number < 1 or number in used_numbers:
            raise ValueError("ring_r_groups number must be a unique positive integer, distinct from fixed R groups")
        used_numbers.add(number)
        if (not isinstance(atoms, list) or len(atoms) < 2
                or any(type(i) is not int or i < 0 or i >= mol.GetNumAtoms() for i in atoms)
                or len(atoms) != len(set(atoms))):
            raise ValueError("ring_r_groups atoms must be two or more distinct zero-based atom indices")
        matching = [ring for ring in rings if set(atoms).issubset(ring)]
        if len(matching) != 1 or any(not mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms):
            raise ValueError("ring_r_groups atoms must identify eligible atoms on a single ring")
        yield number, atoms, matching[0]


def _side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def _meets(a, b, c, d):
    # The bounding boxes also distinguish disjoint, collinear segments.
    return (all(max(min(a[k], b[k]), min(c[k], d[k])) <= min(max(a[k], b[k]), max(c[k], d[k]))
                for k in (0, 1))
            and _side(a, b, c) * _side(a, b, d) <= 0
            and _side(c, d, a) * _side(c, d, b) <= 0)


def ring_r_group_layout(mol, points, sites_by_ring, bond_length):
    """Return (number, sites, centre, label) in the same coordinate space as points."""
    if not sites_by_ring:
        return []
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    placements = []
    for number, sites, ring in sites_by_ring:
        center = tuple(sum(points[i][k] for i in ring) / len(ring) for k in (0, 1))
        candidates = []
        for i, j in zip(ring, ring[1:] + ring[:1]):
            a, b = points[i], points[j]
            mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            distance = math.dist(mid, center)
            if distance < 0.01 or abs(_side(a, b, center)) < 0.01:
                continue
            reach = 1 + 0.6 * bond_length / distance
            label = tuple(center[k] + reach * (mid[k] - center[k]) for k in (0, 1))
            if any(_meets(center, label, points[u], points[v])
                   for u, v in bonds if (u, v) not in ((i, j), (j, i))):
                continue
            clearance = min(math.dist(label, p) for p in points)
            if any(math.dist(label, other) < bond_length or _meets(mid, label, start, other)
                   for _, _, start, other in placements):
                continue
            # Prefer a free edge next to eligible sites when several directions have room.
            horizontal = abs(mid[0] - center[0]) / distance
            score = (round(clearance / bond_length, 2), int(i in sites) + int(j in sites),
                     round(horizontal, 2), -i, -j)
            candidates.append((score, label))
        if not candidates:
            raise ValueError(f"ring_r_groups R{number}: no clear bond from the ring center to its exterior")
        _, label = max(candidates)
        placements.append((number, sites, center, label))
    return placements
