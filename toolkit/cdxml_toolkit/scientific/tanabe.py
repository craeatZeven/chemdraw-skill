"""Tanabe-Sugano diagrams computed from the full d^n ligand-field Hamiltonian.

For d^n (n = 2 to 8) in an octahedral field the Hamiltonian is built over every Slater
determinant of the ten d spin-orbitals: electron repulsion from the Racah parameters B and C
(Slater-Condon rules with Gaunt coefficients) plus the one-electron octahedral potential
(t2g at -4 Dq, eg at +6 Dq). It is diagonalised along Delta/B; each level is labelled by its
spin (from <S^2>), and its O_h symmetry from its degeneracy and the character of a C4 rotation,
so A1/A2 and T1/T2 are told apart by symmetry, not by guesswork. Free-ion terms come from
<L^2> at Delta = 0. Energies are plotted relative to the ground state, as Tanabe and Sugano did.
"""
from __future__ import annotations

import math
from functools import lru_cache
from itertools import combinations
from typing import Any, Mapping

from .axes import Plot
from .canvas import BLUE, Canvas
from .common import number
from .energy import _result

# C/B ratios Tanabe and Sugano used for each configuration.
DEFAULT_C_OVER_B = {2: 4.42, 3: 4.50, 4: 4.61, 5: 4.48, 6: 4.80, 7: 4.63, 8: 4.71}
L_LETTERS = 'SPDFGHIKLMN'
_LABELS = {(1, 1): 'A{sub:1g}', (1, -1): 'A{sub:2g}', (2, 0): 'E{sub:g}', (3, 1): 'T{sub:1g}', (3, -1): 'T{sub:2g}'}
DEGENERACY_TOL = 1e-7


# ---------------------------------------------------------------------------- angular algebra
def _fact(n):
    return math.factorial(int(round(n)))


def wigner_3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (integer arguments), Racah formula."""
    if m1 + m2 + m3 != 0 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    pre = (-1) ** (j1 - j2 - m3) * math.sqrt(
        _fact(j1 + j2 - j3) * _fact(j1 - j2 + j3) * _fact(-j1 + j2 + j3) / _fact(j1 + j2 + j3 + 1)
        * _fact(j1 + m1) * _fact(j1 - m1) * _fact(j2 + m2) * _fact(j2 - m2) * _fact(j3 + m3) * _fact(j3 - m3))
    total = 0.0
    for t in range(0, j1 + j2 + j3 + 1):
        terms = (t, j3 - j2 + t + m1, j3 - j1 + t - m2, j1 + j2 - j3 - t, j1 - t - m1, j2 - t + m2)
        if min(terms) < 0:
            continue
        total += (-1) ** t / math.prod(_fact(v) for v in terms)
    return pre * total


@lru_cache(maxsize=None)
def gaunt(k, m1, m2, l=2):
    """Condon-Shortley c^k(l m1, l m2)."""
    return (-1) ** m1 * (2 * l + 1) * wigner_3j(l, k, l, 0, 0, 0) * wigner_3j(l, k, l, -m1, m1 - m2, m2)


@lru_cache(maxsize=None)
def coulomb(m1, m2, m3, m4, f2, f4):
    """<m1 m2 | 1/r12 | m3 m4> for d electrons, F^0 set to zero (it shifts every level alike)."""
    if m1 + m2 != m3 + m4:
        return 0.0
    return sum(gaunt(k, m1, m3) * gaunt(k, m4, m2) * fk for k, fk in ((2, f2), (4, f4)))


# ---------------------------------------------------------------------------- determinants
def _orbital(index):
    """Spin-orbital index 0..9 -> (m_l, spin) with spin 0 = alpha, 1 = beta."""
    return index // 2 - 2, index % 2


def _sign_below(mask, p):
    return -1 if bin(mask & ((1 << p) - 1)).count('1') % 2 else 1


def _annihilate(mask, p):
    if not mask >> p & 1:
        return None, 0
    return mask & ~(1 << p), _sign_below(mask, p)


def _create(mask, p):
    if mask >> p & 1:
        return None, 0
    return mask | (1 << p), _sign_below(mask, p)


def _block(n, n_alpha):
    alpha = [i for i in range(10) if i % 2 == 0]
    beta = [i for i in range(10) if i % 2 == 1]
    dets = []
    for a in combinations(alpha, n_alpha):
        for b in combinations(beta, n - n_alpha):
            dets.append(sum(1 << i for i in a + b))
    return sorted(dets)


# Octahedral potential in the complex d basis, units of Dq.
_V_OCT = {(0, 0): 6.0, (1, 1): -4.0, (-1, -1): -4.0, (2, 2): 1.0, (-2, -2): 1.0, (2, -2): 5.0, (-2, 2): 5.0}


def hamiltonian(n, c_over_b, n_alpha):
    """(H_ee in units of B, V in units of Dq, S+ to the next block, L+, M_L list, M_S) for one M_S block."""
    import numpy as np
    f4 = c_over_b / 35.0 * 441.0          # F^4 = 441 F4, F4 = C/35 (B = 1)
    f2 = (1.0 + 5.0 * c_over_b / 35.0) * 49.0   # F^2 = 49 F2, F2 = B + 5 F4
    dets = _block(n, n_alpha)
    where = {d: i for i, d in enumerate(dets)}
    size = len(dets)
    h = np.zeros((size, size))
    v = np.zeros((size, size))
    ll = np.zeros((size, size))
    for col, det in enumerate(dets):
        occupied = [i for i in range(10) if det >> i & 1]
        # one-electron octahedral term  sum V_pq a+_p a_q
        for q in occupied:
            mq, sq = _orbital(q)
            for (ma, mb), value in _V_OCT.items():
                if mb != mq:
                    continue
                p = (ma + 2) * 2 + sq
                m1, s1 = _annihilate(det, q)
                m2, s2 = _create(m1, p)
                if m2 is not None:
                    v[where[m2], col] += value * s1 * s2
        # two-electron term  1/2 sum <pq|rs> a+_p a+_q a_s a_r
        for r in occupied:
            for s in occupied:
                if r == s:
                    continue
                mr, sr = _orbital(r)
                ms, ss = _orbital(s)
                m1, s1 = _annihilate(det, r)
                m2, s2 = _annihilate(m1, s)
                for p in range(10):
                    mp, sp = _orbital(p)
                    if sp != sr:
                        continue
                    mq = mr + ms - mp
                    if abs(mq) > 2:
                        continue
                    q = (mq + 2) * 2 + ss
                    value = coulomb(mp, mq, mr, ms, f2, f4)
                    if not value:
                        continue
                    m3, s3 = _create(m2, q)
                    if m3 is None:
                        continue
                    m4, s4 = _create(m3, p)
                    if m4 is None:
                        continue
                    h[where[m4], col] += 0.5 * value * s1 * s2 * s3 * s4
        # L+ = sum sqrt(6 - m(m+1)) a+_{m+1} a_m
        for q in occupied:
            mq, sq = _orbital(q)
            if mq == 2:
                continue
            p = (mq + 3) * 2 + sq
            m1, s1 = _annihilate(det, q)
            m2, s2 = _create(m1, p)
            if m2 is not None:
                ll[where[m2], col] += math.sqrt(6 - mq * (mq + 1)) * s1 * s2
    # S+ = sum_m a+_{m alpha} a_{m beta}, into the block with one more alpha electron
    upper = _block(n, n_alpha + 1) if n_alpha < min(n, 5) else []
    where_up = {d: i for i, d in enumerate(upper)}
    splus = np.zeros((len(upper), size))
    for col, det in enumerate(dets):
        for m in range(-2, 3):
            beta, alpha = (m + 2) * 2 + 1, (m + 2) * 2
            m1, s1 = _annihilate(det, beta)
            if m1 is None:
                continue
            m2, s2 = _create(m1, alpha)
            if m2 is not None:
                splus[where_up[m2], col] += s1 * s2
    ml = [sum(_orbital(i)[0] for i in range(10) if d >> i & 1) for d in dets]
    msz = (2 * n_alpha - n) / 2
    return h, v, splus, ll, ml, msz


@lru_cache(maxsize=None)
def _matrices(n, c_over_b):
    n_alpha = (n + 1) // 2           # the M_S = 0 or 1/2 block holds every spin multiplet once
    import numpy as np
    h, v, splus, ll, ml, msz = hamiltonian(n, c_over_b, n_alpha)
    s2 = splus.T @ splus + msz * (msz + 1) * np.eye(len(ml))
    l2 = ll.T @ ll + np.diag([m * m + m for m in ml])
    phase = np.array([math.cos(m * math.pi / 2) for m in ml])      # C4 about z acts as exp(-i M_L pi/2)
    return h, v, s2, l2, phase


def levels(n, delta_over_b, c_over_b):
    """Labelled levels at one Delta/B: ``[(energy/B, 2S+1, label)]``, lowest first; ``None`` labels
    mark a sample where two different terms coincide (an accidental degeneracy)."""
    import numpy as np
    h, v, s2, l2, phase = _matrices(n, c_over_b)
    dq = delta_over_b / 10.0
    w, u = np.linalg.eigh(h + dq * v)
    out = []
    k = 0
    while k < len(w):
        j = k
        while j + 1 < len(w) and w[j + 1] - w[k] < DEGENERACY_TOL * max(1.0, abs(w[k])) + 1e-9:
            j += 1
        group = u[:, k:j + 1]
        spin = float(np.mean(np.einsum('ij,ik,kj->j', group, s2, group)))
        mult = int(round(2 * (math.sqrt(spin + 0.25) - 0.5) + 1))
        if delta_over_b == 0:
            lsq = float(np.mean(np.einsum('ij,ik,kj->j', group, l2, group)))
            big_l = int(round(math.sqrt(lsq + 0.25) - 0.5))
            label = L_LETTERS[big_l] if j - k + 1 == 2 * big_l + 1 else None
        else:
            chi = float(np.sum(group ** 2 * phase[:, None]))
            key = (j - k + 1, int(round(chi)))
            label = _LABELS.get(key) if abs(chi - round(chi)) < 0.05 else None
        out.append((float(w[k]), mult, label))
        k = j + 1
    return out


def tanabe_sugano(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Tanabe-Sugano diagram of d^n in an octahedral field.

    ``{"d_electrons": 3}`` (2 to 8); ``c_over_b`` (default: Tanabe and Sugano's value for each
    n); ``delta_max`` (Delta/B, default 40) and ``energy_max`` (E/B, default 70); ``bands``:
    ``[nu1, nu2]`` (cm^-1, the first two spin-allowed absorptions) to find Delta/B, B and Delta;
    ``at``: a Delta/B whose transitions are listed. Spin-allowed terms (the weak-field ground
    state's multiplicity) are drawn in blue.
    """
    n = int(spec.get('d_electrons', spec.get('n', 0)))
    if not 2 <= n <= 8:
        raise ValueError('Tanabe-Sugano diagrams cover d2 to d8 (d1 and d9 have one excited level)')
    c_over_b = number(spec.get('c_over_b', DEFAULT_C_OVER_B[n]), 'c_over_b', 2.0, 8.0)
    delta_max = number(spec.get('delta_max', 40.0), 'delta_max', 5.0, 80.0)
    energy_max = number(spec.get('energy_max', 70.0), 'energy_max', 10.0, 200.0)
    steps = int(spec.get('samples', 160))
    grid = [1e-4] + [delta_max * k / steps for k in range(1, steps + 1)]

    free = levels(n, 0.0, c_over_b)
    ground_free = free[0][0]
    free_terms = [{'term': f'{m}{label}', 'multiplicity': m, 'L': label, 'energy': round(e - ground_free, 4)}
                  for e, m, label in free]

    samples = []
    for d in grid:
        lv = levels(n, d, c_over_b)
        if any(label is None for _e, _m, label in lv):
            continue                                   # an accidental degeneracy: skip this sample
        e0 = lv[0][0]
        samples.append((d, [(e - e0, m, label) for e, m, label in lv], lv[0]))
    if len(samples) < 10:
        raise ValueError('Too few samples could be labelled; change "samples" or "c_over_b"')

    # Curves: the k-th level of each (multiplicity, label), in energy order (same symmetry never crosses).
    curves: dict[tuple, list] = {}
    for d, lv, _ground in samples:
        seen: dict[tuple, int] = {}
        for e, m, label in lv:
            key = (m, label)
            index = seen.get(key, 0)
            seen[key] = index + 1
            curves.setdefault((m, label, index), []).append((d, e))
    weak_ground = samples[0][2]
    allowed_mult = weak_ground[1]

    # Ground-state changes (spin crossover) located by bisection between samples.
    crossovers = []
    for (d1, _l1, g1), (d2, _l2, g2) in zip(samples, samples[1:]):
        if (g1[1], g1[2]) != (g2[1], g2[2]):
            lo, hi = d1, d2
            for _ in range(40):
                mid = (lo + hi) / 2
                g = levels(n, mid, c_over_b)[0]
                if (g[1], g[2]) == (g1[1], g1[2]):
                    lo = mid
                else:
                    hi = mid
            crossovers.append({'delta_over_b': round((lo + hi) / 2, 3), 'from': f'{g1[1]}{g1[2]}', 'to': f'{g2[1]}{g2[2]}'})

    c = Canvas()
    plot = Plot(c, (0.0, delta_max), (0.0, energy_max), float(spec.get('width', 230.0)), float(spec.get('height', 230.0)))
    plot.frame('Δ/B', 'E/B')
    ends = []
    for (m, label, index), pts in sorted(curves.items(), key=lambda kv: (kv[0][0] != allowed_mult, kv[0])):
        runs = plot.polyline(pts, color=BLUE if m == allowed_mult else None)
        if runs and pts[-1][1] <= energy_max:
            ends.append((pts[-1][1], f'{{sup:{m}}}{label}', m == allowed_mult))
    warnings = []
    unlabelled = 0
    for e, text, allowed in sorted(ends):
        px, py = plot.at(delta_max, e)
        _box, clear = c.place([((px + 3.0, py + dy), 'Left', 'middle') for dy in (0, -6, 6, -12, 12, -18, 18)],
                              text, size=6, color=BLUE if allowed else None)
        if not clear:
            unlabelled += 1
    shown = set()
    for term in free_terms:
        if term['energy'] > energy_max or term['term'] in shown:
            continue
        shown.add(term['term'])
        _px, py = plot.at(0.0, term['energy'])
        c.place([((-18.0, py + dy), 'Right', 'middle') for dy in (0, -6, 6, -12, 12)],
                f'{{sup:{term["multiplicity"]}}}{term["L"]}', size=6,
                color=BLUE if term['multiplicity'] == allowed_mult else None)
    if unlabelled:
        warnings.append(f'{unlabelled} curve label(s) overlap their neighbours where curves bunch at the right edge')
    title = spec.get('title', f'd{{sup:{n}}}, C/B = {c_over_b:g}')
    if title:
        x1, _y1, x2, y2 = c.extent
        c.label(((x1 + x2) / 2, y2 + 16.0), str(title), align='Center')

    def spin_allowed_at(d):
        lv = levels(n, d, c_over_b)
        g = lv[0]
        return [(e - g[0], f'{m}{label}') for e, m, label in lv[1:] if m == g[1]], f'{g[1]}{g[2]}'

    report: dict[str, Any] = {
        'diagram': 'tanabe_sugano', 'd_electrons': n, 'c_over_b': c_over_b,
        'free_ion_terms': free_terms, 'weak_field_ground': f'{weak_ground[1]}{weak_ground[2]}',
        'strong_field_ground': f'{samples[-1][2][1]}{samples[-1][2][2]}', 'spin_crossover': crossovers,
        'terms': sorted({f'{m}{label}' for m, label, _i in curves}), 'warnings': warnings,
    }
    if spec.get('at') is not None:
        d = number(spec['at'], 'at', 0.0, delta_max)
        allowed, ground = spin_allowed_at(d)
        report['at'] = {'delta_over_b': d, 'ground': ground,
                        'spin_allowed': [{'term': t, 'energy_over_b': round(e, 3)} for e, t in allowed]}
    if spec.get('bands') is not None:
        bands = [number(v, 'band', 100.0, 1e6) for v in spec['bands']]
        if len(bands) != 2 or bands[1] <= bands[0]:
            raise ValueError('"bands" is [nu1, nu2] with nu2 > nu1 (cm-1)')
        target = bands[1] / bands[0]
        fits = []

        def ratio(d):
            allowed, _g = spin_allowed_at(d)
            return allowed[1][0] / allowed[0][0] if len(allowed) > 1 and allowed[0][0] > 1e-9 else None

        previous = None
        for d in [delta_max * k / 400 for k in range(1, 401)]:
            r = ratio(d)
            if r is not None and previous is not None and (previous[1] - target) * (r - target) <= 0:
                lo, hi = previous[0], d
                for _ in range(50):
                    mid = (lo + hi) / 2
                    rm = ratio(mid)
                    if rm is None:
                        break
                    if (rm - target) * (previous[1] - target) > 0:
                        lo = mid
                    else:
                        hi = mid
                d_fit = (lo + hi) / 2
                allowed, ground = spin_allowed_at(d_fit)
                b = bands[0] / allowed[0][0]
                fits.append({'delta_over_b': round(d_fit, 3), 'B': round(b, 1), 'Delta': round(d_fit * b, 0),
                             'ground': ground, 'assignment': [allowed[0][1], allowed[1][1]]})
            previous = (d, r) if r is not None else None
        report['band_fit'] = fits or 'no Delta/B in range reproduces the band ratio'
    return _result(c, output, report)
