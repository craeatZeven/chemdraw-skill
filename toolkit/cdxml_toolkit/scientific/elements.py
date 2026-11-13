"""Element data the diagram generators compute from: symbols, periodic-table positions,
ground-state electron configurations, valence-electron counts and d-electron counts."""
from __future__ import annotations

import re

SYMBOLS = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
    'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
    'Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr '
    'Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'
).split()
ATOMIC_NUMBER = {symbol: z for z, symbol in enumerate(SYMBOLS, start=1)}
NOBLE_GASES = (2, 10, 18, 36, 54, 86, 118)
L_LETTERS = 'spdf'
CAPACITY = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
# Madelung (n + l, then n) filling order.
MADELUNG = sorted(((n, l) for n in range(1, 8) for l in range(0, min(n, 4)) if n + l <= 8),
                  key=lambda nl: (nl[0] + nl[1], nl[0]))
# Ground states that depart from the Madelung order: Z -> {subshell: electrons}.
_EXCEPTIONS = {
    24: {'3d': 5, '4s': 1}, 29: {'3d': 10, '4s': 1},
    41: {'4d': 4, '5s': 1}, 42: {'4d': 5, '5s': 1}, 44: {'4d': 7, '5s': 1}, 45: {'4d': 8, '5s': 1},
    46: {'4d': 10, '5s': 0}, 47: {'4d': 10, '5s': 1},
    57: {'4f': 0, '5d': 1}, 58: {'4f': 1, '5d': 1}, 64: {'4f': 7, '5d': 1},
    78: {'5d': 9, '6s': 1}, 79: {'5d': 10, '6s': 1},
    89: {'5f': 0, '6d': 1}, 90: {'5f': 0, '6d': 2}, 91: {'5f': 2, '6d': 1}, 92: {'5f': 3, '6d': 1},
    93: {'5f': 4, '6d': 1}, 96: {'5f': 7, '6d': 1}, 103: {'6d': 0, '7p': 1},
}


def atomic_number(element):
    """Atomic number of a symbol (``'Fe'``) or of an integer given as such."""
    if isinstance(element, int) and not isinstance(element, bool):
        if not 1 <= element <= len(SYMBOLS):
            raise ValueError(f'There is no element with atomic number {element}')
        return element
    symbol = str(element).strip()
    if symbol not in ATOMIC_NUMBER:
        raise ValueError(f'Unknown element symbol {element!r}')
    return ATOMIC_NUMBER[symbol]


def subshell_name(n, l):
    return f'{n}{L_LETTERS[l]}'


def _madelung_fill(electrons):
    config = {}
    for n, l in MADELUNG:
        if electrons <= 0:
            break
        take = min(electrons, CAPACITY[L_LETTERS[l]])
        config[subshell_name(n, l)] = take
        electrons -= take
    if electrons > 0:
        raise ValueError('Too many electrons for the known subshells')
    return config


def electron_configuration(element, charge=0):
    """Ground-state configuration as ``{'1s': 2, ...}`` ordered by shell, then subshell.

    Cations lose electrons from the outermost shell first (4s before 3d); anions fill on
    in Madelung order. Known anomalous ground states (Cr, Cu, Pd, ...) are built in.
    """
    z = atomic_number(element)
    charge = int(charge)
    if charge >= z:
        raise ValueError(f'{SYMBOLS[z - 1]} cannot lose {charge} electrons')
    config = _madelung_fill(z)
    for name, count in _EXCEPTIONS.get(z, {}).items():
        config[name] = count
    config = {k: v for k, v in config.items() if v}
    if charge > 0:
        for _ in range(charge):
            outer = max(config, key=lambda name: (int(name[0]), L_LETTERS.index(name[1])))
            config[outer] -= 1
            if not config[outer]:
                del config[outer]
    elif charge < 0:
        for _ in range(-charge):
            for n, l in MADELUNG:
                name = subshell_name(n, l)
                if config.get(name, 0) < CAPACITY[L_LETTERS[l]]:
                    config[name] = config.get(name, 0) + 1
                    break
            else:
                raise ValueError('Too many electrons for the known subshells')
    return dict(sorted(config.items(), key=lambda kv: (int(kv[0][0]), L_LETTERS.index(kv[0][1]))))


def noble_gas_core(element, charge=0):
    """``(symbol, core_subshells)`` of the largest noble-gas core inside the species."""
    electrons = atomic_number(element) - int(charge)
    inside = [z for z in NOBLE_GASES if z < electrons]
    if not inside:
        return None, {}
    z = inside[-1]
    return SYMBOLS[z - 1], _madelung_fill(z)


def group_period(element):
    """``(group, period)`` in the 18-column table; lanthanoids and actinoids have group ``None``."""
    z = atomic_number(element)
    period = next(i + 1 for i, last in enumerate(NOBLE_GASES) if z <= last)
    start = NOBLE_GASES[period - 2] + 1 if period > 1 else 1
    offset = z - start
    if period == 1:
        return (1 if z == 1 else 18), 1
    if period in (2, 3):
        return (offset + 1 if offset < 2 else offset + 11), period
    if period in (4, 5):
        return offset + 1, period
    if offset < 2:
        return offset + 1, period
    if offset < 17:
        return (None if offset < 16 else 3), period   # f block; Lu and Lr sit in group 3
    return offset - 13, period


def valence_electrons(element):
    """Valence electrons of a main-group element (Lewis structures)."""
    group, period = group_period(element)
    z = atomic_number(element)
    if z == 2:
        return 2
    if group is None or 3 <= group <= 12:
        raise ValueError(f'{SYMBOLS[z - 1]} is not a main-group element; give "valence" explicitly')
    return group if group <= 2 else group - 10


_ROMAN = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, '0': 0}


def parse_ion(text):
    """``'Fe3+'``, ``'Fe(III)'``, ``'Cl-'`` or ``'Cu'`` -> ``(symbol, charge)``."""
    value = str(text).strip().replace('−', '-')
    match = re.fullmatch(r'([A-Z][a-z]?)\s*\(([IV0]+)\)', value)
    if match and match.group(2) in _ROMAN:
        symbol, charge = match.group(1), _ROMAN[match.group(2)]
    else:
        match = re.fullmatch(r'([A-Z][a-z]?)\s*\^?\s*(\d*)([+-])?', value)
        if not match or (match.group(2) and not match.group(3)):
            raise ValueError(f'Cannot read the ion {text!r}; write it like "Fe3+" or "Fe(III)"')
        symbol = match.group(1)
        size = int(match.group(2) or 1) if match.group(3) else 0
        charge = size if match.group(3) == '+' else -size
    atomic_number(symbol)
    return symbol, charge


def d_electron_count(element, oxidation_state=0):
    """d-electron count of a transition-metal centre: group number minus oxidation state."""
    group, _period = group_period(element)
    symbol = SYMBOLS[atomic_number(element) - 1]
    if group is None or not 3 <= group <= 12:
        raise ValueError(f'{symbol} is not a d-block element; give "d_electrons" explicitly')
    count = group - int(oxidation_state)
    if not 0 <= count <= 10:
        raise ValueError(f'{symbol} in oxidation state {oxidation_state:+d} would have {count} d electrons')
    return count
