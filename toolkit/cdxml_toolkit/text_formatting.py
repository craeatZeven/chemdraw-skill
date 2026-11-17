"""text_formatting.py — Shared chemical text formatting for ChemDraw CDXML.

Provides functions for building properly formatted <s> (styled text run)
elements in CDXML, handling two chemistry-specific typographic conventions:

1. **Subscript digits in chemical formulas.**
   In chemical notation, digits that follow letters are molecular counts and
   must be rendered as subscripts: "CH3OH" → "CH₃OH", "Pd2(dba)3" → "Pd₂(dba)₃".
   Plain numbers (temperatures "80 °C", durations "2 h", percentages "95%")
   are left as normal text.

2. **Italic prefixes in IUPAC / organic nomenclature.**
   Stereochemical descriptors, positional locants, and heteroatom locants at
   the start of a reagent name are italicised per IUPAC convention:
   "n-BuLi" → "*n*-BuLi", "tert-BuOH" → "*tert*-BuOH", "N-Boc" → "*N*-Boc".

ChemDraw CDXML face codes used:
  - face="96"  (0x60 = Formula)   — normal reagent text
  - face="32"  (0x20 = Subscript) — subscript digits
  - face="2"   (0x02 = Italic)    — italic prefix runs

Previously duplicated across scheme_polisher.py and reaction_from_image.py.
Consolidated here for v0.3.
"""

from __future__ import annotations

import re
from typing import Tuple
from xml.sax.saxutils import escape as xml_escape

# ---------------------------------------------------------------------------
# Regex: letter (or closing paren) followed by one or more digits.
# Matches subscriptable digit groups in chemical formulas.
# Examples:  CH3 → ("H", "3"),  Pd2 → ("d", "2"),  (dba)3 → (")", "3")
# ---------------------------------------------------------------------------
SUBSCRIPT_RE = re.compile(r'([A-Za-z)])(\d+)')

# Keep underscore-prefixed alias for backward compatibility with callers that
# import the private name directly.
_SUBSCRIPT_RE = SUBSCRIPT_RE

# One-letter element symbols: a lone "H2", "N2", "O2" or "I2" is a formula.
_SINGLE_LETTER_ELEMENTS = frozenset("HBCNOFPSKVYIWU")


def _is_compound_label(text: str, m: re.Match) -> bool:
    """True when a match is a compound or ligand label such as "L1", "P1" or "A2".

    The token must be one capital letter plus digits standing alone (bounded by
    start/end, space or punctuation), and either not an element symbol or
    followed by "1", which formulas never write.
    """
    letter, digits = m.group(1), m.group(2)
    if not letter.isupper():
        return False
    start, end = m.start(1), m.end(2)
    if start > 0 and (text[start - 1].isalnum() or text[start - 1] in ")]"):
        return False
    if end < len(text) and text[end].isalnum():
        return False
    return letter not in _SINGLE_LETTER_ELEMENTS or digits == "1"


# Mechanism designators, where the trailing digit is a step count and stays on the baseline:
# "Sn2" is written S_N2, never as a subscripted sulfur count. "S2" (disulfide) is a formula and
# still subscripts, so only these whole tokens are exempt.
_MECHANISM_TOKEN = re.compile(r'(?<![\w])Sn[12](?![\w])')


def subscript_matches(text: str) -> list[re.Match]:
    """Formula digit groups in *text* that should be subscripted."""
    exempt = [(m.start(), m.end()) for m in _MECHANISM_TOKEN.finditer(text)]
    return [m for m in SUBSCRIPT_RE.finditer(text)
            if not _is_compound_label(text, m) and not any(a <= m.start(1) < b for a, b in exempt)]

# ---------------------------------------------------------------------------
# Italic prefixes recognised in organic chemistry nomenclature.
# Matched at the start of the display name, case-sensitive.
# Longer forms come first so "tert-" is tried before "t-".
# ---------------------------------------------------------------------------
ITALIC_PREFIXES: list[str] = [
    "tert-", "sec-", "iso-",          # long forms first
    "n-", "t-", "s-", "i-",           # single-letter alkyl descriptors
    "o-", "m-", "p-",                 # arene positional (ortho/meta/para)
    "cis-", "trans-",
    "rac-", "meso-",
    "R-", "S-",
    "syn-", "anti-",
    "exo-", "endo-",
    "E-", "Z-",
    "D-", "L-",
    "N-", "O-", "S-", "C-", "P-",    # heteroatom locants (N-Boc, O-alkyl …)
]

_ITALIC_PREFIXES = ITALIC_PREFIXES  # backward-compat alias

# ---------------------------------------------------------------------------
# Notation used in reaction conditions: temperature, room temperature, times.
# The number must not be glued to a preceding word character, so formulas such
# as "H2C" or "C10H21Br" are never read as temperatures.
# ---------------------------------------------------------------------------
_TEMPERATURE_RE = re.compile(
    r'(?<![\w.])([-−]?\d+(?:\.\d+)?)\s*(?:°|deg(?:rees?)?\.?)?\s*C(?![\w-])'
)
_ROOM_TEMP_RE = re.compile(
    r'(?<!\w)(?:r\.t\.|RT|room\s+temp(?:erature)?\.?)(?!\w)'
)
_HOURS_RE = re.compile(r'(?<![\w.])(\d+(?:\.\d+)?)\s*(?:hrs?|hours?)(?!\w)')
_MINUTES_RE = re.compile(r'(?<![\w.])(\d+(?:\.\d+)?)\s*(?:mins|minutes?)(?!\w)')


def normalize_condition_text(text: str) -> str:
    """Write temperatures, room temperature and times the way papers do.

    "80 C" / "80C" / "80 deg C" -> "80 °C"; "RT" / "r.t." / "room temperature"
    -> "rt"; "12 hr" / "12 hours" -> "12 h"; "30 mins" -> "30 min". Anything
    else, including formulas and reagent names, is returned unchanged, and the
    result is stable when normalised again.
    """
    text = _TEMPERATURE_RE.sub(lambda m: f"{m.group(1)} °C", text)
    text = _ROOM_TEMP_RE.sub("rt", text)
    text = _HOURS_RE.sub(lambda m: f"{m.group(1)} h", text)
    return _MINUTES_RE.sub(lambda m: f"{m.group(1)} min", text)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def needs_subscript(text: str) -> bool:
    """Determine whether *text* contains chemical-formula digits that should
    be rendered as subscripts in ChemDraw.

    Returns ``True`` for reagent formulas like ``"CH3OH"``, ``"Cs2CO3"``,
    ``"Pd2(dba)3"`` where trailing digits represent atom counts.

    Returns ``False`` for non-formula text that happens to contain digits:

    * Temperatures — ``"80 °C"``
    * Durations — ``"2 h"``, ``"30 min"``
    * Percentages — ``"95%"``
    * Pure-numeric / unit-only strings — ``"120 °C, 2 h"``

    Examples::

        >>> needs_subscript("Et3N")
        True
        >>> needs_subscript("DMF")
        False
        >>> needs_subscript("80 °C")
        False
    """
    # Temperature (digits before °)
    if re.search(r'\d+\s*°', text):
        return False
    # Duration (digits before h/m at word boundary)
    if re.search(r'\d+\s*[hm](?:\s|$|,)', text):
        return False
    # Percentage
    if re.search(r'\d+\s*%', text):
        return False
    # Pure numeric / unit strings like "reflux", "rt", "120 °C, 2 h"
    if re.fullmatch(r'[\d\s.,°ChmsMinHr/]+', text, re.IGNORECASE):
        return False
    return bool(subscript_matches(text))


_ELEMENTS = frozenset("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd
Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf
Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og""".split())  # no D or T: TS2 is a number
# Groups that a formula label writes as one token (PhB(OH)2, Pd(OAc)2, Pd2(dba)3, n-Bu3SnH).
_FORMULA_GROUPS = frozenset("""Ph Me Et Pr Bu Ar Ac Bn Bz Ts Ms Tf Tol Mes Cy Cp Py py Boc Cbz Piv Fmoc Tr
Trt TMS TBS TIPS TES OAc OTf OMe OEt dba cod acac pin Bpin dppf dppe dppp bipy phen Ad Np Xyl""".split())
_FORMULA_PUNCTUATION = set("()[]{}\u00b7.+-\u2212\u2013\u2014=")


def is_formula_label(text: str) -> bool:
    """Whether a label is a formula, to be set with subscripts like condition text.

    Conditions are formulas almost always; labels are compound numbers almost always, and the
    condition rule would print TS2, INT3 or IM1 with subscripts. So a label counts as a formula
    only when, past an italic prefix (``n-``), it is made entirely of element symbols, common
    groups (Ph, Et, OAc, dba, ...), brackets and counts, and has a count to subscript: PhB(OH)2,
    K2CO3 and Et3N do; TS2, IM1, 3a and L1 do not. A lone element with a count other than 1
    (S2, N2) reads as a formula, as it does in conditions.
    """
    rest = split_italic_prefix(text)[1]
    if not rest or not subscript_matches(rest):
        return False
    i = 0
    while i < len(rest):
        ch = rest[i]
        if ch.isdigit() or ch in _FORMULA_PUNCTUATION:
            i += 1
            continue
        for size in (4, 3, 2, 1):
            token = rest[i:i + size]
            if len(token) == size and (token in _FORMULA_GROUPS or token in _ELEMENTS):
                i += size
                break
        else:
            return False
    return True


# Private-name alias for callers that import ``_needs_subscript``.
_needs_subscript = needs_subscript


def split_italic_prefix(text: str) -> Tuple[str, str]:
    """Split *text* into ``(italic_prefix, remainder)`` if it starts with a
    recognised chemistry italic prefix (see :data:`ITALIC_PREFIXES`).

    Returns ``("", text)`` when no prefix matches.

    Examples::

        >>> split_italic_prefix("n-BuLi")
        ('n-', 'BuLi')
        >>> split_italic_prefix("tert-BuOH")
        ('tert-', 'BuOH')
        >>> split_italic_prefix("Cs2CO3")
        ('', 'Cs2CO3')
    """
    for prefix in ITALIC_PREFIXES:
        if text.startswith(prefix):
            return prefix, text[len(prefix):]
    return "", text


_split_italic_prefix = split_italic_prefix  # backward-compat alias


def build_formatted_s_xml(
    text: str,
    font: str = "3",
    size: str = "10",
    color: str = "0",
    italic_font: str | None = None,
) -> str:
    """Build one or more CDXML ``<s>`` elements with correct chemical styling.

    This is the primary text-formatting entry point. It handles:

    1. **Italic prefix** (``n-``, ``tert-``, ``sec-``, ``N-``, …) rendered
       with ``face="2"`` (Italic).
    2. **Subscript digits** after letters/closing-parens rendered with
       ``face="32"`` (Subscript).
    3. **Normal formula text** rendered with ``face="96"`` (Formula).

    Parameters
    ----------
    text : str
        The display text for a reagent or chemical name (e.g. ``"n-BuLi"``,
        ``"Cs2CO3"``, ``"Pd2(dba)3"``).
    font : str
        CDXML font id for normal + subscript runs (default ``"3"`` = Arial).
    size : str
        Font size in points (default ``"10"``).
    color : str
        CDXML color id (default ``"0"`` = black).
    italic_font : str or None
        If given, use this font id for the italic prefix run instead of
        *font*. Useful when the italic style lives in a separate font entry.

    Returns
    -------
    str
        Raw XML string of ``<s>`` elements ready to embed inside a ``<t>``
        element.  Example for ``"n-BuLi"``::

            <s font="3" size="10" color="0" face="2">n-</s>
            <s font="3" size="10" color="0" face="96">BuLi</s>

    Notes
    -----
    The function is XML-safe: all text content is escaped via
    ``xml.sax.saxutils.escape``.
    """
    italic_prefix, rest = split_italic_prefix(text)
    ifont = italic_font if italic_font is not None else font

    parts: list[str] = []

    # ---- italic prefix ----
    if italic_prefix:
        parts.append(
            f'<s font="{ifont}" size="{size}" color="{color}" '
            f'face="2">{xml_escape(italic_prefix)}</s>'
        )

    # ---- remainder with subscript handling ----
    # Formula face (96) makes ChemDraw subscript any digit that follows a
    # letter on its own, so a label such as "L1" or "P1" must leave the
    # formula face: it is written as a plain-face (0) run instead.
    if rest:
        def run(chunk: str, face: str) -> None:
            if chunk:
                parts.append(
                    f'<s font="{font}" size="{size}" color="{color}" '
                    f'face="{face}">{xml_escape(chunk)}</s>'
                )

        pos = 0
        for m in SUBSCRIPT_RE.finditer(rest):
            if _is_compound_label(rest, m):
                run(rest[pos:m.start()], "96")
                run(m.group(0), "0")
            else:
                run(rest[pos:m.start(2)], "96")
                run(m.group(2), "32")
            pos = m.end()
        run(rest[pos:], "96")

    return "".join(parts)


# Backward-compatible aliases (used by scheme_polisher and reaction_from_image).
_build_formatted_s_xml = build_formatted_s_xml
build_subscripted_s_xml = build_formatted_s_xml
_build_subscripted_s_xml = build_formatted_s_xml


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    passed = 0
    failed = 0

    def check(label: str, got, expected):
        global passed, failed
        if got == expected:
            print(f"  PASS  {label}")
            passed += 1
        else:
            print(f"  FAIL  {label}")
            print(f"        expected: {expected!r}")
            print(f"        got:      {got!r}")
            failed += 1

    print("text_formatting.py self-test")
    print("=" * 50)

    # --- needs_subscript ---
    check("needs_subscript('CH3OH')", needs_subscript("CH3OH"), True)
    check("needs_subscript('DMF')", needs_subscript("DMF"), False)
    check("needs_subscript('Et3N')", needs_subscript("Et3N"), True)
    check("needs_subscript('Cs2CO3')", needs_subscript("Cs2CO3"), True)
    check("needs_subscript('80 °C')", needs_subscript("80 °C"), False)
    check("needs_subscript('2 h')", needs_subscript("2 h"), False)
    check("needs_subscript('95%')", needs_subscript("95%"), False)

    # --- split_italic_prefix ---
    check("split_italic_prefix('n-BuLi')", split_italic_prefix("n-BuLi"), ("n-", "BuLi"))
    check("split_italic_prefix('Cs2CO3')", split_italic_prefix("Cs2CO3"), ("", "Cs2CO3"))
    check("split_italic_prefix('tert-BuOH')", split_italic_prefix("tert-BuOH"), ("tert-", "BuOH"))
    check("split_italic_prefix('N-Boc')", split_italic_prefix("N-Boc"), ("N-", "Boc"))

    # --- build_formatted_s_xml ---
    xml_et3n = build_formatted_s_xml("Et3N")
    check("build_formatted_s_xml('Et3N') contains <s>",
          "<s " in xml_et3n, True)
    check("build_formatted_s_xml('Et3N') has subscript face",
          'face="32"' in xml_et3n, True)
    check("build_formatted_s_xml('Et3N') has formula face",
          'face="96"' in xml_et3n, True)

    xml_nbuli = build_formatted_s_xml("n-BuLi")
    check("build_formatted_s_xml('n-BuLi') has italic face",
          'face="2"' in xml_nbuli, True)
    check("build_formatted_s_xml('n-BuLi') italic run contains 'n-'",
          'face="2">n-</s>' in xml_nbuli, True)

    xml_dmf = build_formatted_s_xml("DMF")
    check("build_formatted_s_xml('DMF') — no subscript for plain text",
          'face="32"' in xml_dmf, False)

    # --- aliases ---
    check("build_subscripted_s_xml is build_formatted_s_xml",
          build_subscripted_s_xml is build_formatted_s_xml, True)
    check("_build_formatted_s_xml is build_formatted_s_xml",
          _build_formatted_s_xml is build_formatted_s_xml, True)

    print("=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    if failed:
        raise SystemExit(1)
    print("All tests passed.")
