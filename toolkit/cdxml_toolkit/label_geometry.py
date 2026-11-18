"""Where ChemDraw sets the glyphs of an atom label, and of free text, at 10 pt Arial.

Measured on ChemDraw's own saves (assets/coverage/*.native.cdxml): an atom label's glyphs run
from 4.6 pt above the atom to 4.0 pt below it (subscripts reach about 2 pt lower); a
left-justified label centres its first character on the atom and runs right, a right-justified
one centres its last character on the atom and runs left. Widths are Arial's (per 1000 em),
subscript digits at 0.75 of that. Placement code and checks share this one model, so a check
measures what the placer meant.
"""

from __future__ import annotations

SIZE = 10.0
SUBSCRIPT = 0.75
ABOVE = 4.62     # glyph top above an atom
BELOW = 4.0      # glyph bottom below an atom
SUB_DROP = 2.1   # how much lower subscripts reach

_ARIAL = {
    **dict(zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
               [667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778,
                722, 667, 611, 722, 667, 944, 667, 667, 611])),
    **dict(zip("abcdefghijklmnopqrstuvwxyz",
               [556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556,
                333, 500, 278, 556, 500, 722, 500, 500, 500])),
    **{d: 556 for d in "0123456789"},
    "[": 278, "]": 278, "(": 333, ")": 333, "+": 584, "-": 333, "−": 584, " ": 278,
    ",": 278, ".": 278, "=": 584, "#": 556, "*": 389, "/": 278, "'": 191,
}


def char_width(ch: str, size: float = SIZE, subscript: bool = False) -> float:
    return _ARIAL.get(ch, 556) * size / 1000.0 * (SUBSCRIPT if subscript else 1.0)


_SIGNS = "+-\u2212"


def _widths(text: str, size: float) -> list[float]:
    """Widths of a label's characters: digits after a letter or ')' are subscripts (CH2OH) and
    a trailing charge sign is a superscript (O-, N+), both at 0.75."""
    out = []
    for i, ch in enumerate(text):
        small = i > 0 and ((ch.isdigit() and (text[i - 1].isalpha() or text[i - 1] in ")]" or text[i - 1].isdigit()))
                           or (ch in _SIGNS and i == len(text) - 1))
        out.append(char_width(ch, size, small))
    return out


def label_box(text: str, ax: float, ay: float, justification: str = "Left",
              size: float = SIZE) -> tuple[float, float, float, float]:
    """The glyph box of an atom label at (ax, ay)."""
    widths = _widths(text, size) or [char_width("C", size)]
    total = sum(widths)
    if justification == "Right":
        # The last character that is not a trailing charge sign sits on the atom (O in "O-"),
        # and everything else of the label runs to its left (natively -6.39..3.89 for O-).
        last = len(text) - 1 if text and text[-1] not in _SIGNS else max(0, len(text) - 2)
        x2 = ax + widths[last] / 2.0
        x1 = x2 - total
    else:
        x1 = ax - widths[0] / 2.0
        x2 = x1 + total
    scale = size / SIZE
    # Subscripts (CH2OH: natively 6.1 pt below the atom), brackets and descenders reach lower.
    drops = any(ch.isdigit() for ch in text[1:]) or any(ch in "[]()gjpqy" for ch in text)
    return (x1, ay - ABOVE * scale, x2, ay + (BELOW + (SUB_DROP if drops else 0.0)) * scale)


def text_box(text: str, px: float, py: float, justification: str = "Left",
             size: float = SIZE) -> tuple[float, float, float, float]:
    """The glyph box of free text whose baseline starts (or is centred, or ends) at (px, py)."""
    lines = text.split("\n") or [text]
    width = max(sum(char_width(c, size) for c in line) for line in lines)
    if justification in ("Center", "Centre"):
        x1 = px - width / 2.0
    elif justification == "Right":
        x1 = px - width
    else:
        x1 = px
    advance = size * 1.15
    return (x1, py - size * 0.76, x1 + width, py + size * 0.22 + advance * (len(lines) - 1))
