# -*- coding: utf-8 -*-
"""Text: markup parsing, Latin/CJK font split, width and ascent/descent estimates (Arial / SimSun 10 pt).

The markup is the same as acs_builder's: {i:..} {b:..} {bi:..} {sub:..} {sup:..}, not nested.
CJK text (ideographs, CJK punctuation, full-width forms) is set in SimSun (font 5), the rest in Arial (font 3),
as acs_builder.runs() does.

The width table holds Helvetica/Arial AFM advances (1/1000 em). Measured in native ChemDraw renders:
  bold 'TS1' inks 16.4 pt (estimate 18.3 with side bearings), four CJK characters 39.1 pt (estimate 40), '-12.3' 24.2 pt (estimate 25.3);
  above the baseline: Latin capitals 7.2-7.4 pt, CJK 8.2-8.6 pt, superscript top 9.3 pt; below it: CJK 0.8-1.2 pt.
Estimates only set spacing and overlap warnings; the pixel check of the render is what counts.
"""
import re

FONT_LATIN, FONT_CJK = 3, 5
SIZE = 10.0
LINE_PITCH = 12.0          # ACS 1996 caption line height, 12 pt

TOK = re.compile(r'\{(i|b|bi|sub|sup):([^{}]*)\}')
FACE = {'': 0, 'b': 1, 'i': 2, 'bi': 3, 'sub': 32, 'sup': 64}
CJK_RANGE = '\u3000-\u303f\u3400-\u9fff\uff00-\uffef'
_CJK_CHAR = re.compile('[%s]' % CJK_RANGE)
_SPLIT = re.compile('[%s]+|[^%s]+' % (CJK_RANGE, CJK_RANGE))

_W = {}
for _c in '0123456789#$?_':
    _W[_c] = 556
_W.update({' ': 278, '!': 278, '"': 355, '%': 889, '&': 667, "'": 191, '(': 333, ')': 333,
           '*': 389, '+': 584, ',': 278, '.': 278, '-': 333, '/': 278, ':': 278, ';': 278,
           '<': 584, '=': 584, '>': 584, '@': 1015, '[': 278, ']': 278})
_W.update(dict(A=667, B=667, C=722, D=722, E=667, F=611, G=778, H=722, I=278, J=500, K=667,
               L=556, M=833, N=722, O=778, P=667, Q=778, R=722, S=667, T=611, U=722, V=667,
               W=944, X=667, Y=667, Z=611))
_W.update(dict(a=556, b=556, c=500, d=556, e=556, f=278, g=556, h=556, i=222, j=222, k=500,
               l=222, m=833, n=556, o=556, p=556, q=556, r=333, s=500, t=278, u=556, v=500,
               w=722, x=500, y=500, z=500))
_W.update({'−': 584, '–': 556, '—': 1000, '·': 278, '•': 350,
           '°': 400, '‡': 556, '†': 556, '′': 191, '″': 355,
           'Δ': 612, 'α': 580, 'β': 575, 'γ': 500, 'δ': 556,
           'ε': 450, 'φ': 640, 'ν': 500, '±': 584, '×': 584,
           '→': 1000, '←': 1000, '≤': 584, '≥': 584, '≈': 584})
_DESC = set('gjpqy(),;[]_')


def is_cjk(ch):
    return bool(_CJK_CHAR.match(ch))


def parse(markup, base=0):
    """'{i:G}{sub:a}' -> [('G', 2), ('a', 32)]; base is the face of the whole string (e.g. 1 = bold)."""
    if not isinstance(markup, str):
        raise ValueError('text must be a string: %r' % (markup,))
    out, pos = [], 0
    for m in TOK.finditer(markup):
        if m.start() > pos:
            out.append((markup[pos:m.start()], base))
        out.append((m.group(2), FACE[m.group(1)] | base))
        pos = m.end()
    if pos < len(markup):
        out.append((markup[pos:], base))
    for t, _ in out:
        if '{' in t or '}' in t:
            raise ValueError('bad markup (only {i:} {b:} {bi:} {sub:} {sup:}, not nested): %r' % markup)
    return [(t, f) for t, f in out if t]


def plain(runs):
    return ''.join(t for t, _ in runs)


def metrics(runs, size=SIZE):
    """Returns (advance width, height above the baseline, depth below it) in pt."""
    k = size / 10.0
    width, asc, desc = 0.0, 0.0, 0.3
    for t, f in runs:
        script = f & 96
        shrink = 0.75 if script in (32, 64) else 1.0
        bold = 1.04 if f & 1 else 1.0
        for ch in t:
            if is_cjk(ch):
                cw, a, d = 1000, 8.6, 1.3
            else:
                cw = _W.get(ch, 600)
                a = 7.4
                d = 2.1 if ch in _DESC else 0.3
            if script == 64:
                a, d = 9.4, 0.0
            elif script == 32:
                a, d = 5.5, 3.0
            width += cw * shrink * bold * size / 1000.0
            asc = max(asc, a * k)
            desc = max(desc, d * k)
    return width, asc, desc


def _esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')


def runs_xml(runs, color_idx=None, size=SIZE):
    """Split each run into CJK (SimSun, 5) and other (Arial, 3) pieces."""
    c = ' color="%s"' % color_idx if color_idx else ''
    out = []
    for t, f in runs:
        for piece in _SPLIT.findall(t):
            font = FONT_CJK if is_cjk(piece[0]) else FONT_LATIN
            out.append('<s font="%d" size="%g" face="%d"%s>%s</s>' % (font, size, f, c, _esc(piece)))
    return ''.join(out)
