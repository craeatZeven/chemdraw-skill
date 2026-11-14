# -*- coding: utf-8 -*-
"""文字：标记解析、中英文分字体、宽度/上下伸估计（Arial / SimSun 10 pt）。

标记语法与 chemdraw skill 的 acs_builder 相同：{i:..} {b:..} {bi:..} {sub:..} {sup:..}，不可嵌套。
中文（CJK 表意字、CJK 标点、全角字符）写成 SimSun（font 5），其余写成 Arial（font 3），
规则同 figure-guide 1.3 / acs_builder.runs()。

宽度表是 Helvetica/Arial 的 AFM 字宽（1/1000 em）。2026-11-09 用 ChemDraw 原生渲染实测：
  'TS1' 粗体墨迹宽 16.4 pt（估 18.3，含字侧空），'反应进程' 39.1 pt（估 40），'−12.3' 24.2 pt（估 25.3）；
  基线以上：拉丁大写 7.2–7.4 pt，中文 8.2–8.6 pt，上标顶 9.3 pt；基线以下：中文 0.8–1.2 pt。
估计值只用于排版间距和重叠预警，最终以渲染后的像素核对为准。
"""
import re

FONT_LATIN, FONT_CJK = 3, 5
SIZE = 10.0
LINE_PITCH = 12.0          # figure-guide 1.1：Captions line height 12 pt

TOK = re.compile(r'\{(i|b|bi|sub|sup):([^{}]*)\}')
FACE = {'': 0, 'b': 1, 'i': 2, 'bi': 3, 'sub': 32, 'sup': 64}
CJK_RANGE = '　-〿㐀-鿿＀-￯'
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
    """'{i:G}{sub:a}' -> [('G', 2), ('a', 32)]；base 是整串的底色字形（如 1=粗体）。"""
    if not isinstance(markup, str):
        raise ValueError('文字必须是字符串：%r' % (markup,))
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
            raise ValueError('标记写法不对（只支持 {i:} {b:} {bi:} {sub:} {sup:}，不能嵌套）：%r' % markup)
    return [(t, f) for t, f in out if t]


def plain(runs):
    return ''.join(t for t, _ in runs)


def metrics(runs, size=SIZE):
    """返回 (前进宽度, 基线以上高度, 基线以下深度)，单位 pt。"""
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
    """每个 run 再按中英文切开：中文 SimSun(5)，其余 Arial(3)。"""
    c = ' color="%s"' % color_idx if color_idx else ''
    out = []
    for t, f in runs:
        for piece in _SPLIT.findall(t):
            font = FONT_CJK if is_cjk(piece[0]) else FONT_LATIN
            out.append('<s font="%d" size="%g" face="%d"%s>%s</s>' % (font, size, f, c, _esc(piece)))
    return ''.join(out)
