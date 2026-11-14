# -*- coding: utf-8 -*-
"""势能图生成器：JSON → 可编辑 CDXML（ACS Document 1996）→ ChemDraw 原生渲染 PNG → 数值检查。

    D:/python/python.exe energy_profile.py spec.json [-o out.cdxml] [--no-render] [--force] [--keep-work DIR]

产物（与 out.cdxml 同目录同名）：out.cdxml、out.png、out.report.json（全部检查结果与映射参数）。
退出码：0 全部通过；1 有检查失败；2 输入或环境错误。

分工（chemdraw-workbench 规矩）：
  * 平滑曲线：chemdraw skill 的 scientific_workflow(operation=plot) 生成原生 <curve>（ep_plot.py）；
  * 坐标轴、能级短横线、连接线、能垒箭头、文字：本脚本按 figure-guide 第 1 节规格写原生对象；
  * 验收：cdx.py render 原生渲染 + 读回 CDXML 自检 + check_cdxml.py + PNG 像素统计（ep_verify.py）。
JSON 格式见同目录 README.md。
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ep_text as T          # noqa: E402
import ep_plot               # noqa: E402
import ep_verify             # noqa: E402
import ep_structure          # noqa: E402

VERSION = '1.0'
# Resolve the skill this script ships inside first, so a git checkout works without
# installing the Skill anywhere; the home installs stay as fallbacks.
SKILL_DIRS = [os.path.dirname(os.path.dirname(HERE)),
              os.path.expanduser('~/.claude/skills/chemdraw'),
              os.path.expanduser('~/.agents/skills/chemdraw')]
# cdx.py re-launches inside the toolkit runtime, so the interpreter running this
# script is already the right one. Never hard-code a machine's Python path.
CDX_PY = sys.executable
# Optional external ACS-spec checker. Absent on most machines: ep_verify skips it
# and says so in the report rather than failing the run.
CHECKER = next((c for c in (
    os.environ.get('CHEMDRAW_CHECK_CDXML'),
    r'D:\Work\tools\chemdraw-standard\check_cdxml.py',
) if c and os.path.isfile(c)), None)

# figure-guide 1.1 / acs_builder.to_cdxml() 的根属性，原样照抄
ROOT_ATTRS = ('BondLength="14.4" LabelFont="3" LabelSize="10" LabelFace="96" CaptionFont="3" CaptionSize="10" '
              'LineWidth="0.6" BoldWidth="2" HashSpacing="2.5" MarginWidth="1.6" BondSpacing="18" ChainAngle="120" '
              'LabelJustification="Auto" CaptionJustification="Left" color="0" bgcolor="1"')
FONTTABLE = ('<fonttable><font id="3" charset="iso-8859-1" name="Arial"/>'
             '<font id="5" charset="Unicode" name="SimSun"/></fonttable>')
PALETTE = ['#000000', '#CC0000', '#808080', '#336699', '#FF6600']   # figure-guide 1.2
STRUCTURE_GAP = 10.0        # pt between a stationary point and the structure above it
LW, BOLD_W = 0.6, 2.0
GAP = 3.0                 # 文字与线条的净空
MARGIN = 10.0
HEAD_AXIS, HEAD_ANN = 1000, 800
UNITS = {'kj/mol': 'kJ mol{sup:−1}', 'kj mol-1': 'kJ mol{sup:−1}', 'kcal/mol': 'kcal mol{sup:−1}',
         'kcal mol-1': 'kcal mol{sup:−1}', 'ev': 'eV', 'j/mol': 'J mol{sup:−1}', 'kj': 'kJ'}

TOP_KEYS = {'style', 'unit', 'unit_markup', 'energy_symbol', 'y_label', 'y_label_style', 'x_label',
            'reference', 'decimals', 'show_values', 'value_brackets', 'show_plus', 'label_layout',
            'label_face', 'ts_marker', 'unit_width', 'energy_height', 'tail', 'bars', 'bar_width',
            'curve_width', 'y_ticks', 'y_tick_decimals', 'x_ticks', 'annotation_color', 'species', 'pathways',
            'annotations', 'samples_per_segment'}
SP_KEYS = {'id', 'label', 'energy', 'ts', 'x', 'label_position', 'dx', 'dy', 'color', 'value', 'show_value',
           'structure'}
PW_KEYS = {'id', 'nodes', 'color', 'dashed', 'width', 'label', 'label_x', 'label_position', 'label_dy'}
AN_KEYS = {'type', 'from', 'to', 'label', 'x', 'side', 'at', 'heads', 'color', 'ref', 'label_dx', 'label_dy'}


class SpecError(ValueError):
    pass


def _keys(obj, allowed, where):
    if not isinstance(obj, dict):
        raise SpecError('%s 必须是对象' % where)
    bad = [k for k in obj if k not in allowed and not k.startswith('_')]
    if bad:
        raise SpecError('%s 有不认识的键 %s（允许：%s；以 _ 开头的键当注释忽略）' % (where, bad, sorted(allowed)))


def _color(v, where):
    if not (isinstance(v, str) and len(v) == 7 and v[0] == '#'):
        raise SpecError('%s 的颜色要写成 #RRGGBB：%r' % (where, v))
    int(v[1:], 16)
    return v.upper()


def _num(v, where):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise SpecError('%s 必须是数值：%r' % (where, v))
    return float(v)


# ---------------------------------------------------------------- 规格
def normalize(spec):
    _keys(spec, TOP_KEYS, '顶层')
    o = {
        'style': spec.get('style', 'smooth'), 'unit': spec.get('unit'), 'unit_markup': spec.get('unit_markup'),
        'energy_symbol': spec.get('energy_symbol', 'E'), 'y_label': spec.get('y_label'),
        'y_label_style': spec.get('y_label_style', 'paren'),
        'x_label': spec.get('x_label', 'Reaction coordinate'), 'reference': spec.get('reference'),
        'decimals': int(spec.get('decimals', 1)), 'show_values': bool(spec.get('show_values', True)),
        'value_brackets': bool(spec.get('value_brackets', False)), 'show_plus': bool(spec.get('show_plus', False)),
        'label_layout': spec.get('label_layout', 'stack'), 'label_face': spec.get('label_face', 'bold'),
        'ts_marker': spec.get('ts_marker', 'none'),
        'unit_width': _num(spec.get('unit_width', 52.0), 'unit_width'),
        'energy_height': _num(spec.get('energy_height', 150.0), 'energy_height'),
        'curve_width': _num(spec.get('curve_width', 1.2), 'curve_width'),
        'y_ticks': spec.get('y_ticks'), 'y_tick_decimals': int(spec.get('y_tick_decimals', 0)),
        'x_ticks': bool(spec.get('x_ticks', False)),
        'annotation_color': _color(spec.get('annotation_color', '#336699'), 'annotation_color'),
        'samples_per_segment': int(spec.get('samples_per_segment', 60)),
    }
    if o['style'] not in ('smooth', 'levels'):
        raise SpecError('style 只能是 smooth 或 levels')
    o['tail'] = _num(spec.get('tail', 0.35 if o['style'] == 'smooth' else 0.0), 'tail')
    o['bars'] = bool(spec.get('bars', o['style'] == 'levels'))
    o['bar_width'] = _num(spec.get('bar_width', 20.0 if o['style'] == 'levels' else 14.0), 'bar_width')
    for k, allowed in (('label_layout', ('stack', 'inline', 'split')), ('label_face', ('bold', 'plain', 'italic')),
                       ('ts_marker', ('none', 'dagger')), ('y_label_style', ('paren', 'slash'))):
        if o[k] not in allowed:
            raise SpecError('%s 只能是 %s' % (k, '/'.join(allowed)))

    raw = spec.get('species')
    if not isinstance(raw, list) or len(raw) < 2:
        raise SpecError('species 至少要有 2 个物种')
    species, order = {}, []
    for i, s in enumerate(raw):
        where = 'species[%d]' % i
        _keys(s, SP_KEYS, where)
        sid = s.get('id')
        if not isinstance(sid, str) or not sid or sid in species:
            raise SpecError('%s 的 id 缺失或重复' % where)
        lp = s.get('label_position', 'auto')
        if lp not in ('auto', 'above', 'below', 'none'):
            raise SpecError('%s.label_position 只能是 auto/above/below/none' % where)
        species[sid] = {'id': sid, 'label': s.get('label', sid), 'E_in': _num(s.get('energy'), where + '.energy'),
                        'ts': bool(s.get('ts', False)), 'x': _num(s['x'], where + '.x') if 'x' in s else float(i),
                        'pos': lp, 'dx': _num(s.get('dx', 0), where + '.dx'), 'dy': _num(s.get('dy', 0), where + '.dy'),
                        'color': _color(s['color'], where) if 'color' in s else None,
                        'value': s.get('value'), 'show_value': bool(s.get('show_value', o['show_values'])),
                        'structure': s.get('structure')}
        T.parse(species[sid]['label'])
        order.append(sid)
    ref = o['reference']
    if ref is not None and ref not in species:
        raise SpecError('reference=%r 不是任何物种的 id' % ref)
    base = species[ref]['E_in'] if ref else 0.0
    for s in species.values():
        s['E'] = s['E_in'] - base

    raw_p = spec.get('pathways') or [{'id': 'main', 'nodes': order}]
    pathways, used = [], set()
    for i, p in enumerate(raw_p):
        where = 'pathways[%d]' % i
        _keys(p, PW_KEYS, where)
        nodes = p.get('nodes')
        if not isinstance(nodes, list) or len(nodes) < 2 or any(n not in species for n in nodes):
            raise SpecError('%s.nodes 必须是 ≥2 个已定义物种 id 的列表' % where)
        xs = [species[n]['x'] for n in nodes]
        if any(b <= a for a, b in zip(xs, xs[1:])):
            raise SpecError('%s 的物种 x 必须严格递增：%s' % (where, list(zip(nodes, xs))))
        if o['bars'] and any((b - a) * o['unit_width'] < o['bar_width'] + 4 for a, b in zip(xs, xs[1:])):
            raise SpecError('%s 相邻物种太近，放不下 %g pt 的短横线（加大 unit_width 或减小 bar_width）'
                            % (where, o['bar_width']))
        pathways.append({'id': str(p.get('id', 'path%d' % (i + 1))), 'nodes': list(nodes),
                         'color': _color(p['color'], where) if 'color' in p else PALETTE[i % len(PALETTE)],
                         'dashed': bool(p.get('dashed', False)),
                         'width': _num(p.get('width', o['curve_width']), where + '.width'),
                         'label': p.get('label'), 'label_x': p.get('label_x'),
                         'label_position': p.get('label_position', 'above'),
                         'label_dy': _num(p.get('label_dy', 0), where + '.label_dy')})
        used.update(nodes)
    if len({p['id'] for p in pathways}) != len(pathways):
        raise SpecError('pathways 的 id 重复')
    unused = [s for s in order if s not in used]
    if unused:
        raise SpecError('这些物种不在任何 pathway 里：%s' % unused)

    ann = []
    for i, a in enumerate(spec.get('annotations') or []):
        where = 'annotations[%d]' % i
        _keys(a, AN_KEYS, where)
        if a.get('type', 'gap') != 'gap':
            raise SpecError('%s.type 目前只支持 gap（两物种间的能量差箭头）' % where)
        if a.get('from') not in species or a.get('to') not in species:
            raise SpecError('%s 的 from/to 必须是物种 id' % where)
        x = a.get('x', a['to'])
        if isinstance(x, str) and x not in species:
            raise SpecError('%s.x 必须是物种 id 或数值' % where)
        ann.append({'from': a['from'], 'to': a['to'], 'label': a.get('label', ''), 'x': x,
                    'side': a.get('side', 'right'), 'at': _num(a.get('at', 0.5), where + '.at'),
                    'heads': a.get('heads', 'both'), 'ref': a.get('ref', 'both'),
                    'color': _color(a['color'], where) if 'color' in a else o['annotation_color'],
                    'dx': _num(a.get('label_dx', 0), where), 'dy': _num(a.get('label_dy', 0), where)})
        if ann[-1]['heads'] not in ('both', 'to', 'from', 'none') or ann[-1]['side'] not in ('left', 'right') \
                or ann[-1]['ref'] not in ('both', 'from', 'to', 'none'):
            raise SpecError('%s：heads=both/to/from/none，side=left/right，ref=both/from/to/none' % where)
    if o['y_ticks'] is not None:
        o['y_ticks'] = [_num(v, 'y_ticks') for v in o['y_ticks']]
    return o, species, order, pathways, ann


def fmt_value(v, o, force_brackets=False):
    v = round(v, o['decimals'])
    if v == 0:
        v = 0.0
    s = '%.*f' % (o['decimals'], abs(v))
    s = ('−' if v < 0 else ('+' if o['show_plus'] and v > 0 else '')) + s
    return '(%s)' % s if (o['value_brackets'] or force_brackets) else s


def y_label_markup(o):
    if o['y_label'] is not None:
        return o['y_label']
    sym = o['energy_symbol']
    if '{' not in sym:
        sym = ('Δ{i:%s}' % sym[1:]) if sym.startswith('Δ') and len(sym) > 1 else '{i:%s}' % sym
    unit = o['unit_markup'] or (UNITS.get(str(o['unit']).lower(), o['unit']) if o['unit'] else None)
    if not unit:
        return sym
    return '%s / (%s)' % (sym, unit) if o['y_label_style'] == 'slash' else '%s (%s)' % (sym, unit)


# ---------------------------------------------------------------- 版面
class Frame:
    def __init__(self, o, species):
        es = [s['E'] for s in species.values()]
        self.emin, self.emax = min(es), max(es)
        if self.emax - self.emin < 1e-9:
            raise SpecError('所有物种能量相同，画不出势能图')
        self.k = o['energy_height'] / (self.emax - self.emin)
        self.u = o['unit_width']
        self.xmin = min(s['x'] for s in species.values())

    def X(self, x):
        return (x - self.xmin) * self.u

    def Y(self, e):
        return (self.emax - e) * self.k


def smooth_samples(p, species, o):
    """余弦缓动：每段 E(t)=E0+(E1-E0)(1-cos πt)/2，节点处斜率为 0，节点本身是采样点。"""
    xs = [species[n]['x'] for n in p['nodes']]
    es = [species[n]['E'] for n in p['nodes']]
    n = o['samples_per_segment']
    tail_n = 8 if o['tail'] > 0 else 0
    while (len(xs) - 1) * n + 2 * tail_n + 1 > ep_plot.MAX_SAMPLES:
        n -= 5
    pts = [(xs[0] - o['tail'] * (1 - k / tail_n), es[0]) for k in range(tail_n)]
    for i in range(len(xs) - 1):
        for k in range(n):
            t = k / n
            pts.append((xs[i] + (xs[i + 1] - xs[i]) * t, es[i] + (es[i + 1] - es[i]) * (1 - math.cos(math.pi * t)) / 2))
    pts.append((xs[-1], es[-1]))
    pts += [(xs[-1] + o['tail'] * k / tail_n, es[-1]) for k in range(1, tail_n + 1)]
    return pts


def poly_range(polys, xa, xb):
    """多条 x 单调折线在 [xa, xb] 内的 (ymin, ymax)，画布坐标（y 向下）。"""
    ys = []
    for pl in polys:
        for (x1, y1), (x2, y2) in zip(pl, pl[1:]):
            lo, hi = max(min(x1, x2), xa), min(max(x1, x2), xb)
            if lo > hi:
                continue
            for x in (lo, hi):
                ys.append(y1 if x2 == x1 else y1 + (y2 - y1) * (x - x1) / (x2 - x1))
    return (min(ys), max(ys)) if ys else None


def text_el(x, y, runs, color, role, owner):
    w, a, d = T.metrics(runs)
    return {'x': x, 'y': y, 'runs': runs, 'color': color, 'role': role, 'owner': owner, 'w': w, 'asc': a, 'desc': d}


def stack(lines, cx, side, env_y, half, color, role, owner, extra_dy=0.0):
    """lines 自上而下；side=above 时最后一行贴线，below 时第一行贴线。"""
    out = []
    if side == 'above':
        y = env_y - half - GAP - T.metrics(lines[-1])[2] - extra_dy
        for runs in reversed(lines):
            out.append(text_el(cx, y, runs, color, role, owner))
            y -= T.LINE_PITCH
    else:
        y = env_y + half + GAP + T.metrics(lines[0])[1] - extra_dy
        for runs in lines:
            out.append(text_el(cx, y, runs, color, role, owner))
            y += T.LINE_PITCH
    return out


def auto_side(sid, species, pathways):
    s = species[sid]
    if s['ts']:
        return 'above'
    for p in pathways:
        if sid in p['nodes']:
            i = p['nodes'].index(sid)
            nb = [species[p['nodes'][j]]['E'] for j in (i - 1, i + 1) if 0 <= j < len(p['nodes'])]
            if all(e > s['E'] for e in nb):
                return 'below'
            if all(e < s['E'] for e in nb):
                return 'above'
            return 'below'
    return 'below'


def layout(o, species, order, pathways, ann, F, curve_segs, curve_attrs):
    L = {'lines': [], 'arrows': [], 'curves': [], 'texts': [], 'structures': [], 'warnings': []}
    polys = {}
    bw = o['bar_width']
    for p in pathways:
        if o['style'] == 'smooth':
            segs = [[(F.X(x), F.Y(e)) for x, e in s] for s in curve_segs[p['id']]]
            L['curves'].append({'segs': segs, 'attrs': curve_attrs[p['id']], 'color': p['color'],
                                'width': p['width'], 'dashed': p['dashed'], 'owner': p['id']})
            pl = list(segs[0])
            for s in segs[1:]:
                pl.extend(s[1:])
            polys[p['id']] = pl
        else:
            pl = []
            for a, b in zip(p['nodes'], p['nodes'][1:]):
                p1 = (F.X(species[a]['x']) + bw / 2, F.Y(species[a]['E']))
                p2 = (F.X(species[b]['x']) - bw / 2, F.Y(species[b]['E']))
                L['lines'].append({'p1': p1, 'p2': p2, 'style': 'dashed' if p['dashed'] else 'solid',
                                   'color': p['color'], 'role': 'connector', 'owner': p['id']})
            for n in p['nodes']:
                X, Y = F.X(species[n]['x']), F.Y(species[n]['E'])
                pl += [(X - bw / 2, Y), (X + bw / 2, Y)]
            polys[p['id']] = pl
    if o['bars']:
        seen = set()
        for p in pathways:
            for n in p['nodes']:
                if n in seen:
                    continue
                seen.add(n)
                X, Y = F.X(species[n]['x']), F.Y(species[n]['E'])
                L['lines'].append({'p1': (X - bw / 2, Y), 'p2': (X + bw / 2, Y), 'style': 'bold',
                                   'color': p['color'], 'role': 'bar', 'owner': n})

    def half_at(sid):
        ws = [p['width'] for p in pathways if sid in p['nodes']] if o['style'] == 'smooth' else [LW]
        return max(max(ws) / 2, BOLD_W / 2 if o['bars'] else 0)

    face = {'bold': 1, 'plain': 0, 'italic': 2}[o['label_face']]
    for sid in order:
        s = species[sid]
        if s['pos'] == 'none':
            continue
        name = T.parse(s['label'], face)
        if s['ts'] and o['ts_marker'] == 'dagger':
            name.append(('‡', face | 64))
        val = None
        if s['show_value']:
            val = T.parse(s['value']) if s['value'] is not None else \
                [(fmt_value(s['E'], o, o['label_layout'] == 'inline'), 0)]
        side = auto_side(sid, species, pathways) if s['pos'] == 'auto' else s['pos']
        cx, Y = F.X(s['x']) + s['dx'], F.Y(s['E'])
        mine = [polys[p['id']] for p in pathways if sid in p['nodes']]
        groups = []
        if val is None:
            groups.append((side, [name]))
        elif o['label_layout'] == 'stack':
            groups.append((side, [name, val]))
        elif o['label_layout'] == 'inline':
            groups.append((side, [name + [(' ', 0)] + val]))
        else:
            groups += [(side, [name]), ('below' if side == 'above' else 'above', [val])]
        for gside, lines in groups:
            w = max(T.metrics(r)[0] for r in lines)
            rng = poly_range(mine, cx - w / 2 - 1, cx + w / 2 + 1) or (Y, Y)
            env = rng[0] if gside == 'above' else rng[1]
            L['texts'] += stack(lines, cx, gside, env, half_at(sid), s['color'], 'species', sid, s['dy'])

        if s.get('structure'):
            sspec = s['structure']
            if not isinstance(sspec, dict):
                raise SpecError('species %r 的 structure 必须是对象' % sid)
            try:
                item = ep_structure.make(sspec, len(L['structures']))
            except ep_structure.StructureError as exc:
                raise SpecError('species %r 的 structure: %s' % (sid, exc))
            extra_dy = _num(sspec.get('dy', 0), 'structure.dy')
            extra_dx = _num(sspec.get('dx', 0), 'structure.dx')
            scx = cx + extra_dx
            # clear the curve as well as the point: a structure sitting on a
            # rising flank collides with it, so take the topmost curve y over
            # the structure's own x span
            inside = poly_range(mine, scx - item['half_w'], scx + item['half_w']) or (Y, Y)
            top = min(Y, inside[0])
            L['structures'].append({
                'owner': sid, 'item': item,
                'cx': scx,
                'cy': top - STRUCTURE_GAP - item['half_h'] + extra_dy,
                'half_w': item['half_w'], 'half_h': item['half_h'],
            })

    for p in pathways:
        if not p['label']:
            continue
        runs = T.parse(p['label'])
        w = T.metrics(runs)[0]
        lx = p['label_x'] if p['label_x'] is not None else \
            (species[p['nodes'][0]]['x'] + species[p['nodes'][1]]['x']) / 2
        cx = F.X(_num(lx, 'label_x'))
        rng = poly_range([polys[p['id']]], cx - w / 2 - 1, cx + w / 2 + 1)
        if rng is None:
            raise SpecError('路径 %s 的 label_x 不在曲线范围内' % p['id'])
        env = rng[0] if p['label_position'] == 'above' else rng[1]
        L['texts'] += stack([runs], cx, p['label_position'], env, p['width'] / 2, p['color'], 'pathway',
                            p['id'], p['label_dy'])

    for a in ann:
        sf, st = species[a['from']], species[a['to']]
        X = F.X(species[a['x']]['x'] if isinstance(a['x'], str) else _num(a['x'], 'annotation x'))
        ends = []
        for which, s in (('from', sf), ('to', st)):
            Xs, Ys = F.X(s['x']), F.Y(s['E'])
            on_point = abs(X - Xs) < 0.5
            if not on_point and a['ref'] in ('both', which) and abs(X - Xs) > 4:
                x0 = Xs + (bw / 2 if o['bars'] else 0) * (1 if X > Xs else -1)
                L['lines'].append({'p1': (x0, Ys), 'p2': (X, Ys), 'style': 'dashed', 'color': a['color'],
                                   'role': 'ref', 'owner': 'ann:%s-%s' % (a['from'], a['to'])})
            ends.append((Ys, half_at(s['id']) if on_point else LW / 2))
        (y1, i1), (y2, i2) = ends
        sgn = 1 if y2 > y1 else -1
        t1, t2 = (X, y1 + sgn * i1), (X, y2 - sgn * i2)
        length = abs(t2[1] - t1[1])
        heads = a['heads']
        need = {'both': 2, 'to': 1, 'from': 1, 'none': 0}[heads] * 6.0 * HEAD_ANN / 1000 + 2
        owner = 'ann:%s-%s' % (a['from'], a['to'])
        if heads != 'none' and length < need:
            L['warnings'].append('能量差箭头 %s 只有 %.1f pt，放不下箭头，改画成无箭头线段' % (owner, length))
            heads = 'none'
        if heads == 'from':
            t1, t2 = t2, t1
        L['arrows'].append({'tail': t1, 'head': t2, 'heads': heads, 'size': HEAD_ANN, 'color': a['color'],
                            'role': 'gap', 'owner': owner})
        lab = a['label'].replace('{value}', fmt_value(st['E'] - sf['E'], dict(o, value_brackets=False)))
        if lab:
            runs = T.parse(lab)
            w, asc, desc = T.metrics(runs)
            ym = y1 + (y2 - y1) * a['at']
            cx = X + (1 if a['side'] == 'right' else -1) * (GAP + w / 2) + a['dx']
            L['texts'].append(text_el(cx, ym + (asc - desc) / 2 - a['dy'], runs, a['color'], 'gap_label', owner))
    return L, polys


def bbox_of(L, with_text=True):
    xs, ys = [], []
    for l in L['lines']:
        w = BOLD_W if l['style'] == 'bold' else LW
        for p in (l['p1'], l['p2']):
            xs += [p[0] - w, p[0] + w]
            ys += [p[1] - w, p[1] + w]
    for a in L['arrows']:
        for p in (a['tail'], a['head']):
            xs += [p[0] - 3, p[0] + 3]
            ys += [p[1] - 3, p[1] + 3]
    for c in L['curves']:
        hw = c['width'] / 2
        for s in c['segs']:
            for x, y in s:
                xs += [x - hw, x + hw]
                ys += [y - hw, y + hw]
    for s in L.get('structures', ()):
        xs += [s['cx'] - s['half_w'], s['cx'] + s['half_w']]
        ys += [s['cy'] - s['half_h'], s['cy'] + s['half_h']]
    if with_text:
        for t in L['texts']:
            xs += [t['x'] - t['w'] / 2, t['x'] + t['w'] / 2]
            ys += [t['y'] - t['asc'], t['y'] + t['desc']]
    return min(xs), min(ys), max(xs), max(ys)


def add_axes(L, o, F, species=None):
    bx0, by0, bx1, by1 = bbox_of(L)
    origin = (bx0 - 6, by1 + 6)
    ytip, xtip = (origin[0], by0 - 8), (bx1 + 8, origin[1])
    L['arrows'] += [{'tail': origin, 'head': ytip, 'heads': 'to', 'size': HEAD_AXIS, 'color': '#000000',
                     'role': 'axis', 'owner': 'y'},
                    {'tail': origin, 'head': xtip, 'heads': 'to', 'size': HEAD_AXIS, 'color': '#000000',
                     'role': 'axis', 'owner': 'x'}]
    if o.get('x_ticks') and species:
        for sid in sorted(species, key=lambda k: species[k]['x']):
            x = F.X(species[sid]['x'])
            if not (origin[0] + 4 < x < xtip[0] - 4):
                continue
            # inward: the plot sits above the x axis (screen y grows downward)
            L['lines'].append({'p1': (x, origin[1]), 'p2': (x, origin[1] - 3), 'style': 'solid',
                               'color': '#000000', 'role': 'tick', 'owner': sid})
    if o['y_ticks']:
        for v in o['y_ticks']:
            y = F.Y(v)
            if not (ytip[1] + 8 < y < origin[1]):
                L['warnings'].append('刻度 %g 超出纵轴范围，已跳过' % v)
                continue
            # inward, matching x_ticks: the plot lies to the right of the y axis
            L['lines'].append({'p1': (origin[0], y), 'p2': (origin[0] + 3, y), 'style': 'solid',
                               'color': '#000000', 'role': 'tick', 'owner': 'y'})
            runs = [(fmt_value(v, dict(o, decimals=o['y_tick_decimals'], value_brackets=False, show_plus=False)), 0)]
            w, asc, desc = T.metrics(runs)
            L['texts'].append(text_el(origin[0] - 5 - w / 2, y + (asc - desc) / 2, runs, None, 'tick', 'y'))
    yl = T.parse(y_label_markup(o))
    _, a, d = T.metrics(yl)
    L['texts'].append(text_el(origin[0], ytip[1] - 4 - d, yl, None, 'axis_label', 'y'))
    if o['x_label']:
        xl = T.parse(o['x_label'])
        _, a, d = T.metrics(xl)
        L['texts'].append(text_el((origin[0] + xtip[0]) / 2, origin[1] + 5 + a, xl, None, 'axis_label', 'x'))
    return {'origin': origin, 'xtip': xtip, 'ytip': ytip}


def check_overlaps(L):
    warn = []

    def rect(t, s=0.5):
        return (t['x'] - t['w'] / 2 + s, t['y'] - t['asc'] + s, t['x'] + t['w'] / 2 - s, t['y'] + t['desc'] - s)

    def name(t):
        return '%s「%s」' % (t['role'], T.plain(t['runs']))
    ts = L['texts']
    for i in range(len(ts)):
        for j in range(i + 1, len(ts)):
            a, b = rect(ts[i]), rect(ts[j])
            if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]:
                warn.append('文字重叠：%s 与 %s' % (name(ts[i]), name(ts[j])))
    shapes = [('%s:%s' % (l['role'], l['owner']), ep_verify.seg_points(l['p1'], l['p2'])) for l in L['lines']]
    shapes += [('%s:%s' % (a['role'], a['owner']), ep_verify.seg_points(a['tail'], a['head'])) for a in L['arrows']]
    for c in L['curves']:
        pts = []
        for s in c['segs']:
            for p, q in zip(s, s[1:]):
                pts += ep_verify.seg_points(p, q, 1.0)
        shapes.append(('curve:%s' % c['owner'], pts))
    for t in ts:
        r = rect(t)
        for nm, pts in shapes:
            if any(r[0] < x < r[2] and r[1] < y < r[3] for x, y in pts):
                warn.append('文字压线：%s 压到 %s' % (name(t), nm))
    return warn


def shift_all(L, axes, dx, dy):
    mv = lambda p: (p[0] + dx, p[1] + dy)  # noqa: E731
    for l in L['lines']:
        l['p1'], l['p2'] = mv(l['p1']), mv(l['p2'])
    for a in L['arrows']:
        a['tail'], a['head'] = mv(a['tail']), mv(a['head'])
    for c in L['curves']:
        c['segs'] = [[mv(p) for p in s] for s in c['segs']]
    for t in L['texts']:
        t['x'], t['y'] = t['x'] + dx, t['y'] + dy
    for s in L.get('structures', ()):
        s['cx'], s['cy'] = s['cx'] + dx, s['cy'] + dy
    for k in axes:
        axes[k] = mv(axes[k])


# ---------------------------------------------------------------- CDXML
class Colors:
    def __init__(self):
        self.table = [(255, 255, 255), (0, 0, 0)]

    def idx(self, h):
        rgb = tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))
        if rgb not in self.table:
            self.table.append(rgb)
        return str(self.table.index(rgb) + 2)

    def xml(self):
        f = lambda c: '0' if c <= 0 else ('1' if c >= 255 else ('%.4f' % (c / 255)).rstrip('0').rstrip('.'))  # noqa
        return '<colortable>%s</colortable>' % ''.join('<color r="%s" g="%s" b="%s"/>' % tuple(f(c) for c in rgb)
                                                     for rgb in self.table)


def serialize(L, W, H):
    C = Colors()
    out, nid = [], [1000]

    def new():
        nid[0] += 1
        return nid[0]
    for c in L['curves']:
        for s in c['segs']:
            i = new()
            pts = ' '.join('%.3f %.3f %.3f %.3f %.3f %.3f' % (x, y, x, y, x, y) for x, y in s)
            a = c['attrs']
            out.append('<curve id="%d" Z="%d" CurvePoints="%s" Closed="%s" FillType="%s" CurveType="%s" color="%s" '
                       'LineWidth="%g"%s/>' % (i, i, pts, a['Closed'], a['FillType'], a['CurveType'],
                                               C.idx(c['color']), c['width'],
                                               ' LineType="Dashed"' if c['dashed'] else ''))
    for l in sorted(L['lines'], key=lambda l: {'ref': 0, 'tick': 1, 'connector': 2, 'bar': 3}[l['role']]):
        i = new()
        (x1, y1), (x2, y2) = l['p1'], l['p2']
        lt = {'solid': '', 'dashed': ' LineType="Dashed"', 'bold': ' LineType="Bold"'}[l['style']]
        out.append('<graphic id="%d" Z="%d" GraphicType="Line" BoundingBox="%.3f %.3f %.3f %.3f"%s color="%s"/>'
                   % (i, i, x1, y1, x2, y2, lt, C.idx(l['color'])))
    for a in L['arrows']:
        i = new()
        (tx, ty), (hx, hy) = a['tail'], a['head']
        if a['heads'] == 'none':
            out.append('<graphic id="%d" Z="%d" GraphicType="Line" BoundingBox="%.3f %.3f %.3f %.3f" color="%s"/>'
                       % (i, i, tx, ty, hx, hy, C.idx(a['color'])))
            continue
        tail = ' ArrowheadTail="Full"' if a['heads'] == 'both' else ''
        h = a['size']
        out.append('<arrow id="%d" Z="%d" FillType="None" ArrowheadHead="Full"%s ArrowheadType="Solid" HeadSize="%d" '
                   'ArrowheadCenterSize="%d" ArrowheadWidth="%d" Head3D="%.3f %.3f 0" Tail3D="%.3f %.3f 0" color="%s"/>'
                   % (i, i, tail, h, round(h * 0.875), round(h * 0.25), hx, hy, tx, ty, C.idx(a['color'])))
    for t in L['texts']:
        i = new()
        ci = C.idx(t['color']) if t['color'] and t['color'] != '#000000' else None
        c_attr = ' color="%s"' % ci if ci else ''
        out.append('<t id="%d" Z="%d" p="%.3f %.3f" Justification="Center" CaptionJustification="Center"%s '
                   'InterpretChemically="no">%s</t>' % (i, i, t['x'], t['y'], c_attr, T.runs_xml(t['runs'], ci)))
    for s in L.get('structures', ()):
        out.append(ep_structure.emit(s['item'], s['cx'], s['cy']))
    wp, hp = max(2, math.ceil(W / 500) + 1), max(2, math.ceil(H / 650) + 1)
    head = ('<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">\n'
            '<CDXML BoundingBox="%.2f %.2f %.2f %.2f" %s>' % (MARGIN, MARGIN, W - MARGIN, H - MARGIN, ROOT_ATTRS)
            + C.xml() + FONTTABLE +
            '<page id="1" BoundingBox="0 0 %.0f %.0f" Width="%.0f" Height="%.0f" WidthPages="%d" HeightPages="%d" '
            'HeaderPosition="0" FooterPosition="0">\n' % (W, H, W, H, wp, hp))
    return head + '\n'.join(out) + '\n</page></CDXML>\n', C


# ---------------------------------------------------------------- 主流程
def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def guard(path, prev, force):
    """只覆盖本脚本上次写的、之后没被改过的文件（或 --force）。"""
    if not os.path.exists(path) or force:
        return
    if path.endswith('.report.json') and str(prev.get('tool', '')).startswith('chemdraw-energy-profile'):
        return
    if prev.get('sha256', {}).get(os.path.basename(path)) == sha256(path):
        return
    raise SpecError('%s 已存在且不是本脚本上次的产物（或已被手改）；换输出名，或加 --force' % path)


def render(cdxml, png, cdx_cmd, env):
    tmp = tempfile.mkdtemp(prefix='ep_render_')
    try:
        staged = os.path.join(tmp, 'ep_%d_%d.cdxml' % (os.getpid(), int(time.time() * 1000)))
        shutil.copyfile(cdxml, staged)
        p = subprocess.run(cdx_cmd + ['render', staged], capture_output=True, env=env, timeout=600)
        out = os.path.splitext(staged)[0] + '.png'
        if p.returncode != 0 or not os.path.isfile(out):
            return False, (p.stdout + p.stderr).decode('utf-8', 'replace')[-800:]
        shutil.copyfile(out, png)
        return True, ''
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def build(spec, cdx_cmd, env, keep_work=None):
    """规格 → 版面 + CDXML 文本（不写盘、不渲染）。返回 dict，供 main 和测试共用。"""
    o, species, order, pathways, ann = normalize(spec)
    F = Frame(o, species)
    segs, attrs, plot_info = {}, {}, None
    if o['style'] == 'smooth':
        work_parent = None
        if keep_work:
            work = os.path.abspath(keep_work)
        else:
            work_parent = tempfile.mkdtemp(prefix='ep_plot_')
            work = os.path.join(work_parent, 'work')
        try:
            series = [(p['id'], smooth_samples(p, species, o)) for p in pathways]
            segs, attrs, plot_info = ep_plot.run_plot(series, work, cdx_cmd, env)
        finally:
            if work_parent:
                shutil.rmtree(work_parent, ignore_errors=True)
                if plot_info:
                    plot_info.pop('plot_cdxml', None)
    L, polys = layout(o, species, order, pathways, ann, F, segs, attrs)
    axes = add_axes(L, o, F, species)
    L['warnings'] += check_overlaps(L)
    bx0, by0, bx1, by1 = bbox_of(L)
    dx, dy = MARGIN - bx0, MARGIN - by0
    shift_all(L, axes, dx, dy)
    W, H = bx1 - bx0 + 2 * MARGIN, by1 - by0 + 2 * MARGIN
    xml, _ = serialize(L, W, H)
    nodes = {sid: (F.X(s['x']) + dx, F.Y(s['E']) + dy) for sid, s in species.items()}
    structure_texts = [t for s in L.get('structures', ()) for t in s['item'].get('texts', ())]
    expect = {'texts': [T.plain(t['runs']) for t in L['texts']] + structure_texts,
              'style': o['style'], 'bars': o['bars'],
              'node_points': [nodes[n] for p in pathways for n in p['nodes']],
              'curve_nodes': sum(len(c['segs']) for c in L['curves']), 'page': (W, H)}
    return {'L': L, 'axes': axes, 'F': F, 'xml': xml, 'nodes': nodes, 'expect': expect, 'dx': dx, 'dy': dy,
            'plot_info': plot_info, 'o': o}


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding='utf-8')
        except Exception:
            pass
    ap = argparse.ArgumentParser(description='势能图 JSON → CDXML + 原生 PNG')
    ap.add_argument('spec')
    ap.add_argument('-o', '--output', help='输出 .cdxml（默认与 spec 同名同目录）')
    ap.add_argument('--no-render', action='store_true')
    ap.add_argument('--force', action='store_true', help='允许覆盖不是本脚本写的同名文件')
    ap.add_argument('--keep-work', help='保留 plot 中间产物到这个新目录（必须不存在）')
    ap.add_argument('--skill', help='chemdraw skill 目录（默认 ~/.claude/skills/chemdraw）')
    args = ap.parse_args(argv)
    try:
        with open(args.spec, encoding='utf-8') as fh:
            spec = json.load(fh)
        out = os.path.abspath(args.output or os.path.splitext(args.spec)[0] + '.cdxml')
        stem = os.path.splitext(out)[0]
        png, rep_path = stem + '.png', stem + '.report.json'
        prev = {}
        if os.path.isfile(rep_path):
            with open(rep_path, encoding='utf-8') as fh:
                prev = json.load(fh)
        for pth in (out, rep_path) + (() if args.no_render else (png,)):
            guard(pth, prev, args.force)
        skill = args.skill or next((d for d in SKILL_DIRS if os.path.isfile(os.path.join(d, 'scripts', 'cdx.py'))), None)
        if not skill:
            raise SpecError('找不到 chemdraw skill（scripts/cdx.py）；用 --skill 指定')
        cdx_cmd = [CDX_PY, os.path.join(skill, 'scripts', 'cdx.py')]
        env = dict(os.environ, PYTHONIOENCODING='utf-8')
        M = build(spec, cdx_cmd, env, args.keep_work)
    except (SpecError, ValueError, OSError, RuntimeError, json.JSONDecodeError) as exc:
        print('错误：%s' % exc)
        return 2
    L, axes, F = M['L'], M['axes'], M['F']
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write(M['xml'])
    report = {'tool': 'chemdraw-energy-profile %s' % VERSION, 'spec': os.path.abspath(args.spec),
              'outputs': {'cdxml': out, 'png': None if args.no_render else png},
              'mapping': {'pt_per_unit_x': F.u, 'pt_per_energy': F.k, 'energy_range': [F.emin, F.emax],
                          'note': 'X = %.3f + %.4f*(x - %g)；Y = %.3f + %.6f*(%g - E)（pt，y 向下）'
                                  % (M['dx'], F.u, F.xmin, M['dy'], F.k, F.emax)},
              'species_points_pt': {k: [round(v[0], 3), round(v[1], 3)] for k, v in M['nodes'].items()},
              'plot': M['plot_info'], 'layout_warnings': L['warnings']}
    report['selfcheck'] = ep_verify.selfcheck(out, M['expect'])
    report['check_cdxml'] = ep_verify.run_checker(out, CHECKER)
    if not args.no_render:
        ok, err = render(out, png, cdx_cmd, env)
        report['render'] = {'ok': ok, 'error': err or None}
        if ok:
            report['png_check'] = ep_verify.verify_png(png, L, axes)
    report['sha256'] = {os.path.basename(p): sha256(p) for p in (out, png) if os.path.isfile(p)
                        and not (p == png and args.no_render)}
    fails = [k for k in ('selfcheck', 'check_cdxml', 'png_check') if report.get(k) and not report[k].get('ok')]
    if not args.no_render and not report['render']['ok']:
        fails.append('render')
    report['ok'] = not fails
    with open(rep_path, 'w', encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    ep_verify.print_summary(report, rep_path)
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
