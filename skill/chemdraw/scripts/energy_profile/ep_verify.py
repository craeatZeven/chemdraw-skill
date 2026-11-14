# -*- coding: utf-8 -*-
"""检查：读回 CDXML 自检、跑 check_cdxml.py、对 ChemDraw 原生渲染的 PNG 做像素统计。

PNG 不进任何人的上下文：这里只算数（覆盖率、颜色中位数、多余墨迹），结论写进 report.json。
像素核对能抓到的事：曲线/线段被 ChemDraw 静默丢弃（覆盖率≈0）、颜色不对、
文字折行（折下来的那行落在预期框外，算"多余墨迹"）、文字没画出来。
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import ep_text as T

HEADSIZES = {800, 1000, 1500, 2000, 2250}


def seg_points(p, q, step=0.5):
    n = max(1, int(math.hypot(q[0] - p[0], q[1] - p[1]) / step))
    return [(p[0] + (q[0] - p[0]) * k / n, p[1] + (q[1] - p[1]) * k / n) for k in range(n + 1)]


def _item(items, name, ok, detail):
    items.append({'name': name, 'ok': bool(ok), 'detail': detail})


# ---------------------------------------------------------------- 读回 CDXML
def selfcheck(path, expect):
    items = []
    root = ET.parse(path).getroot()
    want = {'BondLength': 14.4, 'LineWidth': 0.6, 'BoldWidth': 2.0, 'LabelSize': 10, 'CaptionSize': 10,
            'LabelFont': 3, 'CaptionFont': 3, 'ChainAngle': 120}
    bad = ['%s=%s' % (k, root.get(k)) for k, v in want.items() if root.get(k) is None or abs(float(root.get(k)) - v) > 1e-6]
    _item(items, 'root-style', not bad, '根属性 = ACS 1996（键长 14.4、线宽 0.6、标签/说明 Arial 10）' if not bad else bad)
    fonts = {f.get('id'): f.get('name') for f in root.iter('font')}
    _item(items, 'fonttable', fonts.get('3') == 'Arial' and fonts.get('5') == 'SimSun', fonts)

    runs_bad, sizes, n_runs, n_cjk = [], set(), 0, 0
    texts = []
    for t in root.iter('t'):
        s_all = []
        for s in t.iter('s'):
            n_runs += 1
            txt = s.text or ''
            s_all.append(txt)
            sizes.add(s.get('size'))
            cjk = [T.is_cjk(c) for c in txt]
            if all(cjk) and txt:
                n_cjk += 1
                if s.get('font') != '5':
                    runs_bad.append(('中文不是 SimSun', txt))
            elif any(cjk):
                runs_bad.append(('中英混在一个 run', txt))
            elif s.get('font') != '3':
                runs_bad.append(('拉丁不是 Arial', txt))
        texts.append(''.join(s_all))
    _item(items, 'fonts-per-run', not runs_bad,
          '%d 个 run：中文 %d 个全为 SimSun(5)，其余全为 Arial(3)' % (n_runs, n_cjk) if not runs_bad else runs_bad[:6])
    _item(items, 'font-size', sizes == {'10'}, '字号 %s pt' % sorted(sizes))
    _item(items, 'texts', sorted(texts) == sorted(expect['texts']),
          {'n': len(texts), 'texts': texts} if sorted(texts) == sorted(expect['texts'])
          else {'got': sorted(texts), 'want': sorted(expect['texts'])})
    leftovers = [x for x in texts if '{' in x or '}' in x]
    _item(items, 'markup-resolved', not leftovers, leftovers or '无残留标记 {..}')

    ntab = len(list(root.find('colortable')))
    refs = {int(e.get('color')) for e in root.iter() if e.get('color') and e.tag != 'CDXML'}
    _item(items, 'colors', all(2 <= c <= ntab + 1 for c in refs), '色表 %d 条，引用 %s' % (ntab, sorted(refs)))

    page = root.find('page')
    W, H = float(page.get('Width')), float(page.get('Height'))
    coords = []
    curves = list(root.iter('curve'))
    anchors = []
    cbad = []
    for c in curves:
        v = [float(x) for x in c.get('CurvePoints').split()]
        if len(v) % 6 or any(c.get(k) is None for k in ('Closed', 'FillType', 'CurveType')):
            cbad.append(c.get('id'))
        pts = list(zip(v[0::2], v[1::2]))
        anchors += pts[0::3]
        coords += pts
    _item(items, 'curve-nodes', len(curves) == expect['curve_nodes'] and not cbad,
          '%d 个 <curve>（应 %d），属性齐全，点数为 3 的倍数' % (len(curves), expect['curve_nodes'])
          if not cbad else ('属性缺失/点数不对', cbad))
    if expect['style'] == 'smooth':
        worst = 0.0
        for nx, ny in expect['node_points']:
            d = min(math.hypot(nx - x, ny - y) for x, y in anchors)
            worst = max(worst, d)
        _item(items, 'curve-through-species', worst < 0.3,
              '曲线经过每个物种的(x, E)点，最大偏差 %.3f pt' % worst)
    gl = []
    for g in root.iter('graphic'):
        b = [float(x) for x in g.get('BoundingBox').split()]
        coords += [(b[0], b[1]), (b[2], b[3])]
        if g.get('LineType') == 'Bold':
            gl.append(((b[0] + b[2]) / 2, (b[1] + b[3]) / 2))
    if expect['bars']:
        worst = max(min(math.hypot(nx - x, ny - y) for x, y in gl) for nx, ny in expect['node_points'])
        _item(items, 'bars-at-species', worst < 0.3, '短横线中点与物种(x, E)点最大偏差 %.3f pt' % worst)
    abad = []
    for a in root.iter('arrow'):
        for k in ('Head3D', 'Tail3D'):
            x, y, _ = (float(v) for v in a.get(k).split())
            coords.append((x, y))
        h, c, w = (float(a.get(k)) for k in ('HeadSize', 'ArrowheadCenterSize', 'ArrowheadWidth'))
        if int(h) not in HEADSIZES or abs(c - 0.875 * h) > 1 or abs(w - 0.25 * h) > 1:
            abad.append(a.get('id'))
    _item(items, 'arrowheads', not abad, '箭头头尺寸合规（HeadSize∈%s，0.875/0.25 比例）' % sorted(HEADSIZES)
          if not abad else abad)
    for t in root.iter('t'):
        x, y = (float(v) for v in t.get('p').split())
        coords.append((x, y))
    out = [p for p in coords if not (0 <= p[0] <= W and 0 <= p[1] <= H)]
    _item(items, 'inside-page', not out, '全部坐标在页面 %.0f×%.0f pt 内' % (W, H) if not out else out[:5])
    return {'ok': all(i['ok'] for i in items), 'items': items}


# ---------------------------------------------------------------- check_cdxml.py
def run_checker(path, checker):
    # The checker is an optional external script; None means none was found.
    if not checker or not os.path.isfile(checker):
        return {'ok': True, 'skipped': 'error: no ACS-spec checker found (set CHEMDRAW_CHECK_CDXML)'}
    fd, js = tempfile.mkstemp(suffix='.json')
    os.close(fd)
    try:
        p = subprocess.run([sys.executable, checker, path, '--json', js, '--quiet'], capture_output=True,
                           env=dict(os.environ, PYTHONIOENCODING='utf-8'), timeout=120)
        with open(js, encoding='utf-8') as fh:
            data = json.load(fh)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        return {'ok': False, 'error': str(exc)}
    finally:
        os.unlink(js)
    rows = data.get('rows', [])
    return {'ok': p.returncode == 0, 'exit': p.returncode, 'summary': data.get('summary'),
            'fail': [r for r in rows if r['level'] == 'fail'], 'warn': [r for r in rows if r['level'] == 'warn'],
            'rows': rows}


# ---------------------------------------------------------------- PNG 像素统计
def verify_png(png, L, axes):
    try:
        import numpy as np
        from PIL import Image
    except ImportError as exc:
        return {'ok': True, 'skipped': '没有 numpy/PIL：%s' % exc}
    Image.MAX_IMAGE_PIXELS = None
    a = np.asarray(Image.open(png).convert('RGBA')).astype(np.int16)
    ink = (a[..., 3] > 128) & (a[..., :3].min(-1) < 190)
    Hp, Wp = ink.shape
    items = []

    def longest(mask):
        best = np.zeros(mask.shape[0], int)
        start = np.zeros(mask.shape[0], int)
        for r in range(mask.shape[0]):
            row = mask[r]
            if not row.any():
                continue
            d = np.diff(np.concatenate(([0], row.astype(np.int8), [0])))
            s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
            k = int(np.argmax(e - s))
            best[r], start[r] = e[k] - s[k], s[k]
        return best, start

    bx, sx = longest(ink)
    by, sy = longest(ink.T)
    r0, c0 = int(np.argmax(bx)), int(np.argmax(by))
    rows = [r for r in range(max(0, r0 - 20), min(Hp, r0 + 21)) if bx[r] >= bx[r0] - 3]
    cols = [c for c in range(max(0, c0 - 20), min(Wp, c0 + 21)) if by[c] >= by[c0] - 3]
    ox, oy = axes['origin']
    s_x = bx[r0] / (axes['xtip'][0] - ox + 0.3)
    s_y = by[c0] / (oy - axes['ytip'][1] + 0.3)
    s = (s_x + s_y) / 2
    y0px, x0px = float(np.mean(rows)), float(np.mean(cols))
    _item(items, 'calibration', abs(s_x - s_y) / s < 0.02,
          '坐标轴定标：%.3f / %.3f px/pt（两轴应一致），原点像素 (%.1f, %.1f)' % (s_x, s_y, x0px, y0px))

    def P(x, y):
        return x0px + (x - ox) * s, y0px + (y - oy) * s

    rad = max(3, int(round(0.6 * s)))

    def hit(x, y):
        px, py = P(x, y)
        i, j = int(round(py)), int(round(px))
        if not (0 <= i < Hp and 0 <= j < Wp):
            return False
        return bool(ink[max(0, i - rad):i + rad + 1, max(0, j - rad):j + rad + 1].any())

    def color_at(pts, want):
        sel = []
        for x, y in pts[::3]:
            px, py = P(x, y)
            i, j = int(round(py)), int(round(px))
            win = a[max(0, i - rad):i + rad + 1, max(0, j - rad):j + rad + 1, :3].reshape(-1, 3)
            m = ink[max(0, i - rad):i + rad + 1, max(0, j - rad):j + rad + 1].reshape(-1)
            sel.append(win[m])
        px = np.concatenate(sel) if sel else np.zeros((0, 3))
        if len(px) < 20:
            return None, None
        lum = px @ np.array([0.299, 0.587, 0.114])
        dark = px[lum <= np.percentile(lum, 30)]
        med = np.median(dark, 0)
        tgt = np.array([int(want[i:i + 2], 16) for i in (1, 3, 5)])
        return '#%02X%02X%02X' % tuple(int(v) for v in med), int(np.abs(med - tgt).max())

    for c in L['curves']:
        pts = []
        for sg in c['segs']:
            for p, q in zip(sg, sg[1:]):
                pts += [p, ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)]
        cov = sum(hit(*p) for p in pts) / len(pts)
        thr = 0.35 if c['dashed'] else 0.9
        col, diff = color_at(pts, c['color'])
        ok = cov >= thr and (diff is None or diff <= 90)
        _item(items, 'curve:%s' % c['owner'], ok, '覆盖率 %.3f（阈值 %.2f），墨色 %s（应 %s，差 %s）'
              % (cov, thr, col, c['color'], diff))
    groups = {}
    for l in L['lines']:
        groups.setdefault((l['role'], l['style'] == 'dashed'), []).append(l)
    for (role, dashed), ls in groups.items():
        pts = [p for l in ls for p in seg_points(l['p1'], l['p2'], 1.0)]
        cov = sum(hit(*p) for p in pts) / len(pts)
        thr = 0.35 if dashed else 0.9
        _item(items, 'line:%s%s' % (role, ':dashed' if dashed else ''), cov >= thr,
              '%d 条，覆盖率 %.3f（阈值 %.2f）' % (len(ls), cov, thr))
    for ar in L['arrows']:
        pts = seg_points(ar['tail'], ar['head'], 1.0)[1:-1]
        cov = sum(hit(*p) for p in pts) / max(1, len(pts))
        _item(items, 'arrow:%s:%s' % (ar['role'], ar['owner']), cov >= 0.9, '箭杆覆盖率 %.3f' % cov)
    missing = []
    for t in L['texts']:
        x0, y0 = P(t['x'] - t['w'] / 2, t['y'] - t['asc'])
        x1, y1 = P(t['x'] + t['w'] / 2, t['y'] + t['desc'])
        n = int(ink[max(0, int(y0)):int(y1) + 1, max(0, int(x0)):int(x1) + 1].sum())
        if n < 15:
            missing.append(T.plain(t['runs']))
    _item(items, 'texts-drawn', not missing, '%d 段文字都有墨迹' % len(L['texts']) if not missing
          else ('没画出来', missing))

    allowed = np.zeros_like(ink)

    def stamp(x, y, r_pt):
        px, py = P(x, y)
        r = int(math.ceil(r_pt * s))
        i, j = int(round(py)), int(round(px))
        allowed[max(0, i - r):i + r + 1, max(0, j - r):j + r + 1] = True
    for c in L['curves']:
        for sg in c['segs']:
            for p, q in zip(sg, sg[1:]):
                for x, y in seg_points(p, q, 0.5):
                    stamp(x, y, c['width'] / 2 + 1.0)
    for l in L['lines']:
        for x, y in seg_points(l['p1'], l['p2'], 0.5):
            stamp(x, y, (2.0 if l['style'] == 'bold' else 0.6) / 2 + 1.0)
    for ar in L['arrows']:
        for x, y in seg_points(ar['tail'], ar['head'], 0.5):
            stamp(x, y, 3.5)
    for t in L['texts']:
        x0, y0 = P(t['x'] - t['w'] / 2 - 1.5, t['y'] - t['asc'] - 1.5)
        x1, y1 = P(t['x'] + t['w'] / 2 + 1.5, t['y'] + t['desc'] + 1.5)
        allowed[max(0, int(y0)):int(y1) + 1, max(0, int(x0)):int(x1) + 1] = True
    # conformer structures are opaque XML emitted verbatim, so allow ink inside
    # their bounding box instead of segment-tracing every bond
    for st in L.get('structures', ()):
        x0, y0 = P(st['cx'] - st['half_w'] - 1.5, st['cy'] - st['half_h'] - 1.5)
        x1, y1 = P(st['cx'] + st['half_w'] + 1.5, st['cy'] + st['half_h'] + 1.5)
        allowed[max(0, int(y0)):int(y1) + 1, max(0, int(x0)):int(x1) + 1] = True
    extra = ink & ~allowed
    n_extra = int(extra.sum())
    detail = '多余墨迹 %d px（阈值 300；折行文字、错位对象会落在这里）' % n_extra
    if n_extra:
        ys, xs = np.nonzero(extra)
        detail += '，范围 pt x %.1f–%.1f y %.1f–%.1f' % ((xs.min() - x0px) / s + ox, (xs.max() - x0px) / s + ox,
                                                       (ys.min() - y0px) / s + oy, (ys.max() - y0px) / s + oy)
    _item(items, 'unexplained-ink', n_extra <= 300, detail)
    return {'ok': all(i['ok'] for i in items), 'png_px': [Wp, Hp], 'px_per_pt': round(s, 4), 'items': items}


def print_summary(rep, rep_path):
    o = rep['outputs']
    print('CDXML : %s' % o['cdxml'])
    print('PNG   : %s' % (o['png'] if rep.get('render', {}).get('ok') else '（未渲染）'))
    sc = rep['selfcheck']
    print('自检  : %d/%d 通过' % (sum(i['ok'] for i in sc['items']), len(sc['items'])))
    for i in sc['items']:
        if not i['ok']:
            print('   ✗ %s %s' % (i['name'], i['detail']))
    ck = rep['check_cdxml']
    print('check_cdxml.py : %s' % (ck.get('summary') or ck.get('skipped') or ck.get('error')))
    for r in ck.get('fail', []) + ck.get('warn', []):
        print('   [%s] %s %s' % (r['level'], r['rule'], r['detail']))
    if 'render' in rep and not rep['render']['ok']:
        print('渲染失败：%s' % rep['render']['error'])
    pc = rep.get('png_check')
    if pc:
        if pc.get('skipped'):
            print('像素核对：跳过（%s）' % pc['skipped'])
        else:
            print('像素核对：%d/%d 通过（%.2f px/pt）' % (sum(i['ok'] for i in pc['items']), len(pc['items']),
                                                     pc['px_per_pt']))
            for i in pc['items']:
                if not i['ok']:
                    print('   ✗ %s %s' % (i['name'], i['detail']))
    for w in rep['layout_warnings']:
        print('版面预警：%s' % w)
    print('报告  : %s' % rep_path)
    print('结论  : %s' % ('全部通过' if rep['ok'] else '有失败，见上'))
