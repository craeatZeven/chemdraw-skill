# -*- coding: utf-8 -*-
"""曲线：用 chemdraw skill 的 scientific_workflow(operation=plot) 生成原生 <curve> 节点，再映射到版面。

为什么不自己写 <curve>：chemdraw-workbench 陷阱 1（自算点写的 <curve> 会被 ChemDraw 原生渲染整条丢弃）。
plot 工具的写法（cdxml_toolkit.scientific.plots / common.Drawing.curve，2026-11-09 读源码确认）：
  * 每个样本点写成退化三次贝塞尔的三元组 (p, p, p)，即可编辑折线；
  * 每 251 个点切一段 <curve>（相邻段共用 1 个端点），避开 ChemDraw 单对象控制点上限；
  * 属性 Closed="no" FillType="None" CurveType="0" color=<色表号>。
这里只做仿射坐标变换，点数、分段、属性全部保留；映射前把每个锚点反解回数据值，
与输入样本逐点比对（相对误差 < 1e-6），对不上就报错，不会悄悄画错。
"""
import json
import os
import re
import subprocess
import xml.etree.ElementTree as ET

MAX_SAMPLES = 1500          # plot 在 >1600 点时会做 min/max 抽稀，逐点比对就失效；留余量


def _tag_colors(n):
    """给每条路径一个独一无二的临时颜色，用来在 plot.cdxml 里认出哪段 <curve> 属于哪条路径。"""
    return ['#%02x%02x%02x' % (1 + i, 2 + i, 3 + i) for i in range(n)]


def _hex_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def run_plot(series, work_dir, cdx_cmd, env):
    """series: [(id, [(x, E), ...]), ...]（数据单位）。返回 {id: [[(x, E) 锚点...] 每段一个列表]} 与元信息。"""
    os.makedirs(work_dir, exist_ok=False)
    tags = _tag_colors(len(series))
    entries = []
    for (sid, pts), tag in zip(series, tags):
        if len(pts) > MAX_SAMPLES:
            raise ValueError('路径 %s 的采样点 %d 超过 %d' % (sid, len(pts), MAX_SAMPLES))
        name = 'series_%d.csv' % len(entries)
        with open(os.path.join(work_dir, name), 'w', encoding='utf-8', newline='') as fh:
            fh.write('x,E\n')
            for x, e in pts:
                fh.write('%.10g,%.10g\n' % (x, e))
        entries.append({'input': name, 'x_column': 'x', 'y_column': 'E', 'color': tag})
    spec_path = os.path.join(work_dir, 'plot_spec.json')
    with open(spec_path, 'w', encoding='utf-8') as fh:
        json.dump({'series': entries, 'x_label': 'x', 'y_label': 'E', 'title': ''}, fh)
    out_dir = os.path.join(work_dir, 'plot_out')
    arg = json.dumps({'operation': 'plot', 'spec_path': spec_path, 'output_path': out_dir})
    p = subprocess.run(cdx_cmd + ['call', 'scientific_workflow', '--json', arg],
                       capture_output=True, env=env, timeout=300)
    so = p.stdout.decode('utf-8', 'replace')
    cdxml = os.path.join(out_dir, 'plot.cdxml')
    if p.returncode != 0 or not os.path.isfile(cdxml):
        raise RuntimeError('scientific_workflow plot 失败（退出码 %d）：%s %s'
                           % (p.returncode, so[-600:], p.stderr.decode('utf-8', 'replace')[-600:]))
    with open(os.path.join(out_dir, 'plot.json'), encoding='utf-8') as fh:
        meta = json.load(fh)
    (x0, x1), (y0, y1) = meta['axes']['x_limits'], meta['axes']['y_limits']

    root = ET.parse(cdxml).getroot()
    # 色表：<color> 按出现顺序编号，从 2 起（2=白，3=黑，4=第一个自定义色）
    rgb2idx = {}
    for i, c in enumerate(root.find('colortable')):
        rgb = tuple(int(round(float(c.get(k)) * 255)) for k in ('r', 'g', 'b'))
        rgb2idx[rgb] = str(i + 2)
    idx2sid = {}
    for (sid, _), tag in zip(series, tags):
        idx = rgb2idx.get(_hex_rgb(tag))
        if idx is None or idx in idx2sid:
            raise RuntimeError('plot.cdxml 色表里找不到路径 %s 的标记色 %s' % (sid, tag))
        idx2sid[idx] = sid

    # 坐标轴：最长的竖线 = 左轴（给出 left/top/bottom），最长的横线 = 底轴（给出 right）
    lines = []
    for g in root.iter('graphic'):
        if g.get('GraphicType') == 'Line':
            a = [float(v) for v in g.get('BoundingBox').split()]
            lines.append(a)
    vert = max((l for l in lines if abs(l[0] - l[2]) < 1e-6), key=lambda l: abs(l[3] - l[1]))
    hori = max((l for l in lines if abs(l[1] - l[3]) < 1e-6), key=lambda l: abs(l[2] - l[0]))
    left, top, bottom = vert[0], min(vert[1], vert[3]), max(vert[1], vert[3])
    right = max(hori[0], hori[2])

    def to_data(px, py):
        return (x0 + (px - left) / (right - left) * (x1 - x0),
                y0 + (bottom - py) / (bottom - top) * (y1 - y0))

    segs = {sid: [] for sid, _ in series}
    attrs = {}
    tol = 1e-6 * max(right - left, bottom - top)
    for c in root.iter('curve'):
        sid = idx2sid.get(c.get('color'))
        if sid is None:
            continue
        v = [float(t) for t in c.get('CurvePoints').split()]
        pts = list(zip(v[0::2], v[1::2]))
        if len(pts) % 3:
            raise RuntimeError('plot 的 <curve> 点数不是 3 的倍数（%d），写法变了，停下检查' % len(pts))
        anchors = []
        for k in range(0, len(pts), 3):
            a, b, cc = pts[k:k + 3]
            if max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[0] - cc[0]), abs(a[1] - cc[1])) > tol:
                raise RuntimeError('plot 的 <curve> 不再是退化三元组（折线），写法变了，停下检查')
            anchors.append(to_data(*a))
        segs[sid].append(anchors)
        attrs[sid] = {k: c.get(k) for k in ('Closed', 'FillType', 'CurveType')}

    # 逐点比对：各段去掉共用端点后应正好等于输入样本
    check = {}
    for sid, pts in series:
        if not segs[sid]:
            raise RuntimeError('plot.cdxml 里没有路径 %s 的 <curve>' % sid)
        joined = list(segs[sid][0])
        for s in segs[sid][1:]:
            joined.extend(s[1:])
        if len(joined) != len(pts):
            raise RuntimeError('路径 %s：plot 输出 %d 个锚点，输入 %d 个（被抽稀或改写了）'
                               % (sid, len(joined), len(pts)))
        span_x = (x1 - x0) or 1.0
        span_y = (y1 - y0) or 1.0
        worst = max(max(abs(a[0] - b[0]) / span_x, abs(a[1] - b[1]) / span_y)
                    for a, b in zip(joined, pts))
        if worst > 1e-6:
            raise RuntimeError('路径 %s：反解出的点与输入样本相差 %.2e（相对），坐标映射不可信' % (sid, worst))
        check[sid] = {'curve_nodes': len(segs[sid]), 'anchors': len(joined),
                      'nodes_points': [len(s) for s in segs[sid]], 'max_rel_err': worst}
    return segs, attrs, {'plot_axes': meta['axes'], 'per_pathway': check,
                         'plot_cdxml': cdxml, 'plot_axis_box_px': [left, top, right, bottom]}
