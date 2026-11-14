# -*- coding: utf-8 -*-
"""Curves: native <curve> nodes from scientific_workflow(operation=plot), mapped onto the layout.

Why not write <curve> here: ChemDraw silently drops a <curve> whose CurvePoints it did not produce itself.
How the plot workflow writes curves (cdxml_toolkit.scientific.plots / common.Drawing.curve):
  * every sample is a degenerate cubic Bezier triple (p, p, p): an editable polyline;
  * a new <curve> every 251 points (neighbours share one end point), under ChemDraw's per-object limit;
  * attributes Closed="no" FillType="None" CurveType="0" color=<colour-table index>.
Only an affine transform is applied; point count, segments and attributes are kept. Every anchor
is mapped back to data and compared with the input samples (relative error < 1e-6); a mismatch is an error.
"""
import json
import os
import re
import subprocess
import xml.etree.ElementTree as ET

MAX_SAMPLES = 1500          # plot thins series above 1600 points (min/max), which breaks the point check


def _tag_colors(n):
    """Give each pathway a unique temporary colour so its <curve> nodes can be found in plot.cdxml."""
    return ['#%02x%02x%02x' % (1 + i, 2 + i, 3 + i) for i in range(n)]


def _hex_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def run_plot(series, work_dir, cdx_cmd, env):
    """series: [(id, [(x, E), ...]), ...] in data units. Returns {id: [[(x, E) anchors...] per segment]} and metadata."""
    os.makedirs(work_dir, exist_ok=False)
    tags = _tag_colors(len(series))
    entries = []
    for (sid, pts), tag in zip(series, tags):
        if len(pts) > MAX_SAMPLES:
            raise ValueError('pathway %s has %d samples, more than %d' % (sid, len(pts), MAX_SAMPLES))
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
        raise RuntimeError('scientific_workflow plot failed (exit code %d): %s %s'
                           % (p.returncode, so[-600:], p.stderr.decode('utf-8', 'replace')[-600:]))
    with open(os.path.join(out_dir, 'plot.json'), encoding='utf-8') as fh:
        meta = json.load(fh)
    (x0, x1), (y0, y1) = meta['axes']['x_limits'], meta['axes']['y_limits']

    root = ET.parse(cdxml).getroot()
    # Colour table: <color> entries are numbered from 2 (2 = white, 3 = black, 4 = first custom colour)
    rgb2idx = {}
    for i, c in enumerate(root.find('colortable')):
        rgb = tuple(int(round(float(c.get(k)) * 255)) for k in ('r', 'g', 'b'))
        rgb2idx[rgb] = str(i + 2)
    idx2sid = {}
    for (sid, _), tag in zip(series, tags):
        idx = rgb2idx.get(_hex_rgb(tag))
        if idx is None or idx in idx2sid:
            raise RuntimeError('the marker colour %s of pathway %s is not in the plot.cdxml colour table' % (tag, sid))
        idx2sid[idx] = sid

    # Axes: the longest vertical line is the left axis (left/top/bottom), the longest horizontal one the bottom axis (right)
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
            raise RuntimeError('a plot <curve> has %d points, not a multiple of 3: the plot output changed, stop and check' % len(pts))
        anchors = []
        for k in range(0, len(pts), 3):
            a, b, cc = pts[k:k + 3]
            if max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[0] - cc[0]), abs(a[1] - cc[1])) > tol:
                raise RuntimeError('a plot <curve> is no longer made of degenerate triples (a polyline): the plot output changed, stop and check')
            anchors.append(to_data(*a))
        segs[sid].append(anchors)
        attrs[sid] = {k: c.get(k) for k in ('Closed', 'FillType', 'CurveType')}

    # Point check: the segments without their shared end points must equal the input samples
    check = {}
    for sid, pts in series:
        if not segs[sid]:
            raise RuntimeError('plot.cdxml has no <curve> for pathway %s' % sid)
        joined = list(segs[sid][0])
        for s in segs[sid][1:]:
            joined.extend(s[1:])
        if len(joined) != len(pts):
            raise RuntimeError('pathway %s: plot returned %d anchors for %d inputs (thinned or rewritten)'
                               % (sid, len(joined), len(pts)))
        span_x = (x1 - x0) or 1.0
        span_y = (y1 - y0) or 1.0
        worst = max(max(abs(a[0] - b[0]) / span_x, abs(a[1] - b[1]) / span_y)
                    for a, b in zip(joined, pts))
        if worst > 1e-6:
            raise RuntimeError('pathway %s: mapped-back points differ from the samples by %.2e (relative); the mapping is not trustworthy' % (sid, worst))
        check[sid] = {'curve_nodes': len(segs[sid]), 'anchors': len(joined),
                      'nodes_points': [len(s) for s in segs[sid]], 'max_rel_err': worst}
    return segs, attrs, {'plot_axes': meta['axes'], 'per_pathway': check,
                         'plot_cdxml': cdxml, 'plot_axis_box_px': [left, top, right, bottom]}
