# -*- coding: utf-8 -*-
"""反向对照：故意做坏的 CDXML 交给 ChemDraw 原生渲染，确认像素核对能抓到。

    D:/python/python.exe tests/negative_controls.py [spec.json]

每个变体都应得到"预期结果"一栏写的结论；任何一个不符，退出码 1。
  ok_shifted       整图右移 480 pt，WidthPages 按脚本规则计算 → 应全部通过（验证 WidthPages 修复）
  wrap_onepage     同样右移，但 WidthPages=1 → 右侧文字被 ChemDraw 折行 → 应抓到（多余墨迹）
  no_curvetype     <curve> 去掉 CurveType（陷阱 1）→ ChemDraw 丢弃曲线 → 应抓到（曲线覆盖率）
  long_label       一段文字在 CDXML 里被换成长串 → 应抓到（多余墨迹）
"""
import copy
import json
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import energy_profile as E   # noqa: E402
import ep_verify             # noqa: E402


def main():
    for s in (sys.stdout, sys.stderr):
        s.reconfigure(encoding='utf-8')
    spec_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(HERE), 'examples', 'two_step_cn.json')
    with open(spec_path, encoding='utf-8') as fh:
        spec = json.load(fh)
    skill = next(d for d in E.SKILL_DIRS if os.path.isfile(os.path.join(d, 'scripts', 'cdx.py')))
    cdx = [E.CDX_PY, os.path.join(skill, 'scripts', 'cdx.py')]
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    M = E.build(spec, cdx, env)
    W, H = M['expect']['page']

    L2, ax2 = copy.deepcopy(M['L']), dict(M['axes'])
    E.shift_all(L2, ax2, 480, 0)
    shifted, _ = E.serialize(L2, W + 480, H)
    latin = next(t for t in M['L']['texts'] if t['role'] == 'species' and len(t['runs']) == 1
                 and all(ord(c) < 128 for c in t['runs'][0][0]))['runs'][0][0]
    long_xml = M['xml'].replace('>%s</s>' % latin, '>%s XXXXXXXXXXXXXXXX</s>' % latin, 1)
    assert long_xml != M['xml']
    # 第 4 项 None = 只观察、不计入对错（ChemDraw 行为本身待定）
    cases = {
        'ok_shifted': (shifted, L2, ax2, True),
        'wrap_onepage': (re.sub(r'WidthPages="\d+"', 'WidthPages="1"', shifted), L2, ax2, False),
        'curve_removed': (re.sub(r'<curve\b[^>]*/>', '', M['xml']), M['L'], M['axes'], False),
        'long_label': (long_xml, M['L'], M['axes'], False),
        'no_curvetype': (re.sub(r' CurveType="\d+"', '', M['xml']), M['L'], M['axes'], None),
    }
    tmp = tempfile.mkdtemp(prefix='ep_neg_')
    bad = 0
    try:
        for name, (xml, L, axes, want_ok) in cases.items():
            cd = os.path.join(tmp, name + '.cdxml')
            png = os.path.join(tmp, name + '.png')
            with open(cd, 'w', encoding='utf-8') as fh:
                fh.write(xml)
            ok, err = E.render(cd, png, cdx, env)
            if not ok:
                print('%-13s 渲染失败：%s' % (name, err))
                bad += 1
                continue
            r = ep_verify.verify_png(png, L, axes)
            failed = [i for i in r['items'] if not i['ok']]
            curve = '; '.join(i['detail'] for i in r['items'] if i['name'].startswith('curve:'))
            if want_ok is None:
                verdict, exp = '观察', '—'
            else:
                verdict, exp = ('OK' if r['ok'] == want_ok else '不符'), ('通过' if want_ok else '抓到')
                bad += r['ok'] != want_ok
            print('%-13s 预期 %-4s 实得 %-4s %s  %s | 曲线 %s' % (
                name, exp, '通过' if r['ok'] else '抓到', verdict,
                '; '.join('%s: %s' % (i['name'], i['detail']) for i in failed)[:260], curve[:80]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('反向对照：%s' % ('全部符合预期' if not bad else '%d 项不符' % bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
