# -*- coding: utf-8 -*-
"""Negative controls: deliberately broken CDXML rendered natively, to confirm the pixel check catches it.

    D:/python/python.exe tests/negative_controls.py [spec.json]

Each variant must reach the expected verdict listed here; any mismatch exits with 1.
  ok_shifted       whole figure shifted 480 pt right, WidthPages as the script computes it -> must pass (WidthPages fix)
  wrap_onepage     the same shift with WidthPages=1 -> ChemDraw wraps the right-hand text -> must be caught (stray ink)
  no_curvetype     <curve> without CurveType -> ChemDraw drops the curve -> must be caught (curve coverage)
  long_label       one text replaced by a long string in the CDXML -> must be caught (stray ink)
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
    # A fourth item None = observe only, not scored (ChemDraw's own behaviour is not settled)
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
                print('%-13s render failed: %s' % (name, err))
                bad += 1
                continue
            r = ep_verify.verify_png(png, L, axes)
            failed = [i for i in r['items'] if not i['ok']]
            curve = '; '.join(i['detail'] for i in r['items'] if i['name'].startswith('curve:'))
            if want_ok is None:
                verdict, exp = 'observed', '-'
            else:
                verdict, exp = ('OK' if r['ok'] == want_ok else 'MISMATCH'), ('pass' if want_ok else 'caught')
                bad += r['ok'] != want_ok
            print('%-13s expected %-6s got %-6s %s  %s | curve %s' % (
                name, exp, 'pass' if r['ok'] else 'caught', verdict,
                '; '.join('%s: %s' % (i['name'], i['detail']) for i in failed)[:260], curve[:80]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('negative controls: %s' % ('all as expected' if not bad else '%d mismatch(es)' % bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
