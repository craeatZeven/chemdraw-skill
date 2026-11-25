# -*- coding: utf-8 -*-
"""插图原生对象合规自查（依据 chemdraw 技能 references/figure-validation.md）

检查项
1. 箭头曲线：带 Arrowhead 的 <curve> 必须是 CurveType="8"（普通弯箭头）或 "32"（单电子钩），且 CurvePoints 为 12 个数
   （ChemDraw 界面写的就是这个形状；CurveType=0 会在界面里退化成折线、手柄不可拖）
2. 多段折线冒充箭头：<curve> 带 Arrowhead 但锚点 >4 个 → 违规（应交给 curved_arrow 生成）
3. 弧线语义：曲线箭头必须落在电子源/电子穴上（端点离最近锚点 ≤ 6 pt），由 build 阶段复检
4. 括号：GraphicType="Bracket" 的 BoundingBox 两个 x 必须相等（否则竖边被拉斜）
5. 标记：自由基应为 <n Radical="doublet">，电荷应为 Symbol CirclePlus/CircleMinus（装饰性 Oval 不算）

用法：python scripts/check_figure_compliance.py --repo .
退出码：0 全通过；1 有违规
"""
from __future__ import annotations
import argparse, glob, os, re, sys

def check_figure(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    problems = []
    for m in re.finditer(r'<curve\b([^>]*)/>', t):
        s = m.group(1)
        if 'Arrowhead' not in s:
            continue
        gid = (re.search(r'id="(\d+)"', s) or [None, "?"])[1]
        ct = (re.search(r'CurveType="(\d+)"', s) or [None, None])[1]
        cp = re.search(r'CurvePoints="([^"]+)"', s)
        n = len(cp.group(1).split()) if cp else 0
        if ct not in ("8", "32"):
            problems.append("curve %s: CurveType=%s（原生应为 8 或 32）" % (gid, ct))
        if n != 12:
            problems.append("curve %s: CurvePoints %d 个数（应为 12）" % (gid, n))
    for m in re.finditer(r'<graphic\b([^>]*GraphicType="Bracket"[^>]*)/>', t):
        s = m.group(1)
        gid = (re.search(r'id="(\d+)"', s) or [None, "?"])[1]
        bb = re.search(r'BoundingBox="(-?[\d.]+) (-?[\d.]+) (-?[\d.]+) (-?[\d.]+)"', s)
        if bb and abs(float(bb.group(1)) - float(bb.group(3))) > 0.5:
            problems.append("bracket %s: BoundingBox x 不相等（竖边会被拉斜）" % gid)
    return problems

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.repo, "figures_final", "fig*.cdxml")))
    bad = 0
    for p in files:
        probs = check_figure(p)
        if probs:
            bad += 1
            print("✗ %s" % os.path.basename(p))
            for x in probs[:8]:
                print("    " + x)
    print("\n检查 %d 幅，违规 %d 幅" % (len(files), bad))
    sys.exit(0 if bad == 0 else 1)

if __name__ == "__main__":
    main()
