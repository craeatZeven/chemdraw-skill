"""Render R-group demo panels through the local ChemDraw and stitch before/after strips.

Usage:
    python scripts/rgroup_demo.py --tag before          # render each case, write <tag>/<case>.png
    python scripts/rgroup_demo.py --tag after
    python scripts/rgroup_demo.py --compare             # stitch before|after compare PNGs

The three cases are the ones the round-2 brief asks for: a monosubstituted benzene
carrying one R, a disubstituted benzene carrying two R, and an R on a heteroaromatic
ring. Each is drawn as a scope panel whose header holds the general structure with the
attachment points, which is where a published substrate-scope scheme puts its R marks.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "toolkit"))

DEMO_DIR = ROOT / "output" / "rgroup-demo"

CASES = {
    # case id -> (title, header substrate SMILES, header product SMILES, one concrete scope member)
    "benzene_mono": (
        "Monosubstituted benzene, one R",
        "[*:1]c1ccccc1",
        "[*:1]c1ccccc1C(=O)N",
        "Oc1ccccc1",
    ),
    "benzene_di": (
        "Disubstituted benzene, two R",
        "[*:1]c1ccc([*:2])cc1",
        "[*:1]c1ccc(C(=O)N[*:2])cc1",
        "Oc1ccc(C)cc1",
    ),
    "heterocycle": (
        "R on a heteroaromatic ring",
        "[*:1]c1cccnc1",
        "[*:1]c1cccnc1C(=O)N",
        "Oc1cccnc1",
    ),
}


def _brief(title: str, general: str, product: str, concrete: str, r_attachment: str) -> dict:
    return {
        "kind": "scope",
        "style": "nature_portfolio",
        "width_mm": 90,
        "columns": 1,
        "title": title,
        "r_attachment": r_attachment,
        "header": {
            "substrates": [{"smiles": general, "label": "1"}],
            "products": [{"smiles": product, "label": "2"}],
            "below": ["MeCN, rt, 12 h"],
        },
        "fields": [{"key": "yield", "label": "Yield"}],
        "species": [
            {"id": "2a", "smiles": concrete, "label": "2a", "values": {"yield": "90%"}},
        ],
    }


def _render(tag: str, dpi: int, r_attachment: str) -> int:
    from cdxml_toolkit.mcp_runtime.figure_scope import draw_scope

    out = DEMO_DIR / tag
    pngs = out / "png"
    pngs.mkdir(parents=True, exist_ok=True)
    failures = []
    for case_id, (title, general, product, concrete) in CASES.items():
        target = out / case_id
        if target.exists():
            import shutil
            shutil.rmtree(target)
        report = draw_scope(_brief(title, general, product, concrete, r_attachment), str(target), render=True)
        ok = bool(report.get("ok"))
        native = (report.get("native") or {})
        print(f"[{case_id}] ok={ok} readback={report.get('readback', {}).get('status')} "
              f"native={native.get('status')}")
        if not ok:
            failures.append((case_id, report.get("errors")))
            continue
        produced = native.get("png_paths") or []
        if not produced:
            failures.append((case_id, native.get("reason") or "no png"))
            continue
        import shutil
        for i, p in enumerate(produced):
            dst = pngs / (f"{case_id}.png" if len(produced) == 1 else f"{case_id}_p{i}.png")
            shutil.copyfile(p, dst)
            print(f"    -> {dst}")
        (out / f"{case_id}.report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if failures:
        print("FAILURES:", failures)
    return 1 if failures else 0


def _compare() -> int:
    from PIL import Image, ImageDraw

    comps = DEMO_DIR / "compare"
    comps.mkdir(parents=True, exist_ok=True)
    made = []
    for case_id in CASES:
        left = DEMO_DIR / "before" / "png" / f"{case_id}.png"
        right = DEMO_DIR / "after" / "png" / f"{case_id}.png"
        if not left.is_file() or not right.is_file():
            print(f"[skip] {case_id}: missing {left.name if not left.is_file() else right.name}")
            continue
        l_im = Image.open(left).convert("RGB")
        r_im = Image.open(right).convert("RGB")
        h = max(l_im.height, r_im.height, 700)
        def _scale(im):
            ratio = h / max(im.height, 1)
            return im.resize((max(1, int(im.width * ratio)), h), Image.Resampling.LANCZOS)
        l_im, r_im = _scale(l_im), _scale(r_im)
        pad, banner = 24, 56
        canvas = Image.new("RGB", (l_im.width + r_im.width + pad * 3, h + banner + pad * 2), "white")
        d = ImageDraw.Draw(canvas)
        d.text((pad, 16), f"{case_id}   |   LEFT: before (ring-vertex R bond)   |   "
                          f"RIGHT: after (ring-centre crossing R bond)", fill="black")
        canvas.paste(l_im, (pad, banner))
        canvas.paste(r_im, (l_im.width + pad * 2, banner))
        d.line([(l_im.width + pad * 1.5, banner), (l_im.width + pad * 1.5, banner + h)],
               fill="#CCCCCC", width=2)
        out = comps / f"{case_id}.compare.png"
        canvas.save(out)
        made.append(out)
        print(f"[compare] {out}  {canvas.size}")
    return 0 if made else 1


def _render_single(tag: str, r_attachment: str) -> int:
    """One structure per page, so nothing else in the frame can be mistaken for the R bond.

    Goes through the same compose_chemical_figure writer the scope panels use.
    """
    from cdxml_toolkit.mcp_runtime.figure_tools import compose_chemical_figure
    from cdxml_toolkit.chemdraw.cdxml_to_image import cdxml_to_image

    singles = {
        "benzene_mono": "[*:1]c1ccccc1",
        "benzene_di": "[*:1]c1ccc([*:2])cc1",
        "heterocycle": "[*:1]c1cccnc1",
    }
    for case, smiles in singles.items():
        work = DEMO_DIR / "single" / tag / case
        if work.exists():
            # compose_chemical_figure refuses to overwrite, so a re-render starts clean.
            shutil.rmtree(work)
        work.mkdir(parents=True, exist_ok=True)
        manifest = work / "manifest.json"
        manifest.write_text(json.dumps({
            "version": 1,
            "style": {"preset": "acs1996"},
            "objects": [{
                "type": "molecule", "id": "m1", "smiles": smiles,
                "position": [150, 110], "bond_length": 40.0,
                **({"r_attachment": "centre"} if r_attachment == "centre" else {}),
            }],
        }), encoding="utf-8")
        result = compose_chemical_figure(str(manifest), str(work / "figure.cdxml"))
        cdxml_to_image(str(work / "figure.cdxml"), str(work / "figure.png"), png_dpi=200)
        print(f"[single {tag}/{case}] {smiles}: ok={result.get('ok')} -> {work / 'figure.png'}")
    return 0


def _compare_single() -> int:
    from PIL import Image, ImageDraw

    comps = DEMO_DIR / "compare"
    comps.mkdir(parents=True, exist_ok=True)
    made = []
    for case in CASES:
        left = DEMO_DIR / "single" / "before" / case / "figure.png"
        right = DEMO_DIR / "single" / "after" / case / "figure.png"
        if not left.is_file() or not right.is_file():
            print(f"[skip] {case}: missing single-structure render")
            continue
        l_im = Image.open(left).convert("RGB")
        r_im = Image.open(right).convert("RGB")
        h = max(l_im.height, r_im.height, 420)
        def _scale(im):
            ratio = h / max(im.height, 1)
            return im.resize((max(1, int(im.width * ratio)), h), Image.Resampling.LANCZOS)
        l_im, r_im = _scale(l_im), _scale(r_im)
        pad, banner = 28, 60
        canvas = Image.new("RGB", (l_im.width + r_im.width + pad * 3, h + banner + pad * 2), "white")
        d = ImageDraw.Draw(canvas)
        d.text((pad, 16), f"{case}   |   LEFT: before, R on the ring vertex   |   "
                          f"RIGHT: after, R from the ring centre through a ring bond", fill="black")
        canvas.paste(l_im, (pad, banner))
        canvas.paste(r_im, (l_im.width + pad * 2, banner))
        d.line([(l_im.width + pad * 1.5, banner), (l_im.width + pad * 1.5, banner + h)],
               fill="#CCCCCC", width=2)
        out = comps / f"{case}.compare.png"
        canvas.save(out)
        made.append(out)
        print(f"[compare] {out}  {canvas.size}")
    return 0 if made else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", choices=("before", "after"), help="render the demo panels into output/rgroup-demo/<tag>")
    ap.add_argument("--compare", action="store_true", help="stitch before|after comparison PNGs")
    ap.add_argument("--single", action="store_true", help="render one structure per page (clearest R-bond view)")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--r-attachment", default="vertex", choices=("vertex", "centre"),
                    help="pass --r-attachment centre for the 'after' render (ring-centre crossing bond)")
    args = ap.parse_args()
    if args.compare:
        return _compare_single() if args.single else _compare()
    if args.tag:
        if args.single:
            return _render_single(args.tag, args.r_attachment)
        return _render(args.tag, args.dpi, args.r_attachment)
    ap.error("pass --tag before|after or --compare")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
