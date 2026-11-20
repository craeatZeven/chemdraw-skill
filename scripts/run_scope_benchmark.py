#!/usr/bin/env python3
"""Run the 26-case cross-journal substrate_scope imitation benchmark.

Portable mode (Linux / CI, no ChemDraw required):
    python scripts/run_scope_benchmark.py --no-render

Local ChemDraw / DeepSeek harness mode (Windows / macOS with activated ChemDraw):
    python scripts/run_scope_benchmark.py --render --side-by-side

Generates editable CDXMLs, runs profile-aware ``cdxml_check`` on every page, and
(when ``--render`` is passed on a machine with ChemDraw) renders 300 dpi PNGs and
stitches side-by-side comparison images against the CC BY 4.0 journal originals.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "assets" / "scope-benchmark" / "manifest.json"
CDX = ROOT / "skill" / "chemdraw" / "scripts" / "cdx.py"


def _load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _stitch_side_by_side(orig_png: Path, drawn_png: Path, out_png: Path, title: str) -> bool:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return False
    if not orig_png.is_file() or not drawn_png.is_file():
        return False
    left = Image.open(orig_png).convert("RGB")
    right = Image.open(drawn_png).convert("RGB")
    target_h = max(min(max(left.height, right.height), 1800), 900)
    def _scale(im: Image.Image) -> Image.Image:
        ratio = target_h / max(im.height, 1)
        return im.resize((max(1, int(im.width * ratio)), target_h), Image.Resampling.LANCZOS)
    l_im = _scale(left)
    r_im = _scale(right)
    pad, banner_h = 24, 48
    canvas = Image.new("RGB", (l_im.width + r_im.width + pad * 3, target_h + banner_h + pad * 2), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((pad, 14), f"{title}   |   LEFT: Published Original (CC BY 4.0)   |   RIGHT: ChemDraw Skill Imitation", fill="black")
    canvas.paste(l_im, (pad, banner_h))
    canvas.paste(r_im, (l_im.width + pad * 2, banner_h))
    draw.line([(l_im.width + pad * 1.5, banner_h), (l_im.width + pad * 1.5, banner_h + target_h)], fill="#CCCCCC", width=2)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_png)
    return True


def _run_in_toolkit(args: argparse.Namespace) -> int:
    # Import toolkit runtime directly (this script is invoked via the runtime Python or bootstraps through cdx.py py).
    from cdxml_toolkit.cdxml_check import check_cdxml
    from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure

    os.environ.setdefault("CHEMDRAW_OFFLINE", "1")
    manifest = _load_manifest()
    cases = manifest["cases"]
    if args.only:
        want = set(args.only)
        cases = [c for c in cases if c["slug"] in want]
    if args.limit:
        cases = cases[: args.limit]

    if args.output_dir:
        out_root = Path(args.output_dir).expanduser().resolve()
        out_root.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        out_root = Path(tempfile.mkdtemp(prefix="scope-benchmark-"))
        cleanup = not args.keep_temp

    results = []
    passed = 0
    for idx, case in enumerate(cases, 1):
        slug = case["slug"]
        brief = json.loads((ROOT / case["brief"]).read_text(encoding="utf-8"))
        case_dir = out_root / slug
        if case_dir.exists():
            shutil.rmtree(case_dir)
        report = draw_figure(brief, output_path=str(case_dir), render=bool(args.render))
        entry = {
            "slug": slug,
            "journal": case["journal"],
            "figure": case["figure"],
            "archetype": case["archetype"],
            "check_profile": case["check_profile"],
            "ok": bool(report.get("ok")),
            "checks_passed": bool(report.get("checks_passed")),
            "readback": report.get("readback", {}).get("status"),
            "pages": len(report.get("pages") or []),
            "lint_issues": len((report.get("lint") or {}).get("issues") or []),
            "warnings": report.get("warnings") or [],
            "errors": report.get("errors") or [],
            "native_status": (report.get("native") or {}).get("status"),
            "check_errors": [],
            "check_warnings": [],
            "side_by_side": None,
        }
        if report.get("ok"):
            for page in report.get("pages") or []:
                chk = check_cdxml(
                    page["path"],
                    profile=case["check_profile"],
                    palette=None,
                    allow_sizes=(7.0, 8.0, 10.0),
                )
                entry["check_errors"].extend(chk.get("errors") or [])
                entry["check_warnings"].extend(chk.get("warnings") or [])
            if args.render and args.side_by_side:
                pngs = (report.get("native") or {}).get("png_paths") or []
                if pngs:
                    comp = out_root / "comparisons" / f"{slug}.compare.png"
                    if _stitch_side_by_side(ROOT / case["corpus_png"], Path(pngs[0]), comp,
                                            f"[{case['journal']}] {slug}"):
                        entry["side_by_side"] = str(comp)
        case_ok = (
            entry["ok"]
            and entry["checks_passed"]
            and entry["readback"] == "matched"
            and entry["lint_issues"] == 0
            and not entry["check_errors"]
        )
        entry["passed"] = case_ok
        if case_ok:
            passed += 1
        status = "PASS" if case_ok else "FAIL"
        print(
            f"[{idx:02d}/{len(cases):02d}] {status}  {slug:48s}  pages={entry['pages']}  "
            f"readback={entry['readback']}  chk_err={len(entry['check_errors'])}  "
            f"warn={len(entry['warnings'])}"
        )
        if not case_ok:
            if entry["errors"]:
                print("    errors:", entry["errors"])
            if entry["check_errors"]:
                print("    check_errors:", entry["check_errors"][:3])
            if entry["warnings"]:
                print("    warnings:", entry["warnings"][:3])
        results.append(entry)

    summary = {
        "figure_type": "substrate_scope",
        "total": len(results),
        "passed": passed,
        "failed": len(results) - passed,
        "rendered": bool(args.render),
        "output_dir": str(out_root),
        "cases": results,
    }
    report_path = out_root / "benchmark-report.json"
    report_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nSUMMARY: {passed}/{len(results)} passed (rendered={bool(args.render)}) -> {report_path}")
    if cleanup:
        shutil.rmtree(out_root, ignore_errors=True)
    return 0 if passed == len(results) else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the 26-case substrate_scope imitation benchmark.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--render", action="store_true", help="Render 300 dpi PNGs via local ChemDraw.")
    mode.add_argument("--no-render", action="store_true", help="Portable check only (no ChemDraw needed).")
    parser.add_argument("--side-by-side", action="store_true", help="Stitch original vs imitation PNGs (requires --render).")
    parser.add_argument("--output-dir", help="Directory for generated CDXMLs, PNGs, comparisons, and benchmark-report.json.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary output directory when --output-dir is omitted.")
    parser.add_argument("--only", nargs="+", help="Run only specific case slugs.")
    parser.add_argument("--limit", type=int, help="Run only the first N cases.")
    parser.add_argument("--in-runtime", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if not args.in_runtime:
        try:
            import cdxml_toolkit  # noqa: F401
        except ImportError:
            cmd = [sys.executable, str(CDX), "py", str(Path(__file__).resolve()), "--in-runtime", *sys.argv[1:]]
            return subprocess.call(cmd)
    return _run_in_toolkit(args)


if __name__ == "__main__":
    raise SystemExit(main())
