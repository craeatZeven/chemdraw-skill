#!/usr/bin/env python3
"""Scope rendering benchmarks.

Two corpora share one CLI:

coverage (default) — render every ``assets/coverage/*.cdxml`` through the local
ChemDraw and diff each PNG against the committed baseline of the same name::

    python scripts/run_scope_benchmark.py --render --side-by-side --output-dir output/scope-benchmark

Runs on any Python with Pillow + pywin32 (e.g. D:\\python\\python.exe); no venv
or toolkit install needed. ``--no-render`` audits the corpus (pairs, baseline
sizes) without ChemDraw.

journal — the 26-case cross-journal substrate_scope imitation benchmark
(briefs in ``assets/scope-benchmark/``, toolkit runtime required)::

    python scripts/run_scope_benchmark.py --corpus journal --no-render
    python scripts/run_scope_benchmark.py --corpus journal --render --side-by-side

coverage mode writes ``<output-dir>/renders/<name>.png``,
``<output-dir>/comparisons/<name>.compare.png`` (baseline | fresh render) and
``<output-dir>/benchmark-report.json`` with per-case render time, sizes, RMSE
and differing-pixel ratio. Render failures land in the report's ``failed``
array; the run always finishes the remaining cases.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "assets" / "scope-benchmark" / "manifest.json"
CDX = ROOT / "skill" / "chemdraw" / "scripts" / "cdx.py"
COVERAGE_DIR = ROOT / "assets" / "coverage"

# Grey-level distance above which a pixel counts as "different".
DIFF_THRESHOLD = 16


def _load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _stitch_side_by_side(
    orig_png: Path,
    drawn_png: Path,
    out_png: Path,
    title: str,
    left_label: str = "Published Original (CC BY 4.0)",
    right_label: str = "ChemDraw Skill Imitation",
) -> bool:
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
    draw.text((pad, 14), f"{title}   |   LEFT: {left_label}   |   RIGHT: {right_label}", fill="black")
    canvas.paste(l_im, (pad, banner_h))
    canvas.paste(r_im, (l_im.width + pad * 2, banner_h))
    draw.line([(l_im.width + pad * 1.5, banner_h), (l_im.width + pad * 1.5, banner_h + target_h)], fill="#CCCCCC", width=2)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_png)
    return True


# ---------------------------------------------------------------------------
# coverage corpus
# ---------------------------------------------------------------------------

def _pixel_metrics(baseline: Path, render: Path) -> dict:
    """Size check + RMSE + differing-pixel ratio (pure Pillow, C-speed histogram)."""
    from PIL import Image, ImageChops

    b = Image.open(baseline).convert("RGB")
    r = Image.open(render).convert("RGB")
    out = {
        "baseline_size": [b.width, b.height],
        "render_size": [r.width, r.height],
        "size_match": b.size == r.size,
        "rmse": None,
        "diff_pixel_ratio": None,
    }
    if not out["size_match"]:
        return out
    hist = ImageChops.difference(b, r).convert("L").histogram()
    total = float(sum(hist))
    if total == 0:
        return out
    out["rmse"] = round(math.sqrt(sum(v * d * d for d, v in enumerate(hist)) / total), 3)
    out["diff_pixel_ratio"] = round(sum(hist[DIFF_THRESHOLD + 1:]) / total, 4)
    return out


def _collect_coverage_cases(only: list | None, limit: int | None, renders_dir: Path) -> list[dict]:
    cases = []
    for cdxml in sorted(COVERAGE_DIR.glob("*.cdxml")):
        name = cdxml.stem  # "foo" or "foo.native"
        baseline = COVERAGE_DIR / f"{name}.png"
        cases.append({
            "name": name,
            "cdxml": cdxml,
            "baseline": baseline if baseline.is_file() else None,
            "render_png": renders_dir / f"{name}.png",
        })
    if only:
        want = set(only)
        cases = [c for c in cases if c["name"] in want]
    if limit:
        cases = cases[: limit]
    return cases


def _render_coverage_windows(cases: list[dict], dpi: int) -> None:
    """Open every cdxml in ONE hidden ChemDraw COM session, timing each case.

    A per-case exception is recorded and the loop continues; a dead session
    marks all remaining cases failed rather than aborting the run.
    """
    sys.path.insert(0, str(ROOT / "toolkit"))
    from cdxml_toolkit.chemdraw.cdxml_to_image import _get_chemdraw

    try:
        app, launched = _get_chemdraw()
    except Exception as e:  # cannot even attach/dispatch
        for c in cases:
            c.setdefault("error", f"COM session: {e}")
        return
    was_visible = app.Visible
    try:
        app.Visible = False
        prefs = app.Preferences
        prefs.TransparentPNGs = False
        prefs.PNGResolution = dpi
        for c in cases:
            t0 = time.perf_counter()
            doc = None
            try:
                doc = app.Documents.Open(str(c["cdxml"].resolve()))
                doc.SaveAs(str(c["render_png"].resolve()))
                c["render_ms"] = round((time.perf_counter() - t0) * 1000, 1)
            except Exception as e:
                c["render_ms"] = round((time.perf_counter() - t0) * 1000, 1)
                c["error"] = f"render: {e}"
            finally:
                if doc is not None:
                    try:
                        doc.Close(False)
                    except Exception:
                        pass
    except Exception as e:  # session died mid-run
        for c in cases:
            if "render_ms" not in c and "error" not in c:
                c["error"] = f"COM session lost: {e}"
    finally:
        try:
            if launched:
                app.Quit()
            else:
                app.Visible = was_visible
        except Exception:
            pass


def _run_coverage(args: argparse.Namespace) -> int:
    if args.output_dir:
        out_root = Path(args.output_dir).expanduser().resolve()
    else:
        out_root = ROOT / "output" / "scope-benchmark"
    renders_dir = out_root / "renders"
    comps_dir = out_root / "comparisons"
    out_root.mkdir(parents=True, exist_ok=True)

    cases = _collect_coverage_cases(args.only, args.limit, renders_dir)
    if not cases:
        print(f"no *.cdxml under {COVERAGE_DIR}", file=sys.stderr)
        return 1

    if args.render:
        renders_dir.mkdir(parents=True, exist_ok=True)
        sys.path.insert(0, str(ROOT / "toolkit"))
        from cdxml_toolkit.chemdraw.mac_automation import is_macos
        if is_macos():
            from cdxml_toolkit.chemdraw.cdxml_to_image import cdxml_to_image
            for c in cases:
                t0 = time.perf_counter()
                try:
                    cdxml_to_image(str(c["cdxml"]), str(c["render_png"]), png_dpi=args.dpi)
                    c["render_ms"] = round((time.perf_counter() - t0) * 1000, 1)
                except Exception as e:
                    c["render_ms"] = round((time.perf_counter() - t0) * 1000, 1)
                    c["error"] = f"render: {e}"
        else:
            _render_coverage_windows(cases, args.dpi)

    rows, failed = [], []
    comparisons = 0
    failed_names = set()
    for c in cases:
        row = {
            "name": c["name"],
            "cdxml": str(c["cdxml"]),
            "baseline_png": str(c["baseline"]) if c["baseline"] else None,
            "render_png": str(c["render_png"]) if args.render else None,
            "render_ms": c.get("render_ms"),
        }
        if c.get("error"):
            failed.append({"name": c["name"], "cdxml": str(c["cdxml"]), "error": c["error"]})
            failed_names.add(c["name"])
            rows.append(row)
            continue
        if c["baseline"] is None:
            row["note"] = "no committed baseline PNG for this cdxml"
            rows.append(row)
            continue
        if not args.render:
            from PIL import Image
            with Image.open(c["baseline"]) as im:
                row["baseline_size"] = [im.width, im.height]
            rows.append(row)
            continue
        if not c["render_png"].is_file():
            failed.append({"name": c["name"], "cdxml": str(c["cdxml"]),
                           "error": "render: ChemDraw reported success but no PNG on disk"})
            failed_names.add(c["name"])
            rows.append(row)
            continue
        row.update(_pixel_metrics(c["baseline"], c["render_png"]))
        if args.side_by_side:
            comp = comps_dir / f"{c['name']}.compare.png"
            if _stitch_side_by_side(c["baseline"], c["render_png"], comp, c["name"],
                                    left_label="Committed Baseline",
                                    right_label=f"Fresh ChemDraw Render ({args.dpi} dpi)"):
                row["compare_png"] = str(comp)
                comparisons += 1
            else:
                row["compare_png"] = None
        rows.append(row)

    rmses = sorted(r["rmse"] for r in rows if r.get("rmse") is not None)
    summary = {
        "figure_type": "coverage-render",
        "corpus": str(COVERAGE_DIR),
        "rendered": bool(args.render),
        "dpi": args.dpi if args.render else None,
        "diff_threshold_grey": DIFF_THRESHOLD,
        "total": len(rows),
        "render_failed": len(failed),
        "comparisons": comparisons,
        "size_match": sum(1 for r in rows if r.get("size_match") is True),
        "size_mismatch": sum(1 for r in rows if r.get("size_match") is False),
        "rmse_median": rmses[len(rmses) // 2] if rmses else None,
        "rmse_worst": rmses[-1] if rmses else None,
        "output_dir": str(out_root),
        "cases": rows,
        "failed": failed,
    }
    report_path = out_root / "benchmark-report.json"
    report_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for r in rows:
        tag = "FAIL" if r["name"] in failed_names else ("AUDIT" if not args.render else "OK")
        extra = ""
        if r.get("rmse") is not None:
            extra = f"  rmse={r['rmse']}  diff%={r['diff_pixel_ratio'] * 100:.2f}"
        elif r.get("size_match") is False:
            extra = f"  SIZE MISMATCH {r['baseline_size']} vs {r['render_size']}"
        elif r.get("note"):
            extra = f"  ({r['note']})"
        print(f"[{tag:5s}] {r['name']:44s} ms={r.get('render_ms')}{extra}")
    print(f"\nSUMMARY: total={summary['total']} failed={summary['render_failed']} "
          f"size_match={summary['size_match']}/{summary['size_match'] + summary['size_mismatch']} "
          f"rmse_median={summary['rmse_median']} comparisons={summary['comparisons']} -> {report_path}")
    return 0 if not failed else 1


# ---------------------------------------------------------------------------
# journal corpus (26-case imitation benchmark, toolkit runtime required)
# ---------------------------------------------------------------------------

def _run_in_toolkit(args: argparse.Namespace) -> int:
    # Import toolkit runtime directly (this script is invoked via the runtime Python or bootstraps through cdx.py py).
    from cdxml_toolkit.cdxml_check import check_cdxml
    from cdxml_toolkit.mcp_runtime.figure_draw import draw_figure
    try:
        from cdxml_toolkit.mcp_runtime.figure_scope import draw_scope
    except Exception:  # Older runtimes may only export draw_figure.
        draw_scope = None

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

    def _draw_for_kind(brief, out, render):
        """Route scope briefs to the validator that accepts the scope-specific fields."""
        kind = (brief or {}).get("kind") or "scheme"
        if kind == "scope" and draw_scope is not None:
            return draw_scope(brief, output_path=out, render=render)
        return draw_figure(brief, output_path=out, render=render)

    results = []
    passed = 0
    for idx, case in enumerate(cases, 1):
        slug = case["slug"]
        brief = json.loads((ROOT / case["brief"]).read_text(encoding="utf-8"))
        case_dir = out_root / slug
        if case_dir.exists():
            shutil.rmtree(case_dir)
        report = _draw_for_kind(brief, str(case_dir), bool(args.render))
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
    parser = argparse.ArgumentParser(
        description="Run the scope rendering benchmark (coverage cdxml\u2192PNG, or journal imitation).")
    parser.add_argument("--corpus", choices=("coverage", "journal"), default="coverage",
                        help="coverage (default): assets/coverage cdxml-to-PNG render benchmark; "
                             "journal: 26-case substrate_scope imitation benchmark.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--render", action="store_true", help="Render PNGs via local ChemDraw.")
    mode.add_argument("--no-render", action="store_true", help="Audit/check only (no ChemDraw needed).")
    parser.add_argument("--side-by-side", action="store_true", help="Stitch baseline vs render PNGs (requires --render).")
    parser.add_argument("--output-dir", help="Directory for renders, comparisons, and benchmark-report.json.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution for coverage renders (default 300).")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary output directory when --output-dir is omitted (journal only).")
    parser.add_argument("--only", nargs="+", help="Run only specific case names/slugs.")
    parser.add_argument("--limit", type=int, help="Run only the first N cases.")
    parser.add_argument("--in-runtime", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.side_by_side and not args.render:
        parser.error("--side-by-side requires --render")

    if args.corpus == "coverage":
        return _run_coverage(args)

    if not args.in_runtime:
        try:
            import cdxml_toolkit  # noqa: F401
        except ImportError:
            cmd = [sys.executable, str(CDX), "py", str(Path(__file__).resolve()), "--in-runtime", *sys.argv[1:]]
            return subprocess.call(cmd)
    return _run_in_toolkit(args)


if __name__ == "__main__":
    raise SystemExit(main())
