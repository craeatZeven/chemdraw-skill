"""Round-2 demos for the abbreviation and colour rules.

Each demo renders two panels from briefs that differ in exactly one key and stitches them
side by side, so the image shows the rule rather than asserting it:

  abbrev_repeated     a protecting group repeated across a scope: expanded | abbreviate: ["Boc"]
  abbrev_must_expand  the group that varies is the point of the figure, so it stays expanded
  colour_provenance   two partners keep their own colour into the product (fragment_colors)
  colour_varying      one core elaborated at one position (highlight_varying)
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "toolkit"))

DEMO_DIR = ROOT / "output" / "round2-demos"

RED = "#8B1E24"
BLUE = "#1B4D89"
TEAL = "#008080"


def _scope(species, header=None, **extra):
    brief = {
        "kind": "scope",
        "style": "nature_portfolio",
        "width_mm": 150,
        "columns": 3,
        "fields": [{"key": "yield", "label": "Yield"}],
        "species": species,
    }
    if header:
        brief["header"] = header
    brief.update(extra)
    return brief


def _sp(ident, smiles, label, group=None, **extra):
    row = {"id": ident, "smiles": smiles, "label": label, "values": {"yield": "90%"}}
    if group:
        row["group"] = group
    row.update(extra)
    return row


def cases() -> dict:
    # A Boc that appears on every member of the series: an auxiliary group, repeated, and not
    # the point of the figure, so the rule contracts it.
    repeated = [
        _sp("1a", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1", "1a"),
        _sp("1b", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1F", "1b"),
        _sp("1c", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1Cl", "1c"),
    ]
    # Here the N-protecting group IS the variable under study, so it must stay expanded.
    varying = [
        _sp("2a", "CC(C)(C)OC(=O)Nc1ccccc1", "2a, PG = Boc"),
        _sp("2b", "CC(=O)Nc1ccccc1", "2b, PG = Ac"),
        _sp("2c", "O=C(Nc1ccccc1)OCc1ccccc1", "2c, PG = Cbz"),
    ]
    header = {
        "substrates": [
            {"smiles": "CC(C)(C)OC(=O)Nc1ccccc1", "label": "1", "color": RED},
            {"smiles": "OC(=O)c1ccc(F)cc1", "label": "2", "color": BLUE},
        ],
        "products": [{"smiles": "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1", "label": "3"}],
        "below": ["EDCI, DMAP, DCM, rt"],
    }
    provenance = [
        _sp("3a", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1", "3a"),
        _sp("3b", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1F", "3b"),
        _sp("3c", "CC(C)(C)OC(=O)Nc1ccc(C(=O)O)cc1Cl", "3c"),
    ]
    varying_core = [
        _sp("4a", "Oc1ccccc1C=O", "4a"),
        _sp("4b", "COc1ccc(C=O)cc1", "4b"),
        _sp("4c", "Oc1ccc(C=O)cc1Cl", "4c"),
        _sp("4d", "COc1ccc(C=O)c(Cl)c1", "4d"),
    ]
    return {
        "abbrev_repeated": (
            _scope(repeated, header=header),
            _scope([dict(s, abbreviate=["Boc"]) for s in repeated], header=header),
            "LEFT: Boc written out  |  RIGHT: abbreviate = ['Boc'] (repeated auxiliary group)",
        ),
        "abbrev_must_expand": (
            _scope(varying),
            _scope([dict(s, abbreviate=["Boc"]) for s in varying]),
            "LEFT: protecting groups expanded, as the rule requires  |  "
            "RIGHT: same brief with abbreviate = ['Boc'] (what NOT to do)",
        ),
        "colour_provenance": (
            _scope(provenance, header=header),
            _scope(provenance, header=header,
                   fragment_colors=[{"smarts": "NC(=O)OC(C)(C)C", "color": RED},
                                    {"smarts": "C(=O)O", "color": BLUE}]),
            "LEFT: no colour  |  RIGHT: fragment_colors keeps each partner's colour",
        ),
        "colour_varying": (
            _scope(varying_core),
            _scope(varying_core, highlight_varying=True),
            "LEFT: no colour  |  RIGHT: highlight_varying marks only what changes",
        ),
    }


def render(demo: str) -> int:
    from cdxml_toolkit.mcp_runtime.figure_scope import draw_scope

    out = DEMO_DIR / demo
    for side, brief in zip(("left", "right"), cases()[demo][:2]):
        target = out / side
        if target.exists():
            shutil.rmtree(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        report = draw_scope(brief, str(target), render=True)
        pngs = (report.get("native") or {}).get("png_paths") or []
        print(f"[{demo}/{side}] ok={report.get('ok')} readback={report.get('readback', {}).get('status')} "
              f"pngs={len(pngs)} warnings={report.get('warnings')}")
        if not report.get("ok"):
            print("   errors:", report.get("errors"))
            return 1
        if pngs:
            shutil.copyfile(pngs[0], out / f"{side}.png")
    return 0


def stitch(demo: str) -> int:
    from PIL import Image, ImageDraw

    out = DEMO_DIR / demo
    caption = cases()[demo][2]
    paths = [out / "left.png", out / "right.png"]
    if not all(p.is_file() for p in paths):
        print(f"[skip] {demo}: missing rendered panel")
        return 1
    ims = [Image.open(p).convert("RGB") for p in paths]
    h = max(im.height for im in ims)
    scaled = [im.resize((int(im.width * h / im.height), h), Image.Resampling.LANCZOS) for im in ims]
    pad, banner = 24, 70
    width = sum(im.width for im in scaled) + pad * 3
    canvas = Image.new("RGB", (width, h + banner + pad * 2), "white")
    d = ImageDraw.Draw(canvas)
    d.text((pad, 18), f"{demo}", fill="black")
    d.text((pad, 40), caption, fill="black")
    x = pad
    for im in scaled:
        canvas.paste(im, (x, banner))
        x += im.width + pad
    dest = DEMO_DIR / f"{demo}.compare.png"
    canvas.save(dest)
    print(f"[stitch] {dest} {canvas.size}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("demo", choices=[*cases(), "all"])
    args = ap.parse_args()
    demos = list(cases()) if args.demo == "all" else [args.demo]
    rc = 0
    for demo in demos:
        rc |= render(demo)
        rc |= stitch(demo)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
