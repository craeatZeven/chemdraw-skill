#!/usr/bin/env python3
"""Audit the openly licensed journal figure corpus against all 23 declared figure types.

Standard-library only. Checks CC BY licence discipline (rejecting -NC/-ND), verifies
that every indexed PNG exists on disk, and reports per-type coverage gaps for
cross-paper comparison work.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS_DIR = ROOT / "assets" / "journal-corpus"
INDEX_CSV = CORPUS_DIR / "INDEX.csv"
HOUSE_STYLES_PY = ROOT / "toolkit" / "cdxml_toolkit" / "publication" / "house_styles.py"

FORBIDDEN_LICENCE = re.compile(r"(-nc\b|-nd\b|noncommercial|non-commercial|no\s*derivatives)", re.IGNORECASE)

FIGURE_TYPES = (
    "reaction_scheme", "substrate_scope", "mechanism", "catalytic_cycle", "energy_profile",
    "total_synthesis", "retrosynthesis", "structure_panel", "sar_table",
    "optimization_table", "concept_panel", "graphical_abstract", "stereo_model", "derivatization",
    "biosynthesis", "biomolecule", "coordination_complex", "polymer", "markush", "process_route",
    "photophysical_diagram", "assignment_panel", "setup_diagram",
)


def audit_corpus(corpus_dir: Path = CORPUS_DIR, target_per_type: int = 3) -> dict:
    index_path = corpus_dir / "INDEX.csv"
    if not index_path.is_file():
        raise FileNotFoundError(f"Corpus index not found: {index_path}")

    with index_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    errors: list[str] = []
    by_type: dict[str, list[dict]] = collections.defaultdict(list)
    by_journal: collections.Counter[str] = collections.Counter()

    for line_no, row in enumerate(rows, start=2):
        fname = (row.get("file") or "").strip()
        ftype = (row.get("figure_type") or "").strip()
        lic = (row.get("license") or "").strip()
        journal = (row.get("journal") or "").strip()
        if not fname:
            errors.append(f"Row {line_no}: empty file name")
            continue
        png = corpus_dir / f"{fname}.png"
        if not png.is_file():
            errors.append(f"Row {line_no} ({fname}): missing PNG file {png.name}")
        if not lic.upper().startswith("CC BY") or FORBIDDEN_LICENCE.search(lic):
            errors.append(f"Row {line_no} ({fname}): non-CC-BY licence {lic!r}")
        if ftype not in FIGURE_TYPES:
            errors.append(f"Row {line_no} ({fname}): unknown figure_type {ftype!r}")
        by_type[ftype].append(row)
        by_journal[journal] += 1

    coverage = {}
    uncovered: list[str] = []
    thin: list[str] = []
    for ftype in FIGURE_TYPES:
        items = by_type.get(ftype, [])
        journals = sorted({(r.get("journal") or "").strip() for r in items if r.get("journal")})
        count = len(items)
        coverage[ftype] = {
            "count": count,
            "journals": journals,
            "needed": max(0, target_per_type - count),
        }
        if count == 0:
            uncovered.append(ftype)
        elif count < target_per_type:
            thin.append(ftype)

    return {
        "ok": not errors,
        "total_figures": len(rows),
        "declared_types": len(FIGURE_TYPES),
        "covered_types": len(FIGURE_TYPES) - len(uncovered),
        "target_per_type": target_per_type,
        "by_journal": dict(by_journal.most_common()),
        "uncovered_types": uncovered,
        "thin_types": thin,
        "coverage": coverage,
        "errors": errors,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit the CC BY journal figure corpus.")
    parser.add_argument("--corpus-dir", type=Path, default=CORPUS_DIR)
    parser.add_argument("--target", type=int, default=3, help="Target figures per type (default: 3)")
    parser.add_argument("--markdown", action="store_true", help="Print a Markdown summary table instead of JSON")
    ns = parser.parse_args(argv)
    report = audit_corpus(ns.corpus_dir, target_per_type=ns.target)
    if ns.markdown:
        print(f"# Corpus Audit ({report['total_figures']} figures, {report['covered_types']}/{report['declared_types']} types)\n")
        print("| Figure Type | Count | Journals | Needed (target " + str(report["target_per_type"]) + ") |")
        print("| --- | ---: | --- | ---: |")
        for ftype in FIGURE_TYPES:
            info = report["coverage"][ftype]
            jtext = ", ".join(info["journals"]) or "—"
            print(f"| `{ftype}` | {info['count']} | {jtext} | {info['needed']} |")
    else:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
