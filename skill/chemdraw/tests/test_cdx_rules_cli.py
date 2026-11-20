"""Tests for new cdx.py CLI ergonomics: check --profile, validate, and schema."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CDX = ROOT / "skill" / "chemdraw" / "scripts" / "cdx.py"


def _run_cdx(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CDX), *args],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )


def test_cdx_schema_lists_all_kinds():
    proc = _run_cdx("schema")
    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    for kind in ("scheme", "scope", "mechanism", "panels", "table", "build"):
        assert kind in data["kinds"]


def test_cdx_schema_returns_table_schema_and_example():
    proc = _run_cdx("schema", "table")
    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["kind"] == "table"
    assert "columns" in data["schema"]["properties"]
    assert "rows" in data["schema"]["properties"]
    assert data["example"]["kind"] == "table"


def test_cdx_check_accepts_profile_and_named_palette(tmp_path):
    brief_path = tmp_path / "brief.json"
    out_cdxml = tmp_path / "out.cdxml"
    brief_path.write_text(json.dumps({
        "style": "acs_jacs",
        "species": [
            {"id": "A", "smiles": "CCO", "label": "1"},
            {"id": "B", "smiles": "CC=O", "label": "2"},
        ],
        "steps": [{"substrates": ["A"], "products": ["B"], "above": ["DMP"]}],
    }), encoding="utf-8")
    draw_proc = _run_cdx("draw", str(brief_path), str(out_cdxml))
    assert draw_proc.returncode == 0, draw_proc.stderr + draw_proc.stdout
    chk_proc = _run_cdx("check", str(out_cdxml), "--profile", "acs_jacs")
    assert chk_proc.returncode == 0, chk_proc.stderr + chk_proc.stdout
    report = json.loads(chk_proc.stdout)
    assert report["ok"] is True
    assert report["profile"] == "acs_jacs"


def test_cdx_validate_runs_basic_validation(tmp_path):
    brief_path = tmp_path / "brief.json"
    out_cdxml = tmp_path / "out.cdxml"
    val_dir = tmp_path / "validation"
    brief_path.write_text(json.dumps({
        "species": [
            {"id": "A", "smiles": "c1ccccc1", "label": "1"},
            {"id": "B", "smiles": "c1ccccc1C", "label": "2"},
        ],
        "steps": [{"substrates": ["A"], "products": ["B"], "above": ["MeCl, AlCl3"]}],
    }), encoding="utf-8")
    draw_proc = _run_cdx("draw", str(brief_path), str(out_cdxml))
    assert draw_proc.returncode == 0, draw_proc.stderr + draw_proc.stdout
    val_proc = _run_cdx("validate", str(out_cdxml), "--output-dir", str(val_dir))
    assert val_proc.returncode == 0, val_proc.stderr + val_proc.stdout
    data = json.loads(val_proc.stdout)
    assert data["ok"] is True
    assert Path(data["report"]).exists()


def test_cdx_schema_includes_energy_profile():
    proc = _run_cdx("schema", "energy_profile")
    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["kind"] == "energy_profile"
    assert "pathways" in data["schema"]["properties"]


def test_audit_corpus_script_passes_cc_by_check():
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "audit_corpus.py")],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    report = json.loads(proc.stdout)
    assert report["ok"] is True
    assert report["declared_types"] == 23
    assert report["total_figures"] == 71
