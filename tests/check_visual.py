"""Redraw the R-group cases and compare native ChemDraw exports with their goldens."""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
import json
import math
import os
from pathlib import Path
import tempfile
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = ROOT / "tests" / "golden"
PAPER_BRIEF = ROOT / "assets" / "scope-benchmark" / "briefs" / "bjoc-articles-22-127-scheme2.json"
CASES = {"mono": 1, "di": 2, "hetero": 1, "paper": 2}
DPI = 600
SIZE_PERCENT_LIMIT = 2
HASH_GRID = (32, 32)
PAPER_HEADER_GRID = (128, 32)
HASH_CUTOFFS = (192, 232, 250)
HASH_LIMIT = 16
COORDINATE_LIMIT_PT = 0.10
COORDINATES = ("ring_center_pt", "anchor_pt", "label_endpoint_pt", "crossed_edge_midpoint_pt")


def perceptual_hash(image, grid: tuple[int, int] = HASH_GRID) -> str:
    """Encode coarse RGB ink coverage at three intensities on a white background."""
    from PIL import Image

    if image.mode != "RGB":
        rgba = image.convert("RGBA")
        flattened = Image.new("RGBA", rgba.size, "white")
        flattened.alpha_composite(rgba)
        image = flattened.convert("RGB")
    thumbnail = image.resize(grid, Image.Resampling.BOX)
    bits = "".join(
        "1" if channel < cutoff else "0"
        for channel in thumbnail.tobytes()
        for cutoff in HASH_CUTOFFS
    )
    return f"{int(bits, 2):0{len(bits) // 4}x}"


def image_metrics(name: str, image) -> dict:
    metrics = {"size_px": list(image.size), "perceptual_hash": perceptual_hash(image)}
    if name == "paper":
        # The reaction header occupies the top seventh of the full scope render.
        header = image.crop((0, 0, image.width, max(1, image.height // 7)))
        metrics["header_perceptual_hash"] = perceptual_hash(header, PAPER_HEADER_GRID)
    return metrics


def hamming_distance(expected: str, actual: str) -> int:
    if len(expected) != len(actual):
        raise ValueError("Perceptual hash width changed")
    return (int(expected, 16) ^ int(actual, 16)).bit_count()


def ring_vertices(bonds, crossed):
    """Find the shortest ring containing the crossed edge."""
    start, stop = crossed.attrib["B"], crossed.attrib["E"]
    neighbors = defaultdict(list)
    for bond in bonds:
        if bond is not crossed:
            left, right = bond.attrib["B"], bond.attrib["E"]
            neighbors[left].append(right)
            neighbors[right].append(left)

    paths = deque([(start,)])
    visited = {start}
    while paths:
        path = paths.popleft()
        if path[-1] == stop:
            return path
        for neighbor in neighbors[path[-1]]:
            if neighbor not in visited:
                visited.add(neighbor)
                paths.append((*path, neighbor))
    raise ValueError("Crossed bond does not belong to a ring")


def measure_marks(path: Path) -> list[dict]:
    """Read R-group landmarks in native CDXML page coordinates, in points."""
    marks = []
    for fragment in ET.parse(path).getroot().iter("fragment"):
        nodes = {node.attrib["id"]: node for node in fragment.findall("n")}
        bonds = {bond.attrib["id"]: bond for bond in fragment.findall("b")}
        for anchor_id, anchor in nodes.items():
            if anchor.get("NodeType") != "VariableAttachment":
                continue
            connector, = (bond for bond in bonds.values() if anchor_id in (bond.get("B"), bond.get("E")))
            label_id, = {connector.attrib["B"], connector.attrib["E"]} - {anchor_id}
            crossed_id, = connector.attrib["CrossingBonds"].split()
            crossed = bonds[crossed_id]
            if connector.attrib["id"] not in crossed.get("CrossingBonds", "").split():
                raise ValueError("Variable-attachment connector lost its reciprocal crossing")

            def point(node_id):
                return list(map(float, nodes[node_id].attrib["p"].split()))

            vertices = [point(node_id) for node_id in ring_vertices(bonds.values(), crossed)]
            center = [round(sum(axis) / len(vertices), 3) for axis in zip(*vertices)]
            edges = [point(crossed.attrib[end]) for end in ("B", "E")]
            midpoint = [round((a + b) / 2, 3) for a, b in zip(*edges)]
            label = "".join(nodes[label_id].find("t").itertext()).strip()
            marks.append({
                "label": label,
                "ring_center_pt": center,
                "anchor_pt": point(anchor_id),
                "label_endpoint_pt": point(label_id),
                "crossed_edge_midpoint_pt": midpoint,
            })
    return marks


def compare(golden: dict, actual: dict, hash_limit: int = HASH_LIMIT) -> tuple[list[str], str]:
    errors = []
    size_changes = []
    for dimension, before, after in zip(("width", "height"), golden["size_px"], actual["size_px"], strict=True):
        if before <= 0:
            raise ValueError(f"Invalid golden {dimension}: {before}")
        delta = abs(after - before)
        size_changes.append(f"{dimension} {before}->{after} px ({100 * delta / before:.2f}%)")
        if 100 * delta > SIZE_PERCENT_LIMIT * before:
            errors.append(f"{size_changes[-1]} exceeds {SIZE_PERCENT_LIMIT}%")

    hash_changes = []
    for key, label in (("perceptual_hash", "hash"), ("header_perceptual_hash", "header hash")):
        if key in golden or key in actual:
            distance = hamming_distance(golden[key], actual[key])
            hash_changes.append(f"{label} distance {distance}/{hash_limit} bits")
            if distance > hash_limit:
                errors.append(f"{label} distance {distance} exceeds {hash_limit} bits")

    expected_marks, actual_marks = golden["marks"], actual["marks"]
    if len(expected_marks) != len(actual_marks):
        errors.append(f"R-group marks {len(expected_marks)}->{len(actual_marks)}")
    for index, (before, after) in enumerate(zip(expected_marks, actual_marks), 1):
        if before["label"] != after["label"]:
            errors.append(f"mark {index} label {before['label']!a}->{after['label']!a}")
        for key in COORDINATES:
            shift = math.dist(before[key], after[key])
            if round(shift, 6) > COORDINATE_LIMIT_PT:
                errors.append(f"mark {index} {key} {before[key]}->{after[key]} (shift {shift:.3f} pt)")
    return errors, ", ".join([*size_changes, *hash_changes])


def draw_case(name: str, output_dir: Path) -> Path:
    if name == "paper":
        from cdxml_toolkit.mcp_runtime.figure_scope import draw_scope

        brief = json.loads(PAPER_BRIEF.read_text(encoding="utf-8"))
        # Use the variable-attachment variant of the two published R1 sites.
        header = brief["header"]
        header["substrates"][0].update(smiles="Nc1ccccc1I", ring_r_groups=[{"number": 1, "atoms": [2, 3, 4, 5]}])
        header["products"][0].update(smiles="O=C(N[*:2])c1nc2ccccc2s1",
                                     ring_r_groups=[{"number": 1, "atoms": [7, 8, 9, 10]}])
        report = draw_scope(brief, output_path=str(output_dir / "paper-scope"), render=False)
        if not report["ok"] or (report.get("readback") or {}).get("status") != "matched":
            raise ValueError(f"Scope drawing or readback failed: {report.get('errors') or report.get('readback')}")
        page, = report["pages"]
        return Path(page["path"])

    from cdxml_toolkit.mcp_runtime.figure_tools import compose_chemical_figure

    atoms = [0, 1, 2, 4, 5] if name == "hetero" else list(range(6))
    spec = {"style": "acs1996", "objects": [{
        "type": "molecule",
        "id": name,
        "smiles": "c1ccncc1" if name == "hetero" else "c1ccccc1",
        "position": [100, 100],
        "rotation": 90,
        "ring_r_groups": [
            {"number": index, "atoms": atoms, "color": "#C81E1E"}
            for index in range(1, CASES[name] + 1)
        ],
    }]}
    spec_path = output_dir / f"{name}.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    source = output_dir / f"{name}.cdxml"
    report = compose_chemical_figure(str(spec_path), str(source))
    if not report["ok"]:
        raise ValueError(f"Composition failed: {report.get('errors') or report}")
    return source


def fingerprint(name: str, output_dir: Path, session) -> dict:
    from PIL import Image

    source = draw_case(name, output_dir)
    png = output_dir / f"{name}.png"
    native_cdxml = output_dir / f"{name}.native.cdxml"
    session.render(source, png, dpi=DPI)
    session.render(source, native_cdxml)
    marks = measure_marks(native_cdxml)
    if len(marks) != CASES[name]:
        raise ValueError(f"Expected {CASES[name]} native R-group marks, found {len(marks)}")
    with Image.open(png) as image:
        return {**image_metrics(name, image), "marks": marks}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="New directory for native exports and fingerprints")
    parser.add_argument("--reference-dir", type=Path, help="Directory of baseline PNGs for side-by-side review")
    parser.add_argument("--update-golden", action="store_true", help="Replace golden JSON after all four exports succeed")
    args = parser.parse_args()

    os.environ.setdefault("CHEMDRAW_OFFLINE", "1")
    from cdxml_toolkit.mcp_runtime.native_renderer import NativeRenderSession
    from cdxml_toolkit.mcp_runtime.resource_lock import native_resource_lock

    if args.output_dir:
        output_dir = args.output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=False)
    else:
        output_root = ROOT / "output" / "visual-guard"
        output_root.mkdir(parents=True, exist_ok=True)
        output_dir = Path(tempfile.mkdtemp(prefix="run-", dir=output_root))
    print(f"Native exports and fingerprints: {output_dir}", flush=True)
    measurements = {}
    failed = 0
    with native_resource_lock("chemdraw_com", 60), NativeRenderSession() as session:
        for name in CASES:
            try:
                actual = fingerprint(name, output_dir, session)
                measurements[name] = json.dumps(actual, indent=2) + "\n"
                (output_dir / f"{name}.fingerprint.json").write_text(measurements[name], encoding="utf-8")
                if args.reference_dir:
                    from scripts.run_scope_benchmark import _stitch_side_by_side

                    reference = args.reference_dir / f"{name}.png"
                    comparison = output_dir / f"{name}.compare.png"
                    if not _stitch_side_by_side(reference, output_dir / f"{name}.png", comparison,
                                                name, "Reference ChemDraw export", "Fresh ChemDraw export"):
                        raise ValueError(f"Cannot compare PNGs with {reference}")
                if args.update_golden:
                    print(f"[RECORDED] {name} {actual['size_px']} px, {len(actual['marks'])} marks")
                    continue
                golden = json.loads((GOLDEN_DIR / f"{name}.json").read_text(encoding="utf-8"))
                errors, summary = compare(golden, actual)
                print(f"[{'FAIL' if errors else 'PASS'}] {name}: {summary}")
                for error in errors:
                    print(f"    {error}")
                failed += bool(errors)
            except (ImportError, OSError, ET.ParseError, KeyError, ValueError, RuntimeError) as exc:
                failed += 1
                print(f"[FAIL] {name}: {exc}")
    if args.update_golden and not failed:
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        for name, measurement in measurements.items():
            (GOLDEN_DIR / f"{name}.json").write_text(measurement, encoding="utf-8")
    print(f"{len(CASES) - failed}/{len(CASES)} cases {'recorded' if args.update_golden else 'passed'}")
    return int(bool(failed))


if __name__ == "__main__":
    raise SystemExit(main())
