"""Frost constructions checked from saved CDXML geometry, without a renderer."""
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

from cdxml_toolkit.cdxml_check import check_cdxml
from cdxml_toolkit.scientific.diagrams import build

ROOT = Path(__file__).resolve().parents[2]


def drawing(spec):
    result = build('frost_circle', spec)
    return result, ET.fromstring(result['cdxml'])


def circle_geometry(root):
    circle = next(g for g in root.iter('graphic') if g.get('OvalType') == 'Circle')
    cx, cy, _ = map(float, circle.get('Center3D').split())
    radius = float(circle.get('MajorAxisEnd3D').split()[0]) - cx
    return cx, cy, radius


def level_lines(root):
    cx, _, radius = circle_geometry(root)
    lines = []
    for node in root.iter('graphic'):
        if node.get('GraphicType') != 'Line' or node.get('LineType') == 'Dashed':
            continue
        x1, y1, x2, y2 = map(float, node.get('BoundingBox').split())
        if min(x1, x2) > cx + radius + 20 and abs(y1 - y2) < 0.01:
            lines.append((y1, min(x1, x2), max(x1, x2)))
    return sorted(lines, reverse=True)


def captions(root):
    return [''.join(t.itertext()) for t in root.iter('t')]


@pytest.mark.parametrize('n,expected,degeneracies', [
    (3, [2, -1], [1, 2]),
    (4, [2, 0, -2], [1, 2, 1]),
    (5, [2, (math.sqrt(5) - 1) / 2, -(math.sqrt(5) + 1) / 2], [1, 2, 2]),
    (6, [2, 1, -1, -2], [1, 2, 2, 1]),
])
def test_saved_levels_follow_the_huckel_formula(n, expected, degeneracies):
    result, root = drawing({'centres': n, 'electrons': n, 'show_construction': True})
    cx, cy, radius = circle_geometry(root)
    rows = level_lines(root)
    assert len(rows) == len(expected)
    theory = sorted({round(2 * math.cos(2 * math.pi * k / n), 10) for k in range(n)}, reverse=True)
    assert theory == pytest.approx(expected)
    for (y, _, _), coefficient in zip(rows, theory):
        assert y == pytest.approx(cy + radius * coefficient / 2, abs=0.5)
    assert [level['degeneracy'] for level in result['levels']] == degeneracies
    for coefficient, degeneracy in zip(expected, degeneracies):
        assert any(f'x = {coefficient:.6f}' in text and f'd = {degeneracy}' in text
                   for text in captions(root))
    polygon = next(root.iter('curve'))
    values = list(map(float, polygon.get('CurvePoints').split()))
    vertices = set(zip(values[::2], values[1::2]))
    assert len(vertices) == n
    assert any(abs(x - cx) < 0.02 and abs(y - cy - radius) < 0.02 for x, y in vertices)
    assert all(math.hypot(x - cx, y - cy) == pytest.approx(radius, abs=0.02) for x, y in vertices)
    assert not list(root.iter('fragment'))


@pytest.mark.parametrize('n,m,expected', [
    (3, 2, ['ud', '', '']), (4, 4, ['ud', 'u', 'u', '']),
    (5, 4, ['ud', 'u', 'u', '', '']), (5, 6, ['ud', 'ud', 'ud', '', '']),
    (6, 6, ['ud', 'ud', 'ud', '', '', '']), (6, 8, ['ud', 'ud', 'ud', 'u', 'u', '']),
    (4, 0, ['', '', '', '']), (4, 8, ['ud', 'ud', 'ud', 'ud']),
])
def test_actual_spin_arrows_obey_aufbau_hund_and_pauli(n, m, expected):
    result, root = drawing({'centres': n, 'electrons': m})
    actual = []
    for y, left, right in level_lines(root):
        arrows = []
        for arrow in root.iter('arrow'):
            if arrow.get('ArrowheadHead') != 'HalfLeft':
                continue
            hx, hy, _ = map(float, arrow.get('Head3D').split())
            tx, ty, _ = map(float, arrow.get('Tail3D').split())
            if left < hx < right and abs((hy + ty) / 2 - (y - 2)) < 0.05:
                arrows.append((hx, 'u' if hy < ty else 'd'))
        degeneracy = next(level['degeneracy'] for level in result['levels']
                          if abs(level['x'] - 2 * (y - circle_geometry(root)[1]) / circle_geometry(root)[2]) < 0.001)
        for slot in range(degeneracy):
            start = left + (right - left) * slot / degeneracy
            end = left + (right - left) * (slot + 1) / degeneracy
            actual.append(''.join(spin for x, spin in sorted(arrows) if start < x < end))
    assert actual == expected
    assert sum(map(len, actual)) == m
    assert result['unpaired_electrons'] == expected.count('u')
    expected_rule = '4q+2' if m >= 2 and m % 4 == 2 else '4q' if m > 0 and m % 4 == 0 else 'other'
    assert result['huckel_rule'] == expected_rule


def test_construction_helpers_and_custom_rules_are_saved():
    lines = ['Rule one', 'Rule two', 'Rule three']
    _, root = drawing({'centres': 5, 'electrons': 6, 'lines': lines})
    assert all(line in captions(root) for line in lines)
    assert 'r = 2|\u03b2|' in captions(root)
    assert len([g for g in root.iter('graphic') if g.get('LineType') == 'Dashed']) == 5
    assert any(a.get('ArrowheadHead') == 'Full' and a.get('ArrowheadTail') == 'Full' for a in root.iter('arrow'))
    assert any(g.get('OvalType') == 'Circle Filled' for g in root.iter('graphic'))
    _, plain = drawing({'centres': 5, 'electrons': 6, 'show_construction': False})
    assert not [g for g in plain.iter('graphic') if g.get('LineType') == 'Dashed']
    assert not [g for g in plain.iter('graphic') if g.get('OvalType') == 'Circle Filled']
    assert 'r = 2|\u03b2|' not in captions(plain)


@pytest.mark.parametrize('variant', ['vertex-up', 'flat-bottom'])
def test_incorrect_constructions_are_geometric_and_marked_red(variant):
    result, root = drawing({'centres': 5, 'electrons': 6, 'variant': variant})
    assert not result['valid_construction'] and result['warnings']
    assert result['huckel_rule'] is None
    assert any('Incorrect' in text for text in captions(root))
    colors = {str(i): tuple(float(c.get(key)) for key in ('r', 'g', 'b'))
              for i, c in enumerate(root.find('colortable'), 2)}
    red = [g for g in root.iter('graphic') if colors.get(g.get('color')) == (1, 0, 0)]
    assert len(red) == 2
    assert all(g.get('GraphicType') == 'Line' for g in red)
    _, cy, radius = circle_geometry(root)
    points = list(map(float, next(root.iter('curve')).get('CurvePoints').split()))
    ys = sorted(set(round(y, 2) for y in points[1::2]))
    assert ys[-1] < cy + radius - 0.5
    if variant == 'vertex-up':
        assert ys[0] == pytest.approx(cy - radius, abs=0.02)


@pytest.mark.parametrize('spec', [
    {}, {'centres': 2, 'electrons': 2}, {'centres': 4.5, 'electrons': 4},
    {'centres': True, 'electrons': 2}, {'centres': 4, 'electrons': -1},
    {'centres': 4, 'electrons': 9}, {'centres': 4, 'electrons': 2.2},
    {'centres': 4, 'electrons': False}, {'centres': 4, 'electrons': 4, 'variant': 'other'},
    {'centres': 4, 'electrons': 4, 'show_construction': 'yes'},
    {'centres': 5, 'electrons': 6, 'lines': 'rules'},
    {'centres': 5, 'electrons': 6, 'lines': [1]},
    {'centres': 6, 'electrons': 6, 'variant': 'vertex-up'},
])
def test_invalid_specs_are_rejected(spec):
    with pytest.raises(ValueError):
        build('frost_circle', spec)


def test_cli_mcp_and_redox_remain_distinct(tmp_path):
    spec = tmp_path / 'spec.json'
    spec.write_text(json.dumps({'centres': 6, 'electrons': 6}), encoding='utf-8')
    output = tmp_path / 'circle.cdxml'
    cli = ROOT / 'skill/chemdraw/scripts/cdx.py'
    process = subprocess.run([sys.executable, str(cli), 'scientific', 'frost_circle', str(spec), str(output)],
                             capture_output=True, text=True, encoding='utf-8', timeout=60)
    assert process.returncode == 0, process.stdout + process.stderr
    assert json.loads(process.stdout)['diagram'] == 'frost_circle'
    assert len(level_lines(ET.parse(output).getroot())) == 4
    help_text = subprocess.run([sys.executable, str(cli), 'scientific', '--help'],
                               capture_output=True, text=True, encoding='utf-8', timeout=60).stdout
    assert 'redox' in help_text and 'frost_circle' in help_text
    from cdxml_toolkit.mcp_runtime.scientific_tool import scientific_workflow
    result = scientific_workflow('frost_circle', str(spec), str(tmp_path / 'mcp.cdxml'))
    assert result['ok'] and result['outputs']['result']['diagram'] == 'frost_circle'
    with pytest.raises(ValueError, match='overwrite|exists'):
        scientific_workflow('frost_circle', str(spec), str(output))
    redox = build('frost', {'species': [{'name': 'A', 'oxidation_state': 1},
                                      {'name': 'B', 'oxidation_state': 0}], 'potentials': [1]})
    assert redox['diagram'] == 'frost' and redox['series'][0]['points'][0]['NE'] == 1


def test_bundled_examples_check_and_rebuild_exactly():
    examples = ROOT / 'toolkit/examples/frost_circle'
    specs = sorted(examples.glob('n*.json')) + sorted(examples.glob('incorrect-*.json'))
    assert len(specs) == 6
    for path in specs:
        result = build('frost_circle', json.loads(path.read_text(encoding='utf-8')))
        saved = path.with_suffix('.cdxml').read_text(encoding='utf-8')
        assert saved == result['cdxml']
        check = check_cdxml(saved, palette=None)
        assert check['ok'], check['errors']
        assert not check['warnings'], check['warnings']
    for name in ('comparison', 'construction-rules'):
        saved = (examples / f'{name}.cdxml').read_text(encoding='utf-8')
        check = check_cdxml(saved, palette=None)
        assert check['ok'] and not check['warnings'], check
