"""House-rule check of a saved CDXML, without ChemDraw or RDKit.

Checks what a file can be held to mechanically: document settings (ACS Document 1996 by
default), fonts (Arial; CJK text in SimSun), text sizes, the colour palette, measured bond
lengths, broken id references and constructs ChemDraw is known to drop or render blank.
It does not judge the chemistry: a wrong ring size, a misplaced electron or two isomers
drawn alike pass this check and must be caught by looking at the native render.
"""
from __future__ import annotations

import math
import re
import statistics
from pathlib import Path
from typing import Any, Iterable, Mapping
from xml.etree import ElementTree as ET

from .cdxml_compose import ID_REFERENCES

PROFILES: dict[str, dict[str, Any]] = {
    'acs1996': {
        'settings': {'BondLength': 14.4, 'LineWidth': 0.6, 'BoldWidth': 2.0, 'HashSpacing': 2.5,
                     'MarginWidth': 1.6, 'LabelSize': 10.0, 'CaptionSize': 10.0, 'ChainAngle': 120.0},
        'latin_fonts': ('Arial',), 'cjk_fonts': ('SimSun',), 'text_size': 10.0, 'min_text_size': 6.0,
        'palette': ('#000000', '#FFFFFF', '#FF6600', '#336699'),
    },
}
_CJK = re.compile(r'[\u3000-\u303f\u3400-\u9fff\uff00-\uffef]')
_CJK_ONLY_PUNCT = re.compile(r'[\u3000-\u303f\uff00-\uffef]')
_COLOR_ATTRS = ('color', 'bgcolor', 'fgcolor')
_MARKUP_LEFT = re.compile(r'\{(i|b|bi|sub|sup):')
# Page size used when writing WidthPages/HeightPages (with room to spare), and the width per page
# beyond which a native render has been seen to wrap text (text at x = 537 pt on one page did not).
PAGE_WIDTH, PAGE_HEIGHT = 500.0, 650.0
WRAP_WIDTH = 600.0


def pages_needed(width, height):
    """WidthPages and HeightPages that keep a drawing of this size from wrapping its text."""
    import math
    return max(2, math.ceil(width / PAGE_WIDTH) + 1), max(2, math.ceil(height / PAGE_HEIGHT) + 1)


def _issue(bucket, code, message, **extra):
    bucket.append({'code': code, 'message': message, **extra})


def _hex(rgb):
    return '#%02X%02X%02X' % tuple(max(0, min(255, round(v * 255))) for v in rgb)


def normalize_palette(palette: Iterable[str] | None):
    if palette is None:
        return None
    out = []
    for item in palette:
        text = str(item).strip()
        if re.fullmatch(r'#?[0-9a-fA-F]{6}', text):
            out.append('#' + text.lstrip('#').upper())
        else:
            raise ValueError(f'A palette colour is #RRGGBB, not {item!r}')
    return tuple(out)


def _run_markup(run) -> str:
    face = int(run.get('face', '0') or 0) & 96
    # 96 is ChemDraw's "formula" face (digits lowered), read here as plain text.
    return ('{sub:%s}' if face == 32 else '{sup:%s}' if face == 64 else '%s') % (run.text or '')


def _free_text_boxes(root) -> list[tuple[str, tuple[float, float, float, float]]]:
    """Estimated boxes of captions that are not atom labels (ChemDraw computes the real ones)."""
    from .scientific.canvas import caption_box
    # Atom labels and the text of native objects (TLC plates, stoichiometry grids) place themselves.
    atom_text = {id(t) for owner in ('fragment', 'objecttag') for item in root.iter(owner) for t in item.iter('t')}
    boxes = []
    for t in root.iter('t'):
        if id(t) in atom_text or not t.get('p'):
            continue
        runs = list(t.iter('s'))
        markup = ''.join(_run_markup(r) for r in runs)
        if not markup.strip():
            continue
        if t.get('BoundingBox'):  # ChemDraw wrote its own measurement when it saved the file
            v = [float(c) for c in t.get('BoundingBox').split()[:4]]
            boxes.append((markup, (min(v[0], v[2]), min(v[1], v[3]), max(v[0], v[2]), max(v[1], v[3]))))
            continue
        x, y = (float(v) for v in t.get('p').split()[:2])
        if x == 0 and y == 0:  # a placeholder position (table cells are placed by their table)
            continue
        size = max((float(r.get('size', 10) or 10) for r in runs), default=10.0)
        boxes.append((markup, caption_box(x, y, markup, size, t.get('Justification', 'Left'))))
    return boxes


def _white_indices(root) -> set[str]:
    """Colour indices that are white: 1, and colour-table entries (numbered from 2) of 1 1 1."""
    white = {'1'}
    for k, color in enumerate(root.findall('./colortable/color'), start=2):
        if all(abs(float(color.get(c, 0)) - 1.0) < 1e-6 for c in 'rgb'):
            white.add(str(k))
    return white


def _strokes(root) -> list[tuple[tuple[float, float], tuple[float, float], int]]:
    """Straight pieces of the drawn lines, rectangles, straight arrows and curves, each with
    the position of its element in drawing order."""
    order = {id(el): k for k, el in enumerate(root.iter())}
    pieces = []
    for graphic in root.iter('graphic'):
        if not graphic.get('BoundingBox'):
            continue
        k = order[id(graphic)]
        v = [float(x) for x in graphic.get('BoundingBox').split()[:4]]
        if graphic.get('GraphicType') == 'Line':
            pieces.append(((v[0], v[1]), (v[2], v[3]), k))
        elif graphic.get('GraphicType') == 'Rectangle':
            corners = [(v[0], v[1]), (v[2], v[1]), (v[2], v[3]), (v[0], v[3])]
            pieces += [(a, b, k) for a, b in zip(corners, corners[1:] + corners[:1])]
    for arrow in root.iter('arrow'):
        if arrow.get('Head3D') and arrow.get('Tail3D') and not arrow.get('Center3D'):
            head = tuple(float(x) for x in arrow.get('Head3D').split()[:2])
            tail = tuple(float(x) for x in arrow.get('Tail3D').split()[:2])
            pieces.append((tail, head, order[id(arrow)]))
    for curve in root.iter('curve'):
        k = order[id(curve)]
        v = [float(x) for x in curve.get('CurvePoints', '').split()]
        pts = list(zip(v[0::2], v[1::2]))
        # CurvePoints are (in-handle, anchor, out-handle) triples: one cubic Bezier per pair of anchors.
        triples = [pts[i:i + 3] for i in range(0, len(pts) - 2, 3)]
        for (_, a, out), (into, b, _) in zip(triples, triples[1:]):
            samples = [tuple((1 - u) ** 3 * a[i] + 3 * (1 - u) ** 2 * u * out[i] + 3 * (1 - u) * u * u * into[i]
                             + u ** 3 * b[i] for i in (0, 1)) for u in (n / 8 for n in range(9))]
            pieces += [(p, q, k) for p, q in zip(samples, samples[1:])]
    return pieces


def _knockouts(root) -> list[tuple[tuple[float, float, float, float], int]]:
    """White filled patches (box, drawing order): they hide whatever was drawn before them."""
    order = {id(el): k for k, el in enumerate(root.iter())}
    white = _white_indices(root)
    patches = []
    for el in root.iter():
        if el.tag not in ('curve', 'graphic') or el.get('FillType') != 'Solid' or el.get('color', '0') not in white:
            continue
        if el.tag == 'curve':
            v = [float(x) for x in el.get('CurvePoints', '').split()]
            xs, ys = v[0::2], v[1::2]
            if not xs:
                continue
            box = (min(xs), min(ys), max(xs), max(ys))
        elif el.get('BoundingBox'):
            v = [float(x) for x in el.get('BoundingBox').split()[:4]]
            box = (min(v[0], v[2]), min(v[1], v[3]), max(v[0], v[2]), max(v[1], v[3]))
        else:
            continue
        patches.append((box, order[id(el)]))
    return patches


def _crosses(a, b, box) -> bool:
    """Whether the segment a-b enters ``box`` (Liang-Barsky)."""
    x0, y0, x1, y1 = box
    if x1 <= x0 or y1 <= y0:
        return False
    dx, dy = b[0] - a[0], b[1] - a[1]
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, a[0] - x0), (dx, x1 - a[0]), (-dy, a[1] - y0), (dy, y1 - a[1])):
        if p == 0:
            if q < 0:
                return False
            continue
        r = q / p
        if p < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return False
    return True


def text_over_ink(source, inset: float = 1.0) -> list[str]:
    """Captions (not atom labels) that a line, arrow, curve or another caption runs through.

    A line drawn before a white filled patch that covers the caption does not count. Text
    boxes are estimated from the font metrics ChemDraw uses for Arial and SimSun, so a
    graze within ``inset`` points of a box edge is not counted.
    """
    if hasattr(source, 'iter'):
        root = source
    else:
        text = str(source)
        if not text.lstrip().startswith('<'):
            text = Path(text).read_text(encoding='utf-8-sig')
        body = re.sub(r'^\s*<\?xml[^>]*\?>', '', text, count=1)
        root = ET.fromstring(re.sub(r'<!DOCTYPE[^>]*>', '', body, count=1).encode('utf-8'))
    boxes = _free_text_boxes(root)
    pieces = _strokes(root)
    patches = _knockouts(root)
    found = []
    for markup, (x0, y0, x1, y1) in boxes:
        inner = (x0 + inset, y0 + inset, x1 - inset, y1 - inset)
        # A stroke drawn before a white patch that covers the caption is hidden by it.
        covers = [k for (px0, py0, px1, py1), k in patches
                  if px0 <= x0 + inset and py0 <= y0 + inset and px1 >= x1 - inset and py1 >= y1 - inset]
        hidden = max(covers, default=-1)
        if any(order > hidden and _crosses(a, b, inner) for a, b, order in pieces):
            found.append(f'a line runs through {markup!r}')
    # Two captions' boxes include room for sub- and superscripts that seldom meet the other
    # line's glyphs, so tight leading (lines 1.04 sizes apart) is allowed.
    for i, (ma, a) in enumerate(boxes):
        for mb, b in boxes[i + 1:]:
            if min(a[2], b[2]) - max(a[0], b[0]) > inset and min(a[3], b[3]) - max(a[1], b[1]) > 1.5 * inset:
                found.append(f'{ma!r} overlaps {mb!r}')
    return found


def check_cdxml(source, *, profile: str = 'acs1996', palette: Iterable[str] | None | str = 'profile',
                allow_sizes: Iterable[float] = (), bond_tolerance: float = 0.03,
                ignore: Iterable[str] = ()) -> dict[str, Any]:
    """Check one CDXML file or text. Returns ``{"ok", "errors", "warnings", "summary"}``.

    ``palette`` is ``"profile"`` (the profile's colours), ``None`` (any colour) or a list of
    ``#RRGGBB``. ``allow_sizes`` lists text sizes accepted besides the profile's. ``ignore``
    drops issue codes (``bond_length`` for a Newman projection, whose template draws 24 pt
    bonds on purpose). Errors are things ChemDraw would drop, misdraw or that break the house
    rules; warnings need a look.
    """
    if profile not in PROFILES:
        raise ValueError(f'Unknown profile {profile!r}; available: {", ".join(PROFILES)}')
    rules = PROFILES[profile]
    allowed = rules['palette'] if palette == 'profile' else normalize_palette(palette)
    errors: list[dict] = []
    warnings: list[dict] = []
    text = str(source)
    path = None
    if not text.lstrip().startswith('<'):
        path = Path(text)
        text = path.read_text(encoding='utf-8-sig')
    declaration = re.match(r'\s*<\?xml[^>]*encoding=["\']([^"\']+)["\']', text)
    if declaration and declaration.group(1).lower() == 'utf8':
        _issue(errors, 'xml_encoding', 'The XML declaration says "utf8"; ChemDraw needs "UTF-8" and fails to open the file')
    try:
        body = re.sub(r'^\s*<\?xml[^>]*\?>', '', text, count=1)      # the declaration was read above
        root = ET.fromstring(re.sub(r'<!DOCTYPE[^>]*>', '', body, count=1).encode('utf-8'))
    except ET.ParseError as exc:
        return {'ok': False, 'file': str(path) if path else None, 'profile': profile,
                'errors': [{'code': 'not_xml', 'message': f'Not well-formed XML: {exc}'}], 'warnings': [], 'summary': {}}
    if root.tag != 'CDXML':
        return {'ok': False, 'file': str(path) if path else None, 'profile': profile,
                'errors': [{'code': 'not_cdxml', 'message': f'The root element is <{root.tag}>, not <CDXML>'}],
                'warnings': [], 'summary': {}}

    # ---- document settings
    for name, expected in rules['settings'].items():
        value = root.get(name)
        if value is None:
            _issue(warnings, 'setting_missing', f'{name} is not set (expected {expected:g})', setting=name)
        elif abs(float(value) - expected) > 1e-6:
            _issue(errors, 'setting', f'{name} is {value}, expected {expected:g}', setting=name, found=value)

    # ---- fonts and text
    fonts = {}
    fonttable = root.find('fonttable')
    if fonttable is not None:
        fonts = {font.get('id'): font.get('name') for font in fonttable.findall('font')}
    for attr in ('LabelFont', 'CaptionFont'):
        value = root.get(attr)
        if value is not None and fonts.get(value) not in rules['latin_fonts']:
            _issue(errors, 'font', f'{attr} is {fonts.get(value, "font " + value)}, expected {rules["latin_fonts"][0]}')
    sizes: dict[float, int] = {}
    default_font = root.get('CaptionFont') or root.get('LabelFont')
    for run in root.iter('s'):
        content = run.text or ''
        if not content.strip():
            continue
        font = fonts.get(run.get('font', default_font))
        has_cjk = bool(_CJK.search(content))
        latin = _CJK.sub('', content).strip()
        if _MARKUP_LEFT.search(content):
            _issue(errors, 'unresolved_markup', f'"{content[:20]}" still contains markup such as {{sub:..}}; it would '
                                                'print literally')
        if has_cjk and re.search(r'[A-Za-z0-9]', latin):
            _issue(errors, 'mixed_script_run', f'"{content[:12]}" mixes Chinese and Latin in one run; split it so the '
                                               'Chinese is SimSun and the rest Arial')
        if has_cjk and font not in rules['cjk_fonts']:
            _issue(errors, 'cjk_font', f'Chinese text "{content[:12]}" is in {font or "an undefined font"}; it must be '
                                       f'{rules["cjk_fonts"][0]} or ChemDraw draws boxes')
        elif latin and font in rules['cjk_fonts'] and not has_cjk:
            _issue(warnings, 'latin_in_cjk_font', f'"{content[:12]}" is set in {font}; Latin text belongs in '
                                                  f'{rules["latin_fonts"][0]}')
        elif not has_cjk and font is not None and font not in rules['latin_fonts']:
            _issue(errors, 'font', f'"{content[:12]}" is in {font}, expected {rules["latin_fonts"][0]}')
        elif font is None:
            _issue(errors, 'font_undefined', f'"{content[:12]}" uses a font id that is not in the font table')
        size = float(run.get('size', root.get('CaptionSize', rules['text_size'])))
        sizes[size] = sizes.get(size, 0) + 1
    accepted = {rules['text_size'], *(float(v) for v in allow_sizes)}
    for size, count in sorted(sizes.items()):
        if size < rules['min_text_size']:
            _issue(errors, 'text_too_small', f'{count} text run(s) at {size:g} pt are below the {rules["min_text_size"]:g} pt minimum')
        elif size not in accepted:
            _issue(warnings, 'text_size', f'{count} text run(s) at {size:g} pt (the house size is {rules["text_size"]:g} pt)',
                   size=size)

    # ---- colours
    table = {0: '#000000', 1: '#FFFFFF'}
    colortable = root.find('colortable')
    if colortable is not None:
        for i, color in enumerate(colortable.findall('color'), start=2):
            table[i] = _hex([float(color.get(k, '0')) for k in 'rgb'])
    used: dict[str, int] = {}
    for node in root.iter():
        if node.tag in ('CDXML', 'color'):
            continue
        for attr in _COLOR_ATTRS:
            value = node.get(attr)
            if value is None or not value.isdigit():
                continue
            index = int(value)
            if index not in table:
                _issue(errors, 'color_index', f'<{node.tag}> uses colour index {index}, which is not in the colour table')
                continue
            used[table[index]] = used.get(table[index], 0) + 1
    if allowed is not None:
        for color, count in sorted(used.items()):
            if color not in allowed:
                _issue(errors, 'palette', f'{color} is used {count} time(s) but the palette is {", ".join(allowed)}',
                       color=color)

    # ---- ids and references
    ids: dict[str, int] = {}
    for node in root.iter():
        if node.get('id') is not None and node.tag != 'font':      # font ids are their own namespace
            ids[node.get('id')] = ids.get(node.get('id'), 0) + 1
    for value, count in ids.items():
        if count > 1 and value != '0':
            _issue(errors, 'duplicate_id', f'id {value} is used by {count} objects')
    for node in root.iter():
        for attr in ID_REFERENCES:
            value = node.get(attr)
            if value is None:
                continue
            missing = [token for token in value.split() if token != '0' and token not in ids]
            if missing:
                _issue(errors, 'dangling_reference', f'<{node.tag} id="{node.get("id")}"> {attr} refers to missing '
                                                     f'object(s) {" ".join(missing)}')

    # ---- constructs ChemDraw drops or renders blank
    for curve in root.iter('curve'):
        values = (curve.get('CurvePoints') or '').split()
        if len(values) % 6 or len(values) < 12:
            _issue(errors, 'curve_points', f'<curve id="{curve.get("id")}"> has {len(values) // 2} points; ChemDraw '
                                           'needs whole (handle, point, handle) triples and silently drops other curves')
            continue
        # A curve that states Closed or FillType but not CurveType is dropped from the native render
        # (curves saved by ChemDraw itself carry none of the three and render fine).
        if curve.get('CurveType') is None and (curve.get('Closed') is not None or curve.get('FillType') is not None):
            _issue(errors, 'curve_type_missing', f'<curve id="{curve.get("id")}"> sets Closed/FillType without '
                                                 'CurveType; ChemDraw drops such a curve from the native render')
    for arrow in root.iter('arrow'):
        try:
            head = float(arrow.get('HeadSize', '1000'))
            centre = float(arrow.get('ArrowheadCenterSize', head * 0.875))
            width = float(arrow.get('ArrowheadWidth', head * 0.25))
        except ValueError:
            continue
        if abs(centre - 0.875 * head) > 1.5 or abs(width - 0.25 * head) > 1.5:
            _issue(warnings, 'arrowhead_ratio', f'<arrow id="{arrow.get("id")}"> head {head:g}/{centre:g}/{width:g} '
                                                'is not the house 1 : 0.875 : 0.25 shape')
    for line_width in {n.get('LineWidth') for n in root.iter() if n.tag != 'CDXML' and n.get('LineWidth')}:
        if abs(float(line_width) - rules['settings']['LineWidth']) > 1e-6:
            _issue(warnings, 'line_width', f'Some objects set LineWidth {line_width} (the house width is '
                                           f'{rules["settings"]["LineWidth"]:g})')
    lengths = []
    fragments = 0
    for fragment in root.iter('fragment'):
        fragments += 1
        nodes = {n.get('id'): n for n in fragment.findall('n')}
        bonds = fragment.findall('b')
        positions = {}
        for nid, node in nodes.items():
            if node.get('p'):
                positions[nid] = tuple(float(v) for v in node.get('p').split()[:2])
        for bond in bonds:
            a, b = positions.get(bond.get('B')), positions.get(bond.get('E'))
            if a and b and nodes[bond.get('B')].get('NodeType') is None and nodes[bond.get('E')].get('NodeType') is None:
                lengths.append(math.dist(a, b))
        visible = bool(bonds) or any(n.find('t') is not None or n.get('Element') not in (None, '6')
                                     or n.get('NodeType') for n in nodes.values())
        if nodes and not visible:
            _issue(errors, 'blank_structure', f'<fragment id="{fragment.get("id")}"> is a lone carbon with no label '
                                              'and no bonds: it renders as nothing; write the hydrogens explicitly')
        elif len(nodes) == 2 and len(bonds) == 1:
            bare = [n for n in nodes.values() if n.find('t') is None and n.get('Element') in (None, '6')
                    and not n.get('NodeType')]
            if len(bare) == 1:
                _issue(warnings, 'bare_terminal_carbon',
                       f'<fragment id="{fragment.get("id")}"> is a two-atom molecule whose carbon has no label: it '
                       'renders as a bond end with nothing on it; label the carbon (CH3, H2C)')
    # ---- the page: coordinates outside it, and drawings wider than its pages (text wraps)
    page = root.find('page')
    if page is not None and page.get('Width') and page.get('Height'):
        width, height = float(page.get('Width')), float(page.get('Height'))
        outside = 0
        right = 0.0
        for node in root.iter():
            for attr in ('p', 'Head3D', 'Tail3D'):
                if node.get(attr):
                    x, y = (float(v) for v in node.get(attr).split()[:2])
                    right = max(right, x)
                    if not (-1 <= x <= width + 1 and -1 <= y <= height + 1):
                        outside += 1
        if outside:
            _issue(warnings, 'outside_page', f'{outside} object position(s) lie outside the {width:.0f} x {height:.0f} pt page')
        pages = float(page.get('WidthPages', '1'))
        if right > pages * WRAP_WIDTH:
            need = pages_needed(width, height)
            _issue(warnings, 'text_may_wrap', f'the drawing reaches x = {right:.0f} pt but the page is {pages:g} page(s) '
                                              f'wide; ChemDraw may wrap text there (set WidthPages="{need[0]}")')
    crowded = text_over_ink(root)
    if crowded:
        _issue(warnings, 'text_over_ink', f'{len(crowded)} caption(s) sit on other ink: ' + '; '.join(crowded[:5]),
               cases=crowded)
    summary: dict[str, Any] = {'fragments': fragments, 'colors': used, 'text_sizes': {f'{k:g}': v for k, v in sizes.items()}}
    if lengths:
        median = statistics.median(lengths)
        summary['median_bond_length'] = round(median, 2)
        expected = rules['settings']['BondLength']
        if abs(median - expected) / expected > bond_tolerance:
            _issue(errors, 'bond_length', f'The median drawn bond is {median:.2f} pt, expected {expected:g} pt')
        outliers = sum(1 for length in lengths if abs(length - median) / median > 0.25)
        if outliers:
            _issue(warnings, 'bond_length_spread', f'{outliers} of {len(lengths)} bonds differ from the median by more '
                                                   'than 25% (normal for chairs, cages and Newman projections)')
    skipped = {str(code) for code in ignore}
    errors = [e for e in errors if e['code'] not in skipped]
    warnings = [w for w in warnings if w['code'] not in skipped]
    return {'ok': not errors, 'file': str(path) if path else None, 'profile': profile, 'errors': errors,
            'warnings': warnings, 'summary': summary,
            'not_checked': 'chemical correctness (ring sizes, electron counts, stereo, identity): inspect the native render'}
