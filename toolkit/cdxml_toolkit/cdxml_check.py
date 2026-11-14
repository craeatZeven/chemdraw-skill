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
