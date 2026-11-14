"""Place a conformer structure at an energy-profile stationary point.

`energy_profile.py` draws curves, levels, transition-state marks and Ea/dH
arrows but knows nothing about molecules.  This module adds the missing half:
it turns a compact `structure` block in a species entry into a CDXML fragment
the parent serializer can drop on the page.

Spec (per species):

    "structure": {"template": "chair", "flip": false,
                  "substituents": {"1": {"ax": "CH3"}},
                  "dy": -10}

`template` is a `cdxml_toolkit.templates.conformers` generator: chair, newman,
haworth, fischer, norbornane.  Everything else is that generator's own keyword
arguments, so the parameters and their real names are the toolkit's, not a
second invention here.  `dy` is an extra offset in points (screen y is down, so
a negative dy lifts the structure).

The returned fragment is centred on (0, 0) and its ids are remapped into a
reserved block, so several structures and the parent's own objects cannot
collide.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

# ids below this are the parent's; structures take 1_000_000 upwards
ID_BLOCK = 1_000_000
ID_STRIDE = 10_000

_GENERATORS = {
    'chair': 'chair_cyclohexane',
    'chair_cyclohexane': 'chair_cyclohexane',
    'newman': 'newman_projection',
    'newman_projection': 'newman_projection',
    'haworth': 'haworth_pyranose',
    'haworth_pyranose': 'haworth_pyranose',
    'fischer': 'fischer_projection',
    'fischer_projection': 'fischer_projection',
    'norbornane': 'bicyclic_cage',
    'bicyclo_222': 'bicyclic_cage',
    'adamantane': 'bicyclic_cage',
}

# keys this module consumes itself; everything else goes to the generator
_OWN_KEYS = {'template', 'dx', 'dy', 'scaffold'}


class StructureError(ValueError):
    pass


def _generate(spec: dict) -> dict:
    try:
        from cdxml_toolkit.templates import conformers
    except Exception as exc:  # pragma: no cover - runtime problem, not a spec problem
        raise StructureError(
            'cdxml_toolkit is not importable here (%s); run through the skill runtime' % exc)

    template = str(spec.get('template') or 'chair').strip().lower()
    fn_name = _GENERATORS.get(template)
    if not fn_name:
        raise StructureError('unknown structure template %r (known: %s)'
                             % (template, ', '.join(sorted(set(_GENERATORS)))))
    fn = getattr(conformers, fn_name, None)
    if fn is None:  # pragma: no cover
        raise StructureError('the toolkit has no %s() generator' % fn_name)

    kwargs = {k: v for k, v in spec.items() if k not in _OWN_KEYS and not k.startswith('_')}
    if template in ('norbornane', 'bicyclo_222', 'adamantane') and 'scaffold' not in kwargs:
        kwargs['scaffold'] = 'norbornane' if template == 'norbornane' else template
    kwargs.setdefault('cx', 0.0)
    kwargs.setdefault('cy', 0.0)

    try:
        return fn(**kwargs)
    except TypeError as exc:
        raise StructureError('%s rejected the given parameters: %s' % (fn_name, exc))


def _all_points(root: ET.Element):
    xs, ys = [], []
    for el in root.iter():
        p = el.get('p')
        if p:
            parts = p.split()
            if len(parts) >= 2:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
        bb = el.get('BoundingBox')
        if bb:
            vals = [float(v) for v in bb.split()]
            if len(vals) == 4:
                xs += [vals[0], vals[2]]
                ys += [vals[1], vals[3]]
    return xs, ys


def make(spec: dict, index: int = 0) -> dict:
    """Generate one structure, centred on the origin, with fresh ids.

    Returns {'root': Element, 'half_w': float, 'half_h': float, 'template': str}.
    """
    if not isinstance(spec, dict):
        raise StructureError('structure must be an object')
    res = _generate(spec)
    xml = res.get('fragment_xml')
    if not xml:
        raise StructureError('the toolkit returned no fragment for this structure')

    root = ET.fromstring(xml)
    xs, ys = _all_points(root)
    if not xs:
        raise StructureError('the generated structure has no coordinates')

    cx, cy = (min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0
    half_w, half_h = (max(xs) - min(xs)) / 2.0, (max(ys) - min(ys)) / 2.0

    base = ID_BLOCK + index * ID_STRIDE
    idmap: dict[str, str] = {}
    for el in root.iter():
        old = el.get('id')
        if old and old not in idmap:
            idmap[old] = str(base + len(idmap))

    for el in root.iter():
        p = el.get('p')
        if p:
            parts = p.split()
            el.set('p', '%.3f %.3f' % (float(parts[0]) - cx, float(parts[1]) - cy)
                   + ((' ' + ' '.join(parts[2:])) if len(parts) > 2 else ''))
        bb = el.get('BoundingBox')
        if bb:
            vals = [float(v) for v in bb.split()]
            if len(vals) == 4:
                el.set('BoundingBox', '%.3f %.3f %.3f %.3f'
                       % (vals[0] - cx, vals[1] - cy, vals[2] - cx, vals[3] - cy))
        # color 0 is the document default; keeping it makes the colour-index
        # check see an index outside its 2..ntab+1 window
        if el.get('color') == '0':
            del el.attrib['color']
        old = el.get('id')
        if old and old in idmap:
            el.set('id', idmap[old])
            el.set('Z', idmap[old])
        for ref in ('B', 'E'):
            v = el.get(ref)
            if v and v in idmap:
                el.set(ref, idmap[v])

    texts = [''.join(t.text or '' for t in t_el.iter('s')) for t_el in root.iter('t')]
    return {'root': root, 'half_w': half_w, 'half_h': half_h, 'texts': texts,
            'template': str(spec.get('template') or 'chair').lower()}


def emit(item: dict, cx: float, cy: float) -> str:
    """Serialize the structure so its centre lands on (cx, cy)."""
    root = ET.fromstring(ET.tostring(item['root'], encoding='unicode'))
    for el in root.iter():
        p = el.get('p')
        if p:
            parts = p.split()
            el.set('p', '%.3f %.3f' % (float(parts[0]) + cx, float(parts[1]) + cy)
                   + ((' ' + ' '.join(parts[2:])) if len(parts) > 2 else ''))
        bb = el.get('BoundingBox')
        if bb:
            vals = [float(v) for v in bb.split()]
            if len(vals) == 4:
                el.set('BoundingBox', '%.3f %.3f %.3f %.3f'
                       % (vals[0] + cx, vals[1] + cy, vals[2] + cx, vals[3] + cy))
    return ET.tostring(root, encoding='unicode')
