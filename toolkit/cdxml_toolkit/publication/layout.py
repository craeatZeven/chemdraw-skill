"""Measured page layout and native editable compound groups."""
from collections import Counter
from copy import deepcopy
from pathlib import Path
import math
import re
from xml.etree import ElementTree as ET
from PIL import ImageFont
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS
from ..chemistry_semantics import inventory_keys, read_molecules, semantic_key, document_inventory
from ..mcp_runtime.figure_tools import compose_chemical_figure, _move
from ..text_formatting import _ELEMENTS, is_formula_label, subscript_matches
from .data import field_text, write_json


def style_settings(spec):
    # 10 pt is the ACS 1996 label and caption size: a panel that runs out of room grows or
    # paginates, it never shrinks the type.
    style = {'width_mm':180, 'height_mm':240, 'font_size':10, 'bond_length':14.4,
             'font':'Arial', **spec.get('style',{})}
    if not 60 <= float(style['width_mm']) <= 400 or not 60 <= float(style['height_mm']) <= 500:
        raise ValueError('Figure width/height must be 60..400 / 60..500 mm')
    if not 7 <= float(style['font_size']) <= 24 or not 10 <= float(style['bond_length']) <= 40:
        raise ValueError('Use 7..24 pt text and 10..40 pt bonds')
    return style


def font_for(style, size=None, bold=False):
    size = size or style['font_size']
    name = style.get('font','Arial')
    plain = (name+'.ttf', name.lower()+'.ttf', 'LiberationSans-Regular.ttf', 'DejaVuSans.ttf')
    heavy = (name+' Bold.ttf', name.lower()+'bd.ttf', name+'-Bold.ttf', 'LiberationSans-Bold.ttf', 'DejaVuSans-Bold.ttf')
    for candidate in (heavy + plain if bold else plain):
        try:
            return ImageFont.truetype(candidate, round(float(size)*8))
        except OSError:
            pass
    return ImageFont.load_default()


def text_width(value, style, size=None, bold=False):
    font = font_for(style,size,bold)
    return max((font.getlength(line)/8 for line in str(value).splitlines()), default=0)


def wrap(value, width, style):
    lines=[]
    for paragraph in str(value).split('\n'):
        line=''
        for char in paragraph:
            if line and text_width(line+char,style)>width:
                lines.append(line.rstrip()); line=''
            line+=char
        lines.append(line)
    return '\n'.join(lines)


def fragment_mol(fragment):
    molecules=read_molecules('<CDXML><page>'+ET.tostring(fragment,encoding='unicode')+'</page></CDXML>')
    mol=molecules[0]
    for other in molecules[1:]: mol=Chem.CombineMols(mol,other)
    return mol


def node_text(node):
    return ''.join(s.text or '' for s in node.iter('s')) if node is not None else ''


def molecule_coordinates(mol, bond_length, reference=None):
    mol=Chem.Mol(mol)
    warning=None
    if reference is not None:
        mcs=rdFMCS.FindMCS([reference,mol],timeout=3,ringMatchesRingOnly=True,completeRingsOnly=True)
        if mcs.canceled or mcs.numAtoms<3:
            rdDepictor.Compute2DCoords(mol)
            warning='No usable common core; independently depicted'
        else:
            query=Chem.MolFromSmarts(mcs.smartsString)
            pairs=list(zip(reference.GetSubstructMatch(query),mol.GetSubstructMatch(query)))
            rdDepictor.GenerateDepictionMatching2DStructure(mol,reference,atomMap=pairs)
    else:
        rdDepictor.Compute2DCoords(mol)
    scale=float(bond_length)/1.5
    def points():
        conf=mol.GetConformer()
        return [(conf.GetAtomPosition(i).x*scale,-conf.GetAtomPosition(i).y*scale) for i in range(mol.GetNumAtoms())]
    xy=points()
    if reference is not None:
        # Aligning to a reference can fold two atoms onto one point when the common core is small:
        # a biaryl matched by one ring puts the second ring on top of the first. Independent
        # depiction is honest about the shape, so fall back to it rather than drawing a knot.
        floor=0.35*float(bond_length)
        bonded={(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()))
                for b in mol.GetBonds()}
        n=len(xy)
        closest=min([math.hypot(xy[i][0]-xy[j][0],xy[i][1]-xy[j][1])
                     for i in range(n) for j in range(i+1,n) if (i,j) not in bonded] or [floor])
        if closest < floor:
            # Independent depiction, then a rigid turn onto the reference's core: a rotation plus a
            # translation cannot fold two atoms together, so the table keeps one orientation without
            # the knot the alignment itself produced.
            rdDepictor.Compute2DCoords(mol)
            xy=points()
            query=Chem.MolFromSmarts(mcs.smartsString)
            here=mol.GetSubstructMatch(query); there=reference.GetSubstructMatch(query)
            if len(here)>=3 and len(here)==len(there):
                rconf=reference.GetConformer()
                rpts=[(rconf.GetAtomPosition(k).x*scale,-rconf.GetAtomPosition(k).y*scale) for k in there]
                ax,ay=xy[here[1]][0]-xy[here[0]][0], xy[here[1]][1]-xy[here[0]][1]
                bx,by=rpts[1][0]-rpts[0][0], rpts[1][1]-rpts[0][1]
                base_ang=math.atan2(by,bx)-math.atan2(ay,ax)
                cx=sum(x for x,_ in xy)/len(xy); cy=sum(y for _,y in xy)/len(xy)
                rcx=sum(x for x,_ in rpts)/len(rpts); rcy=sum(y for _,y in rpts)/len(rpts)
                # The common core can be symmetric, so the first matched pair may sit a symmetry
                # step away from the reference's; keep the turn that lands the core atoms closest.
                best=None
                for extra in [0,30,60,90,120,150,180,210,240,270,300,330]:
                    ang=base_ang+math.radians(extra)
                    c,s=math.cos(ang),math.sin(ang)
                    moved=[((x-cx)*c-(y-cy)*s+rcx, (x-cx)*s+(y-cy)*c+rcy) for x,y in xy]
                    err=sum(math.hypot(moved[here[k]][0]-rpts[k][0], moved[here[k]][1]-rpts[k][1])
                            for k in range(len(there)))
                    if best is None or err<best[0]:
                        best=(err,moved)
                xy=best[1]
            warning=f'Alignment folded two atoms {closest:.1f} pt apart; depicted independently and turned onto the core'
    return mol,xy,warning


# Cell geometry, in multiples of the text size (10 pt by default). Journal scope tables hold
# 18 to 22 compounds on a page at body-text size, so a cell spends its height on the structure
# and its lines rather than on padding. Each value is the least that keeps check_pages' boxes
# apart: it pads a structure by 1.2 sizes and draws a caption from 1 size above its baseline to
# 0.3 below.
CELL = {
    'pad': 1.2,               # around the structure's atoms
    'caption_baseline': 2.3,  # from the lowest atom to the first caption baseline
    'line_advance': 1.25,     # from one caption line to the next, plus line_gap points
    'line_gap': 2.0,
    'wrapped_line': 1.5,      # extra advance for each wrapped line inside one caption
    'foot': 0.3,              # below the last baseline, plus foot_gap points
    'foot_gap': 4.0,
}


# A secondary product shares its main product's cell, "+" between them, its structure grey and
# its number and yield black (Nature Communications PMC13594109, fig. 4: 14a + 14b).
SECONDARY_GAP = 18.0
SECONDARY_COLOR = '#808080'
# A group's name sits on its rule, a short piece before it and the rest after it (fig. 4).
BAND_LEAD = 12.0
BAND_GAP = 4.0
LINE_JOIN = ', '   # "3a, 64%": the number and the values that share its line


_FOOTNOTE = re.compile(r'\^([A-Za-z0-9*\u2020\u2021\u00a7]+)')
_ITALIC_SPAN = re.compile(r'\*([^*\s][^*\n]*?[^*\s]|[^*\s])\*')
_FOOTNOTE_LEAD = re.compile(r'^(?:\^([A-Za-z0-9*\u2020\u2021\u00a7]+)|([a-z])(?=\s+[A-Z0-9(\[]))(\s+.*)$')
_LOCANT = re.compile(r'(?:(?<=^)|(?<=[\s(]))(?:(ortho|meta|para|ipso|cis|trans|syn|anti|endo|exo|tert|sec)(?=[-\s])'
                     r'|([NOSCPRZE])(?=-))')
_STEREO_PREFIX = re.compile(r'^(cis|trans|rac|syn|anti|endo|exo|ent|epi|meso|E|Z)(-\s*)(.+)$')


def _italic_and_chemical_runs(piece):
    """Expand *italic* spans and chemical formula subscripts inside a plain caption piece."""
    out = []
    plain_parts = []
    offset = 0
    for m in _ITALIC_SPAN.finditer(piece):
        if m.start() > offset:
            chunk = piece[offset:m.start()]
            plain_parts.append(chunk)
            out += chemical_runs(chunk) or [{'text': chunk}]
        inner = m.group(1)
        plain_parts.append(inner)
        for r in (chemical_runs(inner) or [{'text': inner}]):
            out.append({**r, 'italic': True})
        offset = m.end()
    if offset < len(piece):
        tail = piece[offset:]
        plain_parts.append(tail)
        out += chemical_runs(tail) or [{'text': tail}]
    return ''.join(plain_parts), out


def caption_runs(text, bold=False):
    """(the text as it reads, runs or None) for a caption line.

    ``^b`` marks a footnote letter: it is drawn as an italic superscript and reads back as the
    letter alone ("40%^b" reads "40%b"). ``*text*`` marks an italic span (e.g. "*from thymol*").
    The compound label (``bold``) is drawn bold, with leading stereochemical prefixes
    (``cis-3a``, ``rac-2a``, ``E-1a``) in italic and subscripts only when it is a formula
    (``is_formula_label``: PhB(OH)2, not TS2); other lines get formula subscripts (``chemical_runs``).
    """
    pieces, offset = [], 0
    for m in _FOOTNOTE.finditer(text):
        pieces += [(text[offset:m.start()], False), (m.group(1), True)]
        offset = m.end()
    pieces.append((text[offset:], False))
    runs = []
    plain_parts = []
    for piece, mark in pieces:
        if not piece:
            continue
        if mark:
            plain_parts.append(piece)
            runs.append({'text': piece, 'superscript': True, 'italic': True})
        elif bold and is_formula_label(piece):
            plain_parts.append(piece)
            runs += [{**r, 'bold': True} for r in chemical_runs(piece)]
        elif bold:
            plain_parts.append(piece)
            sm = _STEREO_PREFIX.match(piece)
            if sm:
                runs.append({'text': sm.group(1), 'italic': True})
                runs.append({'text': sm.group(2)})
                runs.append({'text': sm.group(3), 'bold': True})
            else:
                runs.append({'text': piece, 'bold': True})
        else:
            p_text, p_runs = _italic_and_chemical_runs(piece)
            plain_parts.append(p_text)
            runs += p_runs
    plain = ''.join(plain_parts)
    if not bold and all(set(r) == {'text'} for r in runs):
        return plain, None
    return plain, runs


def footnote_runs(footer):
    """Format multi-line footnotes with italic superscript leading markers and chemical subscripts."""
    lines = footer.splitlines()
    all_runs = []
    plain_lines = []
    for idx, line in enumerate(lines):
        if idx:
            all_runs.append({'text': '\n'})
        m = _FOOTNOTE_LEAD.match(line)
        if m:
            marker = m.group(1) or m.group(2)
            rest = m.group(3)
            rest_plain, rest_runs = caption_runs(rest, bold=False)
            plain_lines.append(marker + rest_plain)
            all_runs.append({'text': marker, 'superscript': True, 'italic': True})
            all_runs += rest_runs or [{'text': rest_plain}]
        else:
            p, r = caption_runs(line, bold=False)
            plain_lines.append(p)
            all_runs += r or [{'text': p}]
    plain = '\n'.join(plain_lines)
    if all(set(r) == {'text'} for r in all_runs):
        return plain, None
    return plain, all_runs


def heading_runs(text):
    """A group's name with only its locants in italic: *para* substituents, *N*-heterocycles."""
    runs, offset = [], 0
    for m in _LOCANT.finditer(text):
        if m.start() > offset:
            runs.append({'text': text[offset:m.start()]})
        runs.append({'text': m.group(0), 'italic': True})
        offset = m.end()
    if offset < len(text):
        runs.append({'text': text[offset:]})
    return runs


SCRIPT_SCALE = 0.7   # a sub- or superscript is drawn this much smaller than its line


def _runs_width(runs, style):
    size = style['font_size']
    return sum(text_width(r['text'], style, size * SCRIPT_SCALE if r.get('subscript') or r.get('superscript') else None,
                          bold=bool(r.get('bold'))) for r in runs)


def _line_width(keys, texts, style, runs=None):
    """Width of a shared caption line: each piece as drawn (its subscripts smaller), joined by ", "."""
    runs = runs or {}
    return (sum(_runs_width(runs[k], style) if runs.get(k) else text_width(texts[k], style, bold=k == 'label')
                for k in keys)
            + text_width(LINE_JOIN, style) * (len(keys) - 1))


def chemical_runs(text):
    """Runs for a caption or conditions line whose formula digits are subscripts.

    The same rule as the scheme conditions (``text_formatting.subscript_matches``): digits
    after a letter or ")" are counts ("C7H15", "Pd2(dba)3"); yields, ratios and compound
    labels such as "L1" are left alone. None when nothing needs a subscript.
    """
    matches = subscript_matches(text)
    if not matches:
        return None
    runs, offset = [], 0
    for m in matches:
        start, end = m.start(2), m.end(2)
        if start > offset:
            runs.append({'text': text[offset:start]})
        runs.append({'text': text[start:end], 'subscript': True})
        offset = end
    if offset < len(text):
        runs.append({'text': text[offset:]})
    return runs


_DEFINITION = re.compile(r"^\s*([A-Za-z][A-Za-z']*)\s*=\s*((?:[A-Z][a-z]?\d*)+)\s*$")


def parse_definition(text):
    """(label, {element: count}) for a definition such as "Alk = C7H15", else None."""
    m = _DEFINITION.match(text or '')
    if not m:
        return None
    counts = {}
    for symbol, number in re.findall(r'([A-Z][a-z]?)(\d*)', m.group(2)):
        if symbol not in _ELEMENTS:
            return None
        counts[symbol] = counts.get(symbol, 0) + int(number or 1)
    return m.group(1), counts


def label_group(mol, definition):
    """The substituent a definition names, to be drawn as its label: {label, atoms, outside}.

    The group is the far side of a single, non-ring bond whose atoms (with their hydrogens)
    make exactly the defined formula: "Alk = C7H15" finds the heptyl group. ``atoms`` starts
    with the atom bonded to the rest (``outside``). None when nothing matches, or the match is
    not unique.
    """
    parsed = parse_definition(definition)
    if parsed is None or mol is None:
        return None
    label, want = parsed
    found = []
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        for root, outside in ((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                              (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())):
            side, queue = [root], [root]
            while queue:
                for nb in mol.GetAtomWithIdx(queue.pop()).GetNeighbors():
                    k = nb.GetIdx()
                    if k != outside and k not in side:
                        side.append(k); queue.append(k)
            if len(side) >= mol.GetNumAtoms():
                continue
            counts = {}
            for k in side:
                atom = mol.GetAtomWithIdx(k)
                counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
                if atom.GetTotalNumHs():
                    counts['H'] = counts.get('H', 0) + atom.GetTotalNumHs()
            if counts == want:
                found.append({'label': label, 'atoms': tuple(side), 'outside': outside})
    return found[0] if len(found) == 1 else None


def header_identities(spec):
    """What the header scheme draws, as the page check counts structures (R numbers removed)."""
    header = spec.get('header') or {}
    found = []
    for item in [*header.get('substrates', []), *header.get('products', [])]:
        if not item.get('smiles'):
            continue  # a reagent written as its abbreviation carries no structure
        mol = Chem.MolFromSmiles(item['smiles'])
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        found += [semantic_key(m) for m in Chem.GetMolFrags(mol, asMols=True)]
    return found


def _shorten_ring_r_bonds(xy, mol, bond_length):
    """Keep a short real attachment for each R group bonded to a ring vertex."""
    rings = mol.GetRingInfo().AtomRings()
    points = list(xy)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(neighbours) != 1:
            continue
        anchor = neighbours[0]
        ring = next((set(r) for r in rings if anchor in r), None)
        if not ring:
            continue
        cx = sum(points[k][0] for k in ring) / len(ring)
        cy = sum(points[k][1] for k in ring) / len(ring)
        ax, ay = points[anchor]
        # Out of the ring, the way a general scheme has it: the mark sits outside the vertex, and the
        # short bond between them stays visible.
        dx, dy = ax - cx, ay - cy
        length = math.hypot(dx, dy) or 1.0
        stub = 0.75 * float(bond_length)
        points[atom.GetIdx()] = (ax + dx / length * stub, ay + dy / length * stub)
    return points


def _upright_ring(xy, marks):
    """Turn the drawing the little it takes to put a hexagon edge upright.

    A six-ring's edges sit 60 degrees apart, so making one of them vertical is never more than a 30
    degree turn: the angle is reduced modulo 60 and taken on the short side. Anything bigger (which is
    what turning by an edge's own angle produced) tips the rest of the molecule, and that is why a
    product ended up upside down.
    """
    if not marks:
        return xy
    _nickname, _start, _end, _position, _align, _label, edge = next(iter(marks.values()))
    ax, ay, bx, by = edge[:4]
    angle = math.degrees(math.atan2(by - ay, bx - ax))
    turn = (90.0 - angle) % 60.0
    if turn > 30.0:
        turn -= 60.0
    if abs(turn) < 1e-6:
        return xy
    xs = [p[0] for p in xy]
    ys = [p[1] for p in xy]
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    radians = math.radians(turn)
    c, s = math.cos(radians), math.sin(radians)
    return [((x - cx) * c - (y - cy) * s + cx, (x - cx) * s + (y - cy) * c + cy) for x, y in xy]


def _ring_edge_marks(xy, mol, bond_length, style):
    """Lay out numbered, terminal R marks across an outer aromatic ring edge."""
    from ..render.renderer import _SUPERSCRIPT

    ring_info = mol.GetRingInfo()
    size = style['font_size']
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule)
    marks = {}
    used_edges = set()
    for atom in mol.GetAtoms():
        number = atom.GetAtomMapNum()
        if atom.GetAtomicNum() != 0 or not number or atom.GetDegree() != 1:
            continue
        attachment, = atom.GetBonds()
        anchor = attachment.GetOtherAtom(atom)
        index = anchor.GetIdx()
        rings = [ring for ring in ring_info.AtomRings() if index in ring]
        if attachment.GetBondType() != Chem.BondType.SINGLE or len(rings) != 1 or not anchor.GetIsAromatic():
            continue
        edges = [bond for bond in anchor.GetBonds() if bond.GetIsAromatic()
                 and ring_info.NumBondRings(bond.GetIdx()) == 1 and bond.GetIdx() not in used_edges]
        if not edges:
            continue
        ax, ay = xy[index]
        mx, my = xy[atom.GetIdx()]
        axis = 0 if abs(mx - ax) >= abs(my - ay) else 1
        # Which way the R atom lay from its ring atom, used only to prefer the edge on that side.
        sign = math.copysign(1.0, (mx - ax, my - ay)[axis] or 1.0)
        # Prefer a double bond if the stroke can cross both its lines at a clear angle.
        def edge_score(bond):
            other = xy[bond.GetOtherAtomIdx(index)]
            transverse = abs(other[1 - axis] - xy[index][1 - axis])
            double = kekule.GetBondWithIdx(bond.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
            return (double and transverse >= 0.6 * math.dist(other, xy[index]), sign * other[axis])

        edge = max(edges, key=edge_score)
        used_edges.add(edge.GetIdx())
        bx, by = xy[edge.GetOtherAtomIdx(index)]
        x, y = (ax + bx) / 2, (ay + by) / 2
        # The stroke is the normal of the edge it crosses, so it stays square to that edge whatever
        # angle the ring was depicted at: an axis-aligned stroke would cross a slanted edge askew.
        ex, ey = bx - ax, by - ay
        length = math.hypot(ex, ey) or 1.0
        ux, uy = -ey / length, ex / length
        # Point it out of the ring, towards where the R atom was drawn.
        if ux * (mx - x) + uy * (my - y) < 0:
            ux, uy = -ux, -uy
        outside, inside = bond_length * 0.22, bond_length * 0.32
        start = (x + ux * outside, y + uy * outside)
        end = (x - ux * inside, y - uy * inside)
        position = (start[0] + ux * 1.5, start[1] + uy * size * 0.65 + size * 0.35)
        align = 'Left' if ux > 0.2 else 'Right' if ux < -0.2 else 'Center'
        label = 'R' + _SUPERSCRIPT.get(number, str(number))
        marks[atom.GetIdx()] = (f'R{number}', start, end, position, align, label, (ax, ay, bx, by))
    return marks


def _shift_obj_x(obj, dx):
    """Translate a manifest object horizontally by dx points."""
    if not dx:
        return
    t = obj.get('type')
    if t == 'molecule':
        obj['coordinates'] = [[px + dx, py] for px, py in obj['coordinates']]
    elif t == 'text':
        obj['position'] = [obj['position'][0] + dx, obj['position'][1]]
    elif t in ('arrow', 'line', 'rectangle'):
        obj['start'] = [obj['start'][0] + dx, obj['start'][1]]
        obj['end'] = [obj['end'][0] + dx, obj['end'][1]]


def header_objects(spec, style, left, top):
    """The general reaction drawn over a scope grid (BJOC 22-127 scheme 2): substrates joined
    by "+", the arrow with its conditions, and the products with their labels. Also returns
    the bounds and the real R attachments to hide behind decorative marks."""
    header = spec['header']; size = style['font_size']; gap = size * 1.2
    def shape(item, side='substrates'):
        if item.get('text') and not item.get('smiles'):
            return {'item': item, 'text': True, 'w': text_width(item['text'], style),
                    'h': size * 1.6, 'x0': 0.0, 'y0': -size * 0.8, 'dropped': []}
        # The general reaction is drawn larger than the grid below it, as journal figures do.
        header_bond = float(style['bond_length']) * 1.7
        mol, xy, _ = molecule_coordinates(Chem.MolFromSmiles(item['smiles']), header_bond)
        if spec.get('rmark_in_ring'):
            xy = _shorten_ring_r_bonds(xy, mol, header_bond)
        # A header item may ask for its own turn (a reagent drawn the other way round in the source
        # figure); otherwise the products follow the grid's orientation.
        turn = item.get('rotation')
        if turn is None and side == 'products':
            turn = spec.get('depiction_rotation')
        if float(turn or 0):
            # The general reaction above the grid must face the same way as the cells below it.
            angle=math.radians(float(turn))
            c,sn=math.cos(angle),math.sin(angle)
            hx=sum(x for x,_ in xy)/len(xy); hy=sum(y for _,y in xy)/len(xy)
            xy=[((x-hx)*c-(y-hy)*sn+hx, (x-hx)*sn+(y-hy)*c+hy) for x,y in xy]
        marks = _ring_edge_marks(xy, mol, header_bond, style) if spec.get('rmark_in_ring') else {}
        if marks:
            xy = _upright_ring(xy, marks)
            marks = _ring_edge_marks(xy, mol, header_bond, style)
        xs, ys = [p[0] for p in xy], [p[1] for p in xy]
        for _, start, end, (tx, ty), align, label, _edge in marks.values():
            label_width = text_width(label, style)
            lx = tx - {'Right': label_width, 'Center': label_width / 2, 'Left': 0}[align]
            xs.extend([start[0], end[0], lx, lx + label_width]); ys.extend([start[1], end[1], ty - size, ty])
        return {'item': item, 'xy': xy, 'marks': marks, 'w': max(xs) - min(xs), 'h': max(ys) - min(ys),
                'x0': min(xs), 'y0': min(ys)}
    subs = [shape(i, 'substrates') for i in header.get('substrates', [])]
    prods = [shape(i, 'products') for i in header.get('products', [])]
    above = [caption_runs(t) for t in header.get('above', [])]
    below = [caption_runs(t) for t in header.get('below', [])]
    line = size * 1.2
    above_h = len(above) * line + 4
    tallest = max(s['h'] for s in subs + prods)
    mid = top + max(tallest / 2 + gap, above_h + 4)
    arrow_len = max([text_width(p, style) for p, _ in above + below] + [6 * style['bond_length']]) + 16
    objects = []; x = left; bottom = mid
    hidden_marks = {}
    def molecule(s, key):
        nonlocal x, bottom
        dx = x + gap - s['x0']; dy = mid - s['y0'] - s['h'] / 2
        if s.get('text'):
            plain = s['item']['text']
            runs = chemical_runs(plain)
            objects.append({'type': 'text', 'id': key, 'text': plain, 'align': 'Center',
                            'size': size, 'position': [x + gap + s['w'] / 2, mid + size * 0.35],
                            **({'runs': runs} if runs else {})})
            bottom = max(bottom, mid + size)
            if s['item'].get('label'):
                plain, runs = caption_runs(s['item']['label'], bold=True)
                objects.append({'type': 'text', 'text': plain, 'align': 'Center', 'size': size,
                                'position': [x + gap + s['w'] / 2, mid + s['h'] / 2 + gap + size],
                                **({'runs': runs} if runs else {})})
            x += s['w'] + 2 * gap
            return
        mol_obj = {'type': 'molecule', 'id': key, 'smiles': s['item']['smiles'],
                   'coordinates': [[px + dx, py + dy] for px, py in s['xy']]}
        if s['item'].get('color'):
            mol_obj['color'] = s['item']['color']
        objects.append(mol_obj)
        for index, (nickname, start, end, position, align, label, _edge) in s['marks'].items():
            hidden_marks.setdefault(key, {})[index] = nickname
            objects.append({'type': 'line', 'start': [start[0] + dx, start[1] + dy],
                            'end': [end[0] + dx, end[1] + dy]})
            # The label leaves the ring line interrupted where it sits, the way a scheme drawn by hand
            # does: a white patch under the text, not a stroke running through the R.
            half_w = len(label) * size * 0.34 + size * 0.10
            half_h = size * 0.62
            box_left = position[0] + dx - (half_w if align == 'Center' else 0.0 if align == 'Left' else 2 * half_w)
            objects.append({'type': 'rectangle', 'start': [box_left, position[1] + dy - half_h],
                            'end': [box_left + 2 * half_w, position[1] + dy + half_h * 1.15],
                            'color': '#FFFFFF', 'fill_type': 'solid', 'background': '#FFFFFF'})
            objects.append({'type': 'text', 'text': label, 'align': align, 'size': size,
                            'position': [position[0] + dx, position[1] + dy],
                            'color': spec.get('header_mark_color') or s['item'].get('color')})
        bottom = max(bottom, mid + s['h'] / 2 + gap)
        if s['item'].get('label'):
            plain, runs = caption_runs(s['item']['label'], bold=True)
            objects.append({'type': 'text', 'text': plain, 'align': 'Center', 'size': size,
                            'position': [x + gap + s['w'] / 2, mid + s['h'] / 2 + gap + size],
                            **({'runs': runs} if runs else {})})
            bottom = max(bottom, mid + s['h'] / 2 + gap + size * 1.3)
        x += s['w'] + 2 * gap
    for i, s in enumerate(subs):
        if i:
            objects.append({'type': 'text', 'text': '+', 'align': 'Center', 'size': size,
                            'position': [x, mid + size * 0.35]})
            x += gap
        molecule(s, f'header_{i}')
    start = x; end = x + arrow_len
    objects.append({'type': 'arrow', 'start': [start, mid], 'end': [end, mid]})
    for n, (plain, runs) in enumerate(reversed(above)):
        objects.append({'type': 'text', 'text': plain, 'align': 'Center', 'size': size,
                        'position': [(start + end) / 2, mid - 4 - n * line - size * 0.3], **({'runs': runs} if runs else {})})
    for n, (plain, runs) in enumerate(below):
        objects.append({'type': 'text', 'text': plain, 'align': 'Center', 'size': size,
                        'position': [(start + end) / 2, mid + 4 + size * 0.84 + n * line], **({'runs': runs} if runs else {})})
        bottom = max(bottom, mid + 4 + size + n * line + size * 0.3)
    x = end + gap * 0.5
    for i, s in enumerate(prods):
        molecule(s, f'header_p{i}')
    return objects, bottom, x, hidden_marks


def _hide_header_attachments(page, object_ids, hidden_marks):
    """Hide the real R atom and bond while retaining both in the editable fragment."""
    for key, marks in hidden_marks.items():
        fragment = page.find(f"fragment[@id='{object_ids[key]}']")
        # The renderer writes these explicit-coordinate fragments in input atom order.
        atoms = fragment.findall('n')
        for index, nickname in marks.items():
            atom = atoms[index]
            if atom.get('NodeType') != 'GenericNickname' or atom.get('GenericNickname') != nickname:
                raise ValueError('Could not bind a ring-edge mark to its R atom')
            atom_id = atom.get('id')
            bond, = (b for b in fragment.findall('b') if atom_id in (b.get('B'), b.get('E')))
            label = atom.find('t')
            # A native save interprets the hidden label as chemistry, so give it plain R1.
            label.find('s').text = nickname
            for element in (atom, bond, label):
                element.set('Visible', 'no')


def _caption_advance(value, size):
    return value.count('\n')*size*CELL['wrapped_line']+size*CELL['line_advance']+CELL['line_gap']


def prepare(records,spec):
    style=style_settings(spec); size=style['font_size']
    usable=style['width_mm']*72/25.4-32
    refs={}; prepared=[]; warnings=[]
    by_id={r['compound_id']:r for r in records}
    for record in records:
        group=record['group']; ref_id=spec.get('core_by_group',{}).get(group)
        if group not in refs and ref_id:
            if ref_id not in by_id: raise ValueError(f'Unknown reference compound: {ref_id}')
            ref=Chem.MolFromSmiles(by_id[ref_id]['smiles']); rdDepictor.Compute2DCoords(ref); refs[group]=ref
        mol,xy,warning=molecule_coordinates(Chem.MolFromSmiles(record['smiles']),style['bond_length'],refs.get(group))
        if group not in refs and float(spec.get('depiction_rotation') or 0):
            # The first compound of a table fixes the orientation every later cell is aligned to, so
            # turning it turns the whole table: 180 puts a benzothiazole on the left where the source
            # figure has it, instead of the right where RDKit's own depiction does. Both the points
            # and the conformer are turned, each in its own units.
            angle=math.radians(float(spec.get('depiction_rotation')))
            c,sn=math.cos(angle),math.sin(angle)
            cx=sum(x for x,_ in xy)/len(xy); cy=sum(y for _,y in xy)/len(xy)
            xy=[((x-cx)*c-(y-cy)*sn+cx, (x-cx)*sn+(y-cy)*c+cy) for x,y in xy]
            conf=mol.GetConformer()
            inner=[(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())]
            icx=sum(x for x,_ in inner)/len(inner); icy=sum(y for _,y in inner)/len(inner)
            for i,(x,y) in enumerate(inner):
                conf.SetAtomPosition(i, ((x-icx)*c-(y-icy)*sn+icx, (x-icx)*sn+(y-icy)*c+icy, 0.0))
        refs.setdefault(group,mol)
        if warning: warnings.append({'compound_id':record['compound_id'],'message':warning})
        # A substituent the compound's definition names ("Alk = C7H15") is drawn as its label:
        # for the layout its atoms stand just past the attachment atom, and the label's width
        # counts instead of the chain's.
        named=None; extent=list(xy)
        source=spec.get('label_groups_from')
        if source and record['values'].get(source):
            try:
                named=label_group(mol,record['values'][source])
            except Exception:
                named=None
        if named:
            rx,ry=xy[named['atoms'][0]]; ox,oy=xy[named['outside']]
            d=math.hypot(rx-ox,ry-oy) or 1.0; ux,uy=(rx-ox)/d,(ry-oy)/d
            xy=list(xy)
            for n,k in enumerate(named['atoms'][1:],1):
                xy[k]=(rx+ux*0.6*n,ry+uy*0.6*n)
            reach=text_width(named['label'],style,ACS_LABEL_SIZE_FALLBACK)+2
            extent=[p for k,p in enumerate(xy)]+[(rx+(reach if ux>=0 else -reach),ry)]
        xs,ys=zip(*extent); pad=size*CELL['pad']; molecule_width=max(xs)-min(xs)+2*pad
        if molecule_width>usable:
            raise ValueError(f"Structure {record['compound_id']} exceeds figure width at requested bond length")
        raw={'label':record['label'],**{f['column']:field_text(record,f) for f in spec.get('fields',[])}}
        raw={k:v for k,v in raw.items() if k=='label' or v!=''}  # an empty optional line is not drawn
        plain={k:caption_runs(v,bold=k=='label')[0] for k,v in raw.items()}
        first_runs={k:caption_runs(v,bold=k=='label')[1] for k,v in raw.items()}
        # The number and the values marked inline share one line ("3a, 64%"); the rest, one each.
        inline=['label']+[f['column'] for f in spec.get('fields',[]) if f.get('inline') and f['column'] in raw]
        # The cell is as wide as the structure, or as the widest single value. A shared line may
        # widen it by half again; longer, it breaks after the number instead.
        pieces=[text_width(v,style,bold=k=='label') for k,v in plain.items()]
        width=min(max(100,molecule_width,min(usable,max(pieces)+16)),usable)
        shared=_line_width(inline,plain,style,first_runs)+16
        if len(inline)>1 and shared<=min(usable,1.5*width):
            width=max(width,shared)
        elif len(inline)>1:
            inline=['label']  # too long for one line: each value on its own line
        joined=inline if len(inline)>1 else []
        raw={k:(v if k in joined else wrap(v,width-12,style)) for k,v in raw.items()}
        texts={k:caption_runs(v,bold=k=='label')[0] for k,v in raw.items()}
        runs={k:caption_runs(v,bold=k=='label')[1] for k,v in raw.items()}
        lines=((_caption_advance('',size) if joined else 0)
               +sum(_caption_advance(v,size) for k,v in texts.items() if k not in joined)
               -size*CELL['line_advance']-CELL['line_gap'])
        height=(pad+max(ys)-min(ys)+size*CELL['caption_baseline']+lines
                +size*CELL['foot']+CELL['foot_gap'])
        # The structure stands in the middle of its cell, over its centred caption.
        shift=pad+(width-molecule_width)/2
        coords=[[x-min(xs)+shift,y-min(ys)+pad] for x,y in xy]
        caption=((_caption_advance('',size) if joined else 0)
                 +sum(_caption_advance(v,size) for k,v in texts.items() if k not in joined))
        prepared.append({'input':record,'coords':coords,'texts':texts,'runs':runs,'inline':joined,'named':named,
                         'caption':caption,'caption_width':max([_line_width(joined,texts,style,runs)]+
                                                              [text_width(v,style) for k,v in texts.items() if k not in joined]),
                         'width':width,'height':height,'molecule_bottom':max(y for x,y in coords)})
    return prepared,warnings


def create_pages(records,spec,stage,start_index=0):
    style=style_settings(spec); prepared,warnings=prepare(records,spec)
    width=style['width_mm']*72/25.4; height=style['height_mm']*72/25.4
    margin=16; gap=12; size=style['font_size']; pages=[]; compounds={}
    columns_limit=spec.get('columns')
    row_align=spec.get('row_align','left')
    band_style=spec.get('band_style','rule')
    # Keep first-occurrence group order, and preserve row order within each group.
    groups=list(dict.fromkeys(p['input']['group'] for p in prepared))
    # A secondary product (a minor isomer, a by-product) shares its main product's cell, drawn
    # beside it in grey (Nature Communications PMC13594109, fig. 4: 14b beside 14a).
    minors={}
    for p in prepared:
        if p['input'].get('secondary_of'): minors.setdefault(p['input']['secondary_of'],[]).append(p)
    # Compounds drawn with another's structure (BJOC 22-127: 3q under 3p's NH-Alk) add their
    # caption lines under that structure's and draw nothing of their own.
    family={}
    for p in prepared:
        if p['input'].get('drawn_with'): family.setdefault(p['input']['drawn_with'],[]).append(p)
    for p in prepared:
        extra=family.get(p['input']['compound_id'],[])
        if extra:
            p['height']+=sum(q['caption'] for q in extra)
            p['width']=min(width-2*margin,max([p['width']]+[q['caption_width']+16 for q in extra]))
            p['family']=extra
    arranged=[p for group in groups for p in prepared
              if p['input']['group']==group and not p['input'].get('secondary_of') and not p['input'].get('drawn_with')]
    manifests=[]; current=None; y=0
    footer='\n'.join(spec.get('footnotes',[]))
    footer=wrap(footer,width-2*margin,style) if footer else ''
    bottom=height-margin-((footer.count('\n')+1)*size*1.5+8 if footer else 0)
    def new_page():
        index=start_index+len(manifests)
        raw_title=spec['title'] if 'title' in spec else ('Substrate scope' if spec.get('mode','scope')=='scope' else 'Structure–activity relationships')
        title=wrap(raw_title,width-2*margin,style) if raw_title else ''
        conditions=wrap(spec.get('conditions',''),width-2*margin,style)
        objects=[]
        top=margin
        right=margin
        if title:
            objects.append({'type':'text','id':'title','text':title,'position':[margin,top+size],'size':size})
            top+=(title.count('\n')+1)*size*1.5+8
            right=max(right,margin+max((text_width(line,style) for line in title.splitlines()),default=0))
        if conditions:
            runs=caption_runs(conditions)[1]
            objects.append({'type':'text','id':'conditions','text':conditions,'position':[margin,top+size],'size':size,
                            **({'runs':runs} if runs else {})})
            top+=(conditions.count('\n')+1)*size*1.5+8
            right=max(right,margin+max((text_width(line,style) for line in conditions.splitlines()),default=0))
        rules=[]
        header_info=None
        hidden_marks={}
        if index==0 and spec.get('header'):
            # The general reaction stands over the first page's grid, a rule under it.
            drawn,header_bottom,header_right,hidden_marks=header_objects(spec,style,margin,top)
            objects+=drawn; right=max(right,header_right)
            header_info=(drawn,header_right)
            rule={'type':'line','start':[margin,header_bottom+4],'end':[width-margin,header_bottom+4],'style':'line'}
            objects.append(rule); rules.append(rule)
            top=header_bottom+4+size+4
        current={'index':index,'objects':objects,'entries':[],'top':top,'bottom':top,'rules':rules,
                 'right':right,'header':header_info,'hidden_marks':hidden_marks,'rows':[],'group':None}
        manifests.append(current)
        return current,margin,top
    group_colors = spec.get('group_colors') or {}
    def band(current,y,name):
        """Open a group: its name on a rule (or soft pill) across the panel, only the locants in italic."""
        base_name = name[:-len(' (continued)')] if name.endswith(' (continued)') else name
        bcolor = group_colors.get(base_name) or group_colors.get(name)
        color_opt = {'color': bcolor} if bcolor else {}
        if band_style == 'pill':
            label_x = margin + 6.0
            heading = wrap(name, width - label_x - margin - 12.0, style)
            baseline = y + size + 2.0
            rule_y = baseline - size * 0.35
            tw = text_width(heading.split('\n')[0], style)
            pill_right = label_x + tw + 6.0
            after = pill_right + BAND_GAP
            pill_color = bcolor or '#4A6FA5'
            runs = [{**r, **color_opt} for r in heading_runs(heading)]
            pill = {'type': 'rectangle', 'start': [margin, baseline - size * 1.05],
                    'end': [pill_right, baseline + size * 0.35],
                    'color': pill_color, 'fill_type': 'faded', 'round_edge': True, 'background': True}
            tail = {'type': 'line', 'start': [after, rule_y], 'end': [width - margin, rule_y],
                    'style': 'line', **color_opt}
            current['objects'] += [pill, tail,
                                   {'type': 'text', 'text': heading, 'runs': runs,
                                    'position': [label_x, baseline], 'size': size}]
            current['rules'].append(tail)
            current['right'] = max(current['right'], after)
            return baseline + heading.count('\n') * size * 1.5 + size * 0.45 + 8
        label_x=margin+BAND_LEAD+BAND_GAP
        heading=wrap(name,width-label_x-margin-BAND_GAP,style)
        baseline=y+size; rule_y=baseline-size*0.35   # through the middle of the lower-case letters
        after=label_x+text_width(heading.split('\n')[0],style)+BAND_GAP
        runs = [{**r, **color_opt} for r in heading_runs(heading)]
        tail={'type':'line','start':[after,rule_y],'end':[width-margin,rule_y],'style':'line',**color_opt}
        current['objects'] += [{'type':'line','start':[margin,rule_y],'end':[margin+BAND_LEAD,rule_y],'style':'line',**color_opt},
                               tail,
                               {'type':'text','text':heading,'runs':runs,
                                'position':[label_x,baseline],'size':size}]
        current['rules'].append(tail)
        current['right']=max(current['right'],after)
        return baseline+heading.count('\n')*size*1.5+size*0.3+8
    keys=[]   # one per placed compound; a family member is placed inside its main compound's call
    def place(current,p,x,y,color=None,ty=None,cx=None,row_mol_bottom=None,row_height=None,collector=None):
        key=f'c{len(keys)}'; keys.append(key); cid=p['input']['compound_id']
        text_only=ty is not None   # a compound drawn with another's structure: its caption alone
        dy_mol = 0.0
        if row_mol_bottom is not None and not text_only:
            dy_mol = row_mol_bottom - p['molecule_bottom']
            if not columns:
                dy_mol /= 2
        objects=[] if text_only else [{'type':'molecule','id':key+'_mol','smiles':p['input']['smiles'],
                  'coordinates':[[a+x,b+y+dy_mol] for a,b in p['coords']],**({'color':color} if color else {})}]
        top=ty
        if ty is None:
            mol_bot=row_mol_bottom if row_mol_bottom is not None else p['molecule_bottom']
            ty=y+mol_bot+size*CELL['caption_baseline']
        cx=x+p['width']/2 if cx is None else cx
        def line(field,at,align):
            runs=p['runs'].get(field)
            return {'type':'text','id':key+'_'+field,'text':p['texts'][field],'position':[at,ty],'size':size,
                    'align':align,**({'runs':runs} if runs else {})}
        if p['inline']:
            # One centred line from separate objects, so each value stays bound to its column.
            lx=cx-_line_width(p['inline'],p['texts'],style,p['runs'])/2
            for i,field in enumerate(p['inline']):
                if i:
                    objects.append({'type':'text','id':f'{key}_join{i}','text':LINE_JOIN,'position':[lx,ty],'size':size})
                    lx+=text_width(LINE_JOIN,style)
                objects.append(line(field,lx,'Left'))
                lx+=_line_width([field],p['texts'],style,p['runs'])
            ty+=_caption_advance('',size)
        for field,value in p['texts'].items():
            if field in p['inline']: continue
            objects.append(line(field,cx,'Center'))
            ty+=_caption_advance(value,size)
        current['objects'].extend(objects)
        if collector is not None:
            collector['objects'].extend(objects)
        if text_only:
            half=p['caption_width']/2
            entry={'key':key,'cid':cid,'bbox':[cx-half,top-size,cx+half,ty-size],'text_only':True}
            current['entries'].append(entry)
            compounds[cid]={'input':p['input'],'page':current['index'],'bbox':entry['bbox'],
                            'display':p['texts'],'locked':False}
            if collector is not None:
                collector['entries'].append(entry)
            return ty
        for member in p.get('family',[]):
            ty=place(current,member,x,y,ty=ty,cx=cx,collector=collector)
        cell_h=row_height if row_height is not None else p['height']
        entry={'key':key,'cid':cid,'bbox':[x,y,x+p['width'],y+cell_h],'named':p.get('named')}
        current['entries'].append(entry)
        current['bottom']=max(current['bottom'],y+cell_h); current['right']=max(current['right'],x+p['width'])
        compounds[cid]={'input':p['input'],'page':current['index'],'bbox':entry['bbox'],
                        'display':p['texts'],'locked':False}
        if collector is not None:
            collector['entries'].append(entry)

    cells = []
    for main in arranged:
        cell = [main, *minors.get(main['input']['compound_id'], [])]
        cell_width = sum(item['width'] for item in cell) + SECONDARY_GAP * (len(cell) - 1)
        cells.append((cell, cell_width))

    # Without row_align, columns requests a fixed grid. An explicit alignment uses flow
    # layout, where columns caps the number of cells and a wide cell can occupy more space.
    columns = columns_limit if 'row_align' not in spec else None
    column_width = None
    if columns:
        # A grid of equal columns (BJOC 22-127 scheme 2: four to a row). Cells that cannot fit
        # are an error with the arithmetic, never a quiet fall back to fewer columns.
        usable=width-2*margin
        column_width=(usable-gap*(columns-1))/columns
        widest=max(cells,key=lambda cell: cell[1]) if cells else None
        if widest is not None and widest[1]>column_width+0.01:
            need=widest[1]*columns+gap*(columns-1)
            need_mm=(need+2*margin)/72*25.4
            raise ValueError(
                f"{columns} columns do not fit: the widest cell ({widest[0][0]['input']['compound_id']}) is "
                f"{widest[1]:.1f} pt at {size:g} pt type and {style['bond_length']:g} pt bonds, so "
                f"{columns} cells and {columns-1} gaps of {gap} pt need {need:.1f} pt; the page has "
                f"{usable:.1f} pt ({float(style['width_mm']):g} mm less {margin} pt margins each side, "
                f"{column_width:.1f} pt a column). Set width_mm to at least {math.ceil(need_mm)}: reduced "
                f"to {float(style['width_mm']):g} mm for print, the type would be "
                f"{size*float(style['width_mm'])/need_mm:.1f} pt. Or ask for fewer columns.")

    row_buffer = []
    row_group = None
    remaining_width = width - 2 * margin

    def flush_row():
        nonlocal current, y
        if not row_buffer:
            return
        all_members = [item for cell, _ in row_buffer for item in cell]
        row_mol_bottom = max(item['molecule_bottom'] for item in all_members)
        row_height = max(item['height'] + row_mol_bottom - item['molecule_bottom'] for item in all_members)
        group_top = y
        opened_band = row_group != current['group']
        if opened_band and row_group:
            band_objects = len(current['objects'])
            band_rules = len(current['rules'])
            band_right = current['right']
            group_top = band(current, y, row_group)
        if group_top + row_height > bottom:
            # Remove the band if the first row of its group has to start on the next page.
            if opened_band and row_group:
                del current['objects'][band_objects:]
                del current['rules'][band_rules:]
                current['right'] = band_right
            continued = current['group'] == row_group
            current, _, group_top = new_page()
            if row_group:
                heading = row_group + ' (continued)' if continued else row_group
                group_top = band(current, group_top, heading)
        if group_top + row_height > bottom:
            raise ValueError(f"Compound {all_members[0]['input']['compound_id']} exceeds page height")
        current['group'] = row_group
        y = group_top
        rx = margin
        if columns:
            # Short rows are centered in the fixed grid, even on a page with no full row.
            rx += (columns - len(row_buffer)) * (column_width + gap) / 2
            current['right'] = max(current['right'], width - margin)
        placed_cells = []
        for cell, cell_width in row_buffer:
            slot_width = column_width if columns else cell_width
            collector = {'objects': [], 'entries': []}
            qx = rx + (slot_width - cell_width) / 2
            for index, item in enumerate(cell):
                if index:
                    previous = cell[index - 1]
                    previous_x = qx - SECONDARY_GAP - previous['width']
                    gap_mid = (previous_x + max(x for x, _ in previous['coords']) +
                               qx + min(x for x, _ in item['coords'])) / 2
                    ys = [pos[1] for pos in previous['coords']]
                    dy = row_mol_bottom - previous['molecule_bottom']
                    if not columns:
                        dy /= 2
                    plus = {'type': 'text', 'text': '+', 'align': 'Center', 'size': size,
                            'position': [gap_mid, y + dy + (min(ys) + max(ys)) / 2 + size * 0.35]}
                    current['objects'].append(plus)
                    collector['objects'].append(plus)
                place(current, item, qx, y, None if index == 0 else SECONDARY_COLOR,
                      row_mol_bottom=row_mol_bottom, row_height=row_height, collector=collector)
                qx += item['width'] + SECONDARY_GAP
            placed_cells.append(collector)
            rx += slot_width + gap
        current['rows'].append({'cells': placed_cells, 'width': rx - gap - margin, 'top': y})
        y += row_height + gap
        row_buffer.clear()
    for cell, cell_width in cells:
        group = cell[0]['input']['group']
        if current is None:
            current, _, y = new_page()
        if row_buffer and (group != row_group or len(row_buffer) == columns_limit
                           or (not columns and cell_width > remaining_width)):
            flush_row()
            remaining_width = width - 2 * margin
        row_group = group
        row_buffer.append((cell, cell_width))
        remaining_width -= cell_width + gap
    flush_row()
    # Crop pages to their actual rows, headers and footnotes.
    for page in manifests:
        if columns:
            page['grid']={'columns':columns,'cell_width':round(column_width,2),'gap':gap,
                          'row_pitch':[round(b['top']-a['top'],2) for a,b in zip(page['rows'],page['rows'][1:])]}
        page_bottom=page['bottom']
        if footer:
            foot_plain,foot_runs=footnote_runs(footer)
            page['objects'].append({'type':'text','id':'footnotes','text':foot_plain,
                                    'position':[margin,page_bottom+8+size],'size':size,
                                    **({'runs':foot_runs} if foot_runs else {})})
            page_bottom+=8+(footer.count('\n')+1)*size*1.5
            page['right']=max(page['right'],margin+max(text_width(line,style) for line in foot_plain.splitlines()))
        content_w=page['right']-margin
        # Center the general reaction header over the panel content width.
        if page.get('header'):
            drawn,header_right=page['header']
            dx=(page['right']-header_right)/2.0
            if dx>0.5:
                for obj in drawn:
                    _shift_obj_x(obj,dx)
        # Apply horizontal row alignment (center or justify) across rows on this page.
        if row_align in ('center','justify') and page['rows']:
            max_cols=max(len(r['cells']) for r in page['rows'])
            for row in page['rows']:
                extra=content_w-row['width']
                if extra<=0.5:
                    continue
                ncells=len(row['cells'])
                if row_align=='justify' and ncells>=2 and ncells==max_cols:
                    step=extra/(ncells-1)
                    for idx,cinfo in enumerate(row['cells']):
                        dx=idx*step
                        for obj in cinfo['objects']:
                            _shift_obj_x(obj,dx)
                        for entry in cinfo['entries']:
                            entry['bbox'][0]+=dx; entry['bbox'][2]+=dx
                else:
                    dx=extra/2.0
                    for cinfo in row['cells']:
                        for obj in cinfo['objects']:
                            _shift_obj_x(obj,dx)
                        for entry in cinfo['entries']:
                            entry['bbox'][0]+=dx; entry['bbox'][2]+=dx
        page['width']=math.ceil(page['right']+margin)
        for rule in page['rules']:
            rule['end'][0]=page['right']  # each band's rule spans the panel as laid out
        page['height']=math.ceil(page_bottom+margin)
        if spec.get('border'):
            bcol=spec['border'] if isinstance(spec['border'],str) else '#CCCCCC'
            page['objects'].insert(0,{'type':'rectangle','start':[margin/2,margin/2],
                                      'end':[page['width']-margin/2,page['height']-margin/2],
                                      'color':bcol,'fill_type':'none','round_edge':True,'background':True})
    stage=Path(stage)
    for page in manifests:
        name=f"page-{page['index']+1:03d}.cdxml"
        manifest=stage/f"page-{page['index']+1:03d}.manifest.json"
        write_json(manifest,{'style':{k:style[k] for k in ('font','font_size','bond_length','line_width','bold_width','hash_spacing','margin_width') if k in style},'objects':page['objects']})
        result=compose_chemical_figure(str(manifest),str(stage/name))
        tree=ET.parse(stage/name); root=tree.getroot(); native_page=root.find('page')
        root.set('BoundingBox',f"0 0 {page['width']:g} {page['height']:g}")
        native_page.set('BoundingBox',f"0 0 {page['width']:g} {page['height']:g}")
        native_page.set('WidthPages','1');native_page.set('HeightPages','1')
        ids=result['metadata']['object_ids']
        _hide_header_attachments(native_page, ids, page['hidden_marks'])
        next_id=max(int(n.get('id','0')) for n in root.iter())+1
        for entry in page['entries']:
            record=compounds[entry['cid']]; key=entry['key']
            group=ET.SubElement(native_page,'group',{'id':str(next_id)});next_id+=1
            object_ids={'group':group.get('id'),'fragment':ids.get(key+'_mol'),'label':ids[key+'_label'],
                        'fields':{f['column']:ids[key+'_'+f['column']] for f in spec.get('fields',[])
                                  if key+'_'+f['column'] in ids}}
            joins=[ids[k] for k in ids if k.startswith(key+'_join')]  # the ", " of a shared line
            for oid in [o for o in [object_ids['fragment'],object_ids['label'],*object_ids['fields'].values(),*joins] if o]:
                node=next(n for n in native_page if n.get('id')==oid)
                native_page.remove(node);group.append(node)
            if entry.get('named'):
                next_id=_draw_as_label(group.find('fragment'),record['input']['smiles'],entry['named'],next_id)
            record['objects']=object_ids
            # A compound drawn with another's structure has none of its own to read back.
            record['displayed_identity']=None if entry.get('text_only') else semantic_key(fragment_mol(group.find('fragment')))
        tree.write(stage/name,encoding='utf-8',xml_declaration=True)
        pages.append({'file':name,'width_pt':page['width'],'height_pt':page['height'],
                      **({'grid':page['grid']} if page.get('grid') else {})})
    return pages,compounds,warnings


ACS_LABEL_SIZE_FALLBACK = 10.0


def _draw_as_label(fragment, smiles, named, next_id):
    """Contract a named substituent into one labelled node that keeps it whole (ChemDraw's
    NodeType="Fragment" abbreviation), as house-style abbreviations are drawn."""
    from ..mcp_runtime.figure_house import _atom_nodes, _contract
    from ..render.renderer import _IDGen
    src = Chem.MolFromSmiles(smiles)
    nodes = _atom_nodes(fragment, fragment_mol(fragment), src)
    if nodes is None:
        return next_id
    root = named['atoms'][0]
    expansion = Chem.MolFragmentToSmiles(src, atomsToUse=list(named['atoms']), rootedAtAtom=root, canonical=False)
    ids = _IDGen(next_id)
    _contract(fragment, nodes, named['atoms'], named['outside'], expansion, named['label'], named['label'], ids)
    return max(int(n.get('id')) for n in fragment.iter() if n.get('id') and n.get('id').isdigit()) + 1


def bounds(node, style):
    if node.tag=='t':
        coords=list(map(float,node.get('p','0 0').split()))
        size=max([float(s.get('size',style['font_size'])) for s in node.findall('s')]+[style['font_size']])
        text=node_text(node)
        def run_width(run):
            # As drawn: sub- and superscripts smaller, bold runs wider (face bits 32, 64, 1).
            face=int(run.get('face','0') or 0) if (run.get('face') or '0').isdigit() else 0
            run_size=float(run.get('size',size))*(SCRIPT_SCALE if face&96 and face&96!=96 else 1.0)
            return text_width(run.text or '',style,run_size,bold=bool(face&1))
        runs=node.findall('s')
        width=(max(text_width(line,style,size) for line in text.splitlines()) if '\n' in text or not runs
               else sum(run_width(r) for r in runs))
        left=coords[0]
        if node.get('Justification')=='Center': left-=width/2
        elif node.get('Justification')=='Right':left-=width
        return [left,coords[1]-size,left+width,coords[1]+max(0,text.count('\n'))*size*1.5+size*.3]
    points=[]
    for n in node.findall('n'):   # an abbreviation's inner fragment is not drawn
        if n.get('p'): points.append(list(map(float,n.get('p').split()))[:2])
    if not points: return None
    pad=style['font_size']*1.2
    return [min(p[0] for p in points)-pad,min(p[1] for p in points)-pad,
            max(p[0] for p in points)+pad,max(p[1] for p in points)+pad]


def overlap(a,b):
    return min(a[2],b[2])-max(a[0],b[0])>1 and min(a[3],b[3])-max(a[1],b[1])>1


def check_pages(project, roots, spec):
    style=style_settings(spec); issues=[]; errors=[]
    expected=[Counter() for _ in roots]; boxes=[[] for _ in roots]
    if roots and spec.get('header'):
        expected[0].update(header_identities(spec))  # the general reaction over the first grid
    for cid,record in project['compounds'].items():
        index=record['page']; root=roots[index]; ids=record['objects']
        group=root.find(f".//group[@id='{ids['group']}']")
        if group is None:
            errors.append({'page':index+1,'compound_id':cid,'error':'missing_group'});continue
        if ids.get('fragment') is None and record['input'].get('drawn_with'):
            # Its structure is the one it is drawn with; only its own lines are bound here.
            for field,oid in {'label':ids['label'],**ids['fields']}.items():
                node=group.find(f"t[@id='{oid}']")
                if node is None or node_text(node)!=record['display'][field]:
                    errors.append({'page':index+1,'compound_id':cid,'error':'label_binding_changed','field':field})
            for element in group.findall('t'):
                box=bounds(element,style)
                if box: boxes[index].append((cid,element.get('id'),box))
            continue
        fragment=group.find(f"fragment[@id='{ids['fragment']}']")
        if fragment is None:
            errors.append({'page':index+1,'compound_id':cid,'error':'missing_structure'});continue
        mol=fragment_mol(fragment)
        if semantic_key(mol)!=record['displayed_identity']:
            errors.append({'page':index+1,'compound_id':cid,'error':'structure_changed'})
        expected[index].update(semantic_key(m) for m in Chem.GetMolFrags(Chem.MolFromSmiles(record['displayed_identity']),asMols=True))
        for field,oid in {'label':ids['label'],**ids['fields']}.items():
            node=group.find(f"t[@id='{oid}']")
            if node is None or node_text(node)!=record['display'][field]:
                errors.append({'page':index+1,'compound_id':cid,'error':'label_binding_changed','field':field})
        elements=[fragment,*group.findall('t')]
        for element in elements:
            box=bounds(element,style)
            if box:
                boxes[index].append((cid,element.get('id'),box))
        for run in group.iter('s'):
            if float(run.get('size',style['font_size']))<7:
                issues.append({'page':index+1,'compound_id':cid,'issue':'text_below_7pt'})
    for index,root in enumerate(roots):
        managed={oid for cid,oid,box in boxes[index]}
        parents={child:parent for parent in root.iter() for child in parent}
        for node in root.iter('t'):
            if node.get('id') in managed:continue
            parent=parents.get(node)
            in_structure=False
            while parent is not None:
                if parent.tag in {'fragment','n'}:in_structure=True;break
                parent=parents.get(parent)
            if in_structure:continue
            box=bounds(node,style)
            if box:boxes[index].append(('annotation',node.get('id'),box))
            if any(float(run.get('size',style['font_size']))<7 for run in node.findall('s')):
                issues.append({'page':index+1,'object_id':node.get('id'),'issue':'text_below_7pt'})
        actual=(inventory_keys(read_molecules(ET.tostring(root,encoding='unicode')),expected[index])
                if root.find('.//fragment') is not None else Counter())
        if actual!=expected[index]: errors.append({'page':index+1,'error':'document_inventory_changed'})
        page=project['pages'][index]
        for cid,oid,box in boxes[index]:
            if box[0]<0 or box[1]<0 or box[2]>page['width_pt'] or box[3]>page['height_pt']:
                issues.append({'page':index+1,'compound_id':cid,'object_id':oid,'issue':'outside_page'})
        for i,(cid,oid,box) in enumerate(boxes[index]):
            for other,oid2,box2 in boxes[index][i+1:]:
                if overlap(box,box2): issues.append({'page':index+1,'compound_id':cid,'other_compound_id':other,
                                                   'objects':[oid,oid2],'issue':'overlap'})
    return {'data_binding':{'status':'failed' if errors else 'passed','errors':errors},
            'chemistry':{'status':'failed' if errors else 'passed','scope':'rdkit_readback_consistency'},
            'layout':{'status':'review_required' if issues else 'passed','issues':issues,
                      'measurement':'font_metrics_and_structure_bounds; confirm native pixels'},
            'native_preview':{'status':'not_run'}}
