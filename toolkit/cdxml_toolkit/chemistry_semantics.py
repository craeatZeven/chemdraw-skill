"""Output-derived chemical identity. Never infer success from drawing metadata."""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from pathlib import Path
from xml.etree import ElementTree as ET
from rdkit import Chem


def semantic_key(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    copy = Chem.RemoveHs(copy)
    # Canonicalize groups, including the equivalent explicit/implicit ABS case.
    copy = Chem.CanonicalizeStereoGroups(copy)
    flags = (
        Chem.CXSmilesFields.CX_ENHANCEDSTEREO
        | Chem.CXSmilesFields.CX_COORDINATE_BONDS
        | Chem.CXSmilesFields.CX_BOND_ATROPISOMER
    )
    return Chem.MolToCXSmiles(copy, Chem.SmilesWriteParams(), flags)


# Elements that are not metals; any other element alone and uncharged is the element itself.
_NONMETALS = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86}


def _elemental_metal(atom):
    """A lone, neutral metal atom ("[Mg]" for magnesium turnings): RDKit counts its
    free valences as radical electrons, but it is drawn and read back as the element."""
    return (atom.GetDegree() == 0 and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 0
            and atom.GetAtomicNum() not in _NONMETALS)


def require_supported(
    mol,
    *,
    allow_radicals: bool = False,
    allow_dative: bool = True,
    allow_atropisomers: bool = True,
):
    if mol is None:
        raise ValueError("Invalid chemical structure.")
    allowed = {Chem.ChiralType.CHI_UNSPECIFIED, Chem.ChiralType.CHI_TETRAHEDRAL_CW,
               Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
    if any(a.GetChiralTag() not in allowed for a in mol.GetAtoms()):
        raise ValueError('unsupported_stereochemistry: non-tetrahedral center needs a native template')
    if not allow_radicals and any(a.GetNumRadicalElectrons() and not _elemental_metal(a) for a in mol.GetAtoms()):
        raise ValueError('unsupported_radical: use a verified native template; radical export is not yet validated')
    if not allow_dative and any(
        b.GetBondType() in (
            Chem.BondType.DATIVE,
            Chem.BondType.DATIVEL,
            Chem.BondType.DATIVER,
            Chem.BondType.DATIVEONE,
        )
        for b in mol.GetBonds()
    ):
        raise ValueError('unsupported_coordination_bond: dative bond needs coordination support')
    if not allow_atropisomers and any(
        b.GetStereo() in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW)
        for b in mol.GetBonds()
    ):
        raise ValueError('unsupported_stereochemistry: atropisomer needs a native template')


def _open_double_bonds(mol) -> bool:
    """Whether the molecule has a double bond that could be E or Z but is given as neither."""
    try:
        return any(info.type == Chem.StereoType.Bond_Double and info.specified == Chem.StereoSpecified.Unspecified
                   for info in Chem.FindPotentialStereo(Chem.Mol(mol)))
    except Exception:
        return False


def loose_key(mol):
    """The semantic key without double-bond geometry (E/Z)."""
    copy = Chem.Mol(mol)
    for bond in copy.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
    return semantic_key(copy)


def same_compound(source, drawn):
    """"exact", "drawn_geometry" or None: whether a drawn molecule is the source compound.

    A drawn double bond has a geometry: RDKit reads E/Z off the coordinates of a C=N or C=C,
    so an oxime given as CC(=NO)... comes back as C/C(=N\\O)... When the source leaves some
    double bond's geometry open, a drawing that differs from it only in E/Z is the compound
    asked for, drawn one way ("drawn_geometry"). A geometry the source fixes must come back
    the same.
    """
    if semantic_key(source) == semantic_key(drawn):
        return "exact"
    if _open_double_bonds(source) and loose_key(source) == loose_key(drawn):
        return "drawn_geometry"
    return None


def inventory_keys(mols, like=None):
    """Semantic keys of drawn molecules; with ``like`` (the keys expected), a molecule that is one
    of them drawn with a geometry its source left open counts as that key."""
    found = Counter()
    open_keys = {}
    for key in set(like or ()):
        mol = Chem.MolFromSmiles(key)
        if mol is not None and _open_double_bonds(mol):
            open_keys.setdefault(loose_key(mol), key)
    for mol in mols:
        key = semantic_key(mol)
        if like is not None and key not in like and open_keys:
            key = open_keys.get(loose_key(mol), key)
        found[key] += 1
    return found


def _remove_ring_r_groups(fragment):
    """Remove native variable-position R marks from a readback copy of a fragment."""
    anchors = {n.get('id'): n for n in fragment.findall('n')
               if n.get('NodeType') == 'MultiAttachment' and n.get('Attachments')}
    if not anchors:
        return
    labels = {n.get('id'): n for n in fragment.findall('n')
              if n.get('NodeType') == 'GenericNickname'
              and n.get('GenericNickname', '')[:1] == 'R'
              and n.get('GenericNickname', '')[1:].isdigit()}
    bonds = fragment.findall('b')
    incident = Counter(endpoint for b in bonds for endpoint in (b.get('B'), b.get('E')))
    for bond in bonds:
        anchor, label = bond.get('B'), bond.get('E')
        if anchor not in anchors:
            anchor, label = label, anchor
        if (anchor in anchors and label in labels and incident[anchor] == incident[label] == 1
                and bond.get('Order', '1') == '1'):
            fragment.remove(anchors[anchor])
            fragment.remove(labels[label])
            fragment.remove(bond)


def read_fragments(value):
    """Read top-level CDXML fragments in document order, using the scaffold for variable R groups."""
    root = ET.parse(str(value)).getroot() if isinstance(value, Path) else ET.fromstring(value)
    # Reaction text IDs are not molecules. Parse each top-level fragment separately.
    parents = {child: parent for parent in root.iter() for child in parent}
    results = []
    for fragment in root.iter('fragment'):
        parent = parents.get(fragment)
        nested = False
        while parent is not None:
            if parent.tag == 'fragment':
                nested = True
                break
            parent = parents.get(parent)
        if nested:
            continue
        if fragment.find('n[@NodeType="MultiAttachment"]') is not None:
            fragment = deepcopy(fragment)
            _remove_ring_r_groups(fragment)
        block = '<CDXML><page id="90000000">' + ET.tostring(fragment, encoding='unicode') + '</page></CDXML>'
        parsed = Chem.MolsFromCDXML(block)
        if not parsed:
            raise ValueError('unreadable_chemical_fragment: ' + fragment.get('id', '?'))
        nodes = fragment.findall('n')
        if len(parsed) == 1 and len(nodes) == 1 and parsed[0].GetNumAtoms() == 1 \
                and nodes[0].get('NumHydrogens', '').isdigit():
            # RDKit's reader fills a lone atom's valence with H and ignores NumHydrogens:
            # magnesium metal would come back as MgH2.
            mol = Chem.RWMol(parsed[0])
            atom = mol.GetAtomWithIdx(0)
            atom.SetNumExplicitHs(int(nodes[0].get('NumHydrogens')))
            atom.SetNoImplicit(True)
            rad_attr = nodes[0].get('Radical', '')
            if rad_attr == 'Doublet':
                atom.SetNumRadicalElectrons(1)
            elif rad_attr in ('Triplet', 'Diradical'):
                atom.SetNumRadicalElectrons(2)
            Chem.SanitizeMol(mol)
            if rad_attr == 'Doublet':
                mol.GetAtomWithIdx(0).SetNumRadicalElectrons(1)
            elif rad_attr in ('Triplet', 'Diradical'):
                mol.GetAtomWithIdx(0).SetNumRadicalElectrons(2)
            parsed = [mol.GetMol()]
        results.extend((fragment.get('id'), mol) for mol in parsed)
    if not results:
        raise ValueError('no_chemical_fragments')
    return results


def read_molecules(value):
    return [mol for _, mol in read_fragments(value)]


def validate_fragment(smiles, xml, *, allow_radicals: bool = True):
    source = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else Chem.Mol(smiles)
    if source is None:
        raise ValueError('invalid_source_smiles')
    require_supported(source, allow_radicals=allow_radicals)
    actual = read_molecules('<CDXML><page id="90000000">' + xml + '</page></CDXML>')
    combined = actual[0]
    for mol in actual[1:]:
        combined = Chem.CombineMols(combined, mol)
    expected, observed = semantic_key(source), semantic_key(combined)
    verdict = same_compound(source, combined)
    if verdict is None:
        raise ValueError(f'chemical_semantics_changed: source={expected}; roundtrip={observed}')
    return {'status': 'preserved', 'source_cxsmiles': expected,
            'roundtrip_cxsmiles': observed, 'method': 'final_cdxml_roundtrip',
            'stereo_groups': len(source.GetStereoGroups()),
            **({'double_bond_geometry': 'drawn'} if verdict == 'drawn_geometry' else {})}


def document_inventory(path, like=None):
    """Counter of the semantic keys of a document's molecules (see ``inventory_keys`` for ``like``)."""
    source = Path(path) if not str(path).lstrip().startswith('<') else str(path)
    return inventory_keys(read_molecules(source), like)


def validate_document_preserved(source, output):
    before, after = document_inventory(source), document_inventory(output)
    if before != after:
        raise ValueError('chemical_semantics_changed: layout operation changed the molecule inventory')
    return {'status': 'preserved', 'method': 'final_cdxml_roundtrip', 'molecules': sum(after.values())}
