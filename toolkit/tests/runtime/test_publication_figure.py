"""Publication workflows use the grounded structures from existing figure tests."""
import csv
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

STEREO = 'N[C@@H]1CC[C@H](O)C1'
OTHER = 'CC[C@H](F)[C@H](C)O |&1:2,4|'


def request(tmp_path, rows=None, **options):
    rows = rows or [['a', '1', STEREO, '82%', 'amine'], ['b', '2', OTHER, '<10 nM', 'alcohol']]
    table = tmp_path/'data.csv'
    with table.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['compound_id','label','smiles','value','group'])
        writer.writerows(rows)
    spec = {'table': str(table), 'mode':'scope', 'fields':[{'column':'value','label':'Result'}],
            'group_column':'group', 'native_render':False, **options}
    path = tmp_path/'spec.json'; path.write_text(json.dumps(spec))
    return path


def execute(spec, output, operation='create'):
    from cdxml_toolkit.mcp_runtime.publication_figure import publication_figure
    return publication_figure(operation, str(spec), str(output))


def test_create_binding_pagination_and_reproducibility(tmp_path):
    rows = [[f'id{i}', str(i), STEREO if i % 2 else OTHER,
             '' if i == 2 else f'>{i} nM', 'Series'] for i in range(48)]
    spec = request(tmp_path, rows, mode='sar')
    a,b=tmp_path/'a',tmp_path/'b'
    assert execute(spec,a)['ok']
    assert execute(spec,b)['ok']
    project=json.loads((a/'project.json').read_text(encoding="utf-8"))
    assert len(project['compounds']) == 48
    assert len(project['pages']) > 1
    assert project['compounds']['id2']['input']['values']['value'] == ''
    for page in project['pages']:
        assert (a/page['file']).read_bytes() == (b/page['file']).read_bytes()
    assert json.loads((a/'checks.json').read_text(encoding="utf-8"))['data_binding']['status']=='passed'
    assert not json.loads((a/'checks.json').read_text(encoding="utf-8"))['layout']['issues']


@pytest.mark.parametrize('rows', [
    [['a','1',STEREO,'1',''],['a','2',OTHER,'2','']],
    [['a','1','not a structure','1','']],
])
def test_invalid_data_never_published(tmp_path,rows):
    with pytest.raises(ValueError, match='row'):
        execute(request(tmp_path,rows),tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_metric_update_preserves_manual_layout_and_annotation(tmp_path):
    spec=request(tmp_path)
    base=tmp_path/'base'; execute(spec,base)
    project=json.loads((base/'project.json').read_text(encoding="utf-8")); record=project['compounds']['a']
    page=record['page']; tree=ET.parse(base/project['pages'][page]['file'])
    group=tree.find(f".//group[@id='{record['objects']['group']}']")
    from cdxml_toolkit.mcp_runtime.figure_tools import _move
    _move(group,5,7)
    label=group.find(f"t[@id='{record['objects']['label']}']")
    label.find('s').set('face','1')
    ET.SubElement(tree.find('page'),'t',{'id':'99001','p':'30 25'}).append(ET.Element('s',{'font':'3','size':'8'}))
    tree.find(".//t[@id='99001']/s").text='Keep this annotation'
    edited=tmp_path/'edited.cdxml';tree.write(edited)
    before=ET.tostring(group.find('fragment'))
    spec=request(tmp_path,[['b','2',OTHER,'<10 nM','alcohol'],['a','1',STEREO,'91%','amine']],
                 project_path=str(base/'project.json'),edited_pages={str(page):str(edited)})
    out=tmp_path/'updated'; result=execute(spec,out,'update'); assert result['ok']
    after=ET.parse(out/project['pages'][page]['file'])
    assert ET.tostring(after.find(f".//group[@id='{record['objects']['group']}']/fragment"))==before
    assert '91%' in ''.join(after.getroot().itertext())
    assert 'Keep this annotation' in ''.join(after.getroot().itertext())
    assert after.find(f".//t[@id='{record['objects']['label']}']/s").get('face')=='1'


def test_concurrent_value_edit_is_conflict_draft(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8"));r=p['compounds']['a']
    tree=ET.parse(base/p['pages'][r['page']]['file'])
    tree.find(f".//t[@id='{r['objects']['fields']['value']}']/s").text='Result: 77%'
    edited=tmp_path/'edited.cdxml';tree.write(edited)
    spec=request(tmp_path,[['a','1',STEREO,'91%','amine'],['b','2',OTHER,'<10 nM','alcohol']],
                 project_path=str(base/'project.json'),edited_pages={str(r['page']):str(edited)})
    out=tmp_path/'update'; result=execute(spec,out,'update')
    assert not result['ok'] and result['status']=='conflict'
    assert (out/'conflicts.json').is_file()
    assert not (out/'project.json').exists()
    assert list((out/'draft').glob('*.cdxml'))


def test_renumbered_native_ids_use_unique_identity_and_label(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8")); tree=ET.parse(base/p['pages'][0]['file'])
    mapping={n.get('id'):str(int(n.get('id'))+100000) for n in tree.iter() if n.get('id')}
    for node in tree.iter():
        if node.get('id'):node.set('id',mapping[node.get('id')])
        for attr in ('B','E','BondOrdering'):
            if node.get(attr):node.set(attr,' '.join(mapping.get(v,v) for v in node.get(attr).split()))
    path=tmp_path/'renumbered.cdxml';tree.write(path)
    spec=request(tmp_path,project_path=str(base/'project.json'),edited_pages={'0':str(path)})
    assert execute(spec,tmp_path/'updated','update')['ok']


def test_missing_group_and_ambiguous_duplicate_are_conflicts(tmp_path):
    from copy import deepcopy
    spec=request(tmp_path,[['a','1',STEREO,'1',''],['b','1',STEREO,'1','']])
    base=tmp_path/'base';execute(spec,base);p=json.loads((base/'project.json').read_text(encoding="utf-8"))
    tree=ET.parse(base/p['pages'][0]['file']);page=tree.find('page')
    for group in page.findall('group'):group.set('id',str(int(group.get('id'))+9000))
    edited=tmp_path/'edited.cdxml';tree.write(edited)
    spec=request(tmp_path,[['a','1',STEREO,'1',''],['b','1',STEREO,'1','']],project_path=str(base/'project.json'),edited_pages={'0':str(edited)})
    assert not execute(spec,tmp_path/'conflict','update')['ok']


def test_new_and_deleted_compounds_keep_existing_pages(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8"))
    spec=request(tmp_path,[['a','1',STEREO,'82%','amine'],['c','3',OTHER,'33%','alcohol']],project_path=str(base/'project.json'))
    out=tmp_path/'updated';assert execute(spec,out,'update')['ok']
    q=json.loads((out/'project.json').read_text(encoding="utf-8"))
    assert set(q['compounds'])=={'a','c'}
    assert q['compounds']['a']['page']==p['compounds']['a']['page']
    assert q['compounds']['c']['page']>=len(p['pages'])


def test_mol_and_sdf_identity_binding(tmp_path):
    from rdkit import Chem
    mol=Chem.MolFromSmiles(STEREO)
    Chem.MolToMolFile(mol,str(tmp_path/'a.mol'))
    spec=request(tmp_path,[['a','1','a.mol','82%','']],structure={'format':'mol','column':'smiles'})
    assert execute(spec,tmp_path/'mol')['ok']
    writer=Chem.SDWriter(str(tmp_path/'library.sdf'));mol.SetProp('compound_id','a');writer.write(mol);writer.close()
    spec=request(tmp_path,[['a','1','','82%','']],structure={'format':'sdf','file':'library.sdf'})
    assert execute(spec,tmp_path/'sdf')['ok']


def test_structure_replacement_has_explicit_core_map(tmp_path):
    # Existing aspirin fixture and its resolved salicylic-acid substructure.
    old='CC(=O)Oc1ccccc1C(=O)O'
    new='O=C(O)c1ccccc1O'
    from rdkit import Chem
    a,b=Chem.MolFromSmiles(old),Chem.MolFromSmiles(new)
    match=a.GetSubstructMatch(b)
    spec=request(tmp_path,[['a','1',old,'82%','']]);base=tmp_path/'base';execute(spec,base)
    spec=request(tmp_path,[['a','1',new,'82%','']],project_path=str(base/'project.json'),atom_maps={'a':[[i,j] for j,i in enumerate(match)]})
    out=tmp_path/'updated';result=execute(spec,out,'update')
    assert result['ok'], (out/'conflicts.json').read_text(encoding="utf-8") if (out/'conflicts.json').exists() else result
    q=json.loads((out/'project.json').read_text(encoding="utf-8"))
    from cdxml_toolkit.chemistry_semantics import semantic_key
    assert q['compounds']['a']['displayed_identity']==semantic_key(b)


def test_check_rejects_tampered_measurement(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8"));r=p['compounds']['a']
    tree=ET.parse(base/p['pages'][0]['file']);tree.find(f".//t[@id='{r['objects']['fields']['value']}']/s").text='wrong'
    edited=tmp_path/'edited.cdxml';tree.write(edited)
    spec=request(tmp_path,project_path=str(base/'project.json'),edited_pages={'0':str(edited)})
    assert not execute(spec,tmp_path/'checked','check')['ok']


def test_xlsx_cached_formula_and_missing_cache(tmp_path):
    # Minimal OOXML fixture: explicitly control the formula cache without Excel.
    import zipfile
    from xml.sax.saxutils import escape
    from cdxml_toolkit.publication.data import table_rows
    path=tmp_path/'data.xlsx'
    def workbook(cached):
        with zipfile.ZipFile(path,'w') as z:
            z.writestr('[Content_Types].xml','<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>')
            z.writestr('_rels/.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
            z.writestr('xl/workbook.xml','<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>')
            z.writestr('xl/_rels/workbook.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>')
            z.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>value</t></is></c></row><row r="2"><c r="A2"><f>1+1</f>'+('<v>2</v>' if cached else '')+'</c></row></sheetData></worksheet>')
    pytest.importorskip('openpyxl')
    workbook(False)
    with pytest.raises(ValueError,match='Data!A2'):table_rows(path)
    workbook(True)
    assert table_rows(path)[1][0][1]['value']=='2'


def test_unresolved_chemistry_persists_and_existing_output_is_protected(tmp_path):
    spec=request(tmp_path,unresolved_chemistry=['reader disagreement for compound a'])
    base=tmp_path/'base';execute(spec,base)
    before=(base/'project.json').read_bytes()
    with pytest.raises(ValueError,match='overwrite'):execute(spec,base)
    assert (base/'project.json').read_bytes()==before
    spec=request(tmp_path,project_path=str(base/'project.json'))
    out=tmp_path/'updated';execute(spec,out,'update')
    assert json.loads((out/'checks.json').read_text(encoding="utf-8"))['chemistry']['status']=='unresolved'


def test_locked_overflow_is_reported_without_moving(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8"));tree=ET.parse(base/p['pages'][0]['file'])
    from cdxml_toolkit.mcp_runtime.figure_tools import _move
    _move(tree.find('page/group'),1000,0)
    edited=tmp_path/'edited.cdxml';tree.write(edited)
    spec=request(tmp_path,project_path=str(base/'project.json'),edited_pages={'0':str(edited)})
    out=tmp_path/'updated';result=execute(spec,out,'update')
    assert any(i['issue']=='outside_page' for i in result['checks']['layout']['issues'])


def test_coordinate_bearing_cxsmiles_does_not_reuse_old_wedges(tmp_path):
    from rdkit import Chem
    from cdxml_toolkit.chemistry_semantics import read_molecules
    source=Path(__file__).resolve().parents[2]/'samples/publication-figures/stereo-reaction.cdxml'
    mol=read_molecules(source)[0]
    spec=request(tmp_path,[['a','1',Chem.MolToCXSmiles(mol),'82%','']])
    assert execute(spec,tmp_path/'out')['ok']


def test_xlsx_preserves_decimal_percent_and_explicit_units(tmp_path):
    xl=pytest.importorskip('openpyxl')
    from cdxml_toolkit.publication.data import table_rows
    w=xl.Workbook();s=w.active;s.append(['percent','activity','missing'])
    s.append([0.397,12.3,None]);s['A2'].number_format='0.0%'
    s['B2'].number_format='0.00 "µM"'
    p=tmp_path/'units.xlsx';w.save(p)
    row=table_rows(p)[1][0][1]
    assert row=={'percent':'39.7%','activity':'12.30 µM','missing':''}


def test_project_page_files_cannot_leave_the_project_folder(tmp_path):
    spec=request(tmp_path);base=tmp_path/'base';execute(spec,base)
    p=json.loads((base/'project.json').read_text(encoding="utf-8"))
    (tmp_path/'outside.cdxml').write_bytes((base/p['pages'][0]['file']).read_bytes())
    p['pages'][0]['file']='../outside.cdxml'
    (base/'project.json').write_text(json.dumps(p))
    spec=request(tmp_path,project_path=str(base/'project.json'))
    with pytest.raises(ValueError,match='page file'):
        execute(spec,tmp_path/'checked','check')


def _scope_page(tmp_path, rows, **options):
    out = tmp_path/'scope'
    result = execute(request(tmp_path, rows, **options), out)
    assert result['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    page = project['pages'][0]
    return project, page, ET.parse(out/page['file']).getroot()


def test_a_scope_panel_keeps_10_pt_type_and_ends_a_margin_past_its_content(tmp_path):
    rows = [['a', '3a', STEREO, '92%', ''], ['b', '3b', OTHER, '0%', '']]
    project, page, root = _scope_page(tmp_path, rows, footnotes=['Isolated yields.'])
    assert root.get('LabelSize') == '10' and root.get('CaptionSize') == '10'
    # The page is cropped to what it holds, not left at the default 180 x 240 mm.
    assert page['height_pt'] < 240 * 72 / 25.4 / 2
    ys = [float(n.get('p').split()[1]) for n in root.iter() if n.get('p')]
    assert page['height_pt'] - max(ys) < 40
    # The footnote follows the last row instead of sitting at the foot of a full page.
    texts = {''.join(t.itertext()): float(t.get('p').split()[1]) for t in root.iter('t') if t.get('p')}
    assert texts['Isolated yields.'] == max(texts.values())
    assert 'Result: 0%' in texts  # a failed entry keeps its value


def test_an_optional_note_is_drawn_only_where_it_is_given(tmp_path):
    rows = [['a', '3a', STEREO, '0%', '', 'not formed'], ['b', '3b', OTHER, '63%', '', '']]
    table = tmp_path/'notes.csv'
    with table.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['compound_id', 'label', 'smiles', 'value', 'group', 'note'])
        writer.writerows(rows)
    spec = {'table': str(table), 'mode': 'scope', 'native_render': False,
            'fields': [{'column': 'value', 'label': ''}, {'column': 'note', 'label': '', 'optional': True}]}
    path = tmp_path/'notes.json'; path.write_text(json.dumps(spec))
    out = tmp_path/'notes'
    assert execute(path, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    assert project['compounds']['a']['display'] == {'label': '3a', 'value': '0%', 'note': 'not formed'}
    assert project['compounds']['b']['display'] == {'label': '3b', 'value': '63%'}
    assert set(project['compounds']['b']['objects']['fields']) == {'value'}
    checks = json.loads((out/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed' and not checks['layout']['issues']


def test_one_scope_field_reads_as_the_value_alone():
    from cdxml_toolkit.mcp_runtime.figure_scope import _cell_fields
    yield_field = {'key': 'yield', 'label': 'Yield', 'unit': None, 'show_label': False}
    assert _cell_fields([yield_field], notes=False) == [{'column': 'v1', 'label': '', 'inline': True}]
    assert _cell_fields([{**yield_field, 'show_label': True}], notes=False) == [{'column': 'v1', 'label': 'Yield'}]
    ee = {'key': 'ee', 'label': 'ee', 'unit': None, 'show_label': False}
    assert [c['label'] for c in _cell_fields([yield_field, ee], notes=True)] == ['Yield', 'ee', '']
    assert _cell_fields([yield_field], notes=True)[-1] == {'column': 'note', 'label': '', 'optional': True}


def test_an_empty_optional_line_is_absent_and_an_empty_value_is_a_dash():
    from cdxml_toolkit.publication.data import field_text
    record = {'values': {'v1': '', 'note': ''}}
    assert field_text(record, {'column': 'v1', 'label': ''}) == '—'
    assert field_text(record, {'column': 'note', 'label': '', 'optional': True}) == ''


def test_an_update_that_lengthens_text_grows_the_cropped_page(tmp_path):
    spec = request(tmp_path); base = tmp_path/'base'; execute(spec, base)
    before = json.loads((base/'project.json').read_text(encoding='utf-8'))['pages'][0]['height_pt']
    longer = 'not formed under the standard conditions; the starting material was recovered'
    # The last compound's caption grows past the foot of the cropped page.
    spec = request(tmp_path, [['a', '1', STEREO, '82%', 'amine'], ['b', '2', OTHER, longer, 'alcohol']],
                   project_path=str(base/'project.json'))
    out = tmp_path/'updated'; result = execute(spec, out, 'update')
    assert result['ok']
    page = json.loads((out/'project.json').read_text(encoding='utf-8'))['pages'][0]
    assert page['height_pt'] > before
    assert not [i for i in result['checks']['layout']['issues'] if i['issue'] == 'outside_page']


def test_twelve_benzamides_with_two_fields_fit_one_page_at_10_pt(tmp_path):
    # Journal scope tables hold 18 to 22 compounds on a page at body-text size; a cell spends its
    # height on the structure and its lines, not on padding.
    rows = [[f'4{chr(97 + i)}', f'4{chr(97 + i)}', smiles, f'{90 - i}%', group]
            for i, (smiles, group) in enumerate(
                [(f'{x}c1ccc(cc1)C(=O)N', 'para') for x in ('C', 'CO', 'F', 'Cl', 'Br', 'FC(F)(F)')]
                + [(f'{x}c1cccc(c1)C(=O)N', 'meta') for x in ('C', 'CO', 'F', 'Cl', 'Br', 'FC(F)(F)')])]
    table = tmp_path/'twelve.csv'
    with table.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['compound_id', 'label', 'smiles', 'yield', 'group', 'ee'])
        writer.writerows([row + ['95%'] for row in rows])
    spec = {'table': str(table), 'mode': 'scope', 'native_render': False, 'group_column': 'group',
            'title': 'Substrate scope', 'conditions': 'Conditions as supplied by the user',
            'footnotes': ['a Isolated yield.'],
            'fields': [{'column': 'yield', 'label': 'Yield'}, {'column': 'ee', 'label': 'ee'}]}
    path = tmp_path/'twelve.json'; path.write_text(json.dumps(spec))
    out = tmp_path/'twelve'
    result = execute(path, out)
    assert result['ok'] and result['pages'] == 1
    assert not result['checks']['layout']['issues']
    root = ET.parse(out/'page-001.cdxml').getroot()
    assert root.get('LabelSize') == '10'


def test_cell_geometry_leaves_the_structure_box_clear_of_its_caption():
    from cdxml_toolkit.publication.layout import CELL
    size = 10.0
    # check_pages pads a structure by 1.2 sizes and draws a caption from 1 size above its
    # baseline to 0.3 below; it reports boxes that overlap by more than 1 pt. The first caption
    # must clear the padded structure, and each caption the one above it.
    assert CELL['caption_baseline'] * size - size >= 1.2 * size - 1.0
    assert CELL['line_advance'] * size + CELL['line_gap'] >= 1.3 * size + 1.0


def test_each_scope_group_is_named_on_its_rule(tmp_path):
    # Nature Communications PMC13594109, fig. 4: "— para substituents ———" — the group's name
    # sits on the rule, which stops short of it on both sides, and only the locant is italic.
    project, page, root = _scope_page(tmp_path, [['a', '3a', STEREO, '92%', 'para substituents'],
                                                 ['b', '3b', OTHER, '81%', 'N-heterocycles']])
    rules = [a for a in root.iter('arrow') if a.get('ArrowheadHead') == 'None'
             and a.get('Head3D').split()[1] == a.get('Tail3D').split()[1]]
    assert len(rules) == 4  # a short piece before each name and a long one after it
    labels = {''.join(t.itertext()): t for t in root.iter('t')}
    for name, italic in (('para substituents', 'para'), ('N-heterocycles', 'N')):
        label = labels[name]
        assert [s.text for s in label.iter('s') if int(s.get('face')) & 2] == [italic]
        x, baseline = (float(v) for v in label.get('p').split())
        pieces = [r for r in rules if abs(float(r.get('Tail3D').split()[1]) - (baseline - 3.5)) < 0.5]
        assert len(pieces) == 2
        ends = sorted((min(float(r.get('Tail3D').split()[0]), float(r.get('Head3D').split()[0])),
                       max(float(r.get('Tail3D').split()[0]), float(r.get('Head3D').split()[0]))) for r in pieces)
        assert ends[0][1] < x < ends[1][0]  # the name stands in the gap between the two pieces
    tops = sorted(min(float(n.get('p').split()[1]) for n in f.iter('n'))
                  for f in root.iter('fragment'))
    second_rule = sorted(float(r.get('Tail3D').split()[1]) for r in rules)[-1]
    assert sum(top > second_rule for top in tops) == 1  # one structure in each band
    checks = json.loads((tmp_path/'scope'/'checks.json').read_text(encoding='utf-8'))
    assert not checks['layout']['issues']


def _secondary_request(tmp_path, rows):
    table = tmp_path/'minor.csv'
    with table.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['compound_id', 'label', 'smiles', 'value', 'group', 'minor_of'])
        writer.writerows(rows)
    spec = {'table': str(table), 'mode': 'scope', 'native_render': False, 'secondary_column': 'minor_of',
            'fields': [{'column': 'value', 'label': ''}]}
    path = tmp_path/'minor.json'; path.write_text(json.dumps(spec))
    return path


def test_a_secondary_product_stands_beside_its_main_product_in_grey(tmp_path):
    # Nature Communications PMC13594109, fig. 4: 14b drawn in grey next to 14a, with its own yield.
    path = _secondary_request(tmp_path, [['a', '14a', STEREO, '61%', '', ''], ['b', '14b', OTHER, '22%', '', 'a'],
                                         ['c', '15', STEREO, '70%', '', '']])
    out = tmp_path/'minor'
    assert execute(path, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    a, b, c = (project['compounds'][k]['bbox'] for k in 'abc')
    assert b[1] == a[1] and 0 < b[0] - a[2] <= 24      # on the same row, right of its main product
    assert c[0] >= b[2] + 12                           # the next cell starts after the pair
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    table = [tuple(round(float(x.get(k)) * 255) for k in 'rgb') for x in root.find('colortable').findall('color')]
    grey = str(table.index((0x80, 0x80, 0x80)) + 2) if (0x80, 0x80, 0x80) in table else 'missing'
    group_b = root.find(f".//group[@id='{project['compounds']['b']['objects']['group']}']")
    # Only the structure is grey; its number and yield stay black, as in fig. 4 (14b, 52% yield).
    assert {e.get('color') for e in group_b.find('fragment').iter() if e.get('color')} == {grey}
    assert grey not in {s.get('color') for t in group_b.findall('t') for s in t.iter('s')}
    plus = [t for t in root.iter('t') if ''.join(t.itertext()) == '+']
    assert len(plus) == 1 and a[2] < float(plus[0].get('p').split()[0]) < b[0]
    group_a = root.find(f".//group[@id='{project['compounds']['a']['objects']['group']}']")
    assert grey not in {e.get('color') for e in group_a.iter()}
    checks = json.loads((out/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed' and not checks['layout']['issues']


@pytest.mark.parametrize('rows, message', [
    ([['a', '1', STEREO, '1%', '', 'z']], 'secondary of an unknown compound'),
    ([['a', '1', STEREO, '1%', '', 'a']], 'its own secondary'),
    ([['a', '1', STEREO, '1%', '', ''], ['b', '2', OTHER, '2%', '', 'a'], ['c', '3', STEREO, '3%', '', 'b']],
     'itself a secondary'),
    ([['a', '1', STEREO, '1%', 'g1', ''], ['b', '2', OTHER, '2%', 'g2', 'a']], 'same group'),
])
def test_a_secondary_product_must_name_a_main_product(tmp_path, rows, message):
    path = _secondary_request(tmp_path, rows)
    spec = json.loads(path.read_text()); spec['group_column'] = 'group'; path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match=message):
        execute(path, tmp_path/'bad')
    assert not (tmp_path/'bad').exists()


def test_scope_brief_fields_can_drop_their_name_and_their_dash():
    from cdxml_toolkit.mcp_runtime.figure_scope import _cell_fields, _fields
    fields, errors = _fields({'species': [], 'fields': [
        'yield', {'key': 'group', 'show_label': False, 'optional': True}]})
    assert errors == []
    assert _cell_fields(fields, notes=False) == [
        {'column': 'v1', 'label': 'yield'}, {'column': 'v2', 'label': '', 'optional': True}]  # two fields: own lines
    _, errors = _fields({'species': [], 'fields': [{'key': 'group', 'optional': 'yes'}]})
    assert any('optional must be true or false' in e for e in errors)


def test_a_definition_line_follows_the_number_without_a_field_name():
    from cdxml_toolkit.mcp_runtime.figure_scope import _cell_fields
    yield_field = {'key': 'yield', 'label': 'Yield', 'unit': None, 'show_label': False}
    assert _cell_fields([yield_field], notes=False, definitions=True) == [
        {'column': 'definition', 'label': '', 'optional': True, 'inline': True},
        {'column': 'v1', 'label': '', 'inline': True}]


def test_a_scope_definition_is_checked():
    from cdxml_toolkit.mcp_runtime.figure_scope import _validate
    brief = {'kind': 'scope', 'fields': ['yield'], 'species': [
        {'id': '3p', 'smiles': 'CCCCCCCc1ccccc1', 'definition': 'Alk = C7H15', 'values': {'yield': '63%'}},
        {'id': '3q', 'smiles': 'CCc1ccccc1', 'definition': 5, 'values': {'yield': '60%'}}]}
    assert _validate(brief) == ['Species 3q: definition must be text (a line such as "Alk = C7H15").']
    clash = {'kind': 'scope', 'fields': ['definition'], 'species': [{'id': 'a', 'smiles': 'C'}]}
    assert any('reserved' in e for e in _validate(clash))


def test_formula_digits_in_a_caption_line_are_subscripts():
    from cdxml_toolkit.publication.layout import chemical_runs
    assert chemical_runs('Alk = C7H15') == [
        {'text': 'Alk = C'}, {'text': '7', 'subscript': True}, {'text': 'H'}, {'text': '15', 'subscript': True}]
    assert chemical_runs('92%') is None and chemical_runs('dr 3:1') is None and chemical_runs('L1') is None


def test_a_family_caption_is_drawn_as_its_own_line_and_clears_the_next_row(tmp_path):
    # Witness: the scope figure whose 3p/3q cells read "Alk = C7H15" and "Alk = C18H37" under one structure.
    rows = [['p', '3p', 'CCCCCCCc1ccccc1', '63%', '', 'Alk = C7H15'],
            ['q', '3q', 'CCCCCCCCCCCCCCCCCCc1ccccc1', '60%', '', 'Alk = C18H37'],
            ['r', '3r', 'Cc1ccccc1', '71%', '', '']]
    table = tmp_path/'family.csv'
    with table.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['compound_id', 'label', 'smiles', 'value', 'group', 'definition'])
        writer.writerows(rows)
    spec = {'table': str(table), 'mode': 'scope', 'native_render': False, 'width_mm': 60,
            'fields': [{'column': 'definition', 'label': '', 'optional': True, 'inline': True},
                       {'column': 'value', 'label': '', 'inline': True}]}
    path = tmp_path/'family.json'; path.write_text(json.dumps(spec))
    out = tmp_path/'family'
    assert execute(path, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    assert project['compounds']['p']['display'] == {'label': '3p', 'definition': 'Alk = C7H15', 'value': '63%'}
    assert project['compounds']['r']['display'] == {'label': '3r', 'value': '71%'}
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    line = root.find(f".//t[@id='{project['compounds']['p']['objects']['fields']['definition']}']")
    runs = [(s.text, int(s.get('face')) & 32) for s in line.findall('s')]
    assert runs == [('Alk = C', 0), ('7', 32), ('H', 0), ('15', 32)]
    label = root.find(f".//t[@id='{project['compounds']['p']['objects']['label']}']")
    value = root.find(f".//t[@id='{project['compounds']['p']['objects']['fields']['value']}']")
    # BJOC 22-127 scheme 2: "3p, Alk = C7H15, 63%" is one line under the structure.
    assert label.get('p').split()[1] == line.get('p').split()[1] == value.get('p').split()[1]
    assert float(label.get('p').split()[0]) < float(line.get('p').split()[0]) < float(value.get('p').split()[0])
    texts = [t for t in root.iter('t') if ''.join(s.text or '' for s in t.iter('s')) not in ('', )]
    assert not any(''.join(s.text or '' for s in t.iter('s')).startswith(('definition:', '—')) for t in texts)
    boxes = sorted((c['bbox'] for c in project['compounds'].values()), key=lambda b: (b[1], b[0]))
    for upper in boxes:
        for lower in boxes:
            if lower[1] > upper[1]:
                assert lower[1] >= upper[3]          # a later row starts below the whole caption
    checks = json.loads((out/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed' and not checks['layout']['issues']


def _caption_line(root, project, cid):
    ids = project['compounds'][cid]['objects']
    return [root.find(f".//t[@id='{oid}']") for oid in (ids['label'], *ids['fields'].values())]


def test_a_number_and_its_yield_are_one_centred_line_in_bold_then_plain(tmp_path):
    # BJOC 22-127 scheme 2: "3a, 64%" under the structure, the number bold, the line centred.
    rows = [['a', '3a', STEREO, '64%', '']]
    out = tmp_path/'one'
    spec = request(tmp_path, rows, fields=[{'column': 'value', 'label': '', 'inline': True}])
    assert execute(spec, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    label, value = _caption_line(root, project, 'a')
    assert label.get('p').split()[1] == value.get('p').split()[1]
    assert {s.get('face') for s in label.iter('s')} == {'1'} and {s.get('face') for s in value.iter('s')} == {'0'}
    between = [t for t in root.iter('t') if ''.join(t.itertext()) == ', ']
    assert len(between) == 1 and float(label.get('p').split()[0]) < float(between[0].get('p').split()[0]) \
        < float(value.get('p').split()[0])
    from cdxml_toolkit.publication.layout import style_settings, text_width
    style = style_settings({})
    left = float(label.get('p').split()[0])
    right = float(value.get('p').split()[0]) + text_width('64%', style)
    xs = [float(n.get('p').split()[0]) for n in root.iter('n') if n.get('p')]
    assert abs((left + right) / 2 - (min(xs) + max(xs)) / 2) < 3  # the line is centred under the structure
    checks = json.loads((out/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed' and not checks['layout']['issues']


def test_with_several_values_each_has_its_own_centred_line(tmp_path):
    # Nature Communications PMC13594109, fig. 4: "10, X = F" / "64% yield" / ">99% e.e.".
    rows = [['a', '10', STEREO, '64% yield', '']]
    out = tmp_path/'many'
    spec = request(tmp_path, rows, fields=[{'column': 'value', 'label': ''}])
    assert execute(spec, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    label, value = _caption_line(root, project, 'a')
    assert float(value.get('p').split()[1]) > float(label.get('p').split()[1])
    assert label.get('Justification') == value.get('Justification') == 'Center'
    assert abs(float(label.get('p').split()[0]) - float(value.get('p').split()[0])) < 0.01


def test_a_line_too_long_for_the_cell_breaks_after_the_number(tmp_path):
    # A value that already fills the widest cell: the number cannot share its line.
    rows = [['a', '3a', 'CC', 'obtained as an inseparable 1:1 mixture of diastereomers after two '
                              'chromatographic purifications on silica gel, 64%', '']]
    out = tmp_path/'long'
    spec = request(tmp_path, rows, fields=[{'column': 'value', 'label': '', 'inline': True}])
    assert execute(spec, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    label, value = _caption_line(root, project, 'a')
    assert float(value.get('p').split()[1]) > float(label.get('p').split()[1])
    assert not [t for t in root.iter('t') if ''.join(t.itertext()) == ', ']


def test_a_caret_marks_a_footnote_letter_as_a_superscript():
    from cdxml_toolkit.publication.layout import caption_runs
    assert caption_runs('9%^b') == ('9%b', [{'text': '9%'}, {'text': 'b', 'superscript': True, 'italic': True}])
    assert caption_runs('12^b', bold=True) == ('12b', [{'text': '12', 'bold': True},
                                                       {'text': 'b', 'superscript': True, 'italic': True}])
    assert caption_runs('Alk = C7H15^a')[1][-1] == {'text': 'a', 'superscript': True, 'italic': True}
    assert caption_runs('92%') == ('92%', None) and caption_runs('3a', bold=True) == ('3a', [{'text': '3a', 'bold': True}])


def test_a_footnote_mark_reads_back_as_the_letter(tmp_path):
    rows = [['a', '12^b', STEREO, '40%^c', '']]
    project, page, root = _scope_page(tmp_path, rows)
    assert project['compounds']['a']['display']['label'] == '12b'
    checks = json.loads((tmp_path/'scope'/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed'


def test_a_structure_is_centred_in_a_cell_widened_by_its_caption(tmp_path):
    rows = [['a', '3a', 'CC', 'a long line that widens the cell well past the structure', '']]
    project, page, root = _scope_page(tmp_path, rows)
    box = project['compounds']['a']['bbox']
    xs = [float(n.get('p').split()[0]) for n in root.iter('n') if n.get('p')]
    assert abs((min(xs) + max(xs)) / 2 - (box[0] + box[2]) / 2) < 1


def test_a_formula_label_gets_its_subscripts_and_stays_bold():
    from cdxml_toolkit.publication.layout import caption_runs
    assert caption_runs('PhB(OH)2', bold=True) == ('PhB(OH)2', [{'text': 'PhB(OH)', 'bold': True},
                                                              {'text': '2', 'subscript': True, 'bold': True}])
    assert caption_runs('TS2', bold=True) == ('TS2', [{'text': 'TS2', 'bold': True}])


HEADER = {'substrates': [{'smiles': 'Ic1ccc([*:1])cc1N'}, {'smiles': 'ClCC(=O)N[*:2]'}],
          'products': [{'smiles': 'O=C(N[*:2])c1nc2cc([*:1])ccc2s1', 'label': '3a–t'}],
          'above': ['Cu(OAc)2 (20 mol %)', 'S8'], 'below': ['H2O, 100 °C, 24 h']}


def test_a_header_scheme_stands_above_the_grid(tmp_path):
    # BJOC 22-127 scheme 2 opens with the general reaction, R1 and R2 marked, over a rule.
    rows = [['a', '3a', STEREO, '64%', ''], ['b', '3b', OTHER, '51%', '']]
    out = tmp_path/'header'
    spec = request(tmp_path, rows, header=HEADER, fields=[{'column': 'value', 'label': '', 'inline': True}])
    result = execute(spec, out)
    assert result['ok'], (result, (out/'checks.json').read_text(encoding='utf-8') if (out/'checks.json').exists() else '')
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    page = root.find('page')
    grouped = {id(f) for g in page.iter('group') for f in g.iter('fragment')}
    header_frags = [f for f in page.iter('fragment') if id(f) not in grouped]
    assert len(header_frags) == 3
    header_arrow = [a for a in page.iter('arrow') if a.get('ArrowheadHead') == 'Full']
    assert len(header_arrow) == 1
    arrow_y = float(header_arrow[0].get('Head3D').split()[1])
    cell_tops = [project['compounds'][c]['bbox'][1] for c in 'ab']
    rules = [a for a in page.iter('arrow') if a.get('ArrowheadHead') == 'None']
    assert rules and all(arrow_y < float(r.get('Tail3D').split()[1]) < min(cell_tops) for r in rules)
    texts = {''.join(t.itertext()).strip() for t in page.iter('t')}
    assert 'S8' in texts and 'Cu(OAc)2 (20 mol %)' in texts
    checks = json.loads((out/'checks.json').read_text(encoding='utf-8'))
    assert checks['data_binding']['status'] == 'passed' and not checks['layout']['issues']


def test_a_scope_brief_header_is_checked():
    from cdxml_toolkit.mcp_runtime.figure_scope import _validate
    base = {'kind': 'scope', 'fields': ['yield'], 'species': [{'id': 'a', 'smiles': 'CCO', 'values': {'yield': '9%'}}]}
    assert _validate({**base, 'header': {'substrates': [{'smiles': 'CC'}], 'products': [{'smiles': 'CCO'}],
                                         'above': ['O2']}}) == []
    for header, message in (({'products': [{'smiles': 'C'}]}, 'substrates'),
                            ({'substrates': [{'smiles': 'C'}], 'products': [{}]}, 'smiles or a name'),
                            ({'substrates': [{'smiles': 'C'}], 'products': [{'smiles': 'C'}], 'arrow': 'x'}, 'unknown keys')):
        assert any(message in e for e in _validate({**base, 'header': header})), header


def test_a_shared_line_measures_subscripts_at_their_drawn_size():
    # "Alk = C7H15" was measured with full-size digits, so the ", " after it stood off by a gap.
    from cdxml_toolkit.publication.layout import _line_width, caption_runs, style_settings, text_width
    style = style_settings({})
    plain, runs = caption_runs('Alk = C7H15')
    full = text_width(plain, style)
    assert _line_width(['d'], {'d': plain}, style, {'d': runs}) < full
    assert _line_width(['d'], {'d': plain}, style, {}) == full


def test_the_plus_stands_midway_between_the_two_structures(tmp_path):
    path = _secondary_request(tmp_path, [['a', '14a', STEREO, '61%', '', ''], ['b', '14b', OTHER, '22%', '', 'a']])
    out = tmp_path/'plus'
    assert execute(path, out)['ok']
    project = json.loads((out/'project.json').read_text(encoding='utf-8'))
    root = ET.parse(out/project['pages'][0]['file']).getroot()
    groups = {c: root.find(f".//group[@id='{project['compounds'][c]['objects']['group']}']") for c in 'ab'}
    right = max(float(n.get('p').split()[0]) for n in groups['a'].find('fragment').iter('n'))
    left = min(float(n.get('p').split()[0]) for n in groups['b'].find('fragment').iter('n'))
    plus = next(t for t in root.iter('t') if ''.join(t.itertext()) == '+')
    assert abs(float(plus.get('p').split()[0]) - (right + left) / 2) < 0.5


@pytest.mark.parametrize('text, parsed', [
    ('Alk = C7H15', ('Alk', {'C': 7, 'H': 15})), ("R' = C18H37", ("R'", {'C': 18, 'H': 37})),
    ('Ar = C6H4OMe', None), ('not formed', None), ('Alk=CH3', ('Alk', {'C': 1, 'H': 3})),
])
def test_a_definition_names_a_label_and_a_formula(text, parsed):
    from cdxml_toolkit.publication.layout import parse_definition
    assert parse_definition(text) == parsed
