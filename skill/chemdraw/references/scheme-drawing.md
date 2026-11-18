# Drawing a new reaction scheme

Use this when the structures are known and the figure does not exist yet: a reaction, a multi-step route, conditions and yields. For a figure to copy from an image, see [image reconstruction](image-visual-review.md).

## Which tool

| Need | Use |
| --- | --- |
| A reaction, route, substrate-scope panel or mechanism from a request, drawn and checked in one call | `draw_figure` / `cdx.py draw` (next sections) |
| Reactants, reagents, products, conditions and yields laid out by hand-written compact text | `render_scheme` (below) |
| Exam or lecture drawing with abbreviations (C10H21), Chinese notes, boxed compound codes, hand-placed objects | [ACS 1996 drawings](acs1996-drawings.md), `cdx.py build` |
| Exact coordinates, hand-placed curved arrows, grids | [Figure manifest](publication-figures.md), `compose_chemical_figure` |
| Many compounds already in a CSV/XLSX table (scope, SAR), or an edit of such a figure | [Data-driven figures](data-driven-figures.md) |

## One call: draw_figure

Turn the request into a brief, then run `cdx.py draw brief.json out.cdxml --render` (or `draw_figure(brief=..., output_path=..., render=true)`). You decide what the user asked for; the tool grounds the structures, checks the chemistry and the geometry, and reports what it verified. It draws one reaction, a linear route, a convergent route (arms that merge into one step) or one starting material with several products. A substrate-scope panel and a mechanism each have their own brief, described after the report table below.

```json
{
  "species": [
    {"id": "SA", "smiles": "OC(=O)c1ccccc1O", "label": "1"},
    {"id": "Ac2O", "name": "Ac2O"},
    {"id": "ASP", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "2"}
  ],
  "steps": [
    {"substrates": ["SA", "Ac2O"], "products": ["ASP"],
     "above": ["H2SO4 (cat.)"], "below": ["rt, 2 h"], "yield": "85%"}
  ],
  "max_width": 504
}
```

- `species` lists every structure with an `id` and either `smiles` (used as given) or `name` (curated reagent table, then PubChem unless `CHEMDRAW_OFFLINE=1`). A name that cannot be resolved stops the run and nothing is written; give a SMILES instead. Do not guess a SMILES to get past this.
  - **Double-bond geometry:** a C=N or C=C whose geometry the SMILES leaves open (an oxime given as `CC(=NO)c1ccccc1`) is drawn one way, and the checks accept the drawing as that compound; the readback notes it. Write `/` and `\` in the SMILES (`C/C(=N\O)c1ccccc1`) when the geometry matters: a fixed one must come back the same.
  - **Coordination/dative bonds (`Order="dative"`):** transition-metal complexes, borane-amine adducts, and Lewis acid-base pairs support native dative bonds (`Order="dative"` in CDXML) via RDKit dative SMILES (`->`, for example `[NH3]->[BH3]`) or explicit `dative_bonds: [[donor_atom_idx, acceptor_atom_idx], ...]` on a species record.
  - **Radicals:** open-shell radical species in SMILES (for example `[CH2]c1ccccc1`, `[O]C(C)(C)C`, or carbenes) preserve their radical state (`Radical="Doublet"`, `"Singlet"`, or `"Triplet"` on the CDXML `<n>` node) and pair with visible radical dot symbols and 1-electron fishhook arrows (`electrons: 1`).
- `steps` take `substrates`, `products`, `above` / `below` (text lines), `yield` and `arrow` (`solid`, `dashed`, `failed`, `reversible` for an equilibrium, `resonance` for a double-headed arrow). A product of one step may be the substrate of the next; list the steps in any order and the tool works out the route:
  - **Convergent branch-merge arrows (`layout: "convergent"`):** two synthesis arms that meet in one step can either be drawn with true Y-junction / T-junction convergent branch-merge arrows (upper and lower branch sequences merging into a vertical connector that feeds the shared coupling arrow) or as stacked rows with the merge row repeating the joined intermediates.
  - **Divergent routes (`layout: "divergent"`):** one substrate with several products is drawn once on the left, each product on its own row with a horizontal arrow, the arrows joined by a vertical line so every product reads as made from that substrate. A branch's other substrates (the nucleophile of a Michael addition, say) and any reagent structures stand above its arrow, the co-substrates counted as reactants of that step; a step that does not use the shared substrate gets its own row underneath. With `"layout": "fan"` the substrate stands in the middle instead and an arrow leaves it for each product, left, right, then down and up (two to four products, each made from it in one step), each branch with its conditions and its co-substrates over its arrow, as a paper draws one compound taken two ways. One compound with two products is drawn as a fan without asking; give `"layout": "divergent"` to keep it stacked. A fan wider than `max_width` is drawn stacked, and the report's `layout` names what was drawn.
  - **Catalytic cycles (`layout: "catalytic-cycle"`):** when `"layout": "catalytic-cycle"` is set (or when the steps form a closed catalytic loop), intermediates are arranged around a circle of radius max(90, 32 × steps) pt (widened until no structures overlap), each step drawn as an arc of that circle from one intermediate to the next (ChemDraw's own arc arrow, so the cycle reads as a ring), with entering reactants / departing products on short radial arrows at each step and the conditions outside the ring. An equilibrium step (`"arrow": "reversible"`) is a pair of arcs on that circle, 3 pt apart, the forward one just outside and the backward one just inside, both over the same span and each with a half head whose barb stands outside the pair, as the two shafts of a straight equilibrium arrow. A failed step, or two intermediates too close for an arc, keeps a straight arrow along the chord.
  - Put in `above` / `below` only what the user gave: a step without conditions is drawn without them and the report says so.
- `max_width` is the row budget in points (default 504). A longer route wraps onto more rows, as many steps to a row as fit: by default (`"wrap": "continue"`) the next row starts at the margin with the arrow that leaves the structure the row above ended with, so every structure is drawn once; `"repeat"` draws that structure again at the start of the row, and `"serpentine"` turns the route back on itself with a vertical arrow between rows.
- `style` is a journal preset from ChemDraw's stock style sheets (default `acs1996`; also `wiley`, `adv_synth_catal`, `rsc`, `rsc_1col`, `nature`, `synthesis` and others; `cdx.py styles` lists them, `cdx.py draw --style PRESET` sets it). The figure is laid out at ACS geometry and rescaled to the preset's bond length, fonts and line widths; `max_width` still applies to the styled figure, and the report's `style` field names the preset and bond length. Scope and mechanism briefs take the same key; `cdx.py build` specs take a top-level `style` too.
- `style` can also name a **house style**: how a journal, publisher or research group actually draws schemes, measured from openly licensed published figures. A measurement preset fixes only bond length, line widths and fonts; a house style names one of those presets as its base and adds drawing conventions. `cdx.py styles` lists both kinds (house styles in a second section with their kind, base preset and number of evidence figures), `cdx.py styles NAME` prints a card with its evidence, and `chemdraw://style-guide/{name}` returns it. The conventions a scheme takes:
  - `compound_label_face` (`bold` or `plain`) and `compound_label_color` for the compound numbers;
  - `yield_placement`: `under_arrow` (the default), `under_product_label` (a separate line under the product's number), `after_product_label` (`3, 85%` on the number's line, joined by `yield_separator`, default `", "`) or `after_product_label_parenthesised` (`3 (85%)`); `yield_format` is a template such as `"{value}% yield"`;
  - `substrate_layout`: `inline_plus` (A + B) or `stacked_plus` (A above B with the "+" between, centred on the arrow);
  - `abbreviate_groups`: groups drawn as abbreviation labels (CO2Me, CO2Et, OMe, OCF3, OAc, NO2, CN, Boc, NHBoc, Cbz, Ph, Bn, OBn, Bz, OBz, Piv, Ts, OTs, Ms, OMs, Tf, OTf, OTBS, OTMS, TMS, TIPS, OTIPS, Bpin, N3, NMe2, CF3, tBu, nBu, iPr, Et, Me, Ac, PhthN); they stay chemically complete, so the saved file reads back as the same molecule. A group on the left of its parent is placed correctly but desktop ChemDraw spells it right-flowing (CO2Me rather than MeO2C) until Clean Up Structure is run on the file; see the [style gallery](../assets/house-styles/gallery/README.md) and the per-style evidence folders beside it;
  - `new_bond_color`: bonds formed in the step, found by atom mapping; when the mapping is incomplete or ambiguous the bonds stay black and the report says so;
  - `substrate_colors`: the i-th substrate of a step drawn in the i-th colour; `highlight_fill`: a pastel fill behind the largest ring of the final product; `reagent_text_color`; `arrow_style` `bold`; `show_plus_sign`; `conditions_split` `all_above` or `all_below`.

  The card's `measured` ratios set the rest of the geometry after the base preset: bond line width and wedge width as fractions of the bond length, the atom-label size from the bond length over the capital-letter height, and the shortest reaction arrow in bond lengths (labels kept between 6 and 14 pt, arrows between 3 and 16 bonds). The report lists them as `measured.*` under `applied`.

  A null convention leaves today's drawing. The report's `style` gives `house`, `preset`, `applied` and `skipped` (with the reason); `compound_label_italic_suffix`, `aromatic_rings: "circle"` and `conditions_split: "reagents_above_solvent_below"` (the brief's own above/below split is kept) are not applied yet.

  A card may also carry `figure_types`: one section per kind of figure the journal or group draws (`reaction_scheme`, `substrate_scope`, `mechanism`, `catalytic_cycle`, `energy_profile`, `total_synthesis`, `retrosynthesis`, `structure_panel`, `sar_table`, `optimization_table`, `concept_panel`, `graphical_abstract`, `stereo_model`, `derivatization`, `biosynthesis`, `biomolecule`, `coordination_complex`, `polymer`, `markush`, `process_route`, `photophysical_diagram`, `assignment_panel`, `setup_diagram`), each with its own conventions, signature features and the evidence crops it was read from (`cdx.py styles NAME` prints them). The tools apply `mechanism.curved_arrow_color` and `mechanism.curved_arrow_weight` (`thin`, `bond` or `bold`, relative to the preset's bond width) to the curved arrows of a `kind: "mechanism"` brief, and `catalytic_cycle.cycle_arrow_color` to the arrows of a catalytic-cycle layout. The other figure-type conventions (lone pairs, circled charges, intermediate numbering, energy units, disconnection marks, SAR table rules and so on) are listed under `skipped` as guidance: follow them when finishing the figure by hand or with `cdx.py build`. A `kind: "scope"` panel takes the base preset, `compound_label_face`, `compound_label_color`, `abbreviate_groups` and the measured geometry; its species accept `color`, `atom_colors`, `bond_colors` and `abbreviate` as in a scheme, and the `substrate_scope` conventions are listed under `skipped` as guidance.

- The brief itself can style what a figure needs beyond the house style, with or without one:
  - a species' `color` paints the whole structure; `atom_colors` `[{"atoms": [0, 1], "color": "#0072B2"}]` paints atoms (a bond between two listed atoms too) and `bond_colors` `[{"bonds": [[7, 8]], "color": "#D55E00"}]` paints bonds, with atoms numbered from 0 in SMILES order;
  - a species' `abbreviate` lists the groups to contract in that structure; it replaces the house style's `abbreviate_groups` for it, and `[]` keeps it fully drawn;
  - a species' `orientation` (`flip_horizontal`, `flip_vertical` or `rotate_180`) mirrors or turns the structure about its centre; wedges and hashes are swapped on a mirror image, so the drawn compound keeps its configuration;
  - `title` and `panel_label` (`"a"`, `"(b)"`) are drawn above the figure at its left edge, the label in bold;
  - a species' `label` is drawn as written when it is a compound number (`3a`, `TS2`, `INT3`, `L1`) and with subscripts when it is a formula (`PhB(OH)2`, `K2CO3`, `Et3N`): a label counts as a formula only when it is made entirely of element symbols, common groups (Ph, Et, OAc, dba, ...), brackets and counts, so compound numbers are never set as formulas. The same holds for a scope's compound labels;
  - a species' `caption` is a line under its compound number (`"gram scale"`, `"98:2 er"`);
  - a species' `box` frames the structure and its number: `true` for a thin black frame, or `{"color": "#1F4E9E", "fill": "#EAF2FB", "style": "dashed", "rounded": true}`; the fill is drawn behind the structure;
  - a species' `generic_rings` draws a ring of unspecified size fused at one of its bonds, as journals draw a substrate class: `[{"atoms": [3, 4], "members": 6}]` names the bonded pair (numbered from 0 in SMILES order) and draws a dashed arc from one atom the long way round to the other, on the circle of a `members`-atom ring (default 6) through that bond, on the far side from the atoms beside them; the bond itself stays solid and nothing is drawn inside. The arc stops just short of a labelled atom (X, N). The structure is laid out as if the ring were there (a ring of `members` atoms is added for the 2D layout and dropped again), so an open chain such as 4-pentynoic acid is drawn with room for the ring instead of zigzag across it, and a compound number the arc would enclose is moved under it. The arc goes on whichever side of the bond keeps further from every drawn bond (by measured distance, its two ends aside); if neither keeps 3 pt, a smaller ring is tried down to four atoms, and `warnings` says when it had to shrink or still meets a bond. Give the SMILES of the drawn part with the two atoms bonded; the readback compares that part. A catalytic cycle keeps its arrows and labels off the arc; in other layouts an arrow or caption that meets it is named in `warnings`. Mechanism species take the same key;
  - top-level `overlays` mark one fragment in several structures at once, for the ring that stays unchanged through the intermediates of a cycle: `[{"fragment": "c1ccccc1", "species": ["I", "II", "III"], "shape": "circle"}]` draws a dashed outline around the fragment (a SMILES, or a SMARTS) in each listed species, or in every drawn structure that has it when `species` is left out. The outline is built from that structure's drawn atoms, so it sits at the same place on the fragment wherever the structure was placed. `circle` (the default) is centred on the fragment, half a bond length outside its farthest atom; `outline` follows the fragment's outer atoms at that distance, for a chain or an unsymmetrical group. A listed species without the fragment stops the run; a fragment found more than once is marked at its first match and `warnings` says so. In a catalytic cycle the chord arrows, the arrows of entering and leaving species and the step labels are laid out clear of each outline, and a compound number goes under it when it reaches below the structure; in other layouts nothing is moved, and an arrow within 3 pt of an outline or a caption it runs through is named in `warnings` (`cdx.py check` reports the same as `overlay_collides`). It is written as a `<curve CurveType="3" Closed="yes" LineType="Dashed">` with no fill whose points end on the first anchor again, never as an oval, which loses its dash on import: ChemDraw draws a curve only from anchor to anchor, so a `CurveType="2"` curve (dashed, not closed) is left with a gap;
  - a condition line in `above` or `below` may be a run `{"text": "Pd(OAc)2", "color": "#DC3220", "bold": true, "italic": false}` or a list of runs that make one line (`[{"text": "cat.", "italic": true}, {"text": " Pd(OAc)2"}]`); `title` and `caption` take the same forms.

  These are applied after the house style, so a colour the brief asks for wins over the card's; the report's `style.applied` lists them and `warnings` names any atom, bond or run that could not be found.

Several schemes set one above the other, as journal figures do ("(a) ... (b) previous work (c) this work"), are one brief of `"kind": "panels"`: `{"kind": "panels", "panels": [{"title": "(a) previous work", "brief": {...a scheme brief...}}, {"title": "(b) this work", "brief": {...}}], "legend": [[{"text": "R1", "color": "#1F4E9E"}, {"text": " = H, MeO, Cl"}]]}`. Each panel is drawn and checked as its own scheme or mechanism brief (`style` and `max_width` given at the top apply to all), under a boxed one-line title as wide as the widest panel; the legend lines, plain text or coloured runs, go under the last panel. Colour a class of substrates with the species' own `color` or `atom_colors` and say what the colour stands for in the legend. A panel can hold, instead of one `brief`:

- `"gallery": [{"id": "v", "smiles": "...", "caption": ["vonoprazan", "antiulcer agent"], "atom_colors": {...}}, ...]`: structures with their names and no arrows, as a figure opens with the known compounds that share a motif. Cells are filled into rows up to `max_width` (480 pt by default) and each row of two or more is spread to the widest row; caption lines are centred under their structure, formula digits as subscripts. Each structure is read back and coloured like a species.
- `"sections": [{"title": "L. Grimaud 2009 [35]", "brief": {...}}, ...]`: several schemes in one panel, each opened by its title in a box on a dashed rule that runs across the figure, as reviews list earlier reports one under another.

The report lists each panel's layout (`gallery` and `sections` for these), readback and warnings, and `checks_passed` holds only when every panel passed.

Read the report before you answer:

| Field | Meaning |
| --- | --- |
| `species[].source` | `user` for a SMILES you were given, otherwise the table that resolved the name |
| `steps[].balanced`, `missing` | atom and charge balance of each step; a missing by-product (for example acetic acid) is named, not drawn |
| `layout` | `linear`, `sequential`, `convergent`, `divergent` or `catalytic-cycle`: how the steps were arranged |
| `lint` | overlaps, text under 6 pt, drawing wider than `max_width`; `lint.relayout` appears when a failed check triggered the one automatic re-layout, with what changed and the original problems |
| `readback` | the saved file's structures compared with the brief (`matched` or `mismatch`) |
| `native` | `rendered` with `png_path`, or `not_run` with the reason; `not_run` means unrendered draft |
| `notes`, `warnings` | steps without conditions, unused species, balance and layout problems that remain |

`ok` means a file was written; `checks_passed` means the layout and readback checks both passed. Open the PNG and look at it even then. The output is never written over an existing file: without an output path a fresh name is made in the temp folder.

### Substrate-scope panel

`"kind": "scope"` draws many compounds with their labels and the numbers the user gave, packed onto as many editable pages as needed:

```json
{
  "kind": "scope",
  "title": "Substrate scope",
  "conditions": "As supplied by the user",
  "footnotes": ["a Isolated yield."],
  "fields": [{"key": "yield", "label": "Yield"}, {"key": "ee", "label": "ee"}],
  "species": [
    {"id": "3a", "smiles": "Cc1ccc(cc1)C(=O)N", "label": "3a", "group": "para", "values": {"yield": "88%", "ee": "95%"}},
    {"id": "3b", "name": "benzamide", "label": "3b", "values": {"yield": "91%"}}
  ]
}
```

- `header` draws the general reaction over the grid, as journal scope tables open: `{"substrates": [{"smiles": "Ic1ccc([*:1])cc1N"}, {"smiles": "ClCC(=O)N[*:2]"}], "products": [{"smiles": "O=C(N[*:2])c1nc2cc([*:1])ccc2s1", "label": "3a–t"}], "above": ["Cu(OAc)2 (20 mol %)", "S8"], "below": ["H2O, 100 °C, 24 h"]}` puts the substrates (joined by "+"), the arrow with its conditions and the products in a row at the top of the first page, a rule under it. Write the varying positions as numbered R groups (`[*:1]`, drawn R¹); species may be given by `name` too. The header's structures are checked like the rest of the page.
- `fields` is optional; without it the keys used in `values` become the columns, in first-seen order. A `unit` is added after the value, so do not repeat it in the text.
- Captions are centred under their structure, the compound number in bold. With one field each compound reads as journal scope tables do, the value on the number's line (**3a**, 92%); say what the number is in a footnote or the conditions, or set `"show_label": true` on the field. With several fields the number has a line to itself and every field the next lines, each naming its field (`ee: 95%`), unless the field sets `"show_label": false`; a field with `"optional": true` is left out where its value is empty instead of drawn as a dash. A shared line that would not fit the cell breaks after the number.
- `definition` names what a structure stands for, and joins the number's line: `"Alk = C7H15"` reads **3p**, Alk = C7H15, 63%, only for the compounds that have one. When the structure has exactly one substituent of that formula (on a single, non-ring bond), it is drawn as the label, NH–Alk, as the journal does, so a long chain does not widen the cell; the label is a ChemDraw abbreviation that keeps the group whole, so the file reads back as the full compound, and `highlight_varying` colours it. Give the full SMILES; a definition whose formula matches no single group leaves the structure drawn in full.
- `^` marks a footnote letter in a label or a value: `"12^b"` and `"40%^c"` draw the letter as an italic superscript, and the report and the bound text read it as `12b`.
- `note` adds a line under one compound and nowhere else: `"not formed"` beside a `0%`.
- Formula digits in the conditions and in every line under a structure (definitions, values, notes) are subscripts, by the same rule as scheme conditions: `C7H15`, `Pd2(dba)3`; yields, ratios and labels such as `L1` stay as written. A failed entry keeps its value; give `0%` rather than leaving it out, which shows a dash.
- `"highlight_varying": true` draws what changes from compound to compound in a second colour, as journal scope tables do: the core every compound shares (the maximum common substructure, whole rings only) stays black and the rest of each structure, with the bond that joins it to the core, turns orange `#FF6600`, a colour the default check palette accepts. Give a hex colour instead of `true` for another. A compound with its own `color`, `atom_colors` or `bond_colors` keeps them. The report's `highlighted` gives how many atoms each compound had coloured; when the compounds share no core of three atoms or more nothing is coloured and a warning says so. With a `header`, its R marks (R¹, R²) take the same colour, so the general reaction points at what varies below.
- `"secondary_of": "14a"` on a species draws it beside that main product in the same cell, a "+" between them, the structure in grey (`#808080`, in the default check palette) and its own number and values in black: a minor isomer or by-product with its yield. It belongs to the main product's group, and a secondary product cannot have one of its own.
- `"drawn_with": "3p"` lists a species under 3p's structure instead of drawing its own, as journals do when one generic drawing (NH–Alk) stands for several compounds: its caption line, **3q**, Alk = C₁₈H₃₇, 60%, follows 3p's in the same cell. Give each its `definition`, so the shared drawing reads as the group's label. It names a species of the same group with a structure of its own (not a secondary product, not another `drawn_with`); the readback counts the structure once, for 3p.
- `group` opens a band: the group's name stands on a rule across the panel (a short piece before it, the rest after it), only its locants in italic (*para* substituents, *N*-heterocycles), then its compounds, each band starting a new row. Groups and compounds keep first-seen order. A value you were not given is left out: the panel shows a dash, never zero, and the report names the compound.
- `width_mm` and `height_mm` set the largest page (default 180 x 240 mm); more compounds continue on further pages, labels and values staying with their structure. Text stays at 10 pt, and each page is cropped to its content, the footnotes following the last row.
- `output` is a new folder: `page-001.cdxml`, ... and, with `--render`, `png/page-001.png`, ... The report lists `pages`, `native.png_paths`, and the same `species`, `lint`, `readback`, `checks_passed`, `notes` and `warnings` as above, plus `checks.data_binding` (each label and value found next to its structure). Pages are packed from measured sizes, so there is no automatic re-layout.
- To change a panel after the user edits it in ChemDraw, use `publication_figure update` ([data-driven figures](data-driven-figures.md)).

### Mechanism

`"kind": "mechanism"` draws the structures of a mechanism and one curved arrow for every electron move you list. You decide which electrons move; the tool places the structures, turns each move into an arrow from the drawn coordinates and checks the result. It never works out a mechanism: a step without `moves` is refused.

```json
{
  "kind": "mechanism",
  "species": [
    {"id": "ald", "smiles": "CC=O", "label": "1"},
    {"id": "hyd", "smiles": "[OH-]"},
    {"id": "alk", "smiles": "CC([O-])O", "label": "2"}
  ],
  "steps": [
    {"substrates": ["ald", "hyd"], "products": ["alk"], "above": ["H2O"],
     "moves": [
       {"from": {"species": "hyd", "atom": 0}, "to": {"species": "ald", "atom": 1}},
       {"from": {"species": "ald", "bond": [1, 2]}, "to": {"species": "ald", "atom": 2}}
     ]}
  ]
}
```

- Steps are drawn in the order given. A step with one substrate and one product continues the same row, the intermediate drawn once; any other step starts a new row and draws its substrates again (the report notes it). A row that would pass `max_width` wraps the same way.
- A species may be a Newman projection instead of a structure: `{"id": "1-anti", "newman": {"front": ["Me", "H", "H"], "rear": ["Br", "H", "Me"], "conformation": "staggered"}, "label": "1-anti"}` The front substituents sit at the top, lower right and lower left, the rear ones (staggered) at the upper right, bottom and upper left, each in the order given, so the conformer is chosen by the order of `rear`: with `front` `["Me", "H", "H"]`, `rear` `["H", "Me", "H"]` is anti and `["Me", "H", "H"]` gauche. `eclipsed` turns the rear set 14 degrees off the front one so both stay visible; `dihedral_offset` sets that angle directly (60 is staggered, it is not added to it), and an angle within 10 degrees of a multiple of 120 is refused because every rear label would print over a front one. It is laid out in the row like any structure, its spokes 1.4 bond lengths long so it sits in scale with the structures beside it, and is read back as drawn. When the front carbon is bonded to the rest of the molecule, as in projections of an intermediate that carry its ring system, write that front substituent `"*"` and give the rest as `"attach": {"smiles": "O=C([*])C1CCCCC1"}`, with `[*]` where it meets the front carbon: the group is drawn at ordinary bond length, running outward along that spoke, and read back with the projection as one molecule. Attached groups use the projection's labels (C, H, N, O, F, Si, P, S, Cl, Br, I, B, with OH, NH, NH2 or SH); one that needs another label, a radical or a charge beyond ±1 is refused, so draw that species as an ordinary structure. A charge of ±1 in the attached SMILES (`[*]C1=[O+]CCC1`) and one given on a substituent with `"charges": {"rear": [-1, 0, 0]}` (in the order of `rear`; the minus beside [Rh] of an ylide) are drawn as ChemDraw draws them: `Charge` on the atom and a circled ⊕ or ⊖ beside it, upper right where that keeps 3 pt from every bond, 1.5 pt from every label and clear of the disc, otherwise in the freest direction. An abbreviation (CO2Me) cannot carry a charge. A projection has no numbered atoms, so no move may start or end on it, and a step that shows one needs no moves (a rotation between conformers, `arrow: reversible`); its balance is not checked.
- `"clashes": [{"between": [1, 10]}]` on a species marks a steric clash, as a small red double zigzag (two short zigzags 2.4 pt apart, 9 pt long) across the line between the two named atoms, in the middle of the clear gap between them; `color` changes the red. Atoms are the species' SMILES atom numbers, or on a Newman projection `"front:0"`..`"front:2"`, `"rear:0"`..`"rear:2"` in the order of `front` and `rear` (a front `"*"` is the attached group's first atom) and `"attach:<n>"` for atom n of the attach SMILES. Where the gap is under 6.4 pt the mark touches a group and the report warns; the report counts the marks in `clash_marks`.
- A move goes `from` / `to` an `atom` (a lone pair, a charge, a radical) or a `bond` `[i, j]` (a pi or sigma bond), in species that are substrates of that step. `electrons` is 2 (full head, the default) or 1 (half head). `color` gives that arrow its own hex colour (`"#2E7D32"` for one pathway, `"#C62828"` for the other), over any house colour. The step `arrow` may be `solid`, `reversible` or `resonance`; `above`, `below` and `yield` work as in schemes.
- Atoms are numbered from 0 in the order of the SMILES you gave, hydrogens not counted. The report lists them per species (`species[].atoms`, for example `0 C, 1 C, 2 O`). Set `show_atom_numbers: true` to print them on a draft and check the moves by eye, then draw again without it.
- Each arrow is a ChemDraw arc arrow (an `<arrow>` with `AngularSize`, editable as an arrow in ChemDraw), bulging to whichever side and height keeps it clear of the atoms and bonds it does not start or end at, and of the arrows already drawn: arrows are routed in the order given, two heads stay 6 pt apart and two bodies 3 pt (arrows may share a tail). A one-electron arrow ends in a half head at 0.6 of the full head's size, as ChemDraw draws its own fishhooks.
- The report adds `moves` (what each arrow connects and its `checks`) and `arrows` (`crossings`, the atoms, bonds and other arrows an arrow runs into, and `relayout` when one automatic retry changed the gaps and loop height). A two-electron arrow that starts on an atom with no lone pair or negative charge is drawn and flagged in `moves[].checks` and `warnings`; atom and charge balance is checked per step as for schemes.
- Limits: arrows are drawn as given; the tool does not judge whether the mechanism is chemically right. A long arrow can still run close to a label, so open the PNG.

## render_scheme

Pass exactly one of `compact_text`, `yaml_text` or `json_path` (a `parse_reaction` result). Ground every structure first: give a SMILES from the user, a resolver or a parser. A structure given only as `name "..."` is looked up in the curated reagent table (common reagents, solvents and abbreviations); there is no PubChem or IUPAC lookup on this path, and an unknown name stops with an error asking for a SMILES.

Compact text, one scheme in a few lines:

```
"1"{OC(=O)c1ccccc1O} + Ac2O{CC(=O)OC(C)=O} --> "2"{CC(=O)Oc1ccccc1C(=O)O} (85%)
  above: "H2SO4 (cat.)"
  below: "rt, 2 h"
```

- `ID{SMILES}` defines a structure inline; a quoted ID such as `"1"` also becomes its compound label below the structure.
- `-->` is a solid arrow, `..>` a dashed one, `<=>` an equilibrium arrow, `<->` a resonance arrow, `X>` a failed one. `A --> B --> C` makes two steps and draws the shared intermediate once.
- `(85%)` after the last product is the yield, printed under the arrow.
- `above:` and `below:` take quoted text lines or bare structure IDs (drawn as structures). For several steps use `step 1:`, `step 2:` blocks, or `-->|a|` with an `@conditions` block of lettered conditions.
- `run: "2.15 g" -> "1.60 g, 72% yield"` or `run[1]: 85%` adds per-run scale and yield annotations.
- Directives start with `@`: `@layout linear|sequential|stacked|divergent|convergent|stacked-rows|catalytic-cycle`, `@wrap repeat|serpentine|none`, `@steps_per_row N`, `@max_row_width PT`, `@title "..."`.

### Compact Scheme DSL (`compact_text` — "Mermaid for reaction schemes")

`compact_text` lets you write complete multi-step schemes, divergent routes, and catalytic cycles in a concise declarative syntax:

```text
@layout catalytic-cycle
Cat: {[Pd]} "Pd(0)"
Int1: {c1ccccc1[Pd]Br} "Ar-Pd-Br"
Int2: {c1ccccc1[Pd]c2ccccc2} "Ar-Pd-Ar'"
Cat -->|a| Int1
  above: "ArBr, oxidative addition"
Int1 -->|b| Int2
  above: "Ar'B(OH)2, base"
Int2 -->|c| Cat
  above: "reductive elimination"
  run[1]: 85%
```

### JACS / Angew publication styling & 3D conformers

When a scheme or figure requires top-journal visual highlights, transition states, organometallic hapticity, orbitals, or 3D conformer projections (via `draw_figure`, `build_acs_figure`, or `acs_builder.py`):

- **`ring_fills` & `atom_halos`:** pastel JACS/Angew ring-interior polygon shading (`ring_fills`) and circular pharmacophore atom highlights (`atom_halos`) placed on a background layer behind bonds and labels.
- **`multi_attachments` ($\eta^5$-Cp ferrocene & organometallics):** native CDXML `NodeType="MultiAttachment"` centroid nodes with `Attachments="..."` atom ID lists and coordination bonds (`Display="Dash"` or `"WedgeBegin"`).
- **`ts_bracket` (`[...]‡`):** square transition-state brackets (`BracketType="SquarePair"`) with a superscript double-dagger `‡` (`\u2021`) annotation. In a `draw_figure` brief, set `ts_bracket: true` (or `{superscript, bracket_type, lip_size}`) on the species, or `ts_bracket: "<species id>"` on a step; a step's `ts_bracket: true` is refused because it does not say which species is the transition state.
- **`orbital` & `interaction`:** native `<orbital>` lobes (`s`, `p`, `sp3`, `dz2`, `lobe`, shaded/filled/unfilled) and non-covalent H-bond / $\pi$-stacking dashed interaction lines with distance labels (`"2.8 Å"`).
- **Conformer templates (`chair`, `newman`, `haworth`, `fischer`, `norbornane`):** deterministic 2D coordinate generators for cyclohexane chairs (`1_ax`/`1_eq` through `6_ax`/`6_eq`, optional `flip`), Newman projections (`front`, `back`, `dihedral_deg`), Haworth pyranose/furanose rings, Fischer projections, and bridged bicyclic `[2.2.1]` norbornane skeletons.

A route stays on one row unless told otherwise. Set `@max_row_width 504` (points, the width of a two-column journal page) and the renderer picks the largest number of steps per row that fits; an explicit `@steps_per_row` always wins. In YAML the keys are `max_row_width` and `steps_per_row`.

Condition lines are tidied on every path: `80 C` becomes `80 °C`, `RT` and `room temperature` become `rt`, `12 hrs` becomes `12 h`, formula digits become subscripts and prefixes such as `n-` and `tert-` are italic. Text that is already correct is left alone.

## Checks

1. Read the per-fragment `chemistry_validation` in the result; anything other than `preserved` means a structure changed while drawing.
2. For a balanced-reaction check run `check_reaction_balance` on the same species; it names a missing by-product but never draws one.
3. Render with `cdx.py render FILE.cdxml` and look at the PNG for clipped labels, arrows running into structures and plus signs touching atoms.

Not available in `render_scheme`: curved arrows (use `draw_figure` with `kind: mechanism`). Neither tool chooses conditions for you. State only the conditions the user gave; a missing one stays missing.
