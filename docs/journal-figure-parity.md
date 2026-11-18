# Journal figure parity

Goal: the drawings this Skill produces should pass as publication figures in an organic
chemistry methodology paper, not merely be chemically correct.

Three sets are meant to be read together. The corpus and the coverage renders are
development evidence and live in the repository's `assets/`, outside the deployable Skill;
the briefs ship with the Skill in `skill/chemdraw/assets/preview/`:

| Set | What it is |
| --- | --- |
| `assets/journal-corpus/` | 51 real figures from Beilstein, Communications Chemistry and Nature Communications, each with DOI, licence and caption (`INDEX.csv`) |
| `assets/coverage/` | 40 figures this Skill currently produces, one per supported type, each as `.png`, `.cdxml` and a ChemDraw-re-saved `.native.cdxml` |
| `skill/chemdraw/assets/preview/` | the JSON briefs behind the preview figures, so any of them can be regenerated |

A gap is only worth acting on when it is named in a corpus figure and reproduced in a
coverage figure. Every item below has both.

## How to work

- Reproduce first: regenerate the coverage figure, look at it, then look at the corpus figure
  beside it. Numbers and pixel measurements beat adjectives.
- Add the failing test before the behaviour, as the repository requires. A layout rule that
  cannot be tested numerically still gets a test on the generated CDXML (attribute present,
  element count, ordering).
- After any change to the renderer, regenerate the affected coverage figures and the `.native`
  copies, then `cdx.py check assets/coverage/*.native.cdxml --allow-sizes 8,6 --palette any`.
- Keep `scripts/validate_distribution.py` and the Skill suite green.
- **Re-run `tests/runtime/test_draw_corpus.py` after any layout change.** Its fixtures record
  what the layout did — pages, sizes, positions — not chemical correctness, so a change that
  is an improvement still turns them red. When a fixture moves, record the new state with a
  note saying why, and put the intention in the gaps above rather than in the fixture.
- Do not add a figure type that has no corpus evidence: a convention invented from taste is
  what the corpus exists to prevent.

## Gaps, in the order they matter

Each item states the corpus witness, what to produce, and how to tell it is done.

### P1 — Coloured curved arrows, and arrow routing that avoids the structure

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` (green and red arrows on the same
cycle), `journal-corpus/bjoc-articles-22-130-scheme5.png`.
Today: `coverage/preview_mechanism.png` — one arrow colour, and an arrow that runs into the
structure because placement is taken from coordinates alone.

Produce: per-arrow colour, and placement that arcs around the atoms it does not start or end
on. The CDXML side is already known — an `<arrow>` with `AngularSize` plus `Tail3D`/`Head3D`
and a `Center3D` origin is what ChemDraw itself writes (see
`toolkit/tests/data/native-arc-arrows.cdxml`); direction is which end carries the head.

Done when: a two-pathway mechanism draws both arcs without either crossing a bond, verified
by regenerating and inspecting, and the generated CDXML carries two distinct arrow colours.

### P2 — Equilibrium arrows between species

Witness: same two figures; the cycle uses `⇌` between intermediates rather than a single step
arrow.
Today: a reaction step is one solid arrow.

Produce: a step whose arrows are an equilibrium pair, and a scheme-level way to say so in the
brief JSON.
Done when: a brief marks a step as equilibrium and the output renders a double half-headed
arrow, with the CDXML attributes asserted in a test.

### P3 — Newman projections inside a scheme

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` (four projections carrying the
stereochemical argument), `journal-corpus/bjoc-articles-22-130-scheme25.png`.
Today: Newman projections exist as their own figure type, not as a component placed inside a
scheme at a given position.

Produce: a scheme object that places a Newman projection at a point, sized like the structures
around it.
Done when: a mechanism brief places two projections among its structures and the scheme lays
out without overlap.

### P4 — Scope tables: highlight the varying group

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` — the substituent that changes is
drawn in red in all eighteen products; `journal-corpus/nature-communications-pmc13594109-fig4.png`
— the same idea with coloured substituents inside colour-coded bands.
Today: `coverage/substrate_scope.png` — every product in black, so the reader has to compare
structures to find what varies.

Produce: a per-binding highlight in the scope model that recolours the bond or label of the
varying group.
Done when: the example scope shows the varying group in a second colour, and the check profile
accepts that colour (`--palette` already takes a list; the default palette will need the second
colour added).

### P5 — Keep the failures visible

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` keeps `3j, 0%` in the grid.
Today: an entry without a yield is drawn without one, and an entry at `0%` looks like an
omission.

Produce: render `0%` and a "not formed" note when the brief says so, rather than dropping the
cell.
Done when: a scope brief with a zero entry shows it, and a test asserts the caption text.

### P6 — Variable-group rows

Witness: same figure, `3p, Alk = C7H15, 63%` and `3q, Alk = C18H37, 60%` under one structure.
Today: no way to express "this one structure stands for a family".

Produce: a caption form that carries a group definition beneath the structure.
Done when: a scope brief emits the two-line caption and the layout keeps the second line clear
of the next row.

### P7 — Secondary products beside the main one

Witness: `journal-corpus/nature-communications-pmc13594109-fig4.png` — `14b`, `17b` and `18b`
drawn in grey next to the main product with their own yields.
Today: one structure per cell.

Produce: an optional secondary structure in a scope cell after a "+", the structure grey at full
size and its own number and yield black (as `14b`, `52% yield` reads in the witness).
Done when: a cell with a secondary product renders both without overlap and the secondary is
visually subordinate.

### P8 — Section bands inside a scope

Witness: same Nature figure — `para substituents`, `ortho substituents`, `N-heterocycles`,
`C-nucleophiles` as labelled bands.
Today: a flat grid.

Produce: an optional grouping where each group's name stands on its rule (a short piece of rule
before it, the rest after it), only the locant in italic (*para* substituents, *N*-heterocycles).
Done when: a scope brief with two groups renders two bands, each name on its rule and clear of the
structures below.

### P9 — Generic (unspecified) rings

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` — in every species, X and the carbon
beside it carry a dashed ring with nothing drawn inside it: a ring of any size fused at that bond,
drawn as a dashed arc from one atom round to the other, the solid bond closing it.
Today: nothing. (This item was first read as a dashed overlay around a drawn fragment; that
annotation exists too, as `overlays`, but it is not what the witness shows.)

Produce: a species-level generic ring on a bonded pair of atoms, drawn as a dashed open arc on the
far side of the bond from the rest of the molecule, sized as a ring of `members` atoms.
Done when: a brief with a generic ring on its substrate and product draws one arc on each, clear
of the arrows, and the check reports nothing on it.

### P10 — Fit the page to the content

Witness: every corpus figure is cropped tight to its drawing
(`journal-corpus/bjoc-articles-22-127-scheme2.png`).
Today: `coverage/substrate_scope.png` puts five structures in the top third of a portrait
page and leaves the rest white — the page carries its default extent instead of the drawing's.

Produce: choose the page box from the laid-out content plus a margin, per figure.
Done when: the PNG's content bounds sit within a small margin of the page, measured as
non-white pixel bounds against image size.

### P11 — Never shrink the text to make it fit

Witness: corpus body text sits at 8-10 pt at print size
(`journal-corpus/nature-communications-pmc13594109-fig4.png`).
Today: `coverage/substrate_scope.native.cdxml` writes `LabelSize="8"` and `CaptionSize="8"`
where the ACS 1996 baseline is 10 — the layout shrank the type rather than growing the page.
`cdx.py check` reports both.

Produce: keep the type at the figure's stated size and make the page or the grid grow instead.
Done when: `check` reports no `setting` deviation for the figure and the text still clears its
neighbours.

### P12 — Say the yield, not "Yield:"

Witness: corpus scope entries read **3a**, 64% and **3b**, 64% — the number in bold and the
yield after it on one line, centred under the structure (BJOC 22-127 scheme 2). Where an entry
carries several values the number keeps the first line and each value takes a centred line of its
own (Nature Communications PMC13594109, fig. 4: **10**, X = F / 64% yield / >99% e.e.).
Today: `coverage/substrate_scope.png` reads `3a` / `Yield: 92%` on two lines per cell, which
doubles the vertical space a cell needs and repeats a word eighteen times.

Produce: the compact caption form; keep the long form available for the case where the value
needs a qualifier (`85% (2 steps)`).
Done when: the example scope shows **3a**, 92% on one centred line.

### P13 — Pack a scope grid densely enough to keep one page

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` puts eighteen compounds on one
page; `journal-corpus/nature-communications-pmc13594109-fig4.png` puts twenty-two on one.
Today: with the type at the correct 10 pt, a twelve-compound scope paginates — page 1 holds
ten cells and page 2 the remaining two (`draw_corpus/11-scope-twelve`), so the grid is roughly
one cell short of fitting.

Produce: cell geometry that packs twelve to twenty cells on one page at 10 pt — the space
between a structure and its values, and between cells, is where the slack is.
Done when: `11-scope-twelve` is a single page again with `LabelSize`/`CaptionSize` at 10, and
the eighteen-cell witness can be reproduced in one page.

## Acceptance

A figure is publication-ready when all of the following hold, and the evidence is quoted with
the claim rather than asserted:

1. `cdx.py check` reports no errors on the `.native.cdxml`, with `--allow-sizes` set for the
   sizes the figure deliberately uses and a palette that includes its highlight colour.
2. The gaps that the figure exercises (P1–P9) are closed for it, verified by regenerating and
   comparing against the corpus witness named above.
3. Text does not touch ink: `check` reports no `text_over_ink`.
4. A human has looked at the PNG at print size and accepted it. Machine checks cannot decide
   whether a scheme reads well; they only rule out the failures they know.

Record the per-figure status in [figure coverage](figure-coverage.md) next to the existing notes, so the next
session knows which figures have been brought to this bar and which have not.

## Status when this was written

- `journal-corpus/`: 71 figures over 16 types from three publishers, every row CC BY 4.0
  with DOI and caption.
- `coverage/`: 47 figures, each with its `.cdxml`, the ChemDraw-resaved `.native.cdxml` and a
  PNG. Six types were drawn for this task — `derivatization`, `optimization_table`,
  `process_route`, `sar_table`, `stereo_model`, `total_synthesis` — and their briefs live in
  `skill/chemdraw/assets/preview/` (14 briefs in total), so every figure can be regenerated.
- Every one of the 23 declared types now has a witness on at least one side; the corpus
  covers 16 of them and the coverage set the rest.
- Measured deviations are recorded in [figure coverage](figure-coverage.md) rather than hidden: `sar_table`
  and `substrate_scope` shrink `LabelSize`/`CaptionSize` to 8, and the conformer templates
  ignore a `bond_length` of 14.4 and draw at 12 whatever the brief asks.
- Reproducibility was checked: rendering `skill/chemdraw/assets/preview/total-synthesis.json` reproduces
  `assets/coverage/total_synthesis.png` with 0 differing pixels.

### P14 — One-line cell captions in a scope

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` writes `3a, 64%` — label and value on
one line. `journal-corpus/nature-communications-pmc13594109-fig4.png` writes them on two, but
with the value tight under the label.
Today: a cell is a label line plus one line per field, so three lines for label + yield +
definition (`assets/parity/p4-p7-p8-scope.png`).

Produce: label and first field on one line when they belong together; a definition keeps its own
line.
Done when: a two-field cell is as tall as the witness's and the eighteen-cell panel still fits
one page.

### P15 — Header scheme above a scope table

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` opens with the general reaction —
substrates with R¹/R² marked, arrow, conditions, product with R¹/R². The grid below then needs
no structure repeated.
Today: `assets/parity/p4-p7-p8-scope.png` has a text line reading `Standard conditions`; the
reader has to infer the transformation from the products.

Produce: an optional `header` block in a scope brief that draws one general scheme above the
grid, with the varying positions marked.
Done when: a scope brief with a header renders the scheme and its marks line up with the
highlighted positions below.

### P16 — Curved arrows around a cycle

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` and `...-130-scheme5.png` — a cycle is
a ring of arcs, each arrow following the circle.
Today: `assets/parity/p9-cycle-overlay.png` joins species with straight chords, so the "cycle"
reads as a diamond with stubs, not as a ring.

Produce: the cycle layout draws its steps as arcs on a circle, reusing the arc code P1 added for
mechanisms.
Done when: a four-species cycle renders four arcs on one circle, none of them straight.

### P17 — Chemistry formatting everywhere, not just in conditions

Witness: every figure in the corpus subscripts formula digits wherever they appear.
Today: `conditions` and `definition` are formatted (`H2O`, `K3PO4` come out right) but a
**species label** is not, and a digit right after a closing parenthesis escapes the rule:
`PhB(OH)2` renders flat in the five-step route of this session and in
`assets/parity/p9-cycle-overlay.png`, while `Na2CO3` beside it is correct.

Produce: one formatter for every text the layout emits, with parentheses handled and species
labels included.
Done when: `PhB(OH)2` subscripts and no rendered string contains a bare formula digit that the
corpus would subscript.

### P18 — Wrap a linear route, and draw each intermediate once

Witness: `journal-corpus/bjoc-articles-22-103-figure1.png` — a route moves left to right and
wraps; an intermediate drawn once serves the arrows on both sides of it.
Today: `layout: sequential` puts one step per row and redraws the shared structure: a five-step
route becomes five rows about 107 x 116 mm, with compounds 2, 3, 4 and 5 each drawn twice.

Produce: a layout that fills a row before starting the next, and reuses the structure at each
junction instead of redrawing it.
Done when: the five-step route occupies two rows instead of five and no structure appears twice.

### P19 — Branch from one structure to several products

Witness: `journal-corpus/bjoc-articles-22-103-figure1.png` panel (c): the reagent sits in the
middle and arrows leave it to both sides, to two different products.
Today: a step has substrates on the left and products on the right, always in a line; it cannot
express one substrate giving several products from a shared centre.

Produce: a divergent layout where a single structure has arrows leaving in more than one
direction.
Done when: a brief with one substrate and three products renders as a fan, not as three stacked
rows.

### P20 — Panels, boxed titles, and a colour legend

Witness: same figure — panels (a), (b), (c) each carry a boxed one-line title, substrates are
coloured by class, and the footer carries `R1 = H, MeO, Cl` / `R2 = CF3, CF2CF2H, ...`.
Today: a figure can hold bands with italic titles (P8) but not boxed panel headers, not a colour
scheme with a legend, and not a footer legend of variable groups.

Produce: panel groups with a boxed header, a colour assigned to a set of species, and a legend
block that pairs colours with the groups they stand for.
Done when: a three-panel figure reproduces the witness's structure — headers, coloured classes,
footer legend.

## Readback differences that are not changes

The guard reads every drawn molecule back and compares it with the figure's input. Two
differences are unavoidable, because a drawing necessarily fixes what a SMILES leaves open, and
both are tolerated as a deliberate verdict rather than a mismatch:

- **Drawn double-bond geometry.** `CC(=NO)c1ccccc1` comes back as `C/C(=N\O)c1ccccc1`: a C=N or
  C=C bond drawn on paper has an E/Z geometry, and the reader derives it from the coordinates.
  `chemistry_semantics.same_compound` returns `drawn_geometry` when the source left that bond's
  geometry open and the two agree once it is ignored, and the report says
  `double_bond_geometry: drawn`. A geometry the source *did* fix still has to come back the same,
  or the guard refuses the figure. Fixed in `492ae65`.
- **Kekulé form.** Aromatic rings come back with whichever alternating double bonds the drawing
  uses; the semantic key compares aromatised forms.

This needs writing down because it was mis-reported once: a session that had tested an oxime
before `492ae65` recorded the C=N case as a live blocker, and later sessions repeated it from
the log instead of re-running the brief. At `a9d4462` both `p19-fan.json` and
`p20-panels.json` render with `ok: true`, and the three panels read back `matched`. Before
reporting a guard refusal, re-run the brief against the commit under test.

One related formatting rule, decided here: a mechanism designator keeps its number on the
baseline. `Sn2` is written S<sub>N</sub>2, not a subscripted sulfur count, so `Sn1` and `Sn2`
are exempt from the formula-subscript rule; `S2` (disulfide) still subscripts.


### P22 — Four compounds to a row

Witness: `journal-corpus/bjoc-articles-22-127-scheme2.png` fits eighteen compounds in four
columns, with `3a, 64%` under each.
Today: `assets/parity/p4-p7-p8-scope.png` puts three to a row. The caption style now matches the
witness (one line, bold number, band name on the rule), so what is left is density: our cells are
wider because the type is fixed at 10 pt (P11) and our cell padding and gaps are ours, not theirs.

Produce: a scope that fits four cells in the manuscript column at 10 pt — by tightening the cell
padding, the gap between cells and the row pitch, not by shrinking type.
Done when: an eighteen-compound scope renders in four columns on one page, and a measured note
records the cell width, gap and row pitch that made it fit. If four cannot fit at 10 pt, the note
says what the arithmetic is rather than quietly drawing three.

### P23 — A cycle whose labels are next to their own arcs

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` — each step of a cycle carries its
label beside the arc it belongs to, so the reading order is unambiguous.
Today: the arcs themselves are right (`p9-cycle-overlay` draws three arcs on one circle of radius
96.0, and the four overlays are closed), but the labels are nowhere near them: measured distance
from `oxidative addition` to the nearest arc is **162.1 pt** and from `reductive elimination`
**190.1 pt**, against a circle radius of 96. The arrowheads all converge on the centre where
`Pd / Pd(0)` is written, and the two upper circles come close to that text.

Produce: place each step label against the middle of its own arc, outside the circle, and keep
the arrowheads off the centre label.
Done when: every step label is within 40 pt of its arc's middle, no label or structure touches
another, and the arrows read in one direction around the circle.
## What the layout cannot do yet

The seven items above share one root: the renderer lays a figure out as rows of well-formed
cells, and every corpus figure that looks like a paper is doing something else — wrapping a
route, branching from a centre, or laying panels side by side. Until the layout model gains
those three shapes (P18, P19, P20), a long route will look like a stack and a cycle like a
diamond, however correct its chemistry.

## Renders kept beside the parity briefs

`assets/parity/` holds one brief per parity item that needed a figure to argue about, plus the
render each produced on the machine with ChemDraw and RDKit. The PNGs are the state of the
drawing when the brief was written, not an approved target: two still show a defect recorded
here so a later fix has something to diff against.

Regenerate with `cdx.py draw assets/parity/<brief>.json <new>.cdxml --render`. Scope briefs
write a directory (`<new>.cdxml/png/page-001.png`); the others write a PNG beside the CDXML.
The output path must not already exist.
One brief per parity item that needed a figure to argue about, plus the render each produced
on the machine that has ChemDraw and RDKit. The PNGs are **the state of the drawing when the
brief was written**, not an approved target: two of them still show a defect that is recorded
here so a later fix has something to diff against.

Regenerate any of them with:

```bash
python skill/chemdraw/scripts/cdx.py draw assets/parity/<brief>.json <new>.cdxml --render
```

Scope briefs write a directory (`<new>.cdxml/png/page-001.png`); the others write a single
PNG beside the CDXML. The output path must not already exist — the tool refuses to overwrite.

## `p4-p7-p8-scope.json`

Two banded groups with italic titles, the varying group in `#FF6600`, and `3d′` as a grey
secondary product standing immediately right of `3d` with its own yield. Matches the Nature
Communications witness (`assets/journal-corpus/nature-communications-pmc13594109-fig4.png`).
No known defect.

## `p3-newman-mechanism.json`

Two butane conformers (anti, gauche) joined by a reversible arrow, and the hydrate mechanism
from the preview set below them.

**Fixed.** The collision came from the brief, not the code: `dihedral_offset` is the rear set angle relative to the front one, and the brief asked for 120 — a multiple of 120, so the rear spokes landed on the front ones. Angles within 10 degrees of a multiple of 120 are now rejected with an explanation, the brief carries two staggered conformations, and the render shows thirteen labels with no shared coordinate (measured after the fix).

## `p9-cycle-overlay.json`

A four-species catalytic cycle with the methoxyphenyl fragment circled in every intermediate.
The circles come out identical (74.3 x 74.3 pt, radius 37.16, spread 0.0) and use the spelling
`native-dashed-overlay.cdxml` suggested: a `<curve>` with `CurveType="2"` and
`LineType="Dashed"`, no `FillType`, no `Closed`, never an oval. That spelling is drawn open:
the render (and the sample's own PNG) leaves a gap between the last anchor and the first,
because `CurveType` 2 is the dashed flag alone and ChemDraw draws no segment back to the start.
Overlays are now written with `CurveType="3"` (closed and dashed), `Closed="yes"`, and the first
anchor repeated at the end.

**Both defects are fixed, and the sample itself was wrong.** The circles in
`toolkit/tests/data/native-dashed-overlay.cdxml` were open — first and last anchor 27.56 pt
apart, `CurveType="2"`, no `Closed` — so the reference this task was built on carried the very
bug being fixed. It has been replaced with a closed pair (`CurveType="3"`, `Closed="yes"`, the
first anchor repeated), verified in the render. After the fix: four circles closed, every arrow
at least 9.0 pt outside every circle, `Overlay:` warnings 0 and `overlay_collides` 0.

A separate finding stands: a generic-ring arc of `p9-generic-ring.json` runs **through** the
molecule — the closest sample point sits 0.00 pt from a bond, 13 points within 10 pt of one —
and the compound label sits inside the arc, which the checker does report (`Overlay:`).

## Where each gap stands

Code exists for every item through P21. What follows is what has been checked on a machine with
ChemDraw and RDKit, and how; an item is only "verified" here when a render was inspected and
the numbers quoted were measured. `assets/parity/` holds the briefs and renders for P3, P4,
P7, P8 and P9.

| Item | State | Evidence |
| --- | --- | --- |
| P1 arc arrows | verified | three native `<arrow AngularSize>`, per-move colours land in the colortable, `crossings: []`, clockwise arcs put the head on the tail end |
| P1 half arrow | verified | head size back to ChemDraw's own 600/525/150 (a half arrow is the full set at 0.6, not a different ratio); renders as a single thin barb |
| P2 equilibrium arrows | verified | the Newman figure draws its reversible pair as a double half-headed arrow |
| P3 Newman in a scheme | verified | after the brief fix: no two labels share a coordinate; anti has the rear CH3 straight down, gauche at the upper right. A front substituent written `*` with `attach: {smiles}` bonds the front carbon to the rest of the molecule, as the witness's V-1 to V-4 carry their ring system; `p3-newman-mechanism.json` adds one (V-4). Charges as in the witness's ylides: ±1 in the attached SMILES and `charges` on a substituent ([Rh]) are drawn as ChemDraw's circled ⊕/⊖ symbols that represent the atom's `Charge`; `p3-charged-newman.json` has a zwitterion with its flat form and V-1 ⇌ V-4 with [Rh]⊖, to render |
| P4 highlight the varying group | verified | `#FF6600` on the differing atoms and bonds of every cell in `p4-p7-p8-scope` |
| P5 keep 0 % visible | verified | the publication-figure test renders `Result: 0%` rather than dropping the cell |
| P6 variable-group rows | **changed, to re-render** | verified once as its own line; now, as in the witness, a `definition` joins the number's line: **3p**, Alk = C₇H₁₅, 63%. One structure for two compounds: a species' `drawn_with` lists its line under another's structure (3q under 3p's NH–Alk), to re-render. The structure now draws the defined group as its label (NH–Alk): the one substituent of the defined formula on a single, non-ring bond is folded up for the layout and contracted into a `NodeType="Fragment"` abbreviation labelled Alk, which keeps the group whole for the readback, so an octadecyl chain no longer widens its cell |
| P7 secondary product | **changed, to re-render** | as fig. 4 draws it: a "+" between the products (18 pt gap), the secondary structure grey and its number and yield black. The verified version had the caption grey too and no "+" |
| P8 section bands | **changed, to re-render** | as fig. 4 draws it: the name stands on the rule (a 12 pt piece before it, the rest after it) and only the locant is italic (*para* substituents, *N*-heterocycles); the verified version had the whole name italic under a full rule |
| P9 generic rings | **rebuilt, to render** | the witness's dashed circles are unspecified rings, not overlays: a species' `generic_rings` [{atoms: [i, j], members?}] draws a dashed open arc (CurveType 2) from atom i the long way round to j on the circle of a `members`-ring through that bond, on the far side from the atoms beside them, stopping 5 pt short of a labelled atom; the bond must exist. Cycle layouts keep arrows and labels off its box; other layouts warn. `overlays` (a dashed outline around a fragment) stays as a separate annotation, now written closed (CurveType 3 and the first anchor repeated): ChemDraw drew the earlier CurveType 2 outlines open, as the sample's own PNG shows — and the first render ran the open-chain substrate's arc through its alkyne, the number inside: a structure with generic rings is now depicted with a temporary ring of `members` carbons through the pair (for the 2D layout only), and a number the arc would enclose is moved under it; to re-render. As a last guard the side is chosen by measured clearance from every drawn bond (ends aside), shrinking the ring to five or four atoms when neither side keeps 3 pt, with a warning when it shrinks or still meets a bond |
| P10 page fits content | verified | the scope page crops to its content; `check` passes |
| P11 10 pt type | verified | `LabelSize`/`CaptionSize` are 10, not 8 |
| P12 yield, not "Yield:" | **changed, to re-render** | **3a**, 92% on one centred line, the number bold (BJOC 22-127 scheme 2); the verified version put `92%` under `3a`, which is fig. 4's form only when a cell has several values |
| P13 dense grid | verified | twelve compounds on one page, and an eighteen-compound benzamide panel on one page |

Compared against the witness figures themselves rather than their written summaries, P6, P7, P8,
P9 and P12 changed form and need one more render; P9 in particular was rebuilt, because the
witness's dashed circles are unspecified rings, not outlines around a drawn fragment. The
earlier P9 outlines were also drawn open (see the P9 row). P3 stays verified.

| Item | State | Evidence |
| --- | --- | --- |
| P14 one-line cell captions | **implemented, to render** | the number and a single value, or a definition and the value, share one centred line; with several values each takes a centred line under the number (fig. 4); a shared line that would not fit the cell breaks after the number |
| P15 header scheme above a scope | **implemented, to render** | a scope brief's `header` {substrates, products, above?, below?} draws the general reaction at the top of the first page (substrates joined by "+", the arrow with its conditions, the products with their labels) and a rule under it. Numbered R groups `[*:1]` are drawn R¹ (a `GenericNickname` "R1" with a raised number, for schemes too); the header's structures join the page's expected inventory and the readback. With `highlight_varying` the header's R marks (the R and its number) take the highlight colour, as the witness draws R¹/R² in the red of the varying groups; `p15-header-scope.json` mirrors BJOC 22-127 scheme 2 |
| P16 curved cycle arrows | **implemented, to render** | each step of a catalytic cycle is a native arc arrow (`AngularSize`, centred on the cycle) on the circle the intermediates sit on, clipped where the circle leaves each structure's box, the head on the next intermediate; entering and leaving species and the conditions are placed from the arc. An equilibrium step is a pair of concentric half-headed arcs, 3 pt apart over one span, forward outside and backward inside (`p16-cycle-equilibrium.json`, to render: check the barbs stand outside the pair); failed steps keep their straight chord arrows |
| P17 formatting everywhere | **verified** | conditions, definitions and scope captions were already formatted; species labels now are too, by a stricter rule than conditions (`text_formatting.is_formula_label`): only a label made entirely of element symbols, common groups, brackets and counts is a formula, so `PhB(OH)2`, `K2CO3`, `Et3N` get subscripts and `TS2`, `INT3`, `IM1`, `3a`, `L1` stay as written. Measured in the render: `PhB(OH)` face 96 with the `2` at face 32, `NaBH4`, `Ph3P=CH2` the same; `Sn1`/`Sn2` were the one false positive and are now exempt. A house style's bold or plain label keeps the subscripts (only the bold bit of each run changes) |
| P18 wrap the route, draw an intermediate once | **implemented, to render** | a brief's route now wraps with `wrap: "continue"` by default: rows hold as many steps as fit the width, and a row after the first starts with the arrow that leaves the structure the row above ended with, so nothing is drawn twice (six structures, two rows, at 504 pt in the test route that took seven structures before). `repeat` and `serpentine` stay available. The five-row case came from one step per row fitting the old repeat rule; re-render the five-step route to confirm it now takes two |
| P19 branch from one structure | **implemented, to render** | `layout: "fan"` puts the starting material in the middle and sends an arrow to each of two to four products, left, right, then down and up; each branch keeps its conditions (beside a vertical arrow) and its co-substrates over its arrow as reactants of the step. The stacked divergent layout stays the default; `p19-fan.json` takes acetophenone four ways. Two products from one compound are now a fan by default (`layout: "divergent"` keeps them stacked; a fan wider than `max_width` is drawn stacked and reported as `divergent`) |
| P20 panels, boxed titles, colour legend | **implemented, to render** | `kind: "panels"` stacks scheme or mechanism briefs, each drawn and checked by `draw_figure`, under boxed one-line titles of one width, with a legend of plain or coloured lines under the last; the parts are joined by `cdxml_compose.merge_cdxml`, which keeps each part's colours, fonts and ids. Classes are coloured with the species' own `color`/`atom_colors`. A panel can also be a `gallery` (structures with caption lines and no arrows, rows spread to the widest, as panel (a)) or `sections` (briefs each opened by a boxed title on a dashed rule across the figure, as panel (b)); `p20-panels.json` has all three kinds, to render |
| P21 steric clash marks | **implemented, to render** | BJOC 22-125 scheme 4 marks the crowded pairs of its unfavoured conformers (V-1: Ar and the ring; V-2, V-3: the ring and [Rh]) with a small red double zigzag. A mechanism species' `clashes: [{between: [atom, atom], color?}]` draws one in the middle of the clear gap between the two atoms, across the line joining them, as two open curves of straight segments; atoms are SMILES numbers, or `front:i`/`rear:i`/`attach:n` on a projection. A gap under 6.4 pt is reported. `p-clash-marks.json` has V-4 ⇌ V-3 with the ring–[Rh] mark, to render; marks on ordinary structures are covered by the toolkit test |
| P22 four compounds to a row | **implemented, to render: the arithmetic says four do not fit at 180 mm** | Measured on the witness: bond 30.6 px, caption type 21 px (1.46, the ACS 14.4/10), figure 1380 px = 66 type sizes, so at 10 pt it is about 657 pt (232 mm) wide and a cell about 170 pt; the journal printed it reduced to about 73 % (type about 7.3 pt). Our page is 180 mm = 510.2 pt, 478.2 pt between 16 pt margins; four columns with 12 pt gaps leave 110.5 pt a column. In `p4-p7-p8-scope` 4-methoxybiphenyl alone is about 124 pt wide, so four of them need 496 pt before any padding. `columns: 4` now draws an equal four-column grid when the cells fit and otherwise refuses with that arithmetic and the `width_mm` that holds them, never three quietly. `p22-four-columns.json` (180 mm, expected to refuse) and `p22-four-columns-wide.json` (270 mm, expected four columns) set the witness's eighteen compounds |
| P23 a cycle label beside its own arc | **implemented, to render** | names were placed after the species entering and leaving the cycle and were pushed out past them. Names now go first, against the middle of their own arc and outside the circle, searched outward and along the arc within 40 pt; the species follow, out along the radius through their arrow's end until they clear everything by 3 pt with an arrow of at least 16 pt that crosses no name or structure (the leaving biaryl's arrow had been a stub). With stand-in structures the Suzuki cycle's names sit 4.4–6.7 pt from their arcs' middles and the three arrows are 47–66 pt; `p9-cycle-overlay.json` re-renders the measured case |

P1-P13 are closed as far as this machine can check them, apart from the re-renders above. P14-P20 are not polish: the layout draws
rows of well-formed cells, while every corpus figure that reads like a paper wraps a route,
branches from one centre, or lays panels side by side.
