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

Code exists for all thirteen items. What follows is what has been checked on a machine with
ChemDraw and RDKit, and how; an item is only "verified" here when a render was inspected and
the numbers quoted were measured. `assets/parity/` holds the briefs and renders for P3, P4,
P7, P8 and P9.

| Item | State | Evidence |
| --- | --- | --- |
| P1 arc arrows | verified | three native `<arrow AngularSize>`, per-move colours land in the colortable, `crossings: []`, clockwise arcs put the head on the tail end |
| P1 half arrow | verified | head size back to ChemDraw's own 600/525/150 (a half arrow is the full set at 0.6, not a different ratio); renders as a single thin barb |
| P2 equilibrium arrows | verified | the Newman figure draws its reversible pair as a double half-headed arrow |
| P3 Newman in a scheme | verified | after the brief fix: no two labels share a coordinate; anti has the rear CH3 straight down, gauche at the upper right. The witness also attaches its projections to the rest of the molecule (the front carbon bonded to an acyl ring system); standalone projections are what P3 covers |
| P4 highlight the varying group | verified | `#FF6600` on the differing atoms and bonds of every cell in `p4-p7-p8-scope` |
| P5 keep 0 % visible | verified | the publication-figure test renders `Result: 0%` rather than dropping the cell |
| P6 variable-group rows | **changed, to re-render** | verified once as its own line; now, as in the witness, a `definition` joins the number's line: **3p**, Alk = C₇H₁₅, 63%. Still not done: one structure standing for two compounds (the witness draws "Alk" once with 3p and 3q under it) |
| P7 secondary product | **changed, to re-render** | as fig. 4 draws it: a "+" between the products (18 pt gap), the secondary structure grey and its number and yield black. The verified version had the caption grey too and no "+" |
| P8 section bands | **changed, to re-render** | as fig. 4 draws it: the name stands on the rule (a 12 pt piece before it, the rest after it) and only the locant is italic (*para* substituents, *N*-heterocycles); the verified version had the whole name italic under a full rule |
| P9 generic rings | **rebuilt, to render** | the witness's dashed circles are unspecified rings, not overlays: a species' `generic_rings` [{atoms: [i, j], members?}] draws a dashed open arc (CurveType 2) from atom i the long way round to j on the circle of a `members`-ring through that bond, on the far side from the atoms beside them, stopping 5 pt short of a labelled atom; the bond must exist. Cycle layouts keep arrows and labels off its box; other layouts warn. `overlays` (a dashed outline around a fragment) stays as a separate annotation, now written closed (CurveType 3 and the first anchor repeated): ChemDraw drew the earlier CurveType 2 outlines open, as the sample's own PNG shows — and the first render ran the open-chain substrate's arc through its alkyne, the number inside: a structure with generic rings is now depicted with a temporary ring of `members` carbons through the pair (for the 2D layout only), and a number the arc would enclose is moved under it; to re-render |
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
| P15 header scheme above a scope | **implemented, to render** | a scope brief's `header` {substrates, products, above?, below?} draws the general reaction at the top of the first page (substrates joined by "+", the arrow with its conditions, the products with their labels) and a rule under it. Numbered R groups `[*:1]` are drawn R¹ (a `GenericNickname` "R1" with a raised number, for schemes too); the header's structures join the page's expected inventory and the readback. Not done: lining the R marks up with the highlighted positions below; `p15-header-scope.json` mirrors BJOC 22-127 scheme 2 |
| P16 curved cycle arrows | **implemented, to render** | each step of a catalytic cycle is a native arc arrow (`AngularSize`, centred on the cycle) on the circle the intermediates sit on, clipped where the circle leaves each structure's box, the head on the next intermediate; entering and leaving species and the conditions are placed from the arc. Equilibrium and failed steps keep their straight chord arrows (an arc has no double half-headed form here) |
| P17 formatting everywhere | **implemented, to render** | conditions, definitions and scope captions were already formatted; species labels now are too, by a stricter rule than conditions (`text_formatting.is_formula_label`): only a label made entirely of element symbols, common groups, brackets and counts is a formula, so `PhB(OH)2`, `K2CO3`, `Et3N` get subscripts and `TS2`, `INT3`, `IM1`, `3a`, `L1` stay as written. A house style's bold or plain label keeps the subscripts (only the bold bit of each run changes) |
| P18 wrap the route, draw an intermediate once | **implemented, to render** | a brief's route now wraps with `wrap: "continue"` by default: rows hold as many steps as fit the width, and a row after the first starts with the arrow that leaves the structure the row above ended with, so nothing is drawn twice (six structures, two rows, at 504 pt in the test route that took seven structures before). `repeat` and `serpentine` stay available. The five-row case came from one step per row fitting the old repeat rule; re-render the five-step route to confirm it now takes two |
| P19 branch from one structure | **implemented, to render** | `layout: "fan"` puts the starting material in the middle and sends an arrow to each of two to four products, left, right, then down and up; each branch keeps its conditions (beside a vertical arrow) and its co-substrates over its arrow as reactants of the step. The stacked divergent layout stays the default; `p19-fan.json` takes acetophenone four ways |
| P20 panels, boxed titles, colour legend | not implemented | bands exist (P8); panels and legends do not |

P1-P13 are closed as far as this machine can check them, apart from the re-renders above. P14-P20 are not polish: the layout draws
rows of well-formed cells, while every corpus figure that reads like a paper wraps a route,
branches from one centre, or lays panels side by side.
