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

Produce: an optional secondary structure in a scope cell, drawn smaller and in grey, sharing the
cell's caption block.
Done when: a cell with a secondary product renders both without overlap and the secondary is
visually subordinate.

### P8 — Section bands inside a scope

Witness: same Nature figure — `para substituents`, `ortho substituents`, `N-heterocycles`,
`C-nucleophiles` as labelled bands.
Today: a flat grid.

Produce: an optional grouping where each group gets a rule and an italic label above its rows.
Done when: a scope brief with two groups renders two bands and the label sits inside the band,
not over a structure.

### P9 — Dashed mapping overlays across a cycle

Witness: `journal-corpus/bjoc-articles-22-125-scheme4.png` — the ring that stays unchanged is
circled with a dashed overlay in every intermediate.
Today: nothing.

Produce: an annotation that draws a dashed shape around a set of atoms in several structures at
once.
Done when: a cycle brief marks one fragment and every intermediate carries the overlay in the
same relative position.

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

Witness: corpus scope entries read `3a, 64%`, `3b, 64%` — number only, compound label above.
Today: `coverage/substrate_scope.png` reads `3a` / `Yield: 92%` on two lines per cell, which
doubles the vertical space a cell needs and repeats a word eighteen times.

Produce: the compact caption form; keep the long form available for the case where the value
needs a qualifier (`85% (2 steps)`).
Done when: the example scope shows `3a` with `92%` beneath it.

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

**Known defect:** in the gauche conformation the rear substituents are not rotated relative to
the front ones, so three label pairs land on identical coordinates — `(297.21, 160.23)` H/CH3,
`(279.75, 129.99)` CH3/H, and `(262.29, 160.23)` H/H drawn twice — and the PNG shows a
collided glyph. The anti conformation is clean.

## `p9-cycle-overlay.json`

A four-species catalytic cycle with the methoxyphenyl fragment circled in every intermediate.
The circles come out identical (74.3 x 74.3 pt, radius 37.16, spread 0.0) and use the spelling
`native-dashed-overlay.cdxml` established: a `<curve>` with `CurveType="2"` and
`LineType="Dashed"`, no `FillType`, no `Closed`, never an oval.

**Two known defects:**

- the arrow `(296,171) -> (349,263)` passes 24.0 pt from a circle centre, inside its radius
  of 37.16, so it cuts through the dashed overlay;
- the step label `reductive elimination` overlaps the top-left circle.

Neither is caught by a check yet; `arrows_collide` covers arrow against arrow and bond, not
arrow or label against an overlay.

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
| P3 Newman in a scheme | **fixed, to re-render** | the render was right for the brief, the brief was wrong: `dihedral_offset` is the rear set's angle from the front one (60 = staggered), so 120 put every rear spoke and label on a front one. Such offsets (within 10 degrees of a multiple of 120) are now refused with that explanation, the docs say the conformer is chosen by the order of `rear`, and `p3-newman-mechanism.json` draws anti (rear `H, CH3, H`) and gauche (rear `CH3, H, H`), both staggered; six distinct label places are asserted for every accepted offset |
| P4 highlight the varying group | verified | `#FF6600` on the differing atoms and bonds of every cell in `p4-p7-p8-scope` |
| P5 keep 0 % visible | verified | the publication-figure test renders `Result: 0%` rather than dropping the cell |
| P6 variable-group rows | **implemented, to render** | a species' `definition` (`"Alk = C7H15"`) is an unnamed line straight under the number, drawn only where given; a field may set `show_label: false` (unnamed beside other fields) and `optional: true` (no dash where empty); formula digits in every caption line but the label, and in the conditions, are subscripts by the scheme rule (`text_formatting.subscript_matches`). `p6-family-scope.json` has 3p/3q with definitions among plain cells. Not done: one structure standing for two compounds (the witness draws "Alk" once with 3p and 3q under it) |
| P7 secondary product | verified | `3d′` sits 0-6 pt right of `3d`, grey, with its own `9%b`; the all-grey assertion needed the colour-table dedupe fixed first |
| P8 section bands | verified | two bands, italic titles, full-width rules |
| P9 dashed overlay | **fixed, to re-render** | the cycle layout now treats each marked structure's outline as part of it: chord arrows are clipped at the outline's box, entering/leaving arrows start outside it, step labels are placed clear of it and a number that would sit inside goes under it. Anything left is a warning (`overlay_collides`: arrows within 3 pt, as in `arrows_collide`, and captions the dashed line runs through); other layouts only warn. The routing test fails with routing switched off |
| P10 page fits content | verified | the scope page crops to its content; `check` passes |
| P11 10 pt type | verified | `LabelSize`/`CaptionSize` are 10, not 8 |
| P12 yield, not "Yield:" | verified | cells read `3a` over `92%` |
| P13 dense grid | verified | twelve compounds on one page, and an eighteen-compound benzamide panel on one page |

P3, P6 and P9 are closed too, verified on this machine after the fixes: the P6 definition line
renders without a field name and with subscripts, the P3 labels no longer share a coordinate, and
every arrow in P9 clears every circle by at least 9.0 pt (`Overlay:` warnings 0,
`overlay_collides` 0).

| Item | State | Evidence |
| --- | --- | --- |
| P14 one-line cell captions | not implemented | a cell is a label line plus one line per field |
| P15 header scheme above a scope | not implemented | the scope carries a text line reading `Standard conditions` |
| P16 curved cycle arrows | not implemented | `assets/parity/p9-cycle-overlay.png` joins species with straight chords |
| P17 formatting everywhere | partly | `H2O` and `K3PO4` subscript in conditions and definitions; a species label does not, and `PhB(OH)2` renders flat because the digit follows a closing parenthesis |
| P18 wrap the route, draw an intermediate once | not implemented | five steps become five rows, 107 x 116 mm, with compounds 2-5 each drawn twice |
| P19 branch from one structure | not implemented | a step is always substrates-left, products-right |
| P20 panels, boxed titles, colour legend | not implemented | bands exist (P8); panels and legends do not |

P1-P13 are closed as far as this machine can check them. P14-P20 are not polish: the layout draws
rows of well-formed cells, while every corpus figure that reads like a paper wraps a route,
branches from one centre, or lays panels side by side.
