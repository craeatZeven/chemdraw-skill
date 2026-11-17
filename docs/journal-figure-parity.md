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
