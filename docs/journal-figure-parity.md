# Journal figure parity

Goal: the drawings this Skill produces should pass as publication figures in an organic
chemistry methodology paper, not merely be chemically correct.

Three sets sit in `skill/chemdraw/assets/` and are meant to be read together:

| Set | What it is |
| --- | --- |
| `journal-corpus/` | 51 real figures from Beilstein, Communications Chemistry and Nature Communications, each with DOI, licence and caption (`INDEX.csv`) |
| `coverage/` | 40 figures this Skill currently produces, one per supported type, each as `.png`, `.cdxml` and a ChemDraw-re-saved `.native.cdxml` |
| `preview/` | the JSON briefs behind the preview figures, so any of them can be regenerated |

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

Record the per-figure status in `coverage/README.md` next to the existing notes, so the next
session knows which figures have been brought to this bar and which have not.
