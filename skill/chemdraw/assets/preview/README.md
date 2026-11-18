# Figure previews

One figure per figure type, drawn with this Skill's own tools on Windows with ChemDraw 26.
Each PNG is the native ChemDraw render of the CDXML beside it — nothing here was made in
another program. The JSON next to each PNG is the exact brief or spec that produced it, so
every image can be regenerated.

| Type | Image | Brief | Command |
| --- | --- | --- | --- |
| Reaction scheme | `reaction-scheme.png` | `scheme.json` | `cdx.py draw scheme.json out.cdxml --render` |
| Substrate scope | `substrate-scope.png` | `scope.json` | `cdx.py draw scope.json out/ --render` |
| Mechanism (curved arrows) | `mechanism.png` | `mechanism.json` | `cdx.py draw mechanism.json out.cdxml --render` |
| Catalytic cycle | `catalytic-cycle.png` | `cycle.json` | `cdx.py draw cycle.json out.cdxml --render` |
| ACS 1996 structure | `acs1996-structure.png` | `build.json` | `cdx.py build build.json out.cdxml --render` |
| Energy profile (curve) | — | `energy-profile.spec.json` | `cdx.py scientific energy_profile energy-profile.spec.json out.cdxml --render` |
| Energy profile (levels) | — | `levels.json` | `cdx.py scientific energy_profile levels.json out.cdxml --render` |
| Journal house style | `scheme-jacs-house-style.png` | `jacs2.json` | `cdx.py draw jacs2.json out.cdxml --render --style acs_jacs` |

## What each one exercises

**Reaction scheme** — a Henry reaction: two substrates and one product, reagent above the
arrow, conditions and yield below. Atom-balanced on purpose, so the step reports no missing
by-product; exercises the one-call layout and the `above` / `below` split.

**Substrate scope** — five compounds with a group heading, a yield column and a footnote,
packed onto an A4-ish page. Exercises data binding: every label and value must be found
next to its own structure.

**Mechanism** — acetaldehyde plus hydroxide, with one curved arrow per electron move.
Exercises arrow placement from drawn coordinates and the atom-crossing check.

**Catalytic cycle** — a Suzuki coupling as a three-step loop with the entering aryl halide
and boronic acid and the departing biaryl placed on connectors. Exercises the cycle layout
and, specifically, that species which enter or leave the loop are actually drawn.

**ACS 1996 structure** — a single compound with Chinese annotation, formula and molar mass.
Exercises the exam/lecture style: bond length 14.4, Arial 10 labels, SimSun for the Chinese
run, and text that stays legible at the printed size.

**Energy profile (curve)** — a torsion profile with three stationary points, a y-axis
carrying unit and symbol, and a gap annotation. Exercises the curve written by the toolkit's
own curve writer, a barrier arrow beside the first species with its label in the margin, and
the barrier report.

**Energy profile (levels)** — a two-pathway free-energy diagram in the journal `levels`
style: level bars and connectors, a solid major path and a dashed red competing path,
y-axis ticks and a ΔG axis. Exercises the second style and `pathways`. The red is outside the
house palette, so check the output with `cdx.py check --palette any`.

**Journal house style** — the same Michael addition drawn with `--style acs_jacs`, which
colours the substrates, colours the bond formed in the step and moves the yield onto the
product's label line (`3, 96%`). Exercises house-style conventions on top of a preset.

## Notes for whoever picks this up

- Every image here was inspected visually; the reports for each run named no layout issues
  except where noted below.
- The catalytic cycle preview is the **fixed** behaviour: before the fix the layout drew one
  intermediate per step and silently dropped every other species, so a Suzuki loop came out
  as three palladium complexes with no entering or departing species. `cycle.json` exercises
  exactly that regression.
- `substrate-scope.png` is the first page of a project folder (`out/page-001.cdxml`), not a
  single file — that is how `kind: "scope"` works.
- Both energy profiles are specs for `cdx.py scientific energy_profile`, which runs inside the
  toolkit runtime it selects. Their native PNGs are made with the commands in the table.
- A `layout: "divergent"` route draws each branch's co-substrate above its arrow and records
  it as a reactant of that step. Before this fix the co-substrates were dropped, so a
  three-way Michael addition reported its nucleophiles as missing: the same class of
  problem as the catalytic-cycle bug above.
