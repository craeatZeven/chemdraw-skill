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
| Energy profile | `energy-profile.png` | `energy-profile.spec.json` | `cdx.py py <skill>/scripts/energy_profile/energy_profile.py spec.json -o out.cdxml` |

## What each one exercises

**Reaction scheme** — two substrates and one product with reagents above the arrow,
conditions below and the yield under the arrow. Exercises the one-call layout, the
boundary between `above` and `below`, and Roman-numeral-free compound numbering.

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

**Energy profile** — a torsion profile with three stationary points, a y-axis carrying unit
and symbol, and a gap annotation. Exercises the curve coming from ChemDraw's own plot
output rather than hand-computed points.

## Notes for whoever picks this up

- Every image here was inspected visually; the reports for each run named no layout issues
  except where noted below.
- The catalytic cycle preview is the **fixed** behaviour: before the fix the layout drew one
  intermediate per step and silently dropped every other species, so a Suzuki loop came out
  as three palladium complexes with no entering or departing species. `cycle.json` exercises
  exactly that regression.
- `substrate-scope.png` is the first page of a project folder (`out/page-001.cdxml`), not a
  single file — that is how `kind: "scope"` works.
- The energy profile needs the toolkit runtime: run it through `cdx.py py`, not a bare Python.
