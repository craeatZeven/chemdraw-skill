# Figure coverage

Every figure type this Skill can draw, rendered once on Windows with ChemDraw 26
from the example spec that ships with it. The files are in
[`assets/coverage/`](../assets/coverage/), outside the deployable Skill: they are the evidence behind the
claim "we cover these figure types", not a gallery of hand-picked favourites: the
same script rendered all of them.

Each figure has three files:

- `<name>.png` — native ChemDraw render, 300 dpi
- `<name>.cdxml` — the file the generator wrote
- `<name>.native.cdxml` — the same drawing after **ChemDraw opened it and saved it again**
  (xml → cdx → xml through ChemDraw's own SaveAs). If a generator ever writes something
  ChemDraw cannot round-trip, this file is where it shows.

## Inventory

| Figure | Type | PNG |
| --- | --- | --- |
| `apparatus` | Laboratory apparatus assembled from the Clipware stationery | yes |
| `binary_phase_diagram` | Binary phase diagram, eutectic | yes |
| `binary_phase_diagram_isomorphous` | Binary phase diagram, isomorphous | yes |
| `correlation_diagram` | Correlation diagram between two states | yes |
| `crystal_field` | Crystal-field splitting diagram | yes |
| `derivatization` | Divergent route: one substrate, three nucleophiles | yes |
| `electrochemical_cell` | Electrochemical cell schematic | yes |
| `electron_configuration` | Electron-configuration orbital box diagram | yes |
| `ellingham` | Ellingham diagram (oxide formation energy vs T) | yes |
| `energy_cycle` | Energy cycle (closed loop of states) | yes |
| `energy_profile` | Free-energy profile (smooth curve, TS daggers) | yes |
| `energy_profile_catalysis_cn` | Free-energy profile, catalysed vs uncatalysed | yes |
| `energy_profile_competing` | Free-energy profile, competing pathways | yes |
| `energy_profile_levels` | Free-energy profile, levels style (journal) | yes |
| `energy_profile_two_step_cn` | Free-energy profile, two steps | yes |
| `flowchart` | Flowchart | yes |
| `frost` | Redox Frost diagram (N E° against oxidation state) | yes |
| `huckel_mo` | Huckel MO diagram (energy levels from the adjacency matrix) | yes |
| `huckel_mo_pyridine` | Huckel MO diagram, pyridine (heteroatom parameter) | yes |
| `jablonski` | Jablonski diagram (S0/S1/S2/T1, A/F/P, IC/ISC) | yes |
| `latimer` | Latimer diagram (reduction potentials in series) | yes |
| `lewis` | Lewis structure (CO2, lone pairs and bonds) | yes |
| `lewis_carbonate` | Lewis structure, carbonate with formal charges | yes |
| `mechanism` | Mechanism with curved arrows | yes |
| `mo_diagram` | Molecular-orbital diagram (O2, 2s/2p mixing, spin filling) | yes |
| `nmr` | NMR spectrum panel (ethanol, synthetic JCAMP-DX) | yes |
| `optimization_table` | Condition-screening table composed from text and rules | yes |
| `periodic_table` | Periodic table panel | yes |
| `phase_diagram` | Phase diagram (one component) | yes |
| `potential_surface` | Potential-energy surface (contour) | yes |
| `pourbaix` | Pourbaix diagram (Fe, pH vs E) | yes |
| `preview_build` | ACS 1996 structure with Chinese annotation | yes |
| `preview_cycle` | Catalytic cycle (preview brief) | yes |
| `preview_energy_profile` | Energy profile (preview brief) | yes |
| `preview_jacs_house_style` | Reaction scheme with the acs_jacs house style | yes |
| `preview_levels` | Energy profile, levels style (preview brief) | yes |
| `preview_mechanism` | Mechanism with curved arrows (preview brief) | yes |
| `preview_scheme` | Reaction scheme (Henry reaction) | yes |
| `process_route` | Three-step scalable route with scale annotation | yes |
| `sar_table` | SAR table: three property fields per structure | yes |
| `stereo_model` | Chair, Newman, Haworth and Fischer templates in one figure | yes |
| `substrate_scope` | Substrate scope grid | yes |
| `tanabe_sugano` | Tanabe-Sugano diagram (d3, C/B = 4.5) | yes |
| `tlc` | TLC plate with lanes and Rf spots | yes |
| `total_synthesis` | Four-step linear total synthesis | yes |
| `unit_cell` | Crystal unit cell | yes |
| `walsh_diagram` | Walsh diagram (geometry vs orbital energy) | yes |

## What this run did not cover

- **nmr** — `toolkit/docs/examples/scientific/nmr.json` is the right kind of file (a spec whose
  `input` names a spectrum), but the spectrum it names, `source-nmr/1d_data.ft` from the
  [nmrglue integration example](https://nmrglue.readthedocs.io/en/latest/examples/integrate_1d.html),
  is not bundled. Download it beside the spec, or point `input` at another spectrum (CSV,
  NMRPipe, JCAMP-DX, Bruker, JEOL, Varian). A missing spectrum is now reported as
  `spectrum not found: <path>`.
- **apparatus** — the Clipware libraries (`clipware-1.cdxml`, `clipware-2.cdxml`) are licensed
  with ChemDraw and never shipped. Convert a copy of each CTP file from the ChemDraw installation
  (`cdx.py call convert_cdx_cdxml`) into `native-templates/`; see
  [laboratory drawings](../skill/chemdraw/references/native-laboratory-drawings.md). A missing template is
  now reported by name with that instruction.
- **substrate scope** — `kind: "scope"` writes a project folder rather than one file, so it
  stays in `skill/chemdraw/assets/preview/substrate-scope.png` instead of being duplicated here.

## Renders to redo

`cdx.py check` now reports captions that a line, arrow, curve or other caption runs through
(`text_over_ink`). On these renders it finds six figures with such captions, all since fixed in the
generators:

| Figure | What runs through a caption |
| --- | --- |
| `binary_phase_diagram` | tie lines through two field names; the eutectic temperature under the leader of β |
| `binary_phase_diagram_isomorphous` | the L + α leader through the point label B |
| `mo_diagram` | spin arrows into the labels under the centre-column levels |
| `pourbaix` | the H₂O/H₂ line through its own label; a water line through Fe₂O₃ |
| `potential_surface` | contour lines through the stationary-point labels |
| `tanabe_sugano` | stacked term labels at the right edge touching each other |

Caption widths now follow Arial's metrics, so every diagram's spacing moves slightly: re-render
the whole set rather than only these six.

`cdx.py scientific` no longer stops when the console code page cannot encode a character of
its report (a true minus sign on a GBK console): such characters are printed as `\uXXXX`
escapes, and the PNG is rendered.

## Re-rendered after the caption-width change (2026-11-16)

All 38 figures were regenerated with the caption metrics in place, and every one was
put through ChemDraw again (`xml -> cdx -> xml`). The six in "Renders to redo" above,
and the rest, were rebuilt rather than patched.

`cdx.py check ... --allow-sizes 8,6 --palette any` over the 38 native files:
**no `text_over_ink` warnings left.** What it does still report, and why:

- `mechanism` and `preview_jacs_house_style` fail the default `acs1996` settings
  because they are not drawn at ACS 1996 geometry (28 pt bond length for the
  mechanism, the acs_jacs line width and 9.1 pt labels for the house style). Pass
  `--profile` per figure or read those as expected.
- `tlc` shows float rounding after the native re-save: `HashSpacing` 2.49 vs 2.5,
  `MarginWidth` 1.59 vs 1.6.

### What the re-save caught

`preview_mechanism`: the condition label `H2O` is a plain text object in the generated
file. After ChemDraw saves the drawing back, the label is still there **and** a
free-standing water molecule has appeared at the same place (fragment `1052`, one
oxygen with two hydrogens). The render is unchanged because the molecule sits behind
the text, but the saved file now carries a chemical entity nobody drew.

A condition label that reads as a formula can therefore become a structure on a native
round-trip. Compare fragment counts in `.native.cdxml`, not just the picture.

### Both gaps are now filled

- **apparatus** — rendered from the Clipware stationery that ships with the ChemDraw
  installation (`ChemDraw Items/Clipware, part 1.ctp` and `part 2.ctp`), converted to
  CDXML and read by page. Those templates are ChemDraw-licensed, so they stay out of
  the repository: `native-templates/.gitignore` excludes them and only the render is
  committed.
- **nmr** — the example used to point at an NMRPipe file that was never in the tree.
  It now ships a synthetic JCAMP-DX ¹H spectrum of ethanol (`source-nmr/ethanol-1h.dx`,
  generated by a script, no third-party data), and integration regions that match it.

### What `check` reports across the set, and how to read it

`cdx.py check *_native.cdxml --allow-sizes 8,6 --palette any` over the 40 files:
**34 come back with no errors and no warnings**. The rest report two kinds of thing, and
neither is a defect in the figure.

**`setting`** — the checker compares the six root attributes against the ACS 1996 baseline
(bond length 14.4, labels 10 pt, line width 0.6, `HashSpacing` 2.5, `MarginWidth` 1.6,
`ChainAngle` 120). Three different causes show up here:

| Report | Files | Cause |
| --- | --- | --- |
| `HashSpacing` 2.49 vs 2.5; `MarginWidth` 1.59 vs 1.6 | 4 | ChemDraw rounds on save. 0.01 apart, no visible difference |
| `BondLength` 28 vs 14.4; `LabelSize` 11 vs 10 | `mechanism` | The `scientific mechanism` workflow draws at 28 pt; the same brief through `cdx.py draw` comes out at 14.40. Two pipelines, two defaults |
| `LineWidth` 0.79 vs 0.6; `LabelSize` 9.1 vs 10 | `preview_jacs_house_style` | The `acs_jacs` house style asks for those values |

**`line_width`** — only `apparatus`, six warnings (`0.11` to `0.50`). The drawing is
assembled from ChemDraw's Clipware stationery, whose outlines, fill levels and ground
glass are drawn at their own widths; the figure places those objects and does not restyle
them.

**Run the check per figure type rather than over the whole set at once.** One baseline
cannot be right for a molecule drawing, a house-style scheme and a Clipware apparatus.
`check` offers `--ignore CODE` (repeatable) and no profile flag, so the commands below are
the ones that were actually run and their measured results are stated with them:

```bash
# 1. the ACS 1996 baseline over everything except the three special figures (36 files)
cdx.py check $(ls *.native.cdxml | grep -v -e mechanism -e house_style -e apparatus) \
    --allow-sizes 8,6 --palette any
#    34 pass. nmr and tlc fail on nothing but the 0.01 rounding below.

# 2. a house-style figure sets its own label size and line width
cdx.py check preview_jacs_house_style.native.cdxml --allow-sizes 8,6 --palette any \
    --ignore setting        # passes

# 3. Clipware art carries its own line widths
cdx.py check apparatus.native.cdxml --allow-sizes 8,6 --palette any \
    --ignore setting --ignore line_width    # passes

# 4. the scientific mechanism workflow draws at 28 pt, so its geometry is its own
cdx.py check mechanism.native.cdxml --allow-sizes 8,6 --palette any \
    --ignore setting --ignore bond_length --ignore text_size    # passes
```

`mechanism` used to carry two `bare_terminal_carbon` warnings, for the CH3Br and CH3Cl of its SN2
example. The rule fired when exactly one end of a two-atom fragment was unlabelled — which is the
ordinary skeletal methyl halide — and stayed quiet when *both* ends were unlabelled, a bare line
that identifies nothing. The condition was inverted; it now warns only on the both-unlabelled case,
so a line ending at Br is left alone.

**The rounding is now tolerated, not documented away.** `setting` compares root
attributes with a tolerance of 0.02, so the 0.01 that ChemDraw's SaveAs introduces no
longer reports while a real deviation (28 pt against 14.4) still does. With that in place
`nmr`, `tlc` and `apparatus` pass the baseline and only the figures whose geometry is
deliberately different need `--ignore`.

### `substrate_scope` is the weakest figure in the set

It is kept as the reference for `docs/journal-figure-parity.md` P4-P12: five structures in the
top third of a portrait page, `Yield: 92%` instead of `92%`, no header scheme, no highlighted
varying group, and `LabelSize`/`CaptionSize` at 8 where the baseline is 10. Those are the
things the parity work is meant to close.

### Six types added to close the declared-type gap

`derivatization`, `optimization_table`, `process_route`, `sar_table`, `stereo_model` and
`total_synthesis` had no coverage figure at all: the house-style cards declare 23 types and
these six were unimplemented witnesses. They now exist, and their checks say what they are:

- `total_synthesis`, `process_route`, `optimization_table` and `derivatization` pass the
  ACS 1996 baseline as drawn.
- `sar_table` fails on `LabelSize`/`CaptionSize` 8 against 10, the same defect as
  `substrate_scope`: the scope renderer shrinks the type to fit its grid (parity item P11).
- `stereo_model` fails on a median bond of 12 pt against 14.4, and the conformer templates
  ignore the `bond_length` the brief sets — the median stays 12 whatever is asked for, so the
  figure is drawn at 12 and the check reports the truth. Its three outlier bonds are the
  chair/Newman geometry itself, which the checker already excludes by name.

`optimization_table` is composed from `text` and `box` primitives because no table object
exists; `text` has no width parameter, so anything longer than a short label wraps inside its
own box and cannot be laid out as running text. That is why it carries no footnote.

## Refreshing the native copies

`convert_cdx_cdxml` never writes over an existing file. Given a target that exists it returns
`"ok": false` with `Refusing to overwrite an existing file` and `cdx.py call` exits with 1, but
a batch script that does not look at the result carries on, and the old `.native.cdxml` stays
beside a new PNG and CDXML. This happened once: the native copies were 40 minutes older than
the renders they claimed to match. So, for each figure:

1. Delete `<name>.native.cdxml` (and the intermediate `.cdx`) before converting.
2. Convert `.cdxml` to `.cdx` and back, and stop on any result whose `ok` is false
   (in PowerShell, check `$LASTEXITCODE` after each `cdx.py call`).
3. Compare timestamps: every `.native.cdxml` must be newer than its `.cdxml`.

## Journal parity status

Progress against [journal figure parity](journal-figure-parity.md), per figure.
"Generator done" means the change is in the code with a test; the figure here still shows the
old output until it is rendered again.

| Gap | Figures | State |
| --- | --- | --- |
| P12 value alone under the label (`3a` over `92%`) | `substrate_scope` | generator done; re-render |
| P11 type kept at 10 pt | `substrate_scope`, `sar_table` | generator done; re-render, then `check` should report no `setting` |
| P10 page cropped to its content, footnotes after the last row | `substrate_scope`, `sar_table` | generator done; re-render, then measure the PNG margins |
| P5 a failed entry keeps its value, with an optional note | `substrate_scope` | generator done (`note` on a species); the example brief has no 0% entry yet |
| P6 a group definition under one structure | `substrate_scope` | generator done through the same `note` line |
| P13 twelve to eighteen cells on one page at 10 pt | `substrate_scope`, `sar_table` | generator done: twelve benzamides with two fields, or eighteen ungrouped compounds, now fit one page; each group still starts its own row (see P8) |
| P1 coloured curved arrows routed clear of the structure | `mechanism`, `preview_mechanism` | generator done: ChemDraw arc arrows (`AngularSize`), a hex `color` per move, routing that also avoids bonds; a clockwise arrow renders correctly (head on the tail end). Since then: arrows are routed around one another and `check` reports `arrows_collide` (the red and blue heads now stand apart, confirmed natively); one-electron arrows use ChemDraw's own half head, 600 / 525 / 150, after a 1 : 1 : 0.5 barb rendered as a solid triangle that read as two electrons |
| P2 equilibrium arrows between species | `preview_scheme`, cycles | already in the renderer (`arrow: reversible`, or `<=>`), written as ChemDraw writes it (`tests/runtime/data/tautomer-step.cdxml`); now also tested inside a catalytic cycle |
| P8 section bands | `substrate_scope` | generator done: a rule and an italic name open each group; re-render. Eighteen compounds in three groups of six still take two pages: six band rows at 10 pt need about 717 pt and a page offers about 579, which no padding can close |
| P4 the varying group in a second colour | `substrate_scope` | generator done (`highlight_varying`, the common core by maximum common substructure, default `#FF6600`); its tests need RDKit; re-render with `"highlight_varying": true` in the scope brief |
| P7 a secondary product beside the main one | `substrate_scope` | generator done (`secondary_of` on a species, `secondary_column` in a table): same cell, grey `#808080` (added to the default check palette), its own label and values; same bond length, so the check's bond rule still holds |
| P3 Newman projections inside a scheme | `mechanism` | generator done: a species with `newman` is laid out among the structures (spokes 1.4 bond lengths), re-numbered into the page so its disc still hides the rear spokes, and read back as drawn; the whole-brief test needs RDKit; render a brief with two projections |
| P9 Generic rings | `scheme`, `mechanism` | generator done: a species' `generic_rings` draws a dashed open arc for an unspecified ring fused at a bond; render `assets/parity/p9-generic-ring.json`. The dashed `overlays` around a fragment are kept as a separate annotation and are now written closed |

The briefs to render for the open checks are in `assets/parity/`: `p4-p7-p8-scope.json` (two bands, the varying part highlighted, a secondary product under 3d), `p3-newman-mechanism.json` (two butane conformers beside an addition) `p9-generic-ring.json` (4-pentynoic acid and its methylene lactone, each with a generic ring) and `p9-cycle-overlay.json` (the anisyl ring marked in every species of a Suzuki cycle, with the closed outlines). `p6-family-scope.json` adds two family members with an `Alk = ...` line among plain cells.
