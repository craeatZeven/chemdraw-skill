# When to abbreviate a group, and when it must stay expanded

Round-2 work item 3. Rules for the `abbreviate` species field and the shorthand groups
(`Ph`, `Boc`, `TBS`, …), derived from the 26 journal briefs and the published figures they
imitate.

## What the corpus actually contains (measured)

Source: `output/round2-stats.json`, computed from `assets/scope-benchmark/briefs/*.json`.

* **0 of 26 briefs use `abbreviate`.** The field exists and works, but the benchmark never
  exercises it — so these rules are a new decision procedure, not a description of current
  behaviour.
* **3 briefs use `definition`** instead (`bjoc-articles-22-127-scheme2`,
  `communications-chemistry-s42004-022-00764-7-fig2`,
  `nature-catalysis-s41929-025-01405-7-fig2`). A `definition` line such as `Alk = C7H15` draws
  the group as its label and names it under the compound. That is the corpus's real vehicle for
  "contract this group", and it already carries the definition the reader needs.
* **294 species** in total; **9 briefs hold more than 12 species** (dense grids, where space is
  the binding constraint).
* **13 of 26 briefs contain stereochemistry** (`@`, `/`, `\` in the SMILES).
* How often a shorthand-able group repeats across the 294 species:

  | group | species containing it | briefs |
  | --- | --- | --- |
  | Ph (`c1ccccc1`) | 66 | 16 |
  | Ac (`CC(=O)`) | 38 | 11 |
  | OMe (`COc`) | 36 | 14 |
  | Bn (`Cc1ccccc1`) | 22 | 8 |
  | CF3 (`C(F)(F)F`) | 13 | 6 |
  | iPr (`CC(C)C`) | 11 | 2 |
  | tBu | 4 | 1 |
  | Boc | 2 | 1 |

* The toolkit resolves **38 shorthand labels**: Ac, Bn, Boc, Bpin, Bz, CF3, CN, CO2Et, CO2Me,
  Cbz, Et, Me, Ms, N3, NHBoc, NMe2, NO2, OAc, OBn, OBz, OCF3, OMe, OMs, OTBS, OTIPS, OTMS, OTf,
  OTs, Ph, PhthN, Piv, TIPS, TMS, Tf, Ts, iPr, nBu, tBu (`house_styles.ABBREVIATIONS`).
  Validation rejects anything outside that set, because an unresolvable label would break the
  chemical readback.

## Abbreviate only when ALL of these hold

1. **The label is resolvable.** It is one of the 38 known groups. A shorthand the reader cannot
   expand back to a structure is a defect, not a space saving.
2. **It is an auxiliary, not the subject.** Protecting and masking groups (Boc, Cbz, Bn, Ts,
   Ac, TMS, TBS, TIPS, Piv), a generic aryl (Ph), small esters (CO2Me, CO2Et). The group is
   there to make the chemistry work, not to be read.
3. **It repeats.** In at least **3 species**, or at least **twice in one structure**. The 66
   phenyls and 38 acetyls above are the shape of that case; the single Boc in the whole corpus
   is not.
4. **Space is actually tight.** The panel holds **more than 12 species**, or the drawn group
   would be the widest thing in its cell. Contracting one group in a three-compound scheme buys
   nothing and costs readability.
5. **It is defined once.** Either a `definition` line, a footnote, or the header scheme shows
   the group expanded somewhere the reader can find it. Abbreviation is a cross-reference, so
   there has to be something to refer to.

## Must stay expanded when ANY of these holds

1. **At the R-group definition, and on first appearance.** The place where `R¹`, `PG` or the
   shorthand is introduced is where the reader learns what it means. Never contract the
   defining instance; contract the repetitions after it.
2. **The group is the variable under study.** If the series reads `PG = Boc / Ac / Cbz`, the
   protecting group *is* the result. Contracting it would delete the information the figure
   exists to convey. (Demo `abbrev_must_expand`: the right-hand panel is the mistake, shown on
   purpose.)
3. **The group carries or determines stereochemistry.** A stereocentre inside the group, or a
   wedge whose reading depends on the group's atom order. 13 of the 26 briefs are in this
   position. A shorthand hides the connectivity the stereo descriptor refers to.
4. **The group sits at the reaction centre.** Atoms whose bonds are made or broken are drawn
   out, always.
5. **The label would be ambiguous.** TMS (trimethylsilyl vs tetramethylsilane), Ts/Tf, Ac
   (acetyl vs actinium), `R` used for two different things on one page.
6. **The group is what a highlight is pointing at.** If `fragment_colors` or
   `highlight_varying` marks it, the reader needs to see the atoms, not a label.

## What the toolkit does with it

`abbreviate: ["Boc"]` on a species contracts that group to a `NodeType="Fragment"` node with an
inner fragment, so the label still resolves to a SMILES on readback. Both demo panels report
`readback=matched`, which is the point: contraction is a drawing decision and must not change
the chemistry.

## Demos

Script: `scripts/round2_demos.py` (`abbrev_repeated`, `abbrev_must_expand`).

* `output/round2-demos/abbrev_repeated.compare.png` — same three Boc-protected compounds, left
  with the tert-butoxycarbonyl drawn out, right with `abbreviate = ["Boc"]`. Verified by eye:
  the left crop shows the tert-butyl + O + carbonyl chain; the right crop reads **Boc**.
* `output/round2-demos/abbrev_must_expand.compare.png` — a series whose protecting group is the
  variable. Left is correct (all expanded); right applies `abbreviate` and is shown as the
  counter-example of rule 2.
