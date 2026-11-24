# Publication Figure Rules & Cross-Paper Comparison Methodology

This document codifies the explicit layout, typography, and colour rules derived from comparing published open-access figures of the same figure type across journals, and defines how to expand automated multi-type figure generation systematically.

## 1. Cross-Paper Comparison Workflow

Every automated figure rule should originate from **comparing multiple published papers of the same figure type**, never from ad-hoc taste:

1. **Collect CC BY 4.0 witnesses by figure type** into `assets/journal-corpus/` (`INDEX.csv`, licence strictly `CC BY 4.0`, rejecting `-NC` and `-ND`). Target at least **3–5 papers across ≥2 publishers** per figure type.
2. **Compare same-type figures side by side** along five fixed dimensions:
   - **Geometry & proportion**: bond length, cap height, arrow length / bond length, margins, page cropping, multi-row wrapping or grid density.
   - **Typography & formatting**: bold compound numbers vs plain formulas, italic locants/stereochemical descriptors (`cis`, `trans`, `o`/`m`/`p`, `N`-), subscript rules, superscript footnote letters (`^a`), yield placement.
   - **Semantic colour system**: how colour encodes chemical meaning (synthon provenance, newly formed bonds, varying substituent, per-band theme, subordinate by-product in grey) — never decorative rainbow colouring.
   - **Structural conventions**: abbreviation contractions, generic rings, equilibrium/retrosynthesis arrows, Newman/chair insets, brackets, section bands.
   - **Failure & edge-case display**: `0%` / `trace` / `not formed` entries, secondary by-products, shared structures for homologous families (`Alk = C7H15, C18H37`).
3. **Encode measurable rules in three places at once**:
   - **Generator** (`draw_figure` brief kinds: `scheme`, `scope`, `mechanism`, `panels`, `table`; `scientific` diagram generators; `acs_builder`).
   - **Style & palette cards** (`toolkit/cdxml_toolkit/publication/house_styles/*.json` + `palettes.py`).
   - **Machine checker** (`cdxml_check.py` via `cdx.py check --profile <style>` and `cdx.py validate --native`).
4. **Verify on native ChemDraw renders**: inspect PNGs at print size and verify `.native.cdxml` round-trip (`molecule_count_delta == 0`, zero `text_over_ink`).

---

## 2. Semantic Colour Palettes & Colour Discipline

Colour in publication chemistry figures is **semantic, restrained, and print-safe**:

- **Hard budget**: use at most **3–4 semantic hues** per figure (plus `#000000` ink, `#FFFFFF` background, and `#808080` grey for subordinate by-products or spectators).
- **Synthon provenance rule**: when colouring substrates (e.g. electrophile blue, nucleophile orange/red), carry the exact same hues into the corresponding fragments of the product and across every cell of a scope panel.
- **Per-band consistency rule**: in multi-section substrate scopes (e.g. *para* substituents vs *C*-nucleophiles), the band rule, italic band heading, and highlighted varying substituents inside that band share one accent hue (`group_colors`).
- **Optimal-entry highlight rule**: in optimization tables, highlight the optimal row with bold + the active palette's primary hue (matching the catalyst/parameter highlighted on the header reaction arrow when present).

Built-in semantic palettes (`toolkit/cdxml_toolkit/publication/palettes.py`, accepted by `cdx.py check --palette <name>` and derived automatically by `--profile <house_style>`):

| Palette | Primary | Secondary | Accent / New Bond | Muted | Use when |
| --- | --- | --- | --- | --- | --- |
| `okabe_ito` | `#0072B2` | `#D55E00` | `#009E73` / `#D55E00` | `#808080` | Default colourblind-safe figures & multi-pathway plots |
| `nature` | `#0078B2` | `#D55E00` | `#E67E22` / `#FF0000` | `#808080` | Nature, Nat. Chem., Nat. Commun., Nat. Synth. |
| `acs_jacs` | `#800000` | `#004080` | `#008060` / `#DC3220` | `#808080` | JACS, Org. Lett., J. Org. Chem., JACS Au |
| `wiley_angew` | `#005AB5` | `#DC3220` | `#008060` / `#DC3220` | `#808080` | Angew. Chem., Chem. Eur. J., Adv. Synth. Catal. |
| `beilstein` | `#D62728` | `#1F4E9E` | `#008080` / `#D62728` | `#808080` | Beilstein J. Org. Chem. |
| `rsc` | `#003B73` | `#C83E37` | `#2E7D32` / `#C83E37` | `#808080` | Chem. Sci., Chem. Commun., Org. Biomol. Chem. |
| `acs1996` | `#336699` | `#FF6600` | `#FF6600` | `#808080` | Default ACS Document 1996 checker baseline |

---

## 3. Typography, Formatting & Native-Save Rules

1. **Never shrink body text to fit**: keep labels and captions at the style's body size (10 pt for ACS 1996, 6–7 pt for Nature/RSC). Grow or paginate the canvas instead.
2. **Strict formula subscripting**: digits in chemical formulas (`K2CO3`, `PhB(OH)2`, `OD600`) are subscripted everywhere (conditions, table cells, species formula labels, captions). Compound identifiers (`3a`, `TS2`, `INT3`, `L1`) and mechanism designators (`Sn1`, `Sn2`) stay on the baseline.
3. **Italic locants & footnotes**: locants and stereochemical prefixes (`ortho`, `meta`, `para`, `cis`, `trans`, `N-`, `O-`, `C-`) are italic; footnote markers (`^a`, `^b`) are italic superscripts.
4. **Native-save phantom-fragment guard**: every non-atom `<t>` caption (conditions, labels, table cells, footnotes) must carry `InterpretChemically="no"` so ChemDraw's SaveAs never turns a condition label such as `H2O` or `MeOH` into a hidden molecule fragment. Verify with `cdx.py validate FILE.cdxml --native` (`molecule_count_delta == 0`).
5. **Zero text-over-ink**: rules, arrows, curves, and bounding boxes must keep measurable clearance from every caption box (`cdx.py check --profile <style>` reports `text_over_ink` warnings).

---

## 4. Figure-Type Coverage & Automation Matrix (All 23 Declared Types)

| Figure Type | Primary Automated Entry | Key Layout & Colour Rules Enforced | Corpus Status & Priority Next Steps |
| --- | --- | --- | --- |
| `reaction_scheme` | `cdx.py draw` (`kind: "scheme"`) | Auto linear/continue-wrap/convergent/divergent/fan, 28 house styles, formula subscripts, `InterpretChemically="no"` | Well covered (10+ witnesses) |
| `substrate_scope` | `cdx.py draw` (`kind: "scope"`) | Omitted title (`title: ""`), centred colour-coded general reaction header, row-aligned caption baselines, `columns` + `row_align: "center"|"justify"`, SMARTS `fragment_colors` synthon provenance, `highlight_varying`, `band_style: "rule"|"pill"`, `border`, `0%` kept, `drawn_with` families, grey `secondary_of`, `group_colors`, italic `*note*` and superscript footnotes | **26 CC BY 4.0 witnesses across 6 journals** + 26 imitation briefs in `assets/scope-benchmark/` |
| `optimization_table` | `cdx.py draw` (`kind: "table"`) | Optional general reaction header, booktabs 3-line rules with zero `text_over_ink`, Arial-metric column widths, formula subscripts, `^a` footnotes, optimal-row highlight (`bold`/`color`/`shaded`/`boxed`) | 3 witnesses; embeds inside `kind: "panels"` |
| `mechanism` | `cdx.py draw` (`kind: "mechanism"`) | Native ChemDraw arc arrows, per-arrow pathway colours, collision avoidance, inline Newman projections (`attach`), equilibrium arrows, `generic_rings` | Well covered (10 witnesses); unify 14.4 pt default with `scientific mechanism` |
| `catalytic_cycle` | `cdx.py draw` (`layout: "catalytic-cycle"`) | Native circular arc arrows, radial enter/leave species, closed dashed `overlays`, `generic_rings`, `cycle_arrow_color` | 1 witness — collect ≥2 more true catalytic-cycle witnesses |
| `retrosynthesis` | `cdx.py draw` (`arrow: "retro"`, `disconnections`) | Native open double-shaft arrows (`⇒`), wavy/dashed/coloured bond disconnection cuts, serpentine/convergent wrapping | 1 witness — collect ≥2 more multi-step retrosynthesis witnesses |
| `concept_panel` | `cdx.py draw` (`kind: "panels"`) | `stacked` or `side_by_side` multi-column grids, `box`/`rule`/`plain` titles, `accent` colour (`this_work_color`), embedded `scheme`/`mechanism`/`table`/`energy_profile`/`gallery`/`sections`, coloured footer `legend` | Well covered (10 witnesses) |
| `structure_panel` | `cdx.py draw` (`kind: "panels"`, `gallery`) | Uniformly spread multi-row structure grids with centred multi-line captions and motif highlighting | Covered (6 witnesses) |
| `sar_table` | `cdx.py draw` (`kind: "scope"` / `kind: "table"`) | Multi-field structure cards or tabular SAR with 10 pt type preserved | **0 corpus witnesses** — collect 3–5 open-access SAR tables |
| `total_synthesis` | `cdx.py draw` (`wrap: "continue"`/`"serpentine"`) | Long routes with shared intermediates drawn once, lettered condition keys | **0 corpus witnesses** — collect 3–5 multi-step synthesis routes |
| `process_route` | `cdx.py draw` (`kind: "scheme"`) | Multi-step scalable route with scale annotations & PMI (`route_summary`) | **0 corpus witnesses** — collect open-access process routes |
| `derivatization` | `cdx.py draw` (`layout: "fan"`/`"divergent"`) | Central product branching radially (`fan`) or stacked (`divergent`) | **0 corpus witnesses** — collect 3–5 late-stage functionalization schemes |
| `energy_profile` | `cdx.py draw` (`kind: "energy_profile"`) / `scientific energy_profile` | Quantitative ΔG / ΔG‡ validation, levels & smooth styles, automatic `house_style` `pathway_colors` / `energy_unit` / `level_style`, named `palette`, embeddable inside `panels` | 2 witnesses — collect ≥1 more across another publisher |
| `photophysical_diagram` | `cdx.py scientific jablonski` | S₀/S₁/S₂/T₁ state checks, A/F/P/IC/ISC arrows | **0 corpus witnesses** — collect open-access photophysics figures |
| `stereo_model` | `cdx.py build` (`conformer`) / `draw` (`newman`) | Chair, Newman, Haworth, Fischer, norbornane templates at requested `bond_length` (14.4 pt ACS baseline) | **0 corpus witnesses** — collect open-access stereochemical model figures |
| `coordination_complex` | `cdx.py draw` (`dative_bonds`) / `build` (`multi_attachments`) | Native `Order="dative"` bonds, $\eta^n$ `MultiAttachment`, metal colour | 2 witnesses; wire `coordination_complex.metal_color` from house styles |
| `polymer` | `cdx.py build` (`polymer_bracket`) | Visual square/round repeat brackets with subscript `n` | 2 witnesses; add first-class bracket support in `draw_figure` |
| `markush` | `cdx.py markush` | R-group labels, variable rings, `MultiAttachment`, CXSMILES + V3000 | 2 witnesses |
| `graphical_abstract` | `cdx.py draw` (`kind: "panels"`, `toc: true`) | High-contrast focal scheme + ACS TOC size check (`82.5 × 44.5 mm` = `234 × 126 pt`) | Covered (6 witnesses) |
| `assignment_panel` | `cdx.py scientific nmr` / `compose` | Peak assignment labels, numbered atoms, inset zoom | Covered (4 witnesses) |
| `biomolecule` / `biosynthesis` | `cdx.py draw` / `build` | Haworth sugars, italic enzyme labels above arrows | 2 + 2 witnesses |
| `setup_diagram` | `cdx.py scientific apparatus` / `tlc` | Clipware apparatus templates, native TLC plates with Rf & stain colours | **0 corpus witnesses** |

Run `python scripts/audit_corpus.py --markdown` at any time to verify CC BY 4.0 licence purity and track remaining corpus gaps across all 23 types.

---

## 5. Verification Checklist for Any Automated Figure

```bash
# 1. Inspect schema & minimal example before writing a brief
python skill/chemdraw/scripts/cdx.py schema scope

# 2. Draw & render
python skill/chemdraw/scripts/cdx.py draw brief.json out.cdxml --render

# 3. Check house rules against the exact target style profile
python skill/chemdraw/scripts/cdx.py check out.cdxml --profile nature_portfolio --allow-sizes 8,10

# 4. Validate structural integrity (add --native --cross-reader on machines with ChemDraw)
python skill/chemdraw/scripts/cdx.py validate out.cdxml

# 5. Run the scope rendering benchmarks
#    coverage (default): render every assets/coverage/*.cdxml, diff against the
#    committed baseline PNGs (RMSE + differing-pixel ratio), report per case:
python scripts/run_scope_benchmark.py --render --side-by-side --output-dir output/scope-benchmark
#    coverage audit only (no ChemDraw):
python scripts/run_scope_benchmark.py --no-render
#    26-case cross-journal substrate_scope imitation benchmark (toolkit runtime):
python scripts/run_scope_benchmark.py --corpus journal --no-render
python scripts/run_scope_benchmark.py --corpus journal --render --side-by-side

# 6. Audit corpus licence discipline and coverage across all 23 types
python scripts/audit_corpus.py --markdown
```

The journal runner uses the scope entry point for scope briefs so that options
such as `row_align` and `fragment_colors` are validated as scope settings. Run it
with the runtime built from the same checkout as the corpus. From a checkout,
`scripts/install.py --runtime-only --apply` builds its local `toolkit/`.
Use `cdx.py doctor` to check the installed fingerprint if the runner reports
unknown brief fields after a checkout change.

---

## 6. Hard Rules Distilled from an 83-Figure Textbook Run (2026-11-24)

> Source: a full textbook illustration run (83 figures, 12 iterative rounds, cross-model split:
> one model draws, another reviews). Every rule below is **measured**, not inferred;
> the long form with symptoms/root causes lives in [`lessons-2026-11-24-figures.md`](./lessons-2026-11-24-figures.md).

### 6.1 Round-trip "matched" is not render equivalence

A CDXML can parse back with every atom, bond, text and symbol present, and still **lose an element when
ChemDraw renders it**. In that run a radical dot lived in a `<graphic GraphicType="Oval">` node: readback
said "present", ChemDraw drew nothing, and 12 rounds of review never noticed.

**Check before sign-off**: render the *delivered* CDXML with ChemDraw itself and diff the ink against the
delivered PNG (per-element bounding boxes, not just file existence). If a PNG was produced by a built-in
or patched renderer, tag it (`render_source: chemdraw|builtin|patch`) and never mix sources inside one set.

### 6.2 Decorative `<graphic>` primitives that ChemDraw does not draw

- **Free radical dot** → use the atom property `<n Radical="doublet"/>`; ChemDraw positions the dot itself.
  Do **not** synthesise it from an Oval graphic or a text glyph.
- **Formal charge** → `<graphic GraphicType="Symbol" SymbolType="CirclePlus|CircleMinus">`.
- **Text glyphs as marks are unreliable**: `·` in Arial renders as a tofu square, `·` in SimSun explodes in
  size at larger point sizes; `•` renders but its size/baseline follow the font.
- Adopt a **native-object allowlist** (atoms, bonds, `<arrow>`, `<bracket>`, `<graphic SymbolType=...>`,
  orbitals). Schematics built from bare lines + text look right but are invisible to ChemDraw's chemical
  model and to every readback check.

### 6.3 Export resolution has two different controls

| Path | Behaviour |
| --- | --- |
| `export_image(dpi=…)` (MCP connector) | fixed 820×388 @144 dpi — the `dpi` argument has no effect |
| `app.Preferences.PNGResolution` + `Document.SaveAs("x.png")` (COM) | honoured: 300 → 3417×1613, 600 → 6834×3224 |

Also: `SaveAs` supports `cdxml / svg / emf / tif / eps / png` and **silently produces nothing for `pdf`**.
A ready-made wrapper lives in the authoring repo as `render_cdxml.py` (ChemDraw COM, resolution as a flag).

### 6.4 Marker clearance, not just proximity to the atom

A radical/charge marker must clear **everything**, not merely sit near its atom:

- distance to any straight bond segment ≥ `radius + 0.6` units
- distance to any text anchor ≥ `radius + 1.5` units
- exclude large ovals (orbital lobes) from the check, or every orbital figure becomes a false positive
  (first implementation flagged 18 of 83 figures that way)

Placement recipe that worked: enumerate candidate points in the upper half-plane around the atom, keep those
with ≥ 11 px clearance at 300 dpi, then take the smallest radius (closest to the atom). The conventional
"upper-right superscript" position is often **occupied by the C–H bond** — if so, go directly above the atom.

### 6.5 Schematic arrows: verify direction and landing, not existence

Measured defects from that run:

- an electron-push arrow drawn **bond → lone pair** instead of lone pair → bond (tail was on the C–Cl bond
  midpoint, head pointed at the lone-pair symbol)
- an arrow whose head landed on the ring centre, 5.6 units away from the intended ortho carbon
- rotation arrows faked with `<curve>` objects instead of `<arrow>` (so they carry no arrowhead semantics)

**Checks**: arrow head within **2 units** of the target atom/bond midpoint; tail at the electron source;
conrotatory/disrotatory sense derived from the arc geometry, not from the label text.

### 6.6 Regeneration must be idempotent; PNG edits must be written back

- Any pixel-level repair of a PNG **must** be mirrored into the CDXML, otherwise the next re-render silently
  reverts it. Tag such figures `render_source: patch` and downgrade their audit row.
- Repair only inside the target radius (paint white, redraw). A full-region paste once erased the atom glyph
  and parts of two bonds without any error.
- Re-rendering the same CDXML twice must be byte-comparable; if it is not, the pipeline is nondeterministic.

### 6.7 Audit ledgers need a fixed schema and a monotonicity gate

In that run a later round rewrote `audit/summary.csv` with a **new header**, populated only its own rows,
and the "OK" count fell from 51 to 34 — every later round then rebuilt from the worse baseline.

- Fix the schema (version it); schema changes ship with a migration.
- After every round assert `OK_count >= previous` and `artifact_count >= previous`; on regression, restore the
  missing rows/files from the reference branch and fail the round.
- Track the evidence class per row (`basis: readback+caption | readback-only | carried | unresolved`);
  only `readback+caption` may be reported as OK.

### 6.8 Silent-failure list (exit 0, nothing done)

Observed in one session, all with exit code 0:

1. `worktree add` failed (dir exists / branch in use) and the script continued → every later path was missing
2. `SaveAs` to an unsupported format produced no file and no exception
3. the reviewer had no declared-artifact list, so it scanned `docs/` before/after images and passed a figure
   whose real PNG had never been regenerated
4. logs written as UTF-16 by `Tee-Object` were read as UTF-8, so a completion marker never matched (2 h idle)

**Rule**: after every external call, assert *the artifact exists and is non-empty*; never trust the exit code.
