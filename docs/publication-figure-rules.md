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
