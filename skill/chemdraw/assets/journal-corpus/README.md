# Journal figure corpus

Real figures from methodology papers, kept as a reference for what a journal figure actually
looks like next to what this Skill produces (`../preview/`, `../coverage/`).

## Licence discipline

Only **CC BY** figures are here. That rule is applied twice, because the first pass was not
enough: an early run selected articles whose licence *text* contained "Creative Commons
Attribution" and kept eleven **CC BY-NC-ND** figures from Nature Communications, which are not
redistributable. They were removed; the filter now rejects anything whose licence contains
`-NC`, `-ND` or "noncommercial". Every row in `INDEX.csv` states its licence, and the whole
file asserts `CC BY` on load.

`INDEX.csv` carries, per figure: journal, DOI, title, figure number, caption and licence —
the attribution CC BY requires. Nothing here is a redrawn copy: each PNG is the publisher's
own image file (Beilstein SVG rasterised through PDF; Nature Communications PNG from the
Springer Nature media server).

## What is here — 25 figures

| Journal | Count |
| --- | --- |
| Beilstein J. Org. Chem. | 16 |
| Nature Communications | 9 |

| Type | Count | What it shows |
| --- | --- | --- |
| reaction_scheme | 7 | Multi-step schemes, conditions above the arrow |
| substrate_scope | 7 | Scope tables: yields under every entry, `0%` kept, variable groups, **colour-coded section bands** in the Nature examples |
| mechanism | 5 | Curved arrows, equilibrium arrows, Newman projections, dashed mapping overlays |
| optimization_table | 2 | Condition screening as a table |
| structure_panel | 2 | Ligand and catalyst families |
| concept_panel | 2 | Overview panels |
| retrosynthesis | 1 | Disconnection arrows |
| assignment_panel | 1 | NOESY evidence drawn onto structures |

## Conventions here that this Skill does not produce yet

Reading the corpus against `../coverage/` gives the gap list:

- **Coloured curved arrows.** Mechanisms distinguish competing pathways with green and red
  arrows.
- **Equilibrium arrows between species** inside a scheme, not only as a step arrow.
- **Newman projections placed inside a scheme** as the substituent argument.
- **Highlighting the varying group** in a scope table (red bond/label in Beilstein, colour
  bands and coloured substituents in Nature Communications).
- **`0%` entries and "not formed"** — the renderer omits a value it was not given; the papers
  keep the failure visible.
- **Variable-group rows** such as `Alk = C7H15, C18H37` under one structure.
- **Secondary products in grey**, drawn beside the main product (Nature scope figures).
- **Section bands** dividing a scope into para/ortho/heterocycle/nucleophile groups, each band
  labelled in italic.
- **Dashed mapping overlays** across a catalytic cycle.

## Extending it

Access differs by publisher and was measured, not assumed:

- **Beilstein** — each figure is its own SVG (`/bjoc/content/inline/*.svg`). Strip the external
  DOCTYPE before parsing or the parser tries to fetch the DTD; rasterise through PDF at 2.2x.
- **Springer Nature** (Nature Communications, Communications Chemistry, Scientific Reports) —
  figure files follow `media.springernature.com/full/springer-static/image/art%3A<DOI>/MediaObjects/<j>_<year>_<id>_Fig<N>_HTML.png`.
  Captions and the licence come from the Europe PMC full-text XML for the same article.
- **Blocked** at the time of collection: `europepmc.org/.../bin/...` figures (403), PMC
  `oa.fcgi` tarballs for these IDs (404), `pubs.rsc.org/image/article/...` (403/404). Do not
  spend time on those paths without re-testing.

## Files

- `<journal>-<volume>-<article>-<figure>.png`, `nature-communications-<pmcid>-fig<N>.png`
- `INDEX.csv` — one row per figure
