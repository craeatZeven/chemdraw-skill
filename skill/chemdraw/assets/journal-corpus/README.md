# Journal figure corpus

Real figures from methodology papers, kept as a reference for what a journal figure actually
looks like, next to what this Skill produces (`../preview/` and `../coverage/`).

## Provenance and licence

Every figure comes from **Beilstein Journal of Organic Chemistry**, which is Diamond Open
Access: all articles are **CC BY 4.0**, so the figures may be redistributed with attribution.
`INDEX.csv` carries, for each file, the article DOI, title, figure number, its caption, and
the licence, which is the attribution CC BY requires. Nothing here is a redrawn or edited
copy: each PNG is a rasterisation of the SVG the publisher serves for that figure.

Collected by fetching the article pages, taking the per-figure SVG that Beilstein publishes
(`/bjoc/content/inline/*.svg`, `/bjoc/content/figures/*.svg`), stripping the external
DOCTYPE (the parser otherwise tries to fetch the DTD) and rasterising through PDF at 2.2x.
The figure type in `INDEX.csv` is derived from the caption keywords, so treat it as a
starting point rather than a classification.

## What is here

| Type | Count | What it shows |
| --- | --- | --- |
| mechanism | 4 | Curved arrows on a real catalytic cycle, coloured by mode, with equilibrium arrows and Newman projections |
| substrate_scope | 5 | Scope tables: varying substituents highlighted in red, yields under every entry, `0%` entries kept, variable groups (`Alk = C7H15, C18H37`) |
| structure_panel | 2 | Families of ligands and catalysts laid out as a panel |
| reaction_scheme | 3 | Multi-step schemes with conditions above the arrow |
| retrosynthesis | 1 | Retrosynthetic analysis with disconnection arrows |
| assignment_panel | 1 | NOESY evidence drawn onto the structures |

## Conventions in these figures that this Skill does not produce yet

Reading the corpus against `../coverage/` gives the gap list:

- **Coloured curved arrows.** Mechanisms distinguish competing pathways with green and red
  arrows; the renderer draws one arrow colour per move at best.
- **Equilibrium arrows inside a scheme** (double half-headed) between species, not only as a
  step arrow.
- **Newman projections placed inside a scheme** as a substituent-level argument.
- **Highlighting the varying group in a scope table** (red bond/label), the convention that
  makes a scope readable at a glance.
- **`0%` entries and "not formed" labels** — the renderer omits a value it was not given,
  the papers keep the failure visible.
- **Variable-group rows** such as `Alk = C7H15, C18H37` under one structure.
- **Highlighted rings** marking a notable substrate in a scope grid.
- **Dashed mapping overlays** in mechanism figures, showing which fragment is which across a
  cycle.

Each of those is a candidate for the renderer next to the layout work already done, and each
can be measured from the figures here rather than guessed.

## Files

- `<journal>-<volume>-<article>-<figure>.png` — the figure
- `INDEX.csv` — one row per figure: `file, figure_type, figure, journal, doi, license, title, caption`
