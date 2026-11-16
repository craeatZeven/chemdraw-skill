# Journal figure corpus

Real figures from methodology papers, kept as a reference for what a journal figure actually
looks like next to what this Skill produces (`../preview/`, `../coverage/`).

## Licence discipline

Only **CC BY** figures are here. Two incidents shaped this section, and both are recorded
because the rule is easy to get wrong:

1. The first pass kept a figure when its licence *text* contained "Creative Commons
   Attribution" — which is also true of **CC BY-NC-ND**. Eleven non-redistributable figures
   from Nature Communications were collected that way and later removed.
2. The "final guard" that removed them compared publisher prose against the already
   normalised string `CC BY 4.0`, decided *every* row was non-CC-BY, and deleted the PNG
   files. The corpus was restored from git. The guard now only excludes rows, never deletes
   files; deletion is a separate, explicit step.

`INDEX.csv` states a normalised licence per row (`CC BY 4.0`) and the collector rejects any
licence containing `-NC`, `-ND`, "noncommercial" or "no derivatives" before a file is
written. `file` holds the bare name; the image is `<file>.png`.

## What is here — 71 figures, 14 MB

| Journal | Count |
| --- | --- |
| Beilstein J. Org. Chem. | 51 |
| Communications Chemistry | 11 |
| Nature Communications | 9 |

| Type | Count |
| --- | --- |
| reaction_scheme | 10 |
| mechanism | 10 |
| concept_panel | 10 |
| substrate_scope | 7 |
| structure_panel | 6 |
| graphical_abstract | 6 |
| optimization_table | 4 |
| assignment_panel | 3 |
| polymer | 2 |
| coordination_complex | 2 |
| markush | 2 |
| biomolecule | 2 |
| biosynthesis | 2 |
| energy_profile | 2 |
| catalytic_cycle | 2 |
| retrosynthesis | 1 |

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
- **Secondary products in grey**, drawn beside the main product.
- **Section bands** dividing a scope into para/ortho/heterocycle/nucleophile groups, each band
  labelled in italic.
- **Dashed mapping overlays** across a catalytic cycle.

## Extending it

Access differs by publisher and was measured, not assumed:

- **Beilstein** — each figure is its own SVG (`/bjoc/content/inline/*.svg`, thematic issues at
  `/bjoc/series/<n>`). Strip the external DOCTYPE before parsing or the parser tries to fetch
  the DTD; rasterise through PDF at 2.2x.
- **Springer Nature** (Nature Communications, Communications Chemistry, Scientific Reports) —
  `media.springernature.com/full/springer-static/image/art%3A<DOI>/MediaObjects/<j>_<year>_<id>_Fig<N>_HTML.png`.
  The article number in the filename **drops leading zeros** (`01180` → `1180`); getting that
  wrong returns 404 and silently yields nothing. Captions and licence come from the Europe PMC
  full-text XML for the same article.
- **Blocked** when tested: `europepmc.org/.../bin/...` figures (403), PMC `oa.fcgi` tarballs
  for these IDs (404), `pubs.rsc.org/image/article/...` (403/404). Do not spend time on those
  paths without re-testing.

## Files

- `<journal>-<volume>-<article>-<figure>.png`, `nature-communications-<pmcid>-fig<N>.png`
- `INDEX.csv` — one row per figure: bare `file` name, type, figure label, journal, DOI,
  licence, title, caption
