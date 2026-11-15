# Figure coverage

Every figure type this Skill can draw, rendered once on Windows with ChemDraw 26
from the example spec that ships with it. This directory is the evidence behind the
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
| `binary_phase_diagram` | Binary phase diagram, eutectic | yes |
| `binary_phase_diagram_isomorphous` | Binary phase diagram, isomorphous | yes |
| `correlation_diagram` | Correlation diagram between two states | yes |
| `crystal_field` | Crystal-field splitting diagram | yes |
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
| `frost` | Frost circle (polygon inscribed in a circle) | yes |
| `huckel_mo` | Huckel MO diagram (energy levels from the adjacency matrix) | yes |
| `huckel_mo_pyridine` | Huckel MO diagram, pyridine (heteroatom parameter) | yes |
| `jablonski` | Jablonski diagram (S0/S1/S2/T1, A/F/P, IC/ISC) | yes |
| `latimer` | Latimer diagram (reduction potentials in series) | yes |
| `lewis` | Lewis structure (CO2, lone pairs and bonds) | yes |
| `lewis_carbonate` | Lewis structure, carbonate with formal charges | yes |
| `mechanism` | Mechanism with curved arrows | yes |
| `mo_diagram` | Molecular-orbital diagram (O2, 2s/2p mixing, spin filling) | yes |
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
| `tanabe_sugano` | Tanabe-Sugano diagram (d3, C/B = 4.5) | yes |
| `tlc` | TLC plate with lanes and Rf spots | yes |
| `unit_cell` | Crystal unit cell | yes |
| `walsh_diagram` | Walsh diagram (geometry vs orbital energy) | yes |

## What this run did not cover

- **nmr** — the example shipped for it (`toolkit/docs/examples/scientific/nmr.json`) is a
  parameter file, but `scientific nmr` takes a spectra file (CSV, NMRPipe, JCAMP-DX,
  Bruker, JEOL, Varian). Nothing was rendered for it here.
- **apparatus** — needs `native-templates/clipware-1.cdxml`, which is not in the tree on
  this branch.
- **substrate scope** — `kind: "scope"` writes a project folder rather than one file, so it
  stays in `assets/preview/substrate-scope.png` instead of being duplicated here.

## Note for whoever runs this again

`cdx.py scientific` prints its report to stdout and crashes with `UnicodeEncodeError`
when the run contains a character the active console code page cannot encode (a true
minus sign `−` is enough on a GBK console). The CDXML is already written at that point,
so the failure only costs the PNG: re-render with `cdx.py render <file>.cdxml`.

