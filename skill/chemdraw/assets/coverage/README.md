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

- **nmr** — `toolkit/docs/examples/scientific/nmr.json` is the right kind of file (a spec whose
  `input` names a spectrum), but the spectrum it names, `source-nmr/1d_data.ft` from the
  [nmrglue integration example](https://nmrglue.readthedocs.io/en/latest/examples/integrate_1d.html),
  is not bundled. Download it beside the spec, or point `input` at another spectrum (CSV,
  NMRPipe, JCAMP-DX, Bruker, JEOL, Varian). A missing spectrum is now reported as
  `spectrum not found: <path>`.
- **apparatus** — the Clipware libraries (`clipware-1.cdxml`, `clipware-2.cdxml`) are licensed
  with ChemDraw and never shipped. Convert a copy of each CTP file from the ChemDraw installation
  (`cdx.py call convert_cdx_cdxml`) into `native-templates/`; see
  [laboratory drawings](../../references/native-laboratory-drawings.md). A missing template is
  now reported by name with that instruction.
- **substrate scope** — `kind: "scope"` writes a project folder rather than one file, so it
  stays in `assets/preview/substrate-scope.png` instead of being duplicated here.

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

### The two figures still missing

Both now fail with a message that names the file to supply (`70e003f`):

- **apparatus** needs `toolkit/docs/examples/scientific/native-templates/clipware-1.cdxml`,
  converted from the Clipware CTP in the ChemDraw installation. That template is
  ChemDraw-licensed, so it must not be committed; only the render may be.
- **nmr** needs `toolkit/docs/examples/scientific/source-nmr/1d_data.ft`, the nmrglue
  example spectrum.
