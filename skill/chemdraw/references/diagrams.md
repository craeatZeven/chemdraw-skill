# Chemistry diagrams

Generators for the diagrams that are not reaction schemes: energy profiles, enthalpy cycles,
Jablonski diagrams, MO diagrams, Hückel pi orbitals, orbital boxes, crystal-field splitting,
Lewis structures, cubic and hexagonal unit cells, one-component and binary phase diagrams,
Pourbaix, Latimer, Frost and Ellingham diagrams, Woodward–Hoffmann correlation diagrams, Walsh
diagrams, potential-energy-surface contours, flow charts, electrochemical cells and the
periodic table. Each one reads a JSON
spec, writes one editable ACS 1996 CDXML and prints a report of what the drawing implies.

```bash
cdx.py scientific energy_profile spec.json new.cdxml --render --report new.json
cdx.py check new.cdxml                     # house rules: settings, fonts, palette, references
```

Example specs for every operation ship in `toolkit/docs/examples/scientific/`. Text in every spec
takes `{sub:..}`, `{sup:..}`, `{i:..}`, `{b:..}` and `{bi:..}` markup; Chinese runs are set in
SimSun automatically. Colours are `black`, `blue` (#336699), `orange` (#FF6600) or `#RRGGBB`.

## The chemistry is computed and checked, not drawn by hand

A spec states the facts (energies, the element, bonds, the structure type). The generator places
the electrons, derives the counts and refuses a spec that contradicts itself. Refusals are a
feature: fix the spec rather than looking for a way around them.

| Operation | Computed | Refused |
| --- | --- | --- |
| `energy_profile` | TS/minimum per point, forward and reverse barriers, overall change | a point marked `"kind": "ts"` that is not above both neighbours |
| `energy_cycle` | missing level energies from step values | step values that do not close the cycle (Hess's law) |
| `jablonski` | energy change per transition | fluorescence between singlet and triplet, ISC between same-spin states, absorption downhill |
| `mo_diagram` | MO filling, bond order, unpaired electrons, magnetism | MO electrons that differ from the atomic total |
| `electron_configuration` | ground state (Cr, Cu, Pd ... exceptions), cation 4s-before-3d loss, unpaired electrons | box occupations that break Pauli or Hund (unless `allow_violations`) |
| `crystal_field` | d count from the ion, high/low-spin filling, spin-only moment, CFSE | an ion outside the d block |
| `lewis` | lone pairs from valence − bonds − formal charge; octet per atom | more bonds than valence electrons, charges that do not add up to `charge` |
| `unit_cell` | atoms per cell (corner 1/8, edge 1/4, face 1/2), formula, coordination numbers | the same site listed twice |
| `huckel_mo` | pi energies (alpha + x beta), coefficients, nodes, HOMO/LUMO/SOMO, pi and delocalisation energy, 4n+2 rule | more electrons than orbitals can hold |
| `phase_diagram` | (schematic) | a critical point not above the triple point |
| `binary_phase_diagram` | the phase field of each marked point; tie line, phase compositions and amounts (lever rule) | an eutectic above a melting point, solubility limits on the wrong side of the eutectic, a liquidus below the solidus |
| `pourbaix` | the stable species in every pH–E region, exactly, from Gibbs energies of formation; species that never appear | fewer than two species |
| `latimer`, `frost` | non-adjacent E° (Σ nE° / Σ n), volt-equivalents, species that disproportionate | species not in falling oxidation state, a missing potential |
| `ellingham` | ΔG° per mole of O{sub:2} along T, crossing temperatures (which element reduces which oxide) | a phase change that makes ΔG° jump instead of kink |
| `correlation_diagram` | orbital links by the non-crossing rule, thermally allowed or forbidden | reactant and product with different symmetry counts |
| `walsh_diagram` | total orbital energy along the coordinate, preferred geometry | an energy list that does not match the coordinate |
| `potential_surface` | contours, minima and first-order saddle points (transition states) with energies and barriers | a malformed grid |
| `electrochemical_cell` | polarity, electron flow, ion migration from the cell type | (nothing to contradict: none of it is stated) |

Warnings in the report (an incomplete octet, a label that sits on other ink, levels spread for
legibility, a not-to-scale note) are listed in `warnings`; read them before delivering.

## energy_profile

```json
{"style": "smooth", "unit": "kJ/mol", "show_values": true,
 "pathways": [
   {"label": "uncatalysed", "points": [{"label": "R", "energy": 0}, {"label": "TS", "energy": 120}, {"label": "P", "energy": -45}]},
   {"label": "catalysed", "dashed": true, "points": [{"label": "R", "energy": 0}, {"label": "TS1", "energy": 62},
     {"label": "Int", "energy": 18}, {"label": "TS2", "energy": 48}, {"label": "P", "energy": -45}]}],
 "arrows": [{"from": "R", "to": "TS", "label": "E{sub:a}"}, {"from": "R", "to": "P", "label": "ΔH", "heads": "to"}]}
```

`style` is `smooth` (one curve) or `levels` (a bar per species, dashed connectors). A single
pathway may be given as `"points"`. Per point: `id` (needed when labels repeat), `kind`
(`ts` or `minimum`, inferred when omitted), `x` (0–1 along the coordinate), `label_side`
(`above`, `below`, `left`, `right`), `show_value`, and `structure`: a CDXML file drawn beside
the point (the stationary-point structure; relative to the spec). Arrows: `pathway` (index),
`heads` (`both`, `to`, `none`), `show_value`, `dx`. Also `x_label`, `y_label`, `y_ticks`,
`value_format` (`"{:.1f}"`), `ts_marker` (the ‡), `width`, `height`, `step`, `legend`.

## energy_cycle (Born–Haber, Hess)

`levels`: `{"id", "label", "energy"}`, `energy` optional when the steps fix it. `steps`:
`{"from", "to", "label", "value", "color", "dashed"}`, drawn left to right in order.
`tolerance` (default 0.5) is how far a stated value may differ from the level difference.

## jablonski

`{"preset": "standard"}` draws S0, S1, S2, T1 with absorption, IC, fluorescence, ISC and
phosphorescence. Otherwise `states`: `{"id", "label", "energy", "multiplicity", "vib"}` and
`transitions`: `{"from", "to", "kind", "label", "color"}`, `kind` one of `absorption`,
`fluorescence`, `phosphorescence`, `ic`, `isc`, `vr`. Radiative transitions are straight arrows,
non-radiative ones wavy.

## mo_diagram

`{"preset": "O2", "charge": -1}` covers H2 to Ne2 with the σ/π order switching after N2. For
any other molecule give `left`, `right` and `center`, each `{"title", "levels": [...]}`; a level
is `{"id", "label", "energy", "degeneracy", "electrons"}`, and centre levels add
`"from": [atomic ids]` (the dashed correlation lines) and `"character"` (`bonding`,
`antibonding`, `nonbonding`; needed for the bond order). `"fill": true` fills the MOs from the
atomic electrons. Preset energies are a schematic order, not to scale.

## electron_configuration

`{"element": "Fe", "charge": 3}` or `{"ion": "Fe3+"}`; `core` is `auto`, `full`, `noble`
(`[Ar]` and the rest) or `valence`; `layout` is `row` or `energy` (subshells stacked by Aufbau
order). Explicit occupations: `"subshells": [{"label": "2p", "boxes": ["ud", "u", "u"]}]` or
`{"label": "2p", "electrons": 4}`.

## crystal_field

`{"geometry": "octahedral" | "tetrahedral" | "square_planar", "ion": "Co3+" or "d_electrons": 6,
"spin": "high" | "low"}`. Octahedral d4–d7 without `spin` draws high spin and warns that both
exist. `orbital_labels` names each d orbital (default on for square planar).

## lewis

```json
{"atoms": [{"id": "O1", "element": "O", "at": [0, 0]}, {"id": "C", "element": "C", "at": [1, 0]},
           {"id": "O2", "element": "O", "at": [2, 0]}],
 "bonds": [["O1", "C", 2], ["C", "O2", 2]]}
```

Atoms sit on a square grid (`at`: column, row); bonded atoms are grid neighbours. Per atom:
`formal_charge`, or `lone_pairs`/`radical` to state the non-bonding electrons (checked against
the formal charge), `valence` for a d-block atom. `"species": [...]` draws several species left
to right, each bracketed with its charge when charged (an ionic compound: `Na+ [:Cl:]−`).
`style` is `dots`, `lines` or `dot_cross` (alternate atoms in dots and crosses). `"grid": "hex"`
puts atoms on a hexagonal grid (axial `[column, row]`; neighbours at 0, 60, 120, ... degrees) for
trigonal planar centres such as CO{sub:3}{sup:2−} and zigzag chains; lone pairs take the free
directions that keep all groups furthest apart.

## unit_cell

`{"preset": "NaCl"}`: `sc`, `bcc`, `fcc`, `diamond`, `CsCl`, `NaCl`, `zinc_blende`, `fluorite`,
`perovskite`; rename species with `"elements": {"A": "K", "B": "Br"}`. A custom basis:
`"sites": [{"element": "Ti", "frac": [0.5, 0.5, 0.5]}, ...]`, one entry per site; boundary copies
are added. `bonds` joins nearest unlike neighbours; `style.radius` and `style.fill` per element.
Hexagonal presets (`hcp`, `wurtzite`, `simple_hexagonal`; `c_over_a`, default ideal 1.633) are drawn
as the textbook hexagonal prism of three primitive cells (prism corner 1/6, face 1/2), so hcp
shows 6 atoms. One-element structures also report the packing fraction (fcc and hcp 0.74, bcc
0.68, sc 0.52, diamond 0.34).

## huckel_mo

`{"system": "butadiene"}`: ethene, allyl (cation, radical, anion), butadiene, pentadienyl
(cation, anion), hexatriene, cyclopropenyl cation, cyclobutadiene, cyclopentadienyl (anion,
cation), benzene, tropylium, cyclooctatetraene. Or `{"atoms": 5, "cyclic": true,
"electrons": 6}`, and `"bonds": [[0, 1], ...]` for another topology. Chains are drawn as rows of
p orbitals (size = coefficient, the filled lobe shows the sign, a dot marks a node); rings from
above, with a circle at each atom sized by the coefficient and filled for a positive sign. The
delocalisation energy is measured against electrons/2 isolated ethene bonds.

Heteroatoms use the Streitwieser parameters (alpha{sub:X} = alpha + h beta, beta{sub:CX} = k beta):
`N_pyridine` (h 0.5, k 1.0, one pi electron), `N_pyrrole` (1.5, 0.8, two), `O_carbonyl` (1.0, 1.0, one),
`O_furan` (2.0, 0.8, two). Presets `formaldehyde`, `acrolein`, `pyridine`, `pyrimidine`, `pyrrole` and
`furan`; or add `"heteroatoms": [{"index": 0, "type": "N_pyridine"}]` (or `"symbol", "h", "k"`) to
any chain or ring. The heteroatom is labelled on the lowest orbital.

## phase_diagram

`{"triple": [T, P], "critical": [T, P], "melting_slope": "negative", "marks": [{"P": 1, "label":
"1 atm"}], "supercritical": "supercritical fluid"}`. Schematic: the triple and critical points
sit at fixed places and other values are placed piecewise between them, so a mark keeps its order
against both points but distances are not to scale. `labels` renames the regions.

## binary_phase_diagram

```json
{"type": "eutectic", "components": ["Pb", "Sn"], "melting": [327, 232], "eutectic": [0.619, 183],
 "solubility": [0.183, 0.977], "basis": "mass", "points": [{"x": 0.4, "T": 230, "label": "P"}]}
```

`type` is `eutectic` (two liquidus branches meeting at the eutectic; `solubility` gives the
limits of alpha and beta at the eutectic temperature, omit it for immiscible solids) or
`isomorphous` (a lens: one solid solution). For an isomorphous system give measured
`liquidus` and `solidus` as `[[x, T], ...]` when the numbers matter; otherwise the lens is
schematic (`lens`, 0.05 to 0.9). Each marked point is reported with its field and, in a
two-phase field, the tie line is drawn and the lever rule gives each phase's composition and
fraction. Field labels that do not fit their field go beside the frame with a leader.
Also `unit`, `basis` (`mole` or `mass`), `curvature`, `x_label`, `y_label`, `title`.

## pourbaix

```json
{"element": "Fe", "activity": 1e-6, "species": [
  {"name": "Fe{sup:2+}", "dGf": -78.90, "charge": 2}, {"name": "Fe{sup:3+}", "dGf": -4.7, "charge": 3},
  {"name": "Fe{sub:2}O{sub:3}", "dGf": -742.2, "atoms": 2, "O": 3},
  {"name": "Fe{sub:3}O{sub:4}", "dGf": -1015.4, "atoms": 3, "O": 4}]}
```

Per species: `dGf` (standard Gibbs energy of formation, kJ/mol, per formula unit), `atoms` of the
element, `O`, `H` and `charge` per formula unit, `aqueous` (default: charged species are
dissolved, at `activity`). The element itself is added when missing. Each species' free energy
relative to the element is linear in pH and E, so every field is found exactly as the region where
that species is lowest; boundaries are straight. Also `T`, `pH` and `E` ranges, `water` (the
dashed O{sub:2}/H{sub:2}O and H{sub:2}O/H{sub:2} lines) and `points` (reported with their stable
species). The diagram is only as good as the Gibbs energies: take them from one data table.

## latimer and frost

`{"species": [{"name": "MnO{sub:4}{sup:−}", "oxidation_state": 7}, ...], "potentials": [0.56, ...]}`,
highest oxidation state first, one E° per adjacent couple. `latimer` draws the chain, adds
`derived` couples as arcs with their computed E° and colours species that disproportionate (the
potential to their right exceeds the one to their left). `frost` plots N E° against oxidation
state (relative to oxidation state 0 when it is listed); `series` draws several chains, such as
acid and base, on one plot.

## ellingham

`{"T": [273, 2273], "lines": [{"label": "2C + O{sub:2} → 2CO", "dH": -221.0, "dS": 178.7}, ...]}`
with dH in kJ and dS in J/K per mole of O{sub:2}, constant between phase changes given as
`segments: [{"from_T", "dH", "dS"}]` (a phase change is a kink; a jump is refused). The report
lists every crossing and which element reduces the other oxide above it.

## correlation_diagram (Woodward–Hoffmann)

Presets `butadiene_conrotatory`, `butadiene_disrotatory`, `hexatriene_conrotatory`,
`hexatriene_disrotatory`, `cycloaddition_2+2`, `cycloaddition_4+2`; or `left` and `right` lists of
`{"label", "symmetry"}` from the lowest orbital up with `electrons` and `symmetry_element`.
Orbitals of the same symmetry are joined lowest to lowest (the non-crossing rule); links that take
an occupied orbital to an antibonding one are orange, and the verdict (thermally allowed or
forbidden) is printed under the diagram.

## walsh_diagram

`{"preset": "AH2", "electrons": 8}` (valence electrons: BeH{sub:2} 4 linear, BH{sub:2} 5 to
H{sub:2}O 8 bent) or `x`, `orbitals: [{"left", "right", "energy": [...]}]` and `electrons`. The
preset's energies are qualitative: the report's preferred geometry is the trend, not an angle.

## potential_surface

`{"preset": "mueller_brown"}` or `{"grid": {"x": [...], "y": [...], "z": [[...], ...]}}` (one row per
y value, e.g. from a relaxed scan); `levels` (count or values), `cap`, `path` (an IRC or other
path drawn on top), `x_label`, `y_label`. Minima (dots) and first-order saddle points (crosses,
TS‡) are located and reported with their energies and the barriers from the lowest minimum.

## flowchart (process and separation charts)

`nodes`: `{"id", "text", "shape": "box" | "round" | "ellipse" | "plain", "row", "col"}` (row and
column are inferred from the edges). `edges`: `{"from", "to", "label", "below", "dashed",
"route": "top" | "bottom"}`; an edge back to an earlier step becomes a return loop outside the
chart. `inputs` (`{"to", "text", "side"}`) and `outputs` (`{"from", "text", "side"}`) add reagents
and by-products as short arrows above or below a step. Charts wider than a double column (504 pt)
are reported; move later steps to a second row with `row`/`col`.

## electrochemical_cell

`{"type": "galvanic" | "electrolytic", "anode": {"electrode", "solution", "reaction"},
"cathode": {...}, "anode_side": "left", "bridge": "salt_bridge" | "membrane" | "none",
"meter": "V", "ions": {"anions": "Cl{sup:−}", "cations": "K{sup:+}"}, "labels": {...}}`.
Oxidation is always at the anode; the anode is negative in a galvanic and positive in an
electrolytic cell. Electron and ion arrows follow from that and cannot be set the wrong way.

## periodic_table

`periods`: `[1, 4]`; `f_block`; `show_number`; `highlight`: `{"elements": [...]}`, `{"group": 17}`,
`{"period": 3}` or `{"block": "d"}` with `color` and a legend `label`; `trends`:
`{"direction": "right" | "left" | "up" | "down", "label"}`.

## compose: several drawings on one page

```json
{"parts": [{"file": "profile.cdxml", "at": [0, 0]},
           {"file": "ts.cdxml", "at": [180, 20], "anchor": "bottom"}]}
```

`cdx.py scientific compose spec.json new.cdxml` merges the parts; ids, references, colours and
fonts are renumbered so every part keeps its own. `anchor` (`top-left`, `top`, `center`,
`bottom`, ...) is the point of the part placed at `at`; a part without `at` goes to the right of
the previous one. Overlapping parts are reported.

## Aromatic circles

`cdx.py circles in.cdxml new.cdxml` redraws Kekulé benzenoid rings (and five-membered
heteroaromatics) with a circle; in `cdx.py build` set `"aromatic_circles": true` on a molecule.
The circle is a depiction: a reader that ignores it sees a saturated ring, so keep the Kekulé
source for chemistry and report the circled file as a drawing.

## cdx.py check

Checks a saved CDXML against the house rules of [ACS 1996 drawings](acs1996-drawings.md):
document settings (bond length 14.4, line 0.6, bold 2.0, hash 2.5, margin 1.6, label and caption
10 pt Arial, chain angle 120), Chinese text in SimSun, text sizes (`--allow-sizes 8,6`), the
palette (black, white, #FF6600, #336699; `--palette any` or a list), the median drawn bond,
duplicate ids, references to missing objects, curves ChemDraw would drop, lone unlabelled
carbons that render blank and an XML declaration of `utf8`. `--ignore bond_length` accepts the
24 pt bonds of a Newman projection and the foreshortened bonds of the cage templates.

The check cannot see chemistry: a wrong ring size, a misfilled orbital or two isomers drawn alike
pass it. The diagram generators compute the chemistry they draw; everything else still needs the
native render inspected.
