# ChemDraw regression report — 2026-11-23

The installed toolkit passed the journal runner on **26/26 cases**, all with
`readback=matched`. All **four variable-attachment cases** met the crossing
criteria. The six visible connector lengths were **0.471–0.494 times the ring
radius**. The journal review also found reduced clearance at fixed R-group
labels in three figures. The installed toolkit and this checkout's runtime
pin differ.

## Reproduction

The journal inputs and runner came from `109bc139b1705fd88fe32aa01ec1b91314125ebb`.
The existing regression worktree already pointed at the fetched branch tip
and was reused. The commands below use the installed ChemDraw venv under
`LOCALAPPDATA`. It reported CPython 3.10.10 and imported
`cdxml-toolkit-community` 0.7.0a3 from its own `site-packages`. The installed
`ring_r_groups.py` matches commit `011f8f3` after newline normalization.
The checkout pins toolkit revision `1255a697814d`. The venv was kept in place
to test the installed toolkit.

The local ChemDraw status reported Professional 26.0.0.6141 and zero open
documents before rendering. The render and native save operations used the
regression files.

```powershell
$venvPython = Join-Path $env:LOCALAPPDATA 'chemdraw-skill/venv/Scripts/python.exe'
$env:CHEMDRAW_MCP_PYTHON = $venvPython
$env:PYTHONIOENCODING = 'utf-8'
$out = 'output/regression-2026-11-23'

& $venvPython scripts/run_scope_benchmark.py --in-runtime --corpus journal --render --side-by-side --output-dir "$out/journal"
& $venvPython skill/chemdraw/scripts/cdx.py py "$out/replay_ring_cases.py" "$out/r-groups-native"
& $venvPython "$out/prepare_ring_review.py"

$env:PYTHONUTF8 = '1'
& $venvPython skill/chemdraw/scripts/cdx.py py -m pytest skill/chemdraw/tests -q --basetemp "$out/.pytest-final"
& $venvPython skill/chemdraw/scripts/cdx.py doctor
Remove-Item Env:PYTHONUTF8
$env:OS = 'Windows_NT'
& $venvPython scripts/validate_distribution.py
```

The two R-group helper scripts and generated evidence are local to the
worktree under `output/regression-2026-11-23/`. Use fresh output paths for
new runs. The replay uses `compose_chemical_figure` with the ACS 1996 style
and 90-degree rotation for benzene with one and two variable marks, and
pyridine with one mark. Benzene supplies all six ring atoms as candidate
sites. Pyridine supplies its five carbon atoms. The paper replay uses
`draw_scope` and reads the Scheme 2 brief from `011f8f3`, which supplies
variable R¹ on the substrate and product. The standard journal input at the
checkout tip uses fixed R¹ there. The paper replay is a separate case.

| Command | Exit | Result |
| --- | ---: | --- |
| Journal runner | 0 | 26 passed, 0 failed, 26 one-page native renders, 26 comparisons, 26 matched readbacks |
| R-group replay | 0 | 4 native PNGs and 4 native CDXML saves. Scaffold inventory and attachment facts survived all four saves. The paper scope also reported `readback=matched`. |
| R-group comparison and pixel measurement | 0 | 4 comparisons, 6 measured connectors |
| Skill pytest | 1 | 184 passed, 1 failed, 5 subtests passed. The failure was `test_generated_mcp_signatures_are_current` against the installed runtime. UTF-8 mode and an ASCII pytest scratch path resolved the earlier subprocess decoding failures. |
| Skill doctor | 1 | ChemDraw COM available, `catalog_current=true`, `runtime_matches_skill=false`, `signatures_current=false`. |
| Portable distribution validation | 1 | Two existing failures. The English-only check flags `scripts/run_scope_benchmark.py` and the two files under `docs/upstream-issues/`. The documentation allowlist also flags those two files. Windows prerequisite detection passed when the shell supplied `OS=Windows_NT`, which was absent on the first run. |

Office opening and rendering were outside this drawing regression and were not
checked.

## R-group acceptance rules and review

1. Trace the connector from inside the ring to the label. It must cross the
   geometric midpoint of a double bond. Both strokes of that double bond
   should break around the continuous connector, leaving white clearance.
2. For each mark, verify an editable `VariableAttachment`, reciprocal
   `CrossingBonds` and unchanged candidate atoms on native save. Multiple
   marks should use distinct ring edges. The raised number should follow R.
3. Measure the visible connector along its axis on the unresized native PNG
   after ChemDraw's label clipping. Compare its length with the ring radius
   at the same scale. Keep the full CDXML anchor-to-label distance as a
   separate measurement.

I opened all four comparison PNGs individually. Files in the table are
relative to `output/regression-2026-11-23/r-groups-native/`.

| Case and comparison | Visible length / radius | Full CDXML ratio | What the comparison showed |
| --- | --- | --- | --- |
| Single substitution, `mono.compare.png` | 113/240 px = **0.471** | 0.628 | Horizontal R¹ connector crosses the middle of the left double bond. Aligned gaps remain on both strokes. The fresh raster matches the prior ChemDraw export pixel for pixel. |
| Double substitution, `di.compare.png` | R¹ 113/240 px = **0.471**, R² 118/240 px = **0.492** | 0.628, 0.884 | R¹ crosses at left and R² crosses the lower-right oblique double bond. Each connector stays continuous and each pair of breaks aligns. The fresh raster matches the prior export pixel for pixel. |
| Pyridine, `hetero.compare.png` | 115/240 px = **0.479** | 0.781 | R¹ crosses the right-hand double bond with white clearance on both strokes. The connector clears the ring N. The fresh raster matches the prior export pixel for pixel. |
| Published Scheme 2, `paper.compare.png` | Substrate and product each 89/180 px = **0.494** | 0.641, 0.641 | Both R¹ marks use the paper's left-side midpoint crossing and interrupted double bond. The product crossing is slightly oblique in the native render. |

Black pixels were sampled along each connector direction in a seven-pixel
strip using a threshold of 128 per RGB channel. The radius is the average
centroid-to-ring-vertex distance measured from native-saved CDXML and mapped
to the PNG scale, verified against the raster ring size. The original paper
crop shows about 17 px of connector for a 35 px ring radius for each mark
(about 0.49). The largest midpoint error in native-saved CDXML was 0.0037 pt.
All six crossings referenced a double bond in both directions. The local
`geometry-with-pixels.json` preserves the measurements.

The full header in `paper-header.compare.png` also shows the fixed R² next
to the substrate nitrogen crowding its label. This issue is outside the
variable R¹ crossing and appears in the standard Scheme 2 run as well.
The published witness is the repository's CC BY 4.0 copy of Scheme 2,
Beilstein J. Org. Chem. 2026, 22, 127, DOI `10.3762/bjoc.22.127`.
[Published figure](../assets/journal-corpus/bjoc-articles-22-127-scheme2.png)

## Journal review

All cases finished with zero lint issues and zero check errors. The runner
reported one drawing warning. In `communications-chemistry-s42004-022-00692-6-fig2`,
compound 41 was depicted independently because it had no usable common core.
The separate CDXML checks reported 132 warnings across the set. The warning
counts and types match the previous run.

I opened all 26 published-versus-rendered comparison panels. The directory
is `output/regression-2026-11-23/journal/comparisons/`, with
`<case>.compare.png` for each entry below. Display-sized copies are in
`journal/review/`. The published figures and the imitation briefs can have
different compounds, values and panel layouts. The runner's readback checks
the brief's scaffold inventory.

| Case | Visible comparison |
| --- | --- |
| `bjoc-articles-22-103-scheme1` | Reaction header and a shorter grid with blue heterocycles and teal fluorinated groups. |
| `bjoc-articles-22-125-scheme2` | Magenta three-membered rings, with yield and trans/cis captions on separate lines. |
| `bjoc-articles-22-127-scheme2` | Red variable portions and the 0% entry are visible. The standard header uses fixed R¹ and has a crowded fixed R² beside N. |
| `bjoc-articles-22-128-scheme10` | Four selected compounds with wedges and yield/ee captions. |
| `bjoc-articles-22-128-scheme17` | Five selected examples in two rows, each with yield, ee and dr. |
| `bjoc-articles-22-128-scheme6` | Four selected examples below the reaction scheme. Header R labels carry their superscripts. |
| `communications-chemistry-s42004-022-00692-6-fig2` | Blue aziridination and amidation bands followed by three 0% entries. |
| `communications-chemistry-s42004-022-00764-7-fig2` | Flow conditions appear as arrow text. The first compound has four stacked definitions and yields. |
| `communications-chemistry-s42004-023-01030-0-fig2` | Red and blue peptide fragments, three section bands and a footnote. |
| `nature-catalysis-s41929-023-01081-5-fig2` | Red and blue bands separate the imine and acylsilane scopes. |
| `nature-catalysis-s41929-024-01153-0-fig3` | Four outlined category labels, blue groups and an e.r. caption for each selected entry. |
| `nature-catalysis-s41929-025-01357-y-fig3` | Three colored bands, red aroyl fragments and two 0% entries. |
| `nature-catalysis-s41929-025-01405-7-fig2` | Blue ethyl groups and a four-line definition/yield caption in the grid. |
| `nature-chemistry-s41557-023-01174-5-fig3` | Two helicene bands. The final digit of the energy unit for compound 2d wraps to its own line in both runs. |
| `nature-chemistry-s41557-026-02231-5-fig3` | Gray aryl groups and separate carbonyl and aryl reagent sections. |
| `nature-communications-pmc13518986-fig4` | Blue prenyl groups, a natural-product band and two low-yield limitation entries. |
| `nature-communications-pmc13594109-fig4` | Four color-coded bands and a gray secondary product next to 14a. |
| `nature-communications-s41467-022-29059-0-fig2` | Stereo wedges and multilined yield, d.r. and e.e. captions across a four-column grid. |
| `nature-communications-s41467-022-32614-4-fig3` | Three aldehyde classes with colored outlines and yields, d.r. and ee. |
| `nature-communications-s41467-022-35062-2-fig2` | Red amino ester motifs, blue alkyne motifs and the trace result for 4aa. |
| `nature-communications-s41467-023-40180-6-fig4` | Red acyl fluoride ends across the grid. The fixed R¹ superscript is tight against B in the header. |
| `nature-communications-s41467-025-62466-7-fig4` | Blue substituents and a separate thymol-derived example in the bottom row. |
| `nature-synthesis-s44160-024-00566-w-fig4` | Pyrrolidine and piperidine bands with orange labels on the heteroatom-containing motif. |
| `nature-synthesis-s44160-024-00631-4-fig4` | Orange acidic and blue neutral bands. The fixed R² superscript crowds B in the header. |
| `nature-synthesis-s44160-025-00771-1-fig4` | Three differently colored scope bands and isolated footnotes below the grid. |
| `nature-synthesis-s44160-025-00821-8-fig4` | Blue nucleophiles and two groups of cage structures with yields and source captions. |

Against the prior journal PNGs in `output/journal-full/comparisons/`, 12
comparison panels were pixel-identical. Another five had at most two pixels
with a channel difference greater than 16. In the other nine, I inspected
the changed regions in `journal/deltas/`. The changes centered on fixed
numbered R labels, generally with the superscript moving to the right of R.
The superscript crowds an adjacent atom in the three cases below. The
automated journal pass does not check this clearance.

| Figure | Follow-up |
| --- | --- |
| BJOC 22-127 Scheme 2 | Preserve clearance between fixed R² and the amide N. Verify both directions of the R label in a native render. |
| Nature Communications 10.1038/s41467-023-40180-6, Fig. 4 | Check fixed R¹ on the borate before and after native export. Verify clearance from B and continuity of the short bond. |
| Nature Synthesis 10.1038/s44160-024-00631-4, Fig. 4 | Add a fixed R²–B rendering regression using the acidic/neutral header. Keep the superscript order and make room for both labels and their bond. |

The next code change should cover fixed `GenericNickname` label anchoring
with a failing native-render clearance test. The Scheme 2 generated CDXML
diff points to the added `LabelDisplay="Left"` and
`InterpretChemically="no"` attributes on the fixed R nodes and text. A
separate checkout-pinned runtime can then be checked with doctor and the
Skill suite to distinguish release packaging from the installed renderer.
