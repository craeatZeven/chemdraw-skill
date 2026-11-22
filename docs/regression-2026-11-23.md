# ChemDraw regression report — 2026-11-23

The installed runtime passed **26/26 journal cases** with `readback=matched`
and rendered all 26 in ChemDraw. It also passed all four variable R-group
cases. Three journal headers still need visual fixes in the default workflow.
The adjusted figures are available for review as editable CDXML and native PNG.

## Reproduction and results

The benchmark used the briefs and runner at `109bc139b1705fd88fe32aa01ec1b91314125ebb`
and the installed venv under `LOCALAPPDATA`. Python reported version 3.10.10
and imported `cdxml-toolkit-community` 0.7.0a3 from that venv. Its
`ring_r_groups.py` matches `011f8f3` after newline normalization. The checkout
pins toolkit revision `1255a697814d`. ChemDraw Professional 26.0.0.6141
reported zero documents open before the native runs.

```powershell
$venvPython = Join-Path $env:LOCALAPPDATA 'chemdraw-skill/venv/Scripts/python.exe'
$env:CHEMDRAW_MCP_PYTHON = $venvPython
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$out = 'output/regression-2026-11-23'

& $venvPython -u scripts/run_scope_benchmark.py --in-runtime --corpus journal --render --side-by-side --output-dir "$out/journal"
& $venvPython skill/chemdraw/scripts/cdx.py py "$out/replay_ring_cases.py" "$out/r-groups-native"
& $venvPython "$out/prepare_ring_review.py"
& $venvPython "$out/prepare_journal_review.py"
& $venvPython skill/chemdraw/scripts/cdx.py py "$out/prepare_review_figures.py" "$out/reviewer"
& $venvPython skill/chemdraw/scripts/cdx.py py -m pytest skill/chemdraw/tests -q --basetemp "$out/.pytest"
& $venvPython skill/chemdraw/scripts/cdx.py doctor
Remove-Item Env:PYTHONUTF8
$env:OS = 'Windows_NT'
& $venvPython scripts/validate_distribution.py
```

Use fresh output paths when repeating the runs. The helper scripts, input
manifests, renders, JSON reports and logs are in the local output directory.

| Check | Exit | Result |
| --- | ---: | --- |
| Journal benchmark | 0 | 26 passed, 0 failed, 26 matched readbacks, 26 one-page native renders, 26 comparisons, zero lint issues and zero CDXML check errors |
| Variable R-group replay | 0 | Four native renders and four native CDXML saves with matching scaffold inventories and attachment facts. The paper replay also had `readback=matched`. |
| R-group comparisons and measurement | 0 | Four comparison PNGs and six visible connector measurements |
| Journal review copies | 0 | 26 published-versus-rendered comparison PNGs |
| Adjusted review figures | 0 | Three native renders and native CDXML saves with scaffold inventories matching the replay sources |
| Skill pytest | 1 | 184 passed, 5 subtests passed, 1 failed. `test_generated_mcp_signatures_are_current` flags the installed runtime mismatch. |
| Skill doctor | 1 | ChemDraw COM available, `catalog_current=true`, `runtime_matches_skill=false`, `signatures_current=false` |
| Portable distribution validation | 1 | The baseline already fails two policy checks. The English-only check flags the journal runner and two files in `docs/upstream-issues/`. The documentation allowlist also flags those two files. |

The journal runner emitted one drawing warning for independently depicting
compound 41 without a common core in `communications-chemistry-s42004-022-00692-6-fig2`.
The CDXML profile checks produced 132 warnings. I opened every journal panel
in `journal/review/` and compared all 26 full comparison PNGs with the previous
run. All 26 were byte-identical to that run. The panels compare published
figures with imitation briefs, which may cover different examples, yields, or
layouts. Readback compares the saved structures with the input brief.

## Variable R-group review rules

1. Trace the connector from inside the ring through the geometric midpoint of
   a double bond. Both strokes of the double bond should break to leave white
   clearance around the continuous connector.
2. Confirm a native `VariableAttachment` with the expected candidate atoms and
   reciprocal `CrossingBonds` after native save. Use distinct edges for two
   marks and retain the raised number after R.
3. Measure the visible connector along its axis on the original native PNG and
   compare it with the radius at the same scale. Record the full CDXML
   anchor-to-label distance separately because ChemDraw clips it at the label.
4. For fixed marks, inspect the superscript and neighboring atom together.
   Both glyphs and the bond segment must be legible. In a reaction header,
   verify the substrate order, standalone reagents, and product orientation
   against the reference.

I opened each of the four comparison images below. They are under
`r-groups-native/` in the output directory. The synthetic cases use benzene
and pyridine through the figure compositor. The Scheme 2 variable R¹ case
uses the brief from `011f8f3`. The standard journal brief specifies fixed R¹
at those two positions.

| Comparison | Visible length / radius | Full CDXML ratio | Native image observation |
| --- | --- | --- | --- |
| `mono.compare.png` | 113/240 px = 0.471 | 0.628 | Left horizontal midpoint crossing with aligned gaps in both double-bond strokes |
| `di.compare.png` | R¹ 113/240 px = 0.471, R² 118/240 px = 0.492 | 0.628, 0.884 | Separate left and lower-right crossings. Each connector stays continuous through both gaps. |
| `hetero.compare.png` | 115/240 px = 0.479 | 0.781 | Right horizontal crossing through both strokes, clear of the ring N |
| `paper.compare.png` | Each R¹ 89/180 px = 0.494 | 0.641, 0.641 | Both connectors leave the left edge and interrupt both strokes as in the published crop. The product crossing is slightly oblique. |

These ratios use measurements from native PNGs requested at 600 dpi, with
dark pixels sampled in a seven-pixel strip along each connector. The radius
is the mean centroid-to-vertex distance of the saved ring, mapped at the
observed scale of 16.67 px/pt and checked against the raster ring size.
The largest native CDXML midpoint error was 0.0037 pt. All six crossed
edges had bond order 2 and reciprocal
crossing references. See `r-groups-native/geometry-with-pixels.json`.

## Figures for visual review

The replay reproduces all three reported defects. I adjusted the coordinates
in local copies of the generated figure manifests and opened the resulting
native PNG comparisons. These adjustments are specific to the review figures.
The Nature adjustments add 7 pt to each affected fixed R–B bond. The Scheme 2
panel uses a 180-degree rotation for 2a,c–r, adds 4.5 pt to the fixed N–R² bond,
places S₈ as a third reactant label, and reorients the product amides.

| Figure | Observation in the replay | Observation in the adjusted figure |
| --- | --- | --- |
| Nature Communications Fig. 4, `nature-communications-s41467-023-40180-6-fig4.header.compare.png` | The left R¹ digit touches the borate B and obscures the bond. | R¹ clears B, and the diagonal R¹–B bond is visible. [Published figure](../assets/journal-corpus/nature-communications-s41467-023-40180-6-fig4.png) |
| Nature Synthesis Fig. 4, `nature-synthesis-s44160-024-00631-4-fig4.header.compare.png` | The R² digit crowds the boronic-acid B and obscures the bond. | The R² digit and B are separate, and the bond segment is visible. [Published figure](../assets/journal-corpus/nature-synthesis-s44160-024-00631-4-fig4.png) |
| BJOC Scheme 2, `bjoc-articles-22-127-scheme2.header.compare.png` and `scheme2-full.compare.png` | Fixed R² crowds N. In 2a,c–r, Cl is on the right and the carbonyl points down. S₈ appears in the conditions above the arrow. The product carbonyl points down. | Cl is on the left and R² is on the right of a legible N–R² bond. The reaction has a separate `+ S₈` before the arrow. The 3a–t header carbonyl and the grid amide carbonyls face up, with the benzothiazole at left. Both variable R¹ crossings keep their midpoint gaps. [Published figure](../assets/journal-corpus/bjoc-articles-22-127-scheme2.png) |

The review files are under `reviewer/`. Each header comparison shows the
published reference, installed-runtime replay and adjusted native render.
The Scheme 2 full comparison also shows the 14 structures in the scope
grid and the shared 3q caption. S₈ is editable text, so the scaffold inventory
comparisons cover the drawn structures. The entries and yields come from the
existing imitation brief. The two Scheme 2 variable marks kept the same
native crossing geometry and connector-to-radius ratios after adjustment.

The default rendering follow-up is to add a native-render regression for
left-facing fixed `GenericNickname` labels and account for their glyph
width when placing the terminal bond. The Scheme 2 default brief also
needs a third header reagent and source-oriented coordinates for 2a,c–r and
the product amides. The measured crossings provide a guard for those changes.
The installed runtime mismatch needs its own check with the checkout-pinned
runtime. An independent visual review of the adjusted figures is pending.
