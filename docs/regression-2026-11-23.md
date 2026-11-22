# ChemDraw regression report, 2026-11-23

The installed runtime passed all 26 journal cases and all four variable
R-group cases. All 26 journal readbacks matched their input briefs. The
journal benchmark checks the supplied imitation briefs. The two published
figures discussed below required separate corrections for visual review.

## Commands and results

The runner and briefs were taken from
`109bc139b1705fd88fe32aa01ec1b91314125ebb`. The specified installed
venv ran Python 3.10.10 and imported `cdxml-toolkit-community` 0.7.0a3 from
its own `site-packages`. The native checks used ChemDraw Professional
26.0.0.6141. ChemDraw reported zero open documents before rendering.

From the repository root, with the local review helpers in the output directory:

```powershell
$venvPython = Join-Path $env:LOCALAPPDATA 'chemdraw-skill/venv/Scripts/python.exe'
$env:CHEMDRAW_MCP_PYTHON = $venvPython
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$out = 'output/regression-2026-11-23'
$cdx = 'skill/chemdraw/scripts/cdx.py'

& $venvPython -u scripts/run_scope_benchmark.py --in-runtime --corpus journal --render --side-by-side --output-dir "$out/journal"
& $venvPython $cdx py "$out/prepare_journal_review.py"
& $venvPython $cdx py "$out/replay_ring_cases.py" "$out/r-groups-native"
& $venvPython $cdx py "$out/prepare_ring_review.py"
& $venvPython $cdx py "$out/ground_review_chemistry.py"
& $venvPython $cdx py "$out/build_review_figures.py" "$out/review-figures"
& $venvPython $cdx py -m pytest skill/chemdraw/tests -q --basetemp "$out/.pytest"
& $venvPython $cdx doctor
Remove-Item Env:PYTHONUTF8
$env:OS = 'Windows_NT'
& $venvPython scripts/validate_distribution.py
```

Use fresh output directories for repeat runs. The local review archive
contains the figures, scripts, chemistry lookup results, and supporting checks.
The three synthetic R-group comparisons use prior ChemDraw exports as their
left panels. The paper R-group case uses a published crop and the
variable-attachment Scheme 2 brief from `011f8f3`.

| Check | Exit | Result |
| --- | ---: | --- |
| Journal benchmark | 0 | 26 passed, 0 failed, 26 matched readbacks, 26 one-page native renders, 26 comparisons, 0 lint issues, 0 CDXML check errors |
| Comparison preparation | 0 | 26 review-size comparison images |
| R-group replay | 0 | Four native renders and four native saves with matched scaffold inventories and attachment facts |
| R-group measurement | 0 | Four comparisons and six visible connector measurements |
| Compound name resolution | 0 | 15 successful OPSIN results |
| Corrected figure rendering | 0 | Two native saves and renders, each with a full comparison and header close-up |
| Skill pytest | 1 | 184 passed, 5 subtests passed, 1 failed. `test_generated_mcp_signatures_are_current` reports a signature mismatch. |
| Skill doctor | 1 | ChemDraw COM available and `catalog_current=true`. `runtime_matches_skill=false` and `signatures_current=false`. |
| Portable distribution validation | 1 | The same two policy checks failed in the clean baseline. The English-language check flags the journal runner and two files under `docs/upstream-issues/`. The repository Markdown check also flags those two files. |

The benchmark issued one drawing warning for compound 41 in
`communications-chemistry-s42004-022-00692-6-fig2` and 132 CDXML
profile warnings. All 26 full comparison PNGs were byte-identical to those
from the previous verification run with the same installed venv. I opened the
26 review-sized comparisons. Several figures depict selected examples or
substitute simplified structures, conditions or captions. The unchanged
benchmark comparisons reproduce those existing differences.

## Variable R-group visual rules

1. Follow the connector from the inside of the ring through the midpoint of
   a double bond. The connector stays continuous, with a visible white gap
   on both sides where it crosses each stroke of the ring bond.
2. After native save, verify a `VariableAttachment` with the expected
   candidate atoms and reciprocal `CrossingBonds`. Two marks must cross
   different edges. Keep each numeric label raised after R.
3. Measure the visible connector in the native PNG at its original scale.
   Report that length over the ring radius. Native label clearance shortens
   the visible line, so record the full CDXML anchor-to-label ratio separately.

I opened all four comparisons in `r-groups-native/`.

| Comparison | Visible PNG length / radius | Full CDXML ratio | What the native render shows |
| --- | ---: | ---: | --- |
| `mono.compare.png` | 113/240 px = 0.471 | 0.628 | Left horizontal connector with gaps in both vertical double-bond strokes |
| `di.compare.png` | 113/240 px = 0.471 and 118/240 px = 0.492 | 0.628 and 0.884 | Separate left and lower-right crossings, each with two clear gaps |
| `hetero.compare.png` | 115/240 px = 0.479 | 0.781 | Right horizontal crossing with two gaps and clearance from pyridine N |
| `paper.compare.png` | 89/180 px = 0.494 for both | 0.641 for both | Two left-edge R¹ marks with the same gap pattern as the paper. The product crossing is slightly oblique. |

The visible lengths were measured from ink in a seven-pixel strip along each
connector in the fresh PNGs requested at 600 dpi. The radius is the mean
CDXML centroid-to-vertex distance multiplied by the observed raster scale
of 16.67 px/pt, checked against the raster ring diameter. The largest native
CDXML midpoint error was 0.0037 pt. All six crossed edges had bond order 2 and
reciprocal crossing references. Measurements are in
`r-groups-native/geometry-with-pixels.json`.

## Corrected figures for visual review

The full comparisons and header close-ups are under `review-figures/`.
The close-ups show the published header, the benchmark replay and the corrected
ChemDraw render. I opened both full comparisons and both close-ups.

- **BJOC Scheme 2.** The corrected panel contains all 20 labelled entries
  from 3a through 3t, with 19 structure drawings because 3p and 3q share
  `NH–Alk`. The panel restores 3i (59%), 3l (58%), 3m (61%), 3n (74%) and
  3o (57%). It shows 3d as 3,5-dimethyl (57%), 3f as unsubstituted
  N-phenyl (66%) and 3k as ortho-CF₃ (63%). It also shows the corrected
  3g, 3h and 3r–3t, keeps the 0% entry for 3j, and places S₈ before the
  arrow. The rebuilt brief passed readback, and the edited figure retained
  its inventory through a native save. [Published Scheme 2 and text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13746456/)
- **Nature Synthesis Fig. 4.** The corrected reaction header depicts the
  aryl–TT/TFT⁺ label bonded to the aryl ring with the BF₄⁻ counterion shown
  explicitly. The blue thianthrene core at right has two sulfurs, four R¹
  marks and the definitions `TT, R¹ = H` and `TFT, R¹ = F`. The boronic
  acid label, reaction conditions and aryl product are legible.
  ChemDraw saved the TT/TFT node with charge +1. The generic label text and
  counterion were checked directly in CDXML, and the scaffold inventory
  matched after mapping generic labels to conventional R markers. The
  grid carries the benchmark's 12 selected products out of the 29 published
  entries. [Published Fig. 4](https://www.nature.com/articles/s44160-024-00631-4)

The benchmark corpus remains the replay input. To incorporate these review
corrections into future default outputs, update the Scheme 2 brief with the
complete compound inventory, yields and verified header layout. For Fig. 4,
use a fixture that supports the generic TT/TFT leaving-group label, its charge
and the independent definition panel. The installed toolkit signature mismatch
needs a separate run using the checkout-pinned toolkit revision. Independent
reviewer sign-off is pending.
