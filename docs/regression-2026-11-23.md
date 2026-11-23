# ChemDraw R-group visual regression

The visual gate redraws four variable-attachment fixtures through the toolkit
and exports each with native ChemDraw at 600 dpi. The fixtures cover a single
benzene R group, two benzene R groups, a pyridine R group, and the BJOC Scheme 2
scope with two variable attachments in the reaction header. The paper fixture
uses the [journal brief](../assets/scope-benchmark/briefs/bjoc-articles-22-127-scheme2.json)
with the two header sites expressed as variable attachments.

## Run the gate

Run from the repository root with licensed ChemDraw for Windows and the toolkit
runtime venv with Pillow and `ring_r_groups` support. The goldens were recorded
with ChemDraw 26.0.0.6141, Pillow 12.3.0 and toolkit 0.7.0a3 with the
variable-attachment implementation from
[`011f8f3`](https://github.com/craeatZeven/chemdraw-skill/commit/011f8f3968aa86e945b12ae901302bd8d8c8aeb8).
The checkout's default toolkit pin predates that implementation. To test edits
to the drawing code, install the candidate toolkit into the selected venv
before running the gate.

```powershell
$venvPython = Join-Path $env:LOCALAPPDATA 'chemdraw-skill/venv/Scripts/python.exe'
$env:CHEMDRAW_MCP_PYTHON = $venvPython
& $venvPython skill/chemdraw/scripts/cdx.py py -m tests.check_visual
```

The command prints each case's before and after dimensions and hash distance.
It exits 1 if any case fails. Fresh PNGs, native CDXML and fingerprint JSON
are written to the displayed run directory under `output/visual-guard/`.
Use `--output-dir` to select a new directory. To create side-by-side PNGs,
pass `--reference-dir` pointing to a retained baseline run with
`mono.png`, `di.png`, `hetero.png` and `paper.png`.

## Pass criteria

- Width and height may each change by up to 2% of the corresponding golden
  dimension.
- The RGB perceptual hash may differ by up to 16 bits. Pillow composites
  each image over white and downsamples it to a 32 × 32 grid with BOX
  resampling. Each channel contributes a bit at each of three intensity
  cutoffs (192, 232, 250). The Hamming distance counts changed bits, including
  channel changes that a grayscale hash can miss. The paper case also hashes
  the top seventh of the image, where the reaction header sits, on a 128 × 32
  grid with the same 16-bit limit. This keeps a local change to an R label
  detectable on a full scope page.
- Every native variable attachment must have a reciprocal crossing bond.
  The expected marks and labels must be present. The ring center, the
  anchor, the label-side endpoint, and the crossed edge midpoint may each
  move by up to 0.10 pt in native CDXML page coordinates. The ring center is
  the mean of the vertices in the shortest ring containing the crossed edge. The anchor
  and label-side endpoint are the two ends of the short bond. The renderer
  may shorten the visible line at the label.

The files in [`tests/golden/`](../tests/golden/) store the PNG dimensions, hash
and landmarks. The PNG measurements depend on the native renderer and fonts.
Use the baseline environment for comparison, and inspect the native PNGs when
the gate flags a change.

## Update the goldens

After reviewing an intentional visual change, run the gate with
`--update-golden --output-dir output/visual-guard/new-baseline`. It writes the
four JSON files only after all cases export successfully. Keep that output
directory if you want PNGs for the next comparison. Review the JSON diff and
run the gate again before committing the new goldens.
