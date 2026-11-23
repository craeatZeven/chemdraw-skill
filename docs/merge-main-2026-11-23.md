# Figure coverage and main integration

The first parent is the figure-coverage branch at `adf783d`. The second parent is
`origin/main` at `109bc13`. At the merge base, there were 109 figure-coverage
commits and 12 main commits. Eighteen files changed on both sides, and Git
flagged twelve as unmerged. The merge parents and conflicts can be checked with
`git show --remerge-diff d1c0297`.

## Merge decisions

- The figure dispatcher retains strict species validation and charge styles
  alongside retrosynthesis, optimization tables and energy profiles. Panel
  validation accepts charge styles, side-by-side layout, title styles, accents
  and a table of contents. Disconnection marks use the same atom mapping as the
  species color tools. The renderer and schema keep the native retrosynthesis
  arrow alongside the existing arrow styles.
- Optimization tables use the shared reaction header at the document bond
  length. Scope figures retain their larger header geometry. The generic R
  label test checks the literal superscript that ChemDraw renders after the R.
- Scope briefs combine header text, header colors and rotation with the editable
  R attachments. Scope rows share a caption baseline. A brief with `columns`
  alone gets an equal-width grid and a width diagnostic when a cell cannot fit.
  An explicit `row_align` uses variable-width rows and makes `columns` a per-row
  limit. Use `row_align: left` for the earlier main branch's left-aligned,
  variable-width layout when a brief already supplies `columns`. The layout
  also retains semantic band colors, SMARTS fragment colors, borders, and
  formatted footnotes. Font deduplication remaps the document defaults as well
  as text runs.
- Windows prerequisites retain both COM registration discovery and the
  platform fallback when the shell omits its OS variable.
- The merge recovered `test_cdx_rules_cli.py` byte-for-byte from main. A
  regression test exercises the combined header and fixed-grid path. The
  interface inventories were regenerated from the merged implementation. The
  CLI schema documents the two column modes.
- The distribution test includes the publication figure rules, textbook figure
  requests and this report. A workstation-specific README was removed from
  this branch's root and remains in the first parent's history. The journal
  runner and installation guidance from the incoming issue drafts was
  consolidated in the publication figure rules.

## Runtime pin

| Field | Value |
| --- | --- |
| Merge commit | `d1c0297bf38ba012bbd795f19a1e8db9af23e22a` |
| Toolkit source commit in `runtime.json` | `213a6ae6beb249ba940a553fd3a15925f56b453f` |
| Toolkit version | `0.7.0a3` |
| Toolkit fingerprint | `38982671523f58f6ad0e2c851bd2fd29cedb89208c874e8f2804681f17c91021` |

The runtime wheel was built from this checkout and installed into an ignored
directory. The designated Windows virtual environment's Python 3.10.10 loaded
that wheel for the journal run. The Toolkit test suite imported the checkout
with the same interpreter. The shared virtual environment was preserved.

## Verification

| Check | Result |
| --- | --- |
| `python -m pytest toolkit/tests -q` | 2,829 passed, 1 failed, 108 skipped, 5 warnings, 29 subtests passed |
| `python -m pytest tests -q` from `toolkit/`, serial rerun | 2,830 passed, 108 skipped, 5 warnings, 29 subtests passed |
| `python -m pytest -m "not network" -q` from `toolkit/` | 2,827 passed, 108 skipped, 3 deselected, 5 warnings, 29 subtests passed |
| `python -m pytest skill/chemdraw/tests -q` | 185 passed, 5 subtests passed |
| `python -m pytest tests -q` | 202 passed, 4 skipped, 1,325 subtests passed |
| `python scripts/validate_distribution.py` | Passed, 206 tests run and 4 skipped |
| `python scripts/run_scope_benchmark.py --corpus journal --no-render` | 26/26 passed |
| `python scripts/run_scope_benchmark.py --corpus journal --render --side-by-side` | 26/26 passed, 26 matched readbacks, 30 native pages |
| `cdx.py doctor --native --office` | Runtime and signatures matched. Native PNG and Office OLE probes passed. |
| Office open and render | Word and PowerPoint opened their generated documents and exported a PDF from each. |

The repository suite skipped three font checks because fontTools with Brotli
was unavailable and one browser check because `SITE_BROWSER` was unset.
The toolkit skips include platform-specific tests, optional engines and local
sample files. The pytest logs record each declared skip.
The ChemScript probe reported that its Python interpreter was unavailable.
The native PNG and Office checks used a temporary path with ASCII characters.
The first PPTX smoke attempt returned `Embedding failed validation:
(-2147221008, 'CoInitialize has not been called', None, None)`. The retry
initialized COM, then opened the PPTX and exported it to PDF.

The local `output/merge-r19-failures.txt` records each test failure and the
first error line from the merge and verification runs. Initial failures
covered scope width, default fonts, the shared table header return value,
the superscript assertion, and distribution rules. The first verification
rerun also inherited offline mode and a UTF-8 override. Offline mode blocked
tests that exercise mocked network calls. The UTF-8 override made the
installer test read Windows subprocess output with the wrong encoding.
The corrected runs used the test suite's default network and encoding settings.

One full Toolkit run reported `test_a_lasting_403_still_ends_the_range_read`
with `AssertionError: assert 4 == (8 + 1)`. The local server received four
requests while the test expected nine. The test passed in isolation. The merge
leaves both the test and the HTTP retry implementation from its first parent
intact. The subsequent serial run passed the entire Toolkit suite, including
the three network-marked tests. No test failures remain on the serial run.
The failing and passing logs are `output/merge-r19-verified-toolkit.log` and
`output/merge-r19-serial-toolkit.log`.

The native benchmark reported seven depiction warnings and 157 style warnings,
including one text intersection. The fresh run reproduced the case results
and all 30 native image dimensions. The largest mean channel difference from
the reviewed renders was 0.000964 on a 0–255 scale. The fresh report is
`output/merge-r19-verified-journal/benchmark-report.json`.

## Journal visual review

The runner's comparison PNGs are under
`output/merge-r19-journal/comparisons/`. Append `.compare.png` to a case name for
its first page. The four cases with a continuation also have a
`.page-002.compare.png`. All 30 comparisons were opened and reviewed. The
left side shows the publication and the right side shows the native ChemDraw
render of the supplied brief. A matched readback checks the render against
its brief. The review below records differences against the source publication.

| Case | Visual observation |
| --- | --- |
| `bjoc-articles-22-103-scheme1` | The header and blue and teal highlights are visible. The render shows ten examples and uses different conditions and several different yields. |
| `bjoc-articles-22-125-scheme2` | The render shows one header product and draws the ratio on its own line. The published header shows both cis and trans products, and the illustrated scaffolds differ. The ester abbreviation overlaps the ring oxygen on several rendered entries. |
| `bjoc-articles-22-127-scheme2` | The red R labels and shared 3p/3q captions render. The generated rows have three drawings where the source uses four, and the scaffold orientation is reversed. |
| `bjoc-articles-22-128-scheme10` | Four examples carry yield and ee on one line under a group rule. The source catalyst comparison and part B are absent. |
| `bjoc-articles-22-128-scheme17` | The 3+2 layout and yield, ee and dr captions are readable. The prepared structures show a different scaffold from the source's bridged ring. |
| `bjoc-articles-22-128-scheme6` | Examples 17–20 render in a 3+1 arrangement. The Ph label on 19 touches a ring bond. The source shows all four examples together and also includes a mechanism panel. |
| `communications-chemistry-s42004-022-00692-6-fig2` | The aziridination and amidation rules and motifs are blue. Page 2 contains the three unsuccessful substrates under a grey rule. |
| `communications-chemistry-s42004-022-00764-7-fig2` | The 3a–3d captions share a drawing. The reactor diagram is represented by arrow text, and some source examples are absent. |
| `communications-chemistry-s42004-023-01030-0-fig2` | The red and blue peptide fragments and footnote are visible. The source HOPO inset and fragment-coupling panel are absent. |
| `nature-catalysis-s41929-023-01081-5-fig2` | The imine and acylsilane sections have colored rules and the footnotes render. The source has more examples and a separate application panel. |
| `nature-catalysis-s41929-024-01153-0-fig3` | Four pill labels and the rounded border are visible. The crystal inset and many scope examples are absent. |
| `nature-catalysis-s41929-025-01357-y-fig3` | Three colored scope bands are visible. The photocatalyst drawing, crystal insets and nested bands are absent. |
| `nature-catalysis-s41929-025-01405-7-fig2` | Red and blue fragments and shared captions render. The X-ray inset is absent, and the chart carries a subset of source examples. |
| `nature-chemistry-s41557-023-01174-5-fig3` | Red [4]helicenes and teal [5]helicenes appear on separate pages. A barrier unit wraps, and the source's [6]helicenes and X-ray views are absent. |
| `nature-chemistry-s41557-026-02231-5-fig3` | Carbonyl and aryl-reagent bands render. Page 2 has a continued heading, one 5bq entry and the two footnotes. Source panel b is absent. |
| `nature-communications-pmc13518986-fig4` | Blue side chains and the grey limitations band are visible. The plus sign in the header sits against the phenol OH label. The selectivity model is absent. |
| `nature-communications-pmc13594109-fig4` | Colored group rules, italic locants, and grey secondary product 14b are visible. The source includes an enzyme reaction header, additional examples, and a different scaffold. |
| `nature-communications-s41467-022-29059-0-fig2` | Three rows retain the yield and stereoselectivity lines. The source uses orange and teal markers at the stereocenters, which the render omits. |
| `nature-communications-s41467-022-32614-4-fig3` | The three differently colored pill labels distinguish the aldehyde categories. The source has denser grids and additional structures. |
| `nature-communications-s41467-022-35062-2-fig2` | Red and blue fragments and the border render. The continuation contains 4u and 4aa and repeats the footnotes. |
| `nature-communications-s41467-023-40180-6-fig4` | The red acyl-fluoride groups and paired yield values render. The asterisks around the 2ag note are visible as literal characters. The source has a second reaction arrow and more entries. |
| `nature-communications-s41467-025-62466-7-fig4` | The border and blue substituents render. The source's green ring fills and yellow highlights are absent. |
| `nature-synthesis-s44160-024-00566-w-fig4` | The orange and blue scope headings divide the two rows of pyrrolidines and the piperidine section. The source's X-ray inset and second header are absent. |
| `nature-synthesis-s44160-024-00631-4-fig4` | Acidic and neutral sections render with orange and blue pill labels. The source's separate TT/TFT motif is absent. |
| `nature-synthesis-s44160-025-00771-1-fig4` | The three bands have distinct colors and the two footnotes are legible. The separate cholesteryl coupling scheme is absent. |
| `nature-synthesis-s44160-025-00821-8-fig4` | Blue nucleophile fragments, two labeled bands and both footnotes render. The narrow left column in the source becomes a horizontal row. |

See [Publication figure rules](publication-figure-rules.md) for drawing and
validation guidance.
