# Figure coverage and main integration

The first parent is the figure-coverage branch at `adf783d`. The second parent is
`origin/main` at `109bc13`. At the merge base, there were 109 figure-coverage
commits and 12 main commits. Eighteen files changed on both sides, and Git
flagged twelve as unmerged.

## Merge decisions

- The figure dispatcher retains strict species validation and charge styles
  alongside retrosynthesis, optimization tables and energy profiles. Panel
  validation accepts charge styles, side-by-side layout, title styles, accents
  and a table of contents. Disconnection marks use the same atom mapping as the
  species color tools. The renderer and schema keep the native retrosynthesis
  arrow alongside the existing arrow styles.
- Scope briefs combine header text, header colors and rotation with the editable
  R attachments. Scope rows share a caption baseline. A brief with `columns`
  alone gets an equal-width grid and a width diagnostic when a cell cannot fit.
  An explicit `row_align` uses variable-width rows and makes `columns` a per-row
  limit. The layout also retains semantic band colors, SMARTS fragment colors,
  borders, and formatted footnotes. Font deduplication remaps the document
  defaults as well as text runs.
- Windows prerequisites retain both COM registration discovery and the
  platform fallback when the shell omits its OS variable.
- The merge recovered `test_cdx_rules_cli.py` byte-for-byte from main. A
  regression test exercises the combined header and fixed-grid path. The
  interface inventories were regenerated from the merged implementation.
- The distribution test includes the publication figure rules, textbook figure
  requests and this report. A workstation-specific README was removed from
  this branch's root and remains in the first parent's history.

See [Publication figure rules](publication-figure-rules.md) for drawing and
validation guidance.
