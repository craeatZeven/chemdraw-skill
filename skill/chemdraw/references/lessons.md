# Lessons and fixes

Observed failures and the fix that worked. Read only the relevant row.

## Everyday drawing

| Symptom | Fix |
| --- | --- |
| `No module named cdxml_toolkit.mcp_runtime`, or old toolkit behaviour | A script ran on the wrong interpreter. Start everything through `cdx.py` (it re-launches in the runtime); run other scripts as `cdx.py py SCRIPT`. |
| Re-rendered PNG did not change | The renderer refuses to overwrite. `cdx.py render` replaces a PNG it rendered before, and only after a successful render; a same-named PNG it did not write is refused (render to `--out-dir NEW_DIR`, or pass `--replace` when the user wants it replaced). Other routes must use a new name. |
| Chinese text shows as boxes | The run needs a CJK font entry (SimSun, font 5). `acs_builder` puts CJK characters in SimSun and keeps Latin text in Arial automatically. |
| Abbreviation reads "C10H21" on the wrong side of the bond, or digits not subscripted | Use a `GenericNickname` node with a reversed label (H21C10) when the bond leaves to the right; put digits in subscript runs. `acs_builder` does both. |
| Ylide or resonance form "corrected" by RDKit | Write the nodes explicitly (`coords`/`label_text`, or raw fragment XML) instead of going through SMILES. |
| Condition text collides with a molecule | Leave about 25 pt between the molecule box and the text; long abbreviations widen the box. |
| Figures look inconsistent in a Word document | Insert every 300 dpi PNG at native size (px / 300 × 2.54 cm) so bond length and font size match across figures. |
| A figure "looked fine" in the log but was wrong | Open every rendered PNG and look at it before reporting. |
| An SVG preview shows objects the native render does not | An SVG or other non-native preview draws objects ChemDraw drops (a curve with hand-computed points, an unknown attribute value). Accept a figure on the native `cdx.py render` PNG only. |
| A `<curve>` vanished from the native render | ChemDraw silently drops curves whose `CurvePoints` are not whole (handle, point, handle) triples. Let a tool write them (`curved_arrow` in `cdx.py build`, the diagram generators); `cdx.py check` reports malformed curves. |
| Edited a build spec, the figure did not change | Older `cdx.py build` refused to replace any existing output, so a stale file survived next to the error. It now rebuilds its own unchanged output in place and still refuses a hand-edited or foreign file (`--replace` to force). |
| A spec parameter seemed to do nothing | Compare PNG (or CDXML) hashes before and after the change, not the spec. Unknown keys are now an error naming the keys an object reads; older builds ignored a conformer `kind` and an orbital `angle` without a word. |
| Methane, methanal or bromomethane came out blank or as a bare bond end | A one- or two-atom molecule whose carbon has no label has no skeleton to draw. `acs_builder` labels such carbons (`carbon_labels`, default `auto`); in hand-written CDXML give the carbon a label (CH{sub:3}, H{sub:2}C). `cdx.py check` reports lone unlabelled carbons. |
| A molecule from SMILES has the right connectivity but a poor layout | SMILES fixes topology, not the 2D layout. For a textbook layout give `coords` (one point per heavy atom) or a conformer template. |
| `Invalid class string` when ChemDraw starts | Some installations register only `ChemDraw_x64.Application`. The toolkit tries both ProgIDs; set `CHEMDRAW_PROGID` to force one. |
| A diagram generator refused the spec | It found a contradiction (a transition state below its neighbours, a cycle that breaks Hess's law, electrons that do not add up). Fix the numbers; see [chemistry diagrams](diagrams.md). |
| The house-rule check passed but the figure is chemically wrong | `cdx.py check` sees settings, fonts, colours and references only. Ring sizes, electron counts and isomers need a person to look at the native render. |

## Reconstruction

| Symptom | Fix |
| --- | --- |
| A C–N bridge shows as two segments and is taken for a dashed bond | It is occlusion: keep one single bond, remove `Display="Dash"`, set depth with `crossings`. Never add an intersection atom. |
| One crossing has a gap, the matching one does not | Give explicit front/back pairs to the figure composer; remove old hand-drawn cover lines. `Z` alone or a white overlay did not fix it. |
| DECIMER returns a parseable but wrong bridged structure | Check connectivity, carbonyls, N-ethyl, protecting-group abbreviations and charges before any styling. |
| MCS is null but the change lists are empty | Read the comparison status; an incomplete mapping is not "no change". |
| Native save adds bridgehead configurations | Compare the saved file with `cdx.py chemscript-info`; do not delete stereo or reverse wedges to make readers agree. |
| Native open fails with no output or `NoneType ... SaveAs` | The XML declaration said `utf8`; it must be `UTF-8`. |
| Time lost restyling | Settle structure and depth before typography; reuse unchanged crops. |

## Markush

| Symptom | Fix |
| --- | --- |
| Read-back of a Markush CDXML shows carbon for R groups, or loses the variable attachment | Readers turn generic labels into carbon and ignore `MultiAttachment`. Take chemistry from `markush.json` or the CXSMILES, not from the CDXML. |
| `m:` missing from CXSMILES | RDKit drops `m:` on write. Use `markush_compile`, which writes the CXSMILES itself. |
| NH label shows as "HN" or the H lands on the wrong side | Use an `N` atom with `h: 1`; the writer adds the H as its own caption on the side of the drawing. |
| Label flips to the other side of its bond | Set the label display (Left or Right) from the bond direction; the writer does this from the source layout. |
| `(R^C)_p` or `R^B` sits shifted from its bond end | ChemDraw centres the first character (Left) or the last (Right) on a near-horizontal single bond; `markush_compile` compensates for this. Steep bonds keep the stored position. |
| A bond drawn hugging a ring circle is lost | The detector reads it as part of the circle. Answer the review task with `missing_bonds` for the two marks. |
| A short double bond between two labels reads as single | Check `bond_orders` in review; short pairs near text are easy to miss. |

For perspective-bond disagreement and abbreviation hydrogen loss, see [stereo diagnosis](stereo-reconstruction.md). Workflows: [image reconstruction](image-visual-review.md), [figure manifest](publication-figures.md), [saved-file validation](figure-validation.md).
