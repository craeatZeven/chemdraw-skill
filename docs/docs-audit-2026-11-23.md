# ChemDraw Skill documentation audit, 2026-11-23

The audit covered `skill/chemdraw/SKILL.md` and the 34 top-level Markdown files under `skill/chemdraw/references/`. It checked tool names, arguments, payload fields, CLI commands and flags, and local paths against the bundled toolkit and Skill scripts. The full registry contains 69 tools, all with generated signatures. The four tools that previously lacked a workflow reference now have one. All 28 `cdx.py` subcommands have documentation, including the three maintainer and diagnostic commands.

## Corrections

| Reference | Discrepancy | Correction |
| --- | --- | --- |
| [Scheme drawing](../skill/chemdraw/references/scheme-drawing.md) | Species records do not implement `dative_bonds` or the illustrated `donor_atom_idx` and `acceptor_atom_idx` pair. | Use grounded dative SMILES and check the written CDXML and readback. |
| [Chemistry reports](../skill/chemdraw/references/chemistry-reports.md) | `molecule_report` has no `include_safety` option. | Point to `si_materials_tables` with `reagent_safety` for reagent reminders. |
| [Chemistry reports](../skill/chemdraw/references/chemistry-reports.md) | `correlate_lcms_files` is neither registered nor implemented. `isotope_verification`, `expected_m2_ratio` and `observed_m2_ratio` are absent from the documented output. | Route normalization to `analyze_lcms_series`. The lower-level species match can carry an `isotope_verified` note. `molecule_report` calculates its isotope pattern from the input structure. |
| [Scientific workflows](../skill/chemdraw/references/scientific-workflows.md) | The NMR output uses neither `is_solvent` nor `solvent_assignment` on peaks. | Describe `solvent_peaks`, `reference_ppm` and the `exclude_solvent` option. |
| [Scientific workflows](../skill/chemdraw/references/scientific-workflows.md) | `j_couplings_hz`, `dd`, `dt` and `td` are absent from multiplet analysis. The standalone SVG and PNG plots have no zoom inset. | Document `frequency_mhz`, `j_hz`, the implemented labels and the inset in editable `plot.cdxml`. |
| [Chemistry resolution](../skill/chemdraw/references/toolkit-chemistry-resolution-interfaces.md) | `naming.reagent_database` does not exist. | Name `resolve.reagent_db` and `get_reagent_db`. |
| [Tool decisions](../skill/chemdraw/references/toolkit-tools.md) and [image interfaces](../skill/chemdraw/references/toolkit-perception-image-interfaces.md) | `reaction_image_to_cdxml` appeared as an unregistered API. Its only code mention is a test asserting its absence. | Direct readers to candidate review and reaction-role assignment in the image workflow. |
| [Tool decisions](../skill/chemdraw/references/toolkit-tools.md) and [workflow router](../skill/chemdraw/references/workflow-router.md) | Render guidance conflated the CLI's per-file `render_to_png` loop with the `render_cdxml_files` batch worker. | State which route uses one session, supports SVG and stops on the first failure. |
| [Workflow router](../skill/chemdraw/references/workflow-router.md) | The Office replacement step said to update objects in place. | Specify a new Office copy and replacement paths inside the manifest directory. |
| [Native laboratory drawings](../skill/chemdraw/references/native-laboratory-drawings.md) | `templates/glassware.cdxml` looked like a bundled path. | Identify it as the example path for a converted local template. |
| [Toolkit CLI](../skill/chemdraw/references/toolkit-cli-interfaces.md) | The guide directed agents to run installed console scripts. | Route module and helper commands through `cdx.py py` and point toolkit developers to the console script declarations. |
| [Agent integration](../skill/chemdraw/references/agent-integration.md) and [scheme drawing](../skill/chemdraw/references/scheme-drawing.md) | The style resource used `{name}` where the registered template uses `{preset}`. The reference and template resources and four MCP prompts lacked a Skill reference. | Use the registered URI templates and route readers to the client's prompt list. The registered prompts are `draw_reaction_scheme`, `reconstruct_figure`, `analyze_nmr_and_lcms` and `build_acs1996_figure`. |

## Workflow coverage added

| Tool | Guidance |
| --- | --- |
| `live_chemdraw_session` | Start with `status`. The router now explains selection fallback on Windows, new-document insertion, and the save side effects of reading a live document on macOS. [Workflow router](../skill/chemdraw/references/workflow-router.md#work-with-an-open-chemdraw-window) |
| `get_toolkit_capabilities` | Use the structured runtime and registry inventory in programmatic checks. `detail: "summary"` omits the per-tool listing. [Operations](../skill/chemdraw/references/operations.md#health-check) |
| `list_engines` | Query installed engines through the tool registry. [Prediction engines](../skill/chemdraw/references/engines.md) |
| `compare_nmr_shifts` | Compare measured shifts with a supplied prediction or local candidate groups. Report whether the result comes from a prediction or an RDKit topology estimate. [Prediction engines](../skill/chemdraw/references/engines.md#predictions) |

## Rules for future edits

1. Audit in both directions. Search the Skill entry page and top-level references for each changed name with `rg`. Resolve tool names against the bundled `build_registry("full")` in `cdxml_toolkit.mcp_runtime.tool_registry`. Compare all registered names with the curated guides. Put each missing use case in its owning guide and route it from the entry page or workflow router.
2. Verify parameters at the layer that consumes them. Use `cdx.py sig TOOL` for tool arguments, `cdx.py COMMAND --help` for CLI flags, and the receiving schema or validator for nested JSON fields. For commands that forward arguments, check the delegated parser. Verify output keys where the result is built. Finding a word elsewhere in the repository does not establish that a specific call accepts it.
3. Resolve Markdown links and anchors from their source pages. Check literal repository paths against the tracked tree and packaged paths against the deployable Skill. Mark paths that a user supplies or a workflow creates. Check resource slugs against their registered templates.
4. Keep exact callable signatures in the generated reference. Use the generators when a registry or catalog changes. Run portable distribution validation and the Skill suite with the bundled toolkit importable, and record installed-runtime drift separately if `cdx.py doctor` reports it.
5. When a change needs native visual verification, use a fresh output directory. Open each comparison image and record what is visible. A journal benchmark pass establishes that the emitted file satisfies the checks against its input brief. Fidelity to the published figure requires inspection of the two panels.

## Verification

| Check | Result |
| --- | --- |
| Generated tool signatures, public interface inventory and native interface catalog | Match the bundled toolkit. |
| Skill suite with checkout `toolkit/` on `PYTHONPATH` | 185 passed, with 5 subtests passed. `cdx.py doctor` also passed. |
| Skill suite using the installed toolkit | 184 passed, with 5 subtests passed. The generated-signature check failed. `cdx.py doctor` reports `runtime_matches_skill: false` and `signatures_current: false`. |
| Portable distribution validation | Two policy checks failed in the clean baseline and in the audit checkout. The language check flags the benchmark runner and both notes under `docs/upstream-issues/`. The Markdown allowlist flags those two notes and the requested dated audit report. Updating the test allowlist or the runner requires a code change. |

The macOS live-session and Office replacement guidance was verified against source. The native visual check below covers the Windows journal benchmark.

The journal benchmark ran with the installed ChemDraw Skill venv (Python 3.10.10 and the installed `cdxml_toolkit`) and ChemDraw Professional 26.0.0.6141. ChemDraw reported zero open documents before the native run. The check-only and native runs each passed 26 of 26 cases with matched readback, zero lint issues and zero CDXML check errors. The native run rendered one page per case and reported one drawing warning and 132 profile warnings. Compound 41 in the Communications Chemistry 00692 comparison was depicted independently because it had no usable common core.

```powershell
$venvPython = Join-Path $env:LOCALAPPDATA 'chemdraw-skill/venv/Scripts/python.exe'
& $venvPython scripts/run_scope_benchmark.py --corpus journal --no-render --output-dir output/docs-audit-2026-11-23/check-only
& $venvPython scripts/run_scope_benchmark.py --corpus journal --render --side-by-side --output-dir output/docs-audit-2026-11-23/journal
```

The full comparisons and `benchmark-report.json` are in the local `journal` output directory. Twelve full comparisons match the committed `output/journal-full/comparisons/` files byte for byte. I opened all 26 resized review copies, each named `<case>.compare.png` in the local `output/docs-audit-2026-11-23/review/` directory. The published figure is on the left and the fresh native render is on the right.

| Case | What the comparison shows |
| --- | --- |
| `bjoc-articles-22-103-scheme1` | Ten rendered entries versus twenty published. The header conditions and several caption values differ. |
| `bjoc-articles-22-125-scheme2` | One product in the rendered header versus separate trans and cis products in the paper. Some orange ester labels touch neighboring oxygen labels. |
| `bjoc-articles-22-127-scheme2` | The 3j zero yield and two-line Alk caption are visible. Entries 3i and 3l–3o are absent, and the header puts sulfur in the conditions. |
| `bjoc-articles-22-128-scheme10` | The render has the selected-examples band. The published ligand and quadrant panel is absent, and the reaction conditions and several yields differ. |
| `bjoc-articles-22-128-scheme17` | Both panels show yield, ee and dr. The render uses a monocyclic scaffold where the original shows a bridged ring. |
| `bjoc-articles-22-128-scheme6` | The four selected captions are present. The published initial reaction, second generic product and DFT mechanism are absent. |
| `communications-chemistry-s42004-022-00692-6-fig2` | Blue aziridination and amidation sections and three zero-yield captions are visible. The rendered header combines the two experimental paths. |
| `communications-chemistry-s42004-022-00764-7-fig2` | The four-line 3a–3d caption is legible. The microreactor schematic is reduced to text above the arrow. |
| `communications-chemistry-s42004-023-01030-0-fig2` | Red and blue peptide fragments and three section bands are visible. The catalyst inset and fragment-coupling panel are absent. |
| `nature-catalysis-s41929-023-01081-5-fig2` | Imine and acylsilane scopes are shown in colored bands. The lower antibiotic-scaffold example is absent. |
| `nature-catalysis-s41929-024-01153-0-fig3` | The render has four boxed group labels and a single menthol-derived example in the last group. The published X-ray inset is absent. |
| `nature-catalysis-s41929-025-01357-y-fig3` | Three scope groups and the zero yields for 15 and 58 are visible. The catalyst inset and most subgroups are absent. |
| `nature-catalysis-s41929-025-01405-7-fig2` | A compact grid preserves red and blue distinctions and a grouped X caption. The purple energy-transfer icon and X-ray inset are absent. |
| `nature-chemistry-s41557-023-01174-5-fig3` | The render includes [4]- and [5]helicene bands. The [6]helicenes, X-ray structures and panel b are absent. |
| `nature-chemistry-s41557-026-02231-5-fig3` | Two scope sections and footnotes are visible. The difunctionalized-product panel b is absent. |
| `nature-communications-pmc13518986-fig4` | Blue prenyl fragments and a gray limitations section are visible. The site-selectivity model inset and other alcohol panels are absent. |
| `nature-communications-pmc13594109-fig4` | The four section colors and a gray secondary product are visible. The reaction header is absent, and the depicted entries and yields differ. |
| `nature-communications-s41467-022-29059-0-fig2` | Twelve labeled entries are laid out 4 × 3. The paper's colored stereocenter dots are absent, and the structures are oriented differently. |
| `nature-communications-s41467-022-32614-4-fig3` | Three aldehyde groups and yield/dr/ee captions are visible. The render uses bands in place of the published boxes and includes fewer entries. |
| `nature-communications-s41467-022-35062-2-fig2` | Red amino fragments, blue aryl groups and the 4aa trace caption are visible. The published multi-entry alkyl caption is absent. |
| `nature-communications-s41467-023-40180-6-fig4` | Red acyl fluoride highlights and paired yield values are visible. The follow-on amidation arrow is absent. |
| `nature-communications-s41467-025-62466-7-fig4` | The boxed scope shows blue variable groups and the thymol example. Yellow highlights, the green ring fill and the ligand inset are absent. |
| `nature-synthesis-s44160-024-00566-w-fig4` | Orange and blue bands divide two ring-size groups. The second reaction header and crystal-structure inset are absent. |
| `nature-synthesis-s44160-024-00631-4-fig4` | The acidic and neutral bands are visible. The header renders an aryl iodide labeled as an Ar–TT salt, and the TT/TFT definition panel is absent. |
| `nature-synthesis-s44160-025-00771-1-fig4` | The arylsodium, secondary-halide and primary-bromide groups are visible. The bottom cholesteryl coupling sequence is absent. |
| `nature-synthesis-s44160-025-00821-8-fig4` | Two scope bands are visible. The starting TT/triflate depiction is replaced by generic R labels, and some of the later products are absent. |
