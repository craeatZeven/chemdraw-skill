# chemdraw-energy-profile

Turn a JSON description of stationary points into an editable ACS 1996 energy profile:
a native ChemDraw `<curve>`, energies, transition-state daggers, energy-gap arrows and
an optional structure at each point — plus a native PNG and a check report.

```bash
python <skill>/scripts/cdx.py py <skill>/scripts/energy_profile/energy_profile.py spec.json -o out.cdxml
```

Products: `out.cdxml`, `out.png`, `out.report.json`. Exit code 0 = every check passed,
1 = a check failed, 2 = bad input or environment. One run takes about 20 s.

The curves come from `scientific_workflow(operation=plot)` (`ep_plot.py`); the points are
never hand-computed, because ChemDraw silently drops `<curve>` nodes whose `CurvePoints`
do not match its own plot output.

## Minimal spec

```json
{"unit": "kJ/mol", "species": [
  {"id": "R", "energy": 0},
  {"id": "TS1", "energy": 75, "ts": true},
  {"id": "P", "energy": -40}]}
```

## Top-level keys

| Key | Default | Meaning |
| --- | --- | --- |
| `style` | `smooth` | `smooth` = textbook curve; `levels` = level bars with connectors |
| `species` | required | stationary points, see below |
| `pathways` | one path in listed order | several pathways, e.g. catalysed vs uncatalysed |
| `annotations` | none | energy-gap arrows (Ea, ΔH) |
| `unit` | none | `kJ/mol`, `kcal/mol`, `eV`; written as kJ mol⁻¹ automatically |
| `energy_symbol` | `E` | y-axis symbol, italicised automatically (`ΔG` keeps Δ upright) |
| `y_label` / `x_label` | derived / `Reaction coordinate` | full axis strings, may contain markup and Chinese |
| `y_label_style` | `paren` | `paren` → *E* (kJ mol⁻¹); `slash` → *E* / (kJ mol⁻¹) |
| `reference` | none | species id; every energy is plotted relative to it |
| `decimals`, `show_plus`, `value_brackets` | 1, false, false | number formatting; minus sign is a real − |
| `show_values` | true | print the value under each name |
| `label_layout` | `stack` | `stack` name over value; `inline` one line; `split` two sides |
| `label_face` | `bold` | `bold` / `plain` / `italic` |
| `ts_marker` | `none` | `dagger` adds a ‡ superscript to transition-state names |
| `unit_width` / `energy_height` | 52 / 150 | horizontal pt per x unit; vertical pt for the whole energy range |
| `curve_width` | 1.2 | curve line width pt |
| `bars` / `bar_width` | on for `levels` / 20 | level bars (20 pt, or 14 for `smooth`) |
| `tail` | 0.35 (`smooth`), 0 (`levels`) | flat run at each end, in x units |
| `y_ticks` / `y_tick_decimals` | none / 0 | y-axis tick values; no ticks unless asked |
| `annotation_color` | `#336699` | default colour for gap arrows and their labels |

`species[]`: `id` (required, unique), `energy` (required), `label` (defaults to `id`,
may contain markup and Chinese), `ts` (transition state — label goes above by default),
`x` (explicit position), `label_position`, `dx`/`dy`, `color`, `value`, `show_value`,
and `structure` (see below).

`pathways[]`: `nodes` (species ids, x strictly increasing), `color` (defaults
`#000000`, `#CC0000`, `#808080`…), `dashed`, `width`, `label`, `label_x`,
`label_position`, `label_dy`.

`annotations[]` (`type` is only `gap`): `from`/`to` (species ids), `label` (`{value}`
is replaced by the signed difference), `x` (arrow position: species id or number,
default = the x of `to`; off a species x it draws a dashed reference line across),
`side` (`left`/`right`), `at` (position along the arrow, 0–1), `heads`
(`both`/`to`/`from`/`none`), `ref` (reference line `both`/`from`/`to`/`none`),
`color`, `label_dx`, `label_dy`.

Markup: `{i:..}` `{b:..}` `{bi:..}` `{sub:..}` `{sup:..}`; no nesting. Chinese text is
written in SimSun, everything else in Arial.

## Structures at a stationary point

A species may carry a conformer:

```json
{"id": "A", "energy": 0, "x": 0,
 "structure": {"template": "chair", "flip": false, "substituents": {"1": {"ax": "CH3"}}}}
```

`template` is one of `chair`, `newman`, `haworth`, `fischer`, `norbornane`
(also `bicyclo_222`, `adamantane`); the remaining keys are passed to
`cdxml_toolkit.templates.conformers` unchanged. `dx`/`dy` shift the structure;
it is placed above the curve's highest point plus `STRUCTURE_GAP` (10 pt), so it never
sits on the line.

Conformer templates only — there is no "any molecule" form. For a structure that is not
one of those templates, draw it with `cdx.py build` into its own CDXML and merge it into
the profile afterwards (that is what the 图5-3 SN1 figure did).

## Checks

`out.report.json` records: a self-readback of the saved CDXML (11 checks), an optional
external ACS-spec check when `CHEMDRAW_CHECK_CDXML` (or the bundled path) resolves to a
checker, a PNG pixel check when numpy and Pillow are available, and a layout warning list
(text on a curve, text below 6 pt, drawing wider than `max_width`). Overlaps are
**reported, not moved**.

## Requirements

ChemDraw desktop, a Python with numpy + Pillow (for the pixel check; skipped and noted
when absent), and `cdxml_toolkit`. Run it through `cdx.py py` so the toolkit runtime is
used — invoking `energy_profile.py` with a bare system Python may report
"cdxml_toolkit is not importable".
