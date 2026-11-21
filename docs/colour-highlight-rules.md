# When to colour a structure, and which colour

Round-2 work item 4. Default advice for `fragment_colors`, `highlight_varying` and
`group_colors`, from a census of the 26 published journal figures and the briefs that imitate
them.

## Measured: how much colour, and which colours

Sources: `output/round2-stats.json` (briefs + corpus PNGs) and
`output/round2-accent-colours.json` (strongly saturated pixels only, so anti-aliasing haze
around black text does not dominate).

**Colour is sparse.** Strongly saturated pixels per published figure run from **0 % to 3.7 %**
of the image; the median is about **1 %**. Three of the 26 figures
(`bjoc-articles-22-128-scheme17`, `bjoc-articles-22-128-scheme6`,
`nature-chemistry-s41557-026-02231-5`) use **no colour at all**. In the briefs, 24 of 26 declare
some colour, 14 use `fragment_colors`, 14 `group_colors`, 14 `highlight_varying`.

**Two hues do almost all the work.** Counting strongly saturated pixels across all 26 images:

| hue family | pixels | share |
| --- | --- | --- |
| blue | 70 448 | 43 % |
| red / maroon | 60 953 | 37 % |
| teal / cyan | 12 915 | 8 % |
| orange | 7 519 | 5 % |
| yellow / olive | 4 195 | 3 % |
| green | 2 725 | 2 % |
| magenta / pink, violet | 2 284 | 1 % |

And the hex values the briefs actually declare (count of uses):

`#8B1E24` dark maroon **20**, `#1B4D89` dark blue **15**, `#1A47B8` **9**, `#0B4F8A` **7**,
`#2B7BBA` **7**, `#008080` teal **6**, `#D95F02` orange **6**, `#2B4C8C` **6**,
then single digits: `#C2591A`, `#0072B2`, `#2E7D32`, `#C00000`, `#DC3220`, greys `#7F7F7F` /
`#6E6E6E`.

So the house palette is a **muted dark red and a muted dark blue**, with teal and orange as the
third and fourth accent. Nobody in this corpus uses pure `#FF0000` or `#0000FF`.

## What the colour is *for* (read off the published figures)

* **`nature-communications-s41467-022-35062-2-fig2`** (blue + red, 3.7 % chromatic): the
  **red** parts are the invariant chiral core — `CH3`, `NH2`, `COOtBu` — present unchanged in
  every product. The **blue** parts are the alkyne-derived fragment, *including* the `R` that
  varies across the scope. Colour traces **which starting material each atom came from**.
* **`nature-communications-s41467-022-32614-4-fig3`** (red + green): red is the
  aldehyde-derived aryl group, green the other partner's fragment. Same provenance idea.
* **BJOC 22-127 scheme 2** (already documented in `docs/journal-figure-parity.md`): the
  header's `R¹`/`R²` marks are red **because** the varying groups in the grid below are red.
  Colour marks **what changes**.

Two distinct purposes, and they need different keys.

## Default advice

**0. Default to black and white.** Colour is an information channel. If removing it loses the
   reader nothing, remove it. A third of the decision is whether to colour at all.

**1. Use `fragment_colors` for provenance — the commonest case.**
   When a scheme joins two or more partners, give each partner's fragment one colour and carry
   it into the products. The reader can then trace an atom from starting material to product
   without reading atom maps. This is the `nature-communications-35062` pattern.
   ```jsonc
   "fragment_colors": [
     {"smarts": "NC(=O)OC(C)(C)C", "color": "#8B1E24"},   // the carbamate partner
     {"smarts": "C(=O)O",          "color": "#1B4D89"}    // the acid partner
   ]
   ```
   Set `include_attaching_bond: true` when the bond to the partner should take its colour too.

**2. Use `highlight_varying` when one core is elaborated.**
   A single scaffold, substituted at one or more positions: colour **only** what differs and
   leave the shared core black. That is exactly what the maximum-common-substructure comparison
   computes, and it is the BJOC 22-127 pattern. If you also want the general reaction's R marks
   to match, that already happens (`paint_header_marks`).

**3. Use `group_colors` for named sub-series.**
   When the scope is split into labelled families ("Aromatic aldehydes" / "Aliphatic
   aldehydes"), one colour per family. Prefer `group_colors: "auto"` so the palette comes from
   the house style instead of a hand-picked hex.

**4. One or two accents, rarely three.** The corpus median is ~1 % coloured pixels. If a figure
   needs four colours to be legible, the layout is the problem, not the palette.

**5. Palette defaults.**
   * first accent `#8B1E24` (dark maroon), second `#1B4D89` (dark blue);
   * third `#008080` (teal), fourth `#D95F02` (orange);
   * for a two-colour figure where colour-blindness matters, prefer blue + orange
     (`#0072B2` / `#D95F02`) over red + green;
   * grey `#7F7F7F` for the subordinate case (a minor isomer, a by-product), which is how the
     toolkit already draws `secondary_of`.

**6. Do not colour the thing that carries the geometry.** A wedge, a hashed bond, an arrow or a
   stereodescriptor should stay readable on its own; colour must not compete with it. Colour a
   *fragment*, not a single stereo bond.

**7. Keep the header consistent with the grid.** If the grid's varying groups are coloured, the
   general reaction's R marks take the same colour — never a different accent for the same
   meaning on one page.

## Demos

Script: `scripts/round2_demos.py` (`colour_provenance`, `colour_varying`).

* `output/round2-demos/colour_provenance.compare.png` — left uncoloured, right with the
  `fragment_colors` block above. Measured strongly-chromatic fraction **0.0015 → 0.0098**; by
  eye the tert-butyl / O / N of the carbamate are red, the carbonyl and OH of the acid are blue,
  the ring stays black.
* `output/round2-demos/colour_varying.compare.png` — left uncoloured, right with
  `highlight_varying: true`. Measured **0.0000 → 0.0041**; by eye only the substituents that
  differ (OH, OMe, Cl) are orange, while the rings, the aldehyde and the labels stay black.

Both panels keep `readback=matched`: colour, like abbreviation, is a drawing decision and must
not alter the chemistry.
