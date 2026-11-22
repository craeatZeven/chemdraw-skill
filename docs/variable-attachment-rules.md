# Variable-Attachment Ring R-Group Rules

This document specifies the CDXML structure, geometry, crossing-break, and typography rules for variable-position ring substituents (`R¹`, `R²`, …) drawn with native ChemDraw `VariableAttachment` nodes.

## 1. Native CDXML Representation

- **Anchor node**: Use native ChemDraw `<n NodeType="VariableAttachment" Attachments="..."/>` listing the XML node IDs of the candidate ring atoms (at least two distinct non-dummy atoms on one unambiguous ring). Do **not** use `MultiAttachment` (reserved for $\eta^n$ haptic coordination to a metal center) and do **not** fake a floating graphic line.
- **Crossing bond**: Connect the `VariableAttachment` anchor (`B`) to the terminal R-group label node (`E`) with a single bond (`Order="1"`).
- **Bidirectional crossing break**: Set `CrossingBonds="<crossed_bond_id>"` on the variable-attachment bond and include the variable-attachment bond ID in `CrossingBonds` on the crossed ring double bond so ChemDraw renders a crisp white halo/gap breaking both lines of the double bond at the intersection.

## 2. Geometry & Proportion Bounds

All ratios are measured directly from CDXML node coordinates relative to the perceived ring center $C$ and mean ring vertex radius $R = \frac{1}{|V|}\sum_{v \in V} \|p_v - C\|$:

| Metric | Rule | Target (6-membered ring) | Target (5-membered ring) |
| --- | --- | --- | --- |
| Bond length ratio (`bond_len / ring_radius`) | **$0.30 \le \text{ratio} \le 0.55$** | `0.50` (horizontal), `0.54` (angled) | `0.50` (horizontal), `0.54` (angled) |
| Anchor distance from center (`dist(anchor, C) / ring_radius`) | **$\ge 0.45$** (never at ring center) | `0.60` (horizontal), `0.63` (angled) | `0.52` (horizontal), `0.54` (angled) |
| Crossed ring edge | **Must intersect a ring double bond** (`Order="2"`) | Unobstructed `C=C` edge between candidate sites | Unobstructed `C=C` edge between candidate sites |
| Crossing break | **White gap required** | Bidirectional `CrossingBonds` | Bidirectional `CrossingBonds` |

### Why the anchor sits near the double bond rather than at the ring center

In published organic chemistry schemes (e.g., *Beilstein J. Org. Chem.* 2026, 22, 127, Scheme 2), a variable-position ring substituent is drawn as a **short** bond starting just past the inner line of the ring double bond (inner line offset $\approx 0.18 \times \text{edge length}$ at standard 18% `BondSpacing`) and ending just outside the ring edge. Placing `VariableAttachment` at the ring center produces an overly long bond ($\approx 1\text{–}2\,R$) that visually crowds the ring interior.

## 3. Multi-Substituent & Heterocycle Placement

- **Multiple variable R groups on one ring (`R¹`, `R²`, …)**: Each group claims a distinct ring double bond with at least $0.85 \times \text{bond\_length}$ label separation and zero bond-to-bond intersection.
- **Heterocycles (pyridine, thiophene, furan, etc.)**: Exclude non-substitutable heteroatoms from `atoms` (e.g., `[0, 1, 2, 4, 5]` on pyridine, `[0, 1, 2, 4]` on thiophene). Candidate edge scoring strictly prioritizes double bonds (`Order="2"`) whose endpoints belong to `atoms` (preferring `C=C` over `C=N`) before clearance and horizontal preference.

## 4. Colour & Typography Discipline

- **Red label, black bond**: When highlighting varying substituents in a reaction header or standalone structure (e.g., `#C81E1E`), paint **only** the `R¹` / `R²` label text node. Keep the short crossing bond **black** (matching scaffold ink).
- **Never reverse to `¹R` or `²R`**:
  - ChemDraw reverses multi-run `GenericNickname` labels (`R` + `<s face="64">1</s>`) into `¹R` / `²R` when right-aligned on left-pointing bonds, or anchors `R` over the ring edge when forced to `LabelDisplay="Left"`.
  - For variable ring labels `R¹`–`R³`, emit a single non-chemical ISO-8859-1 superscript run (`R¹`, `R²`, `R³`) on `NodeType="Unspecified"` with `LabelDisplay` / `LabelAlignment` / `LabelJustification` matched to the outward bond direction (`Right` for left-pointing bonds, `Left` for right-pointing bonds, `Center` otherwise) and `InterpretChemically="no"`.
  - For fixed numbered R groups from SMILES (`[*:2]`), keep `NodeType="GenericNickname"` with `LabelDisplay="Left"` and `InterpretChemically="no"` on `<t>` so RDKit chemical readback preserves `[*:2]` while ChemDraw never flips `R²` to `²R`.

## 5. Verification

1. **CDXML coordinate audit**: Measure `bond_len / ring_radius`, `dist(anchor, C) / ring_radius`, double-bond intersection, and bidirectional `CrossingBonds` directly from the CDXML coordinates.
2. **Native ChemDraw render & round-trip**: Render through native ChemDraw COM (`cdxml_to_image`) and verify the crossing break white gap in the rendered PNG pixels and identical label order after `.native.cdxml` round-trip.
3. **Chemical scaffold readback**: `cdxml_toolkit.chemistry_semantics.read_fragments` strips validated terminal variable-attachment R groups cleanly so RDKit verifies the underlying scaffold and any fixed R groups.
