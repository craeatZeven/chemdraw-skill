# Variable-attachment ring R groups

Use `ring_r_groups` on a molecule in a publication figure or in a reaction
header. Each entry supplies a positive `number`, at least two distinct
zero-based `atoms` from one unambiguous ring, and an optional `color` for
the label. The list of candidate atoms remains editable in ChemDraw.

## Midpoint and clearance

- Draw the connector perpendicular to the selected edge and through its
  geometric midpoint. Use the ring centroid to determine which direction
  is outward. Prefer an unobstructed double bond and use different edges
  for multiple marks. Keep the primary ring mark on the left where the
  reaction header can be oriented that way.
- Place the internal endpoint a short distance past the inner stroke of
  a double bond. In the default style, the inset is about 0.33 of the
  edge length for a double bond. The target visible extension outside
  the ring is 0.15 of the nominal bond length. Allow extra space for
  ChemDraw's label clipping so an oblique line still projects past the
  ring edge with a clear gap to the label.
- Create one native `VariableAttachment` node with `Attachments` set to the
  candidate atom IDs. Connect it to the label with a regular single bond
  and link the connector and selected edge with reciprocal `CrossingBonds`.
  For oblique connectors, narrow the crossed edge's margin locally. A
  double bond has one aligned interruption across its two strokes, and
  the connector and neighbouring edges stay continuous.

## Labels and readback

Apply the highlight colour to the label and use the scaffold colour for
the connector. Keep the raised number after **R**. For R¹, R² and R³,
one text run with the matching Latin-1 superscript preserves the order
in native ChemDraw. Fixed numbered R groups are written as
`GenericNickname` nodes.

Check the midpoint, direction and crossing references in CDXML. Render
with ChemDraw and examine the visible connector and clearance at both
horizontal and oblique edges. The scaffold inventory check projects
terminal variable R groups off the scaffold. Read the original CDXML
attachment nodes for the candidate sites.

The horizontal reference is the reaction header in
[Beilstein J. Org. Chem. 2026, 22, 127, Scheme 2](../assets/journal-corpus/bjoc-articles-22-127-scheme2.png).
