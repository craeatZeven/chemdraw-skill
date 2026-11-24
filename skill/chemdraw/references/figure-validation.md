# Saved-file figure validation

Use for the final CDXML, and always for bridges, charge symbols and stereo-sensitive edits. Native saving goes through the locked worker.

```bash
cdx.py validate FILE.cdxml [--output-dir NEW_DIR] [--native] [--cross-reader]
```

`--output-dir` creates a new diagnostic folder holding `report.json` and any saved native files (default: a fresh folder in temp). `--cross-reader` requires `--native`.

- **Basic:** XML references, molecule inventory, connectivity, formal charge, specified stereo.
- **Native:** CDXML → CDX → CDXML, then the inventories are compared (count delta plus `added_molecules` and `removed_molecules` flag phantom fragments spawned on save). The source is untouched.
- **Cross-reader:** saved fragments read directly by ChemScript and RDKit (no MOL step).

`ok=true` only means the report was written. Read `native.status`, `cross_reader.status` and `source_identity`; the values `not_run`, `failed`, `unavailable`, `unresolved`, `changed` and `preserved` mean different things. The report does not cover appearance: render with `cdx.py render` and look at the PNG.

For structural corrections, `modify_molecule.diff` gives `status`, `atom_changes`, `bond_changes` and `stereo_changes`. Only `completed` supports a "no change" statement; `partial`, `not_completed` or a null mapping do not.

## Render path

Round-trip validation proves the node exists in the XML; it does not prove ChemDraw draws it. Keep the two apart:

| Symptom | Cause | Fix |
| --- | --- | --- |
| Element present in the CDXML, absent from the native PNG | Decorative `<graphic GraphicType="Oval">` is not drawn; a text glyph used as a mark (`·` in Arial renders as a tofu square, SimSun's `·` grows without a stable size) | Put the mark on a real object: radical as `<n Radical="doublet">`, charge as `<graphic GraphicType="Symbol" SymbolType="CirclePlus|CircleMinus">` |
| A `<curve>` disappears on save | `CurvePoints` is not whole (handle, point, handle) triples | Let a tool write it (`curved_arrow`, the diagram generators); `cdx.py check` reports malformed curves |
| Print-size PNG cannot be re-exported | The export path and the resolution setting are different controls | Render through ChemDraw with `Preferences.PNGResolution` + `Document.SaveAs("x.png")` (300 / 600 honoured); an MCP `export_image(dpi=…)` may stay at a fixed size. `SaveAs` returns nothing, without error, for unsupported formats such as `pdf` |

Accept a figure only when the delivered PNG is the native render of the delivered CDXML. If a raster was patched by hand, record `render_source: patch` and mirror the repair in the CDXML, otherwise the next render silently reverts it. Never mix `chemdraw`, `patch` and `builtin` sources inside one figure set.

## Mark geometry

Marks and arrows are chemistry, not decoration, so check numbers rather than appearance:

- **Clearance**: a free-standing radical or charge marker keeps ≥ `radius + 0.6` units to every bond segment and ≥ `radius + 1.5` units to every caption box. Exclude large ovals (orbital lobes) or every orbital figure becomes a false positive.
- **Placement**: prefer directly above the atom. The conventional upper-right superscript position is frequently occupied by a C–H bond; if it is, move the mark rather than shortening the bond.
- **Arrow head**: within 2 units of the intended atom or bond midpoint, tail at the electron source (a lone pair or a π bond), every arc carrying `Center3D`.
- **Rotation sense**: derive conrotatory / disrotatory from the arc geometry; a rotation arrow built from `<curve>` carries no arrowhead semantics at all.

Repair rasters only inside the target radius (paint white, then redraw). A full-region paste can erase an atom glyph and two bonds without any error.

For reader disagreement, hydrogen loss on abbreviations and version-specific stereo handling, see [stereo diagnosis](stereo-reconstruction.md).
