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

For reader disagreement, hydrogen loss on abbreviations and version-specific stereo handling, see [stereo diagnosis](stereo-reconstruction.md).
