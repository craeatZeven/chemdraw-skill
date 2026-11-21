# R-group attachment geometry: vertex bond vs ring-centre crossing bond

Round-2 work item 1. This documents how an `R` / `R¹` / `[*:n]` label is bonded to a ring,
why there are two legitimate answers, and the exact geometry the toolkit now writes.

## The finding that decides the design: the corpus uses BOTH conventions

Reading the 26 published figures at 5x magnification shows that "how an R is bonded to a
ring" is not one rule but two, and the difference carries meaning:

| Case in the corpus | How it is drawn | What it means |
| --- | --- | --- |
| BJOC 22-103 scheme 1, **1a–c**, `R¹` on the benzene | bond starts **on a ring vertex** | the position is fixed (para to the CH₂CN); only the group varies (`1a R¹ = H`, `1b R¹ = OMe`, `1c R¹ = Cl`) |
| BJOC 22-103 scheme 1, **2a–l**, `R³` on the chromenone's benzo ring | bond starts **inside the ring**, exits **through the middle of a C=C ring bond** | the position is *not* fixed; R³ may sit at any free position of that ring |
| `assets/readme/markush/markush.cdxml`, `(R^B)_q` | `NodeType="MultiAttachment"` on the ring centroid | same convention, in ChemDraw's own file format |

So the "环心穿键" drawing is the **unspecified-position** convention, not a universal
improvement. Applying it everywhere would misstate a fixed-position scope. It is therefore an
explicit choice in the brief, and the default is unchanged.

## The authentic representation (decoded from `markush.cdxml`)

ChemDraw's own file for this carries:

```
<n id="186" p="109.80 143.15" NodeType="MultiAttachment" Attachments="107 116 108"/>
<n id="187" p="121.04 163.46" NodeType="GenericNickname" GenericNickname="RB" .../>
<b id="188" B="186" E="187"/>
```

Measured against that ring (centroid computed from its six atoms):

* the MultiAttachment node sits **1.80 pt** from the ring centroid — i.e. on the centre,
  within hand-editing tolerance;
* the bond to the R label is **23.21 pt**, which is **1.61 × the 14.4 pt bond length**;
* the bond crosses ring bond **108–116**, which is **`Order="2"` — a double bond** — and passes
  **0.48 pt** from that bond's midpoint, i.e. through its middle;
* `Attachments` lists the ring atoms that carry **no other substituent** (107, 116, 108). The
  other three ring atoms (103, 106, 115) each bear something else and are excluded.

## The rules the toolkit implements

Given a `[*:n]` attached to a ring atom and `r_attachment: "centre"`:

1. **Start point** — the ring's centroid (mean of the ring atoms' coordinates). A
   `MultiAttachment` node is written there, not a carbon.
2. **Which bond is crossed** — a ring bond whose **two atoms are both free** (neither carries
   another substituent). Among those, **a double bond is preferred** (`Order=2`), matching both
   the witness and the published `R³`. If several remain, the one whose outward normal best
   matches the direction of the brief's declared atom wins, so the R leaves by the side the
   reader expects. If no ring bond has two free atoms, the writer leaves the R on the ring atom
   and warns (`ring_attachment_notes`).
3. **Bond length** — the R label sits `1.61 × bond_length` from the centroid along the crossed
   bond's outward normal (`R_ATTACHMENT_BOND_RATIO` in `overlay_geometry.py`, measured above).
   ChemDraw breaks the line where it crosses the ring, which is why the render shows a short
   inner segment and a short outer segment.
4. **Label offset** — the label node itself is moved to that point, and its `<t>` box is
   translated by the same delta, so the text stays centred on the bond end. No extra offset is
   added: the 1.61 ratio already includes the label's clearance.
5. **`Attachments`** — the brief's declared atom **first**, then the other free ring atoms.
   First position is what makes readback stable (below).
6. **Ring traversal** — the ring is walked in cyclic order from its bonds, not from
   `GetRingInfo().AtomRings()` order, because the crossed bond must be a real ring edge.

Note the geometry is **more compact** than a vertex bond: the label ends 1.61 bond lengths from
the centroid instead of roughly 2.0 (vertex at 1.0 plus one bond). The demo renders confirm it —
the after images are narrower (601 px vs 735 px for the mono case). Cells therefore do not grow.

## Readback: why 26/26 still passes

A `MultiAttachment` node is a drawing convention, not an atom. RDKit's CDXML reader treats it
as a carbon, which would turn one R into a stray methyl plus a bare ring. Two places fold it
back onto the **first** listed attachment atom, restoring exactly the brief's structure:

* `rdkit_utils.collapse_multi_attachments` — used by `frag_to_mol` and by
  `chemistry_semantics.read_fragments` (the `document_inventory` path);
* `figure_tools._strip_multi_attachments_for_inventory` — the composer's own inventory check.

Only a MultiAttachment bonded to a **`GenericNickname`** is folded. The same node type also
serves metal–ligand coordination (`multi_attachments` in a figure manifest), and those contacts
keep their previous reading.

## Using it

```jsonc
{
  "kind": "scope",
  "r_attachment": "centre",          // default "vertex": R stays on the atom the SMILES names
  "header": { "substrates": [{ "smiles": "[*:1]c1ccccc1", "label": "1" }], ... },
  "species": [...]
}
```

It applies to the header structures and to the species. Use `centre` when the scope means "R
anywhere on this ring"; leave it at `vertex` when the position is fixed.

## Evidence and reproduction

* Demo script: `scripts/rgroup_demo.py`
  (`--tag before --single`, `--tag after --single --r-attachment centre`, `--compare --single`).
* Before/after strips: `output/rgroup-demo/compare/{benzene_mono,benzene_di,heterocycle}.compare.png`.
* Objective check that does not depend on eyesight: the ring's interior is found by flood
  filling the background, and the ink at its centroid is counted. Before: **0** ink pixels near
  the centroid. After: **333** (mono), **672** (di). A vertex bond leaves the centre empty; a
  centre bond starts there.
* Full scope panels (header + species) also pass with `readback=matched` for all three cases.

## Known limits

* Opt-in only. No heuristic decides "unspecified" from SMILES, because SMILES always names one
  atom; the two corpus examples above show why guessing would be wrong.
* The **scheme** path (`figure_draw`) is not wired: the flag is implemented in the scope /
  `compose_chemical_figure` path, which is where the benchmark's R marks live.
* Pre-existing, unrelated to this change: an R label on the **left** of a structure renders as
  `²R` rather than `R²` (ChemDraw reverses the runs for a right-justified label). It appears in
  the before images too, and is not fixed here.
