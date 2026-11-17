# Parity briefs

One brief per parity item that needed a figure to argue about, plus the render each produced
on the machine that has ChemDraw and RDKit. The PNGs are **the state of the drawing when the
brief was written**, not an approved target: two of them still show a defect that is recorded
here so a later fix has something to diff against.

Regenerate any of them with:

```bash
python skill/chemdraw/scripts/cdx.py draw assets/parity/<brief>.json <new>.cdxml --render
```

Scope briefs write a directory (`<new>.cdxml/png/page-001.png`); the others write a single
PNG beside the CDXML. The output path must not already exist — the tool refuses to overwrite.

## `p4-p7-p8-scope.json`

Two banded groups with italic titles, the varying group in `#FF6600`, and `3d′` as a grey
secondary product standing immediately right of `3d` with its own yield. Matches the Nature
Communications witness (`assets/journal-corpus/nature-communications-pmc13594109-fig4.png`).
No known defect.

## `p3-newman-mechanism.json`

Two butane conformers (anti, gauche) joined by a reversible arrow, and the hydrate mechanism
from the preview set below them.

**Known defect:** in the gauche conformation the rear substituents are not rotated relative to
the front ones, so three label pairs land on identical coordinates — `(297.21, 160.23)` H/CH3,
`(279.75, 129.99)` CH3/H, and `(262.29, 160.23)` H/H drawn twice — and the PNG shows a
collided glyph. The anti conformation is clean.

## `p9-cycle-overlay.json`

A four-species catalytic cycle with the methoxyphenyl fragment circled in every intermediate.
The circles come out identical (74.3 x 74.3 pt, radius 37.16, spread 0.0) and use the spelling
`native-dashed-overlay.cdxml` established: a `<curve>` with `CurveType="2"` and
`LineType="Dashed"`, no `FillType`, no `Closed`, never an oval.

**Two known defects:**

- the arrow `(296,171) -> (349,263)` passes 24.0 pt from a circle centre, inside its radius
  of 37.16, so it cuts through the dashed overlay;
- the step label `reductive elimination` overlaps the top-left circle.

Neither is caught by a check yet; `arrows_collide` covers arrow against arrow and bond, not
arrow or label against an overlay.
