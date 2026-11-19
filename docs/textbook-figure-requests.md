# Textbook figures the toolkit cannot draw yet

Chapter 1 of *Advanced Organic Chemistry* (83 figures across 9 chapters) is being drawn with this
skill, one figure per `assets` style spec. Eleven of its fifteen figures are done. These four are
not: each needs a drawing capability the toolkit does not have. They are listed with what the
textbook figure has to show, because the shape of the demand decides the generator.

All four are **schematic orbital drawings**, not chemical structures: nothing here is a molecule
with atoms and bonds, so `lewis` and the structure pipeline cannot express them.

## 1. `orbital_shapes` (needed for figure 1-6)

Three panels of hybrid orbitals, drawn as lobes with phase shading:

- four sp3 lobes at 109.5 degrees (tetrahedral projection),
- three sp2 lobes at 120 degrees plus the unhybridised p lobe perpendicular to the plane,
- two sp lobes at 180 degrees plus two p lobes.

A fourth panel shows unequal hybridisation: sp3 lobes of unequal size, as in CH3Cl, with the bond
angle marked.

Produce: a generator that places lobe shapes (filled/hollow, size given per lobe) at given angles,
plus axis lines and angle arcs between lobes.
Done when: a spec of `{"hybridisation": "sp3"}` and friends renders the four panels, and a test
asserts each panel's lobe count and its angles to within a degree.

## 2. `orbital_overlap` (needed for figure 1-7)

Two panels: sigma overlap (two lobes pointing at each other head-on, along the internuclear axis)
and pi overlap (two p-orbital pairs side by side, above and below the bond axis), each with the
phase of the combining lobes marked and the internuclear axis drawn.

Produce: the same lobe primitive as above, placed in pairs, with the overlap region shaded and the
axis line between the two nuclei.
Done when: both panels render and a test asserts the phase parity of the four pi lobes
(with-phase/with-phase on one side, opposite on the other).

## 3. `frost_circle` (needed for figure 1-12)

The Frost circle construction for an n-centre cyclic pi system: the circle of radius 2 beta with
the inscribed n-gon, the energy axis through its centre, the n vertices projected onto the axis as
levels, and the electron filling of those levels.

Note the name collision: the existing `frost` operation is the **redox Frost diagram** (Mn species
against oxidation state), which is a different figure.

Produce: a generator taking n and the pi electron count.
Done when: n = 3, 4, 5, 6 render with the levels at 2, 0, -2 (and the degenerate pairs for 4 and
6), and the filling follows Huckel's rule, asserted in a test.

## 4. `fmo_interaction` (needed for figure 1-14)

A frontier-orbital interaction diagram: a filled donor HOMO on the left, an empty acceptor LUMO on
the right, the gap marked dE, and the two new levels the mixing produces, with the stabilisation of
the filled level shaded.

`energy_profile` with `style: levels` is not this figure: it draws a stepped reaction coordinate,
which is a different diagram.

Produce: a generator with two columns of levels and a mixing block between them.
Done when: a spec with one donor and one acceptor level renders, the generated levels straddle the
donor (filled stabilised, empty destabilised), and a test asserts the ordering and that the donor
level carries two electrons.

## What already exists and should be reused

- `scientific/lewis.py` has the electron-dot grammar but no lobes.
- `scientific/huckel.py` draws p lobes internally (`_p_orbital`) for the MO cartoons; that function
  is the natural seed for the lobe primitive, but it is private and only draws one lobe shape.
- `scientific/page.py` composes finished CDXML files, so once a generator exists the multi-panel
  layout, labels and arrows are already solved.

## Verified against the current skill

Every figure listed above was attempted with the current toolkit and refused or drawn as the wrong
diagram: `energy_profile` with `style: levels` produced a stepped reaction coordinate for 1-14, and
`frost` produced a redox diagram for 1-12. The other eleven chapter-1 figures render (see the
parity assets and the delivery folder).
