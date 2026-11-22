# Scientific Workflows Public Symbols

> Generated from cdxml-toolkit 0.7.0a3. Curated guidance: [../scientific-workflows.md](../scientific-workflows.md).

## `scientific.__main__`

- **function**, line 8: `run(operation, spec_path, output)` - No public docstring in the audited version.

## `scientific.apparatus`

- **function**, line 20: `draw_apparatus(spec, output)` - No public docstring in the audited version.
- **function**, line 9: `template_catalog(path)` - No public docstring in the audited version.

## `scientific.axes`

- **function**, line 129: `clip_polygon(polygon, a, b, c)` - Keep the part of a convex polygon where ``a x + b y <= c`` (Sutherland-Hodgman, one edge).
- **function**, line 109: `clip_segment(a, b, box)` - No public docstring in the audited version.
- **function**, line 10: `nice_ticks(low, high, count = 5)` - Round tick values covering ``[low, high]`` (1, 2, 2.5 or 5 times a power of ten).
- **class**, line 36: `Plot` - Map ``(x, y)`` data into a ``width`` by ``height`` frame whose top-left corner is the origin.
- **method**, line 49: `Plot.at(x, y)` - No public docstring in the audited version.
- **method**, line 55: `Plot.frame(x_title = '', y_title = '', ticks = True, size = 8)` - No public docstring in the audited version.
- **method**, line 52: `Plot.inside(x, y)` - No public docstring in the audited version.
- **method**, line 77: `Plot.polyline(points, dashed = False, color = None, clip = True)` - Draw data points as one or more lines, cut where they leave the frame.
- **method**, line 94: `Plot.segment(a, b, dashed = False, color = None)` - A data-space straight line clipped to the frame (Liang-Barsky).
- **function**, line 26: `tick_text(value, step)` - No public docstring in the audited version.

## `scientific.binary`

- **function**, line 41: `binary_phase_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Temperature-composition diagram of components A and B (x is the mole or mass fraction of B).

## `scientific.canvas`

- **function**, line 143: `boxes_overlap(a, b, gap = 0.0)` - No public docstring in the audited version.
- **class**, line 147: `Canvas` - A :class:`Drawing` that tracks its extent and writes rich captions.
- **method**, line 207: `Canvas.arrow(tail, head, color = None, dashed = False, head_style = 'Full', tail_style = None, head_size = 1000, parent = None)` - A straight arrow. ``head_style``/``tail_style``: ``Full``, ``HalfLeft``, ``HalfRight`` or ``None``.
- **method**, line 327: `Canvas.block(box, pad = 1.0)` - Reserve ``box`` so that :meth:`place` keeps captions out of it.
- **method**, line 333: `Canvas.block_path(points, pad = 1.5)` - Reserve the ink of a polyline: small boxes along each segment.
- **method**, line 365: `Canvas.choose(candidates, markup, size = None)` - The first free candidate ``(xy, align, valign)``, or the least crowded one.
- **method**, line 224: `Canvas.circle(center, radius, fill = None, outline = BLACK, dashed = False, parent = None)` - A circle; ``fill`` paints the disc first so it hides what lies behind it.
- **property**, line 173: `Canvas.extent` - No public docstring in the audited version.
- **method**, line 410: `Canvas.finish(margin = 12.0)` - Move the content to ``margin`` from the page origin and size the page around it.
- **method**, line 349: `Canvas.is_free(box) -> bool` - Whether ``box`` meets no reserved ink and no caption drawn so far.
- **method**, line 278: `Canvas.label(xy, markup, size = None, align = 'Left', valign = 'baseline', color = None, bold = False, parent = None)` - Rich caption at ``xy``. Returns its estimated box ``(x1, y1, x2, y2)``.
- **method**, line 316: `Canvas.measure(xy, markup, size = None, align = 'Left', valign = 'baseline')` - The box :meth:`label` would report, without drawing anything.
- **method**, line 401: `Canvas.overlapping_labels(gap = 0.5)` - Pairs of captions whose estimated boxes overlap: a layout warning, not a proof.
- **method**, line 193: `Canvas.path(points: Iterable[Sequence[float]], color = None, dashed = False, closed = False, filled = False, parent = None)` - A polyline through ``points`` (an editable curve with degenerate handles).
- **method**, line 354: `Canvas.path_is_free(points, step = 1.0) -> bool` - Whether a polyline would cross no caption drawn so far (a leader or connector).
- **method**, line 380: `Canvas.place(candidates, markup, size = None, color = None, bold = False, parent = None, knockout = False)` - Draw ``markup`` at the first candidate ``(xy, align, valign)`` that is free.
- **method**, line 239: `Canvas.rect(x1, y1, x2, y2, color = None, dashed = False, rounded = False, fill = None, parent = None)` - A rectangle outline, optionally over a solid fill.
- **method**, line 164: `Canvas.see(x, y, pad = 0.0)` - No public docstring in the audited version.
- **method**, line 169: `Canvas.see_box(box)` - No public docstring in the audited version.
- **method**, line 187: `Canvas.segment(a, b, color = None, dashed = False, parent = None)` - A straight line from ``a`` to ``b``.
- **method**, line 442: `Canvas.to_string()` - No public docstring in the audited version.
- **method**, line 256: `Canvas.wavy(a, b, color = None, amplitude = 2.2, wavelength = 7.0, arrow = True, parent = None)` - A wavy line from ``a`` to ``b`` (a non-radiative transition), ending in an arrowhead.
- **method**, line 447: `Canvas.write(output = None, margin = 12.0)` - Finish the page and save it; returns ``(saved path or None, CDXML text)``.
- **function**, line 122: `caption_box(x, baseline, markup, size = FONT_SIZE, align = 'Left')` - Box of a caption whose first baseline is at ``baseline``, as ChemDraw measures it.
- **function**, line 78: `parse_markup(markup)` - ``'E{sub:a}'`` -> ``[('E', 0), ('a', 32)]`` (text, CDXML face bits).
- **function**, line 92: `plain(markup)` - The text of a markup string without its formatting tokens.
- **function**, line 64: `resolve_color(value, default = BLACK)` - A palette name or ``#RRGGBB`` (or ``#RGB``) as upper-case ``#RRGGBB``.
- **function**, line 139: `text_height(markup, size = FONT_SIZE)` - No public docstring in the audited version.
- **function**, line 97: `text_width(markup, size = FONT_SIZE)` - Approximate width of one line of Arial/SimSun text, in points.

## `scientific.common`

- **class**, line 23: `Drawing` - No public docstring in the audited version.
- **method**, line 75: `Drawing.add(kind, attrs, parent = None)` - No public docstring in the audited version.
- **method**, line 37: `Drawing.color(rgb_hex: str) -> str` - No public docstring in the audited version.
- **method**, line 102: `Drawing.curve(points, parent = None, closed = False, filled = False, color: str | None = None, dashed: bool = False)` - No public docstring in the audited version.
- **method**, line 120: `Drawing.filled_rect(x1, y1, x2, y2, *, fill_color: str, border_color: str | None = None, parent = None) -> ET.Element` - No public docstring in the audited version.
- **method**, line 92: `Drawing.line(a, b, parent = None, color: str | None = None, dashed: bool = False)` - No public docstring in the audited version.
- **method**, line 152: `Drawing.save(output)` - No public docstring in the audited version.
- **method**, line 80: `Drawing.text(xy, value, size = 10, parent = None, color: str | None = None)` - No public docstring in the audited version.
- **function**, line 8: `number(value, name, low = None, high = None)` - No public docstring in the audited version.
- **function**, line 17: `point(value)` - No public docstring in the audited version.

## `scientific.correlation`

- **function**, line 40: `correlation_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Orbital correlation diagram.
- **function**, line 134: `walsh_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Walsh diagram: orbital energies against a geometric coordinate.

## `scientific.crystal`

- **function**, line 74: `unit_cell(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Cubic unit cell.

## `scientific.diagrams`

- **function**, line 75: `build(operation: str, spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Run one diagram workflow; the returned report omits the CDXML text when a file was written.

## `scientific.elements`

- **function**, line 32: `atomic_number(element)` - Atomic number of a symbol (``'Fe'``) or of an integer given as such.
- **function**, line 153: `d_electron_count(element, oxidation_state = 0)` - d-electron count of a transition-metal centre: group number minus oxidation state.
- **function**, line 61: `electron_configuration(element, charge = 0)` - Ground-state configuration as ``{'1s': 2, ...}`` ordered by shell, then subshell.
- **function**, line 103: `group_period(element)` - ``(group, period)`` in the 18-column table; lanthanoids and actinoids have group ``None``.
- **function**, line 93: `noble_gas_core(element, charge = 0)` - ``(symbol, core_subshells)`` of the largest noble-gas core inside the species.
- **function**, line 136: `parse_ion(text)` - ``'Fe3+'``, ``'Fe(III)'``, ``'Cl-'`` or ``'Cu'`` -> ``(symbol, charge)``.
- **function**, line 44: `subshell_name(n, l)` - No public docstring in the audited version.
- **function**, line 122: `valence_electrons(element)` - Valence electrons of a main-group element (Lewis structures).

## `scientific.energy`

- **function**, line 583: `energy_cycle(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Enthalpy-level diagram with stepwise arrows: Born-Haber and Hess cycles.
- **function**, line 196: `energy_profile(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Reaction energy profile (potential energy against reaction coordinate).
- **function**, line 92: `from_species_spec(spec: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]` - Translate the species-list energy-profile format into this builder's pathway format.
- **function**, line 737: `jablonski(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Jablonski diagram: electronic states with vibrational levels and their transitions.

## `scientific.frost_circle`

- **function**, line 24: `frost_circle(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Draw a cyclic pi HMO Frost circle, distinct from the redox ``frost`` operation.

## `scientific.hexagonal`

- **function**, line 52: `hexagonal_cell(key: str, spec: Mapping[str, Any], output = None) -> dict[str, Any]` - No public docstring in the audited version.

## `scientific.huckel`

- **function**, line 69: `huckel_mo(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Hückel pi-MO diagram.
- **function**, line 42: `solve(atoms, cyclic = False, bonds = None, h_values = None, k_values = None)` - Energies ``x`` (E = alpha + x beta, bonding x > 0) and coefficient vectors, highest x first.

## `scientific.lewis`

- **function**, line 24: `grid_xy(col, row, grid)` - Page position of a grid point in bond lengths: square, or axial hexagonal coordinates.
- **function**, line 186: `lewis_structure(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Lewis structure of one or more species.

## `scientific.mechanism`

- **function**, line 36: `compile_mechanism(spec, *, checks = None)` - No public docstring in the audited version.
- **function**, line 8: `validate_mechanism_steps(steps, molecules)` - Check fully specified equations; partial illustrations stay explicitly unchecked.

## `scientific.orbital_lobes`

- **function**, line 290: `draw_orbital_overlap(canvas, spec: Mapping[str, Any]) -> dict[str, Any]` - Draw sigma, pi, or both overlap panels on an existing canvas.
- **function**, line 161: `draw_orbital_shapes(canvas, spec: Mapping[str, Any]) -> dict[str, Any]` - Draw one hybridisation panel on an existing canvas.
- **function**, line 95: `lobe(canvas, at, angle_deg = -90.0, size = 36.0, width = 16.0, coefficient = 1.0, phase = 1, fill = 'filled', color = None, parent = None)` - Place one ellipse lobe, with its narrow end at ``at`` and axis at ``angle_deg``.
- **function**, line 305: `orbital_overlap(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Draw head-on sigma and side-by-side pi overlap with phase-coloured lobes.
- **function**, line 224: `orbital_shapes(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Draw sp3, sp2, sp or unequal hybrid lobes with axes and angle annotations.
- **function**, line 106: `p_orbital(canvas, at, angle_deg = -90.0, size = 36.0, width = 16.0, coefficient = 1.0, phase = 1, fill = 'filled', color = None, parent = None, negative_color = None, axis = 'xy')` - Two opposite lobes with opposite phases; ``axis='z'`` is normal to the page.

## `scientific.orbitals`

- **function**, line 421: `crystal_field(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Crystal-field splitting diagram of the d orbitals with the electrons filled in.
- **function**, line 399: `crystal_field_fill(geometry, d_electrons, spin = None)` - Occupation of each orbital set: ``[['ud', 'u', 'u'], ['u', 'u']]`` from the lowest set up.
- **function**, line 84: `electron_configuration_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Orbital-box diagram of an electron configuration.
- **function**, line 33: `hund_fill(electrons, orbitals)` - Occupation of ``orbitals`` degenerate orbitals: singly with parallel spins first, then paired.
- **function**, line 242: `mo_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Molecular-orbital energy-level diagram: atomic orbitals left and right, MOs between.

## `scientific.page`

- **function**, line 166: `compose_page(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Lay out ``rows`` of CDXML ``parts`` as one page.
- **function**, line 129: `runs(text)` - [(piece, face)]: 32 on an atom count after a letter (CH3), 64 on a charge written ^- or ^2-.

## `scientific.periodic`

- **function**, line 24: `periodic_table(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Periodic table.
- **function**, line 13: `table_position(z)` - ``(row, column, block)`` of element ``z`` in the 18-column table; the f block sits in two rows below.

## `scientific.phase`

- **function**, line 10: `phase_diagram(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Pressure-temperature phase diagram of one substance, schematic.

## `scientific.plots`

- **function**, line 375: `build_lcms_plate_cdxml(source_or_spec: Any, output_path: str | Path | None = None, *, product_mass: float | None = None, sm_mass: float | None = None, plate_format: int | str | None = None, wavelength: int | str = 254, title: str = 'LCMS Plate Purity Heatmap') -> dict[str, Any]` - Build an editable 24-well or 96-well LCMS purity heatmap in CDXML.
- **function**, line 13: `envelope_indices(values, limit = 1600)` - No public docstring in the audited version.
- **function**, line 59: `plot_xy(x = None, y = None, output_dir = None, *, title = 'Scientific plot', xlabel = 'x', ylabel = 'y', reverse_x = False, provenance = None, dpi = 300, series = None)` - No public docstring in the audited version.

## `scientific.process`

- **function**, line 277: `electrochemical_cell(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Galvanic or electrolytic cell.
- **function**, line 20: `flowchart(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Left-to-right flow chart.

## `scientific.redox`

- **function**, line 305: `ellingham(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Ellingham diagram: ``lines``: ``[{"label": "2Mg + O{sub:2} → 2MgO", "dH": -1203.4, "dS": -216.8}]`` (kJ/mol and J/(mol K), both per mole of O2 and taken as constant), optional ``segments`` ``[{"from_T": 1380, "dH": ..., "dS": ...}]`` after a phase change; ``T`` range in K.
- **function**, line 252: `frost(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Frost diagram (N E° against oxidation state N) from the same spec as :func:`latimer`. Several chains go in ``"series": [{"label", "species", "potentials"}]`` (e.g. acid and base). For the cyclic HMO polygon construction, use ``frost_circle`` instead.
- **function**, line 209: `latimer(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Latimer diagram: ``{"species": [{"name": "MnO4{sup:−}", "oxidation_state": 7}, ...], "potentials": [0.56, ...], "medium": "acid"}``, highest oxidation state first. ``derived`` lists couples (``[["MnO4{sup:−}", "MnO2"]]``) to add as arcs with their computed E°.
- **function**, line 27: `pourbaix(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Potential-pH (Pourbaix) diagram of one element.

## `scientific.spectra`

- **function**, line 626: `analyze_spectrum(path, options = None)` - No public docstring in the audited version.
- **function**, line 866: `build_nmr_cdxml(ppm: Any = None, intensity: Any = None, *, peaks: Optional[Sequence[Any]] = None, inset_region: Optional[Union[Tuple[float, float], list[float]]] = None, title: str = '1H NMR Spectrum', nucleus: str = '1H', output_path: Optional[Union[str, Path]] = None, x: Any = None, y: Any = None, spectrum: Any = None) -> str` - Build a ChemDraw CDXML document for an NMR spectrum with optional zoomed ``inset_region``.
- **function**, line 1244: `build_stacked_nmr_cdxml(spectra: Sequence[str | Path | Mapping[str, Any]], output_path: str | Path | None = None, *, labels: Sequence[str] | None = None, ppm_range: tuple[float, float] | None = None, title: str = 'Stacked NMR Comparison', normalize: str = 'each', show_peaks: bool = True, match_shifts: Sequence[Mapping[str, Any]] | None = None) -> dict[str, Any]` - Build a multi-trace stacked NMR comparison CDXML figure (2 to 8 spectra).
- **function**, line 1020: `nmr_to_cdxml(ppm: Any = None, intensity: Any = None, *, peaks: Optional[Sequence[Any]] = None, inset_region: Optional[Union[Tuple[float, float], list[float]]] = None, title: str = '1H NMR Spectrum', nucleus: str = '1H', output_path: Optional[Union[str, Path]] = None, x: Any = None, y: Any = None, spectrum: Any = None) -> str` - Convert NMR spectrum arrays, peak list, or file path into a CDXML string (with optional ``inset_region``).
- **function**, line 1048: `plot_nmr(ppm: Any = None, intensity: Any = None, *, peaks: Optional[Sequence[Any]] = None, inset_region: Optional[Union[Tuple[float, float], list[float]]] = None, title: str = '1H NMR Spectrum', nucleus: str = '1H', output_path: Optional[Union[str, Path]] = None, output_dir: Optional[Union[str, Path]] = None, x: Any = None, y: Any = None, spectrum: Any = None) -> Union[str, dict[str, Any]]` - Plot an NMR spectrum to CDXML (and optional ``output_dir`` bundle) with ``inset_region`` zoom support.
- **function**, line 525: `read_spectrum(path)` - No public docstring in the audited version.

## `scientific.surface`

- **function**, line 33: `contour_segments(xs, ys, z, level)` - Line segments of one contour level (marching squares, saddle cells split by the centre value).
- **function**, line 22: `mueller_brown(x, y)` - The Mueller-Brown potential (Theor. Chim. Acta 1979, 53, 75): three minima, two saddles.
- **function**, line 139: `potential_surface(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Contour plot of a 2D potential-energy surface.

## `scientific.tanabe`

- **function**, line 60: `coulomb(m1, m2, m3, m4, f2, f4)` - <m1 m2 | 1/r12 | m3 m4> for d electrons, F^0 set to zero (it shifts every level alike).
- **function**, line 54: `gaunt(k, m1, m2, l = 2)` - Condon-Shortley c^k(l m1, l m2).
- **function**, line 103: `hamiltonian(n, c_over_b, n_alpha)` - (H_ee in units of B, V in units of Dq, S+ to the next block, L+, M_L list, M_S) for one M_S block.
- **function**, line 193: `levels(n, delta_over_b, c_over_b)` - Labelled levels at one Delta/B: ``[(energy/B, 2S+1, label)]``, lowest first; ``None`` labels mark a sample where two different terms coincide (an accidental degeneracy).
- **function**, line 222: `tanabe_sugano(spec: Mapping[str, Any], output = None) -> dict[str, Any]` - Tanabe-Sugano diagram of d^n in an octahedral field.
- **function**, line 35: `wigner_3j(j1, j2, j3, m1, m2, m3)` - Wigner 3j symbol (integer arguments), Racah formula.

## `scientific.tlc`

- **function**, line 338: `draw_tlc(spec, output = None)` - No public docstring in the audited version.
- **function**, line 49: `measure_rf(origin, solvent_front, spot)` - Project a measured spot along the calibrated solvent travel direction.
