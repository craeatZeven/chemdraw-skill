"""Registry of the diagram workflows: operation name -> builder ``(spec, output) -> report``."""
from __future__ import annotations

from typing import Any, Callable, Mapping


def _compose(spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Merge the CDXML ``parts`` of a spec into one page (see :func:`cdxml_compose.merge_cdxml`)."""
    from ..cdxml_compose import merge_cdxml
    from .energy import _save_text
    parts = spec.get('parts')
    if not isinstance(parts, list) or not parts:
        raise ValueError('compose needs "parts": [{"file": ..., "at": [x, y], "anchor": ...}]')
    result = merge_cdxml(parts, gap=float(spec.get('gap', 20.0)), margin=float(spec.get('margin', 12.0)))
    result['diagram'] = 'compose'
    if output is not None:
        result['output_path'] = _save_text(result['cdxml'], output)
    return result


def _load() -> dict[str, Callable[..., dict[str, Any]]]:
    from .binary import binary_phase_diagram
    from .crystal import unit_cell
    from .energy import energy_cycle, energy_profile, jablonski
    from .huckel import huckel_mo
    from .lewis import lewis_structure
    from .orbitals import crystal_field, electron_configuration_diagram, mo_diagram
    from .periodic import periodic_table
    from .phase import phase_diagram
    from .process import electrochemical_cell, flowchart
    return {
        'energy_profile': energy_profile,
        'energy_cycle': energy_cycle,
        'jablonski': jablonski,
        'mo_diagram': mo_diagram,
        'electron_configuration': electron_configuration_diagram,
        'crystal_field': crystal_field,
        'lewis': lewis_structure,
        'unit_cell': unit_cell,
        'flowchart': flowchart,
        'electrochemical_cell': electrochemical_cell,
        'periodic_table': periodic_table,
        'huckel_mo': huckel_mo,
        'phase_diagram': phase_diagram,
        'binary_phase_diagram': binary_phase_diagram,
        'compose': _compose,
    }


# Names only, so that importing this module stays cheap; ``build`` imports the builder.
DIAGRAM_OPERATIONS = ('energy_profile', 'energy_cycle', 'jablonski', 'mo_diagram', 'electron_configuration',
                      'crystal_field', 'lewis', 'unit_cell', 'flowchart', 'electrochemical_cell', 'periodic_table',
                      'huckel_mo', 'phase_diagram', 'binary_phase_diagram', 'compose')

def build(operation: str, spec: Mapping[str, Any], output=None) -> dict[str, Any]:
    """Run one diagram workflow; the returned report omits the CDXML text when a file was written."""
    builders = _load()
    if operation not in builders:
        raise ValueError(f'Unknown diagram {operation!r}; available: {", ".join(DIAGRAM_OPERATIONS)}')
    result = builders[operation](spec, output)
    if output is not None:
        result.pop('cdxml', None)
    return result
