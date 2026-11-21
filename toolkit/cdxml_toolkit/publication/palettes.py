"""Curated semantic colour palettes for publication-grade chemical figures.

Each named palette pairs high-contrast print hues with explicit semantic roles
(``primary``, ``secondary``, ``accent``, ``new_bond``, ``muted``,
``highlight_fill``, ``band_colors``) and an ``allowed`` set shared between the
figure generators and :mod:`cdxml_toolkit.cdxml_check`.

Values are grounded in the Okabe-Ito colourblind-safe set and in measurements
from the openly licensed journal corpus (Nature Portfolio, ACS/JACS,
Wiley/Angewandte, Beilstein J. Org. Chem., RSC).
"""
from __future__ import annotations

import re
from types import MappingProxyType
from typing import Any, Mapping

_HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")
_NEUTRALS = ("#000000", "#FFFFFF", "#808080")


def _norm(color: str) -> str:
    text = str(color).strip()
    if not _HEX_RE.match(text):
        raise ValueError(f"Expected #RRGGBB colour, got {color!r}")
    return text.upper()


def _dedupe(colors: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    seen: list[str] = []
    for c in colors:
        h = _norm(c)
        if h not in seen:
            seen.append(h)
    return tuple(seen)


def _make(
    *,
    description: str,
    primary: str,
    secondary: str,
    accent: str,
    new_bond: str,
    muted: str = "#808080",
    highlight_fill: str = "#F2F6FA",
    band_colors: tuple[str, ...],
    extra_allowed: tuple[str, ...] = (),
) -> Mapping[str, Any]:
    allowed = _dedupe([
        *_NEUTRALS,
        primary,
        secondary,
        accent,
        new_bond,
        muted,
        highlight_fill,
        *band_colors,
        *extra_allowed,
    ])
    return MappingProxyType({
        "description": description,
        "primary": _norm(primary),
        "secondary": _norm(secondary),
        "accent": _norm(accent),
        "new_bond": _norm(new_bond),
        "muted": _norm(muted),
        "highlight_fill": _norm(highlight_fill),
        "band_colors": _dedupe(band_colors),
        "allowed": allowed,
    })


PALETTES: Mapping[str, Mapping[str, Any]] = MappingProxyType({
    "okabe_ito": _make(
        description="Okabe-Ito colourblind-safe scientific palette",
        primary="#0072B2",
        secondary="#D55E00",
        accent="#009E73",
        new_bond="#D55E00",
        highlight_fill="#EAF2F8",
        band_colors=("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"),
    ),
    "nature": _make(
        description="Nature Portfolio palette (teal, vermilion/orange, red, purple)",
        primary="#0078B2",
        secondary="#D55E00",
        accent="#E67E22",
        new_bond="#FF0000",
        highlight_fill="#EBF5FB",
        band_colors=("#0078B2", "#D55E00", "#009E73", "#6A3D9A", "#E67E22"),
        extra_allowed=("#0000FF", "#FF6600"),
    ),
    "acs1996": _make(
        description="Default ACS Document 1996 checker palette (steel blue, house orange, neutral grey)",
        primary="#336699",
        secondary="#FF6600",
        accent="#FF6600",
        new_bond="#FF6600",
        band_colors=("#336699", "#FF6600", "#808080"),
    ),
    "acs_jacs": _make(
        description="ACS / JACS / Org. Lett. palette (burgundy, navy, teal, purple)",
        primary="#800000",
        secondary="#004080",
        accent="#008060",
        new_bond="#DC3220",
        highlight_fill="#F9EBEB",
        band_colors=("#800000", "#004080", "#008060", "#6A3D9A"),
        extra_allowed=("#005AB5", "#FF6600", "#336699"),
    ),
    "wiley_angew": _make(
        description="Wiley / Angew. Chem. / Chem. Eur. J. palette",
        primary="#005AB5",
        secondary="#DC3220",
        accent="#008060",
        new_bond="#DC3220",
        highlight_fill="#EDF4FC",
        band_colors=("#005AB5", "#DC3220", "#008060", "#D97706", "#6A3D9A"),
    ),
    "beilstein": _make(
        description="Beilstein J. Org. Chem. palette (red, royal blue, teal)",
        primary="#D62728",
        secondary="#1F4E9E",
        accent="#008080",
        new_bond="#D62728",
        highlight_fill="#FDECEC",
        band_colors=("#1F4E9E", "#D62728", "#008080", "#2CA02C"),
        extra_allowed=("#FF6600",),
    ),
    "rsc": _make(
        description="Royal Society of Chemistry palette (navy, vermilion, forest green)",
        primary="#003B73",
        secondary="#C83E37",
        accent="#2E7D32",
        new_bond="#C83E37",
        highlight_fill="#EEF3F8",
        band_colors=("#003B73", "#C83E37", "#2E7D32", "#D97706"),
    ),
})

_PRESET_TO_PALETTE = {
    "acs1996": "acs1996",
    "nature": "nature",
    "wiley": "wiley_angew",
    "adv_synth_catal": "wiley_angew",
    "rsc": "rsc",
    "rsc_1col": "rsc",
}


def _collect_hex_colors(obj: Any, bucket: list[str]) -> None:
    if isinstance(obj, str):
        if _HEX_RE.match(obj.strip()):
            bucket.append(_norm(obj))
    elif isinstance(obj, Mapping):
        for v in obj.values():
            _collect_hex_colors(v, bucket)
    elif isinstance(obj, (list, tuple, set)):
        for item in obj:
            _collect_hex_colors(item, bucket)


def _palette_from_house(house: Any) -> dict[str, Any]:
    base_name = _PRESET_TO_PALETTE.get(getattr(house, "base_preset", ""), "okabe_ito")
    if getattr(house, "slug", "") in PALETTES:
        base_name = house.slug
    base = dict(PALETTES[base_name])

    card_colors: list[str] = []
    _collect_hex_colors(getattr(house, "conventions", {}), card_colors)
    _collect_hex_colors(getattr(house, "figure_types", {}), card_colors)

    subs = [c for c in (house.conventions.get("substrate_colors") or []) if isinstance(c, str) and _HEX_RE.match(c)]
    primary = _norm(subs[0]) if subs else base["primary"]
    secondary = _norm(subs[1]) if len(subs) > 1 else base["secondary"]
    new_bond = house.convention("new_bond_color") or base["new_bond"]
    accent = (
        house.figure_convention("substrate_scope", "variation_color")
        or house.convention("new_bond_color")
        or base["accent"]
    )
    highlight_fill = house.convention("highlight_fill") or base["highlight_fill"]
    band_seed = [c for c in [*subs, accent, *base["band_colors"]] if _norm(c) not in _NEUTRALS]

    return {
        "description": f"House-style palette for {house.display_name}",
        "primary": _norm(primary),
        "secondary": _norm(secondary),
        "accent": _norm(accent),
        "new_bond": _norm(new_bond),
        "muted": base["muted"],
        "highlight_fill": _norm(highlight_fill),
        "band_colors": _dedupe(band_seed),
        "allowed": _dedupe([
            *base["allowed"],
            *card_colors,
            primary,
            secondary,
            accent,
            new_bond,
            highlight_fill,
        ]),
    }


def resolve_palette(name_or_style: Any) -> dict[str, Any]:
    """Return a semantic palette dict for a palette name, preset, or HouseStyle.

    Raises ``ValueError`` for an unknown string name.
    """
    from .house_styles import HouseStyle, house_style, load_house_styles

    if isinstance(name_or_style, HouseStyle):
        return _palette_from_house(name_or_style)
    key = str(name_or_style or "").strip()
    if not key:
        return dict(PALETTES["okabe_ito"])
    lowered = key.lower()
    if lowered in PALETTES:
        return dict(PALETTES[lowered])
    if lowered in _PRESET_TO_PALETTE:
        return dict(PALETTES[_PRESET_TO_PALETTE[lowered]])
    if lowered in load_house_styles():
        return _palette_from_house(house_style(lowered))
    raise ValueError(
        f"Unknown palette {name_or_style!r}; available palettes: {', '.join(PALETTES)}"
    )


__all__ = ["PALETTES", "resolve_palette"]
