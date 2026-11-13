"""Which COM ProgID starts ChemDraw on this machine.

Older ChemDraw registers ``ChemDraw.Application``; some 64-bit installations register only
``ChemDraw_x64.Application``, and dispatching the missing name fails with "Invalid class
string". ``CHEMDRAW_PROGID`` overrides the choice; otherwise the first registered name wins.
"""
from __future__ import annotations

import os
from functools import lru_cache

CHEMDRAW_PROGIDS = ("ChemDraw.Application", "ChemDraw_x64.Application")


def registered(progid: str) -> bool:
    """Whether ``progid`` has a CLSID under HKEY_CLASSES_ROOT (always False off Windows)."""
    try:
        import winreg
    except ImportError:
        return False
    try:
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, progid + "\\CLSID"):
            return True
    except OSError:
        return False


@lru_cache(maxsize=1)
def _detected() -> str:
    for progid in CHEMDRAW_PROGIDS:
        if registered(progid):
            return progid
    return CHEMDRAW_PROGIDS[0]


def chemdraw_progid() -> str:
    """The ProgID to dispatch: ``$CHEMDRAW_PROGID``, else the first registered one."""
    override = os.environ.get("CHEMDRAW_PROGID", "").strip()
    return override or _detected()
