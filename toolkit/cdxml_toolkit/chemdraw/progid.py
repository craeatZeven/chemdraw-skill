"""Read ChemDraw's Windows COM registration without starting the application."""
from __future__ import annotations

from dataclasses import dataclass
import ntpath
import os
from pathlib import Path
import re
import sys

CHEMDRAW_PROGIDS = ("ChemDraw.Application", "ChemDraw_x64.Application")


@dataclass(frozen=True)
class ChemDrawComRegistration:
    """The ProgID, LocalServer32 command, and resolved executable."""

    progid: str
    command: str
    executable: Path


def _registry_value(winreg, subkey: str, view: int) -> str:
    try:
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, subkey, 0, winreg.KEY_READ | view) as key:
            value, _ = winreg.QueryValueEx(key, "")
    except OSError:
        return ""
    return value.strip() if isinstance(value, str) else ""


def _server_executable(command: str) -> Path | None:
    # The default LocalServer32 value can include quotes and /Automation.
    match = re.match(
        r'^(?:"([^"]+\.exe)"|([^"].*?\.exe))(?:\s|$)',
        ntpath.expandvars(command).strip(),
        re.IGNORECASE,
    )
    if not match:
        return None
    path = Path(match.group(1) or match.group(2))
    try:
        return path.resolve() if path.is_file() else None
    except OSError:
        return None


def _find_registration(progids: tuple[str, ...]) -> ChemDrawComRegistration | None:
    try:
        import winreg
    except ImportError:
        return None

    alternate_view = winreg.KEY_WOW64_32KEY if sys.maxsize > 2**32 else winreg.KEY_WOW64_64KEY
    for view in (0, alternate_view):
        for progid in progids:
            clsid = _registry_value(winreg, rf"{progid}\CLSID", view)
            if not clsid:
                continue
            command = _registry_value(winreg, rf"CLSID\{clsid}\LocalServer32", view)
            executable = _server_executable(command)
            if executable is not None:
                return ChemDrawComRegistration(progid, command, executable)
    return None


def registered(progid: str) -> bool:
    """Whether the ProgID resolves to an installed ChemDraw COM server."""
    return _find_registration((progid,)) is not None


def chemdraw_com_registration() -> ChemDrawComRegistration | None:
    """Find a registered server, checking the process registry view first."""
    override = os.environ.get("CHEMDRAW_PROGID", "").strip()
    return _find_registration((override,) if override else CHEMDRAW_PROGIDS)


def chemdraw_progid() -> str:
    """Use the configured ProgID, a registered server, or the legacy default."""
    override = os.environ.get("CHEMDRAW_PROGID", "").strip()
    if override:
        return override
    registration = _find_registration(CHEMDRAW_PROGIDS)
    return registration.progid if registration else CHEMDRAW_PROGIDS[0]
