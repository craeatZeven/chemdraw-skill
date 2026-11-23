"""Exercise Windows ChemDraw discovery without launching COM or writing the registry."""

from contextlib import nullcontext
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from cdxml_toolkit import doctor
from cdxml_toolkit.chemdraw import chemscript_bridge, progid
from cdxml_toolkit.mcp_runtime import runtime_diagnostics, runtime_discovery


LEGACY, X64 = progid.CHEMDRAW_PROGIDS
KEY_READ, VIEW_64, VIEW_32 = 1, 256, 512


@pytest.fixture
def registry(monkeypatch):
    values = {}

    def open_key(hive, name, reserved=0, access=KEY_READ):
        assert hive == "HKCR"
        key = (access & (VIEW_32 | VIEW_64), name)
        if key not in values:
            raise FileNotFoundError(name)
        return nullcontext(values[key])

    fake = SimpleNamespace(
        HKEY_CLASSES_ROOT="HKCR",
        KEY_READ=KEY_READ, KEY_WOW64_64KEY=VIEW_64, KEY_WOW64_32KEY=VIEW_32,
        OpenKey=open_key, QueryValueEx=lambda value, name: (value, 1),
    )
    monkeypatch.setitem(sys.modules, "winreg", fake)
    monkeypatch.setattr(runtime_discovery, "_common_chemdraw_candidates", lambda: [])
    monkeypatch.setattr(chemscript_bridge.mac_automation, "is_macos", lambda: False)
    monkeypatch.delenv("CHEMDRAW_PROGID", raising=False)
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    return values


def register(registry, name, clsid, command, view=0):
    registry[view, rf"{name}\CLSID"] = clsid
    registry[view, rf"CLSID\{clsid}\LocalServer32"] = command


@pytest.mark.parametrize("legacy_command", [None, "", "bad command", "Missing.exe /Automation"])
def test_stale_legacy_registration_falls_back_to_x64(registry, tmp_path, legacy_command):
    exe = tmp_path / "ChemDraw.exe"
    exe.touch()
    registry[0, rf"{LEGACY}\CLSID"] = "{old}"
    if legacy_command is not None:
        register(registry, LEGACY, "{old}", legacy_command)
    register(registry, X64, "{new}", f'"{exe}" /Automation')

    assert progid.chemdraw_progid() == X64
    assert runtime_discovery.find_chemdraw() == runtime_discovery.Discovery(str(exe), "registry-HKCR")
    assert runtime_diagnostics._chemdraw_capabilities()[1]["registered"]
    assert doctor._check_chemdraw_com()[0]


def test_alternate_registry_view_and_flat_install(registry, tmp_path, monkeypatch):
    exe = tmp_path / "Flat Install" / "ChemDraw.exe"
    exe.parent.mkdir()
    exe.touch()
    monkeypatch.setenv("CHEMDRAW_TEST_ROOT", exe.parent.as_posix())
    alternate_view = VIEW_32 if sys.maxsize > 2**32 else VIEW_64
    register(registry, X64, "{new}", '"%CHEMDRAW_TEST_ROOT%/ChemDraw.exe" /Automation', alternate_view)

    assert progid.chemdraw_progid() == X64
    assert Path(runtime_discovery.find_chemdraw().path) == exe
    assert chemscript_bridge._find_chemdraw_root() == str(exe.parent)


def test_legacy_is_preferred_and_prog_id_override_is_strict(registry, tmp_path, monkeypatch):
    exe = tmp_path / "Program Files" / "Vendor" / "ChemDraw" / "ChemDraw.exe"
    exe.parent.mkdir(parents=True)
    exe.touch()
    register(registry, LEGACY, "{old}", f"{exe} /Automation")
    register(registry, X64, "{new}", f'"{exe}" /Automation')

    assert progid.chemdraw_progid() == LEGACY
    assert chemscript_bridge._find_chemdraw_root() == str(exe.parent.parent)
    monkeypatch.setenv("CHEMDRAW_PROGID", f" {X64} ")
    assert progid.chemdraw_progid() == X64
    assert runtime_discovery.find_chemdraw().path == str(exe)

    monkeypatch.setenv("CHEMDRAW_PROGID", "Absent.Application")
    monkeypatch.setenv("CHEMDRAW_EXE", str(exe))
    assert progid.chemdraw_progid() == "Absent.Application"
    executable, com = runtime_diagnostics._chemdraw_capabilities()
    assert executable["source"] == "environment"
    assert com["status"] == "missing"
    assert not doctor._check_chemdraw_com()[0]


def test_executable_fallback_does_not_imply_a_com_registration(registry, tmp_path, monkeypatch):
    exe = tmp_path / "ChemDraw.exe"
    exe.touch()
    register(registry, LEGACY, "{old}", "Missing.exe /Automation")
    monkeypatch.setattr(runtime_discovery, "_common_chemdraw_candidates", lambda: [exe])

    assert progid.chemdraw_progid() == LEGACY
    executable, com = runtime_diagnostics._chemdraw_capabilities()
    assert executable["source"] == "common-path"
    assert com["status"] == "missing"
