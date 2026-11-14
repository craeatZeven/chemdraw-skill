"""The shipped energy-profile generator: portable defaults and the built-in house-rule check."""
import importlib
import re
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parents[1]
EP_DIR = SKILL_ROOT / 'scripts' / 'energy_profile'


def _verify():
    sys.path.insert(0, str(EP_DIR))
    try:
        return importlib.import_module('ep_verify')
    finally:
        sys.path.remove(str(EP_DIR))


def test_checker_falls_back_to_the_toolkit_check(tmp_path):
    from cdxml_toolkit.scientific.diagrams import build
    out = tmp_path / 'p.cdxml'
    out.write_text(build('periodic_table', {'periods': [1, 2]})['cdxml'], encoding='utf-8')
    result = _verify().run_checker(str(out), None)
    assert result['ok'] is True and result['checker'] == 'cdxml_toolkit.cdxml_check'
    bad = tmp_path / 'b.cdxml'
    bad.write_text(out.read_text(encoding='utf-8').replace('BondLength="14.4"', 'BondLength="30"'), encoding='utf-8')
    assert _verify().run_checker(str(bad), None)['ok'] is False


def test_no_machine_paths_or_agent_homes_in_the_generator():
    for path in EP_DIR.rglob('*.py'):
        text = path.read_text(encoding='utf-8')
        assert 'D:\\' not in text and 'C:\\' not in text, path.name
        # No agent-specific home directory as a default (the Skill stays client-neutral).
        assert not re.search(r'~/\.(?!agents/)[a-z]+/skills', text), path.name
