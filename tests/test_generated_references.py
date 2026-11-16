"""The generated references must match the bundled toolkit before anything is pushed.

Both used to be caught only by the Skill suite, after a push. The interface inventory is
read from the toolkit's source, so it is checked here with the standard library. The tool
signatures come from the live tool registry, which imports the toolkit's dependencies, so
that check runs in the toolkit runtime ``cdx.py py`` selects, against this checkout's
toolkit, and is skipped (saying so) on a machine without a runtime.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "skill" / "chemdraw" / "scripts"
TOOLKIT = ROOT / "toolkit"
NO_RUNTIME = "no working cdxml-toolkit runtime"


class GeneratedReferenceTests(unittest.TestCase):
    def test_interface_inventory_matches_the_bundled_toolkit(self) -> None:
        init = (TOOLKIT / "cdxml_toolkit" / "__init__.py").read_text(encoding="utf-8")
        version = re.search(r'__version__\s*=\s*["\']([^"\']+)', init).group(1)
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "audit_toolkit_interfaces.py"),
             "--package-root", str(TOOLKIT / "cdxml_toolkit"),
             "--output-dir", str(ROOT / "skill" / "chemdraw" / "references"),
             "--docs-dir", str(TOOLKIT / "docs"), "--version", version, "--check"],
            capture_output=True, text=True, encoding="utf-8", timeout=600,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_tool_signatures_match_the_bundled_toolkit(self) -> None:
        code = ("import sys; sys.path.insert(0, %r); import cdx; "
                "sys.exit(0 if cdx._signatures_current() else 3)" % str(SCRIPTS))
        # Run from inside toolkit/: `python -c` imports the current directory first, so the
        # runtime reads this checkout's toolkit rather than the copy installed in it.
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "cdx.py"), "py", "-c", code],
            cwd=str(TOOLKIT), capture_output=True, text=True, encoding="utf-8", timeout=600,
        )
        if result.returncode != 0 and NO_RUNTIME in result.stdout + result.stderr:
            self.skipTest("no toolkit runtime on this machine; the Skill suite checks the signatures")
        self.assertNotEqual(result.returncode, 3, "the signatures are stale: run "
                                                  "`python skill/chemdraw/scripts/cdx.py regen-signatures`")
        self.assertEqual(result.returncode, 0, (result.stdout + result.stderr)[-3000:])


if __name__ == "__main__":
    unittest.main()
