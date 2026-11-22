# Command-Line Workflows

## When To Load

Load only when no registered tool covers an operation, or when diagnosing the Python layer. Normal work goes through `cdx.py`.

## Preferred Entry Points

Check the installed runtime with `cdx.py doctor`. Use `cdx.py py -m MODULE` for a documented toolkit module CLI, or `cdx.py py` to run a helper script in the runtime. For toolkit developers, `[project.scripts]` in `toolkit/pyproject.toml` maps standalone `cdxml-*` commands to their implementation modules.

Prefer the equivalent tool through `cdx.py call` during normal work because it adds isolated execution, timeout handling, overwrite refusal, and a stable JSON contract.

## Inputs And Outputs

Pass absolute paths on Windows. Request JSON or structured errors when supported. Direct every modifying command to a new output path and validate the resulting file.

## Failure Modes

Console scripts can resolve to the wrong Python environment, write status text to stdout, exit through `SystemExit`, or lack worker timeouts. Verify the executable path and capture stdout/stderr separately.

## Supporting APIs

CLI implementations use the high-level Python functions described in the domain references. Use each command's `--help` as the authority for command-only flags.

## Do Not Use Directly

Do not assemble commands from the exhaustive symbol inventory. Do not invoke private ChemScript server commands or rely on machine-specific Scripts paths in Skill instructions.
