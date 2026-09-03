#!/usr/bin/env python3
"""Run strict unit diagnostics through one explicitly selected R installation."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from run_driver import _require_executable


def _resolve_rscript(r_executable: str | os.PathLike[str]) -> Path:
    r_binary = _require_executable(r_executable, label="R executable")
    if r_binary.name not in {"R", "R.exe"}:
        raise ValueError("R executable basename must be R or R.exe")
    sibling_name = "Rscript.exe" if r_binary.name == "R.exe" else "Rscript"
    return _require_executable(
        r_binary.with_name(sibling_name), label="Rscript sibling"
    )


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the unit-test reporter without any PATH-based R fallback."
    )
    parser.add_argument("--r-executable", required=True)
    parser.add_argument("--package", required=True)
    parser.add_argument("--load-package", required=True, choices=("source", "installed"))
    parser.add_argument("--context", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--filter")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _argument_parser().parse_args(argv)
    try:
        rscript = _resolve_rscript(arguments.r_executable)
        runner = Path(__file__).with_name("run-unit-tests.R").resolve(strict=True)
    except (OSError, ValueError) as error:
        print(f"unit driver configuration error: {error}", file=sys.stderr)
        return 2

    command = [
        os.fspath(rscript),
        "--vanilla",
        os.fspath(runner),
        f"--package={arguments.package}",
        f"--load-package={arguments.load_package}",
        f"--context={arguments.context}",
        f"--policy={arguments.policy}",
        f"--output={arguments.output}",
    ]
    if arguments.filter is not None:
        command.append(f"--filter={arguments.filter}")
    environment = dict(os.environ)
    environment["R_PROFILE_USER"] = os.devnull
    environment["R_ENVIRON_USER"] = os.devnull
    try:
        completed = subprocess.run(
            command, check=False, shell=False, env=environment
        )
    except OSError as error:
        print(f"unit driver execution error: {error}", file=sys.stderr)
        return 2
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
