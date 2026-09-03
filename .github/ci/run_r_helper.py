#!/usr/bin/env python3
"""Launch a fixed repository R helper through one explicit R installation."""

from __future__ import annotations

import argparse
import os
import stat
import subprocess
import sys
from pathlib import Path

from run_unit_driver import _resolve_rscript


HELPERS = {
    "prepare-rhub-dependencies": "prepare-rhub-dependencies.R",
    "verify-mkl": "verify-mkl.R",
    "parse-check": "run-r-cmd-check.R",
}


def _fixed_helper(helper_id: str) -> Path:
    filename = HELPERS.get(helper_id)
    if filename is None:
        raise ValueError(f"unknown R helper id: {helper_id!r}")
    path = Path(__file__).with_name(filename)
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"R helper is missing: {filename}") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"R helper must be a regular non-symlink file: {filename}")
    return path.resolve(strict=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-executable", required=True)
    parser.add_argument("--helper", required=True, choices=sorted(HELPERS))
    parser.add_argument("helper_argv", nargs=argparse.REMAINDER)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    helper_argv = arguments.helper_argv
    if helper_argv[:1] == ["--"]:
        helper_argv = helper_argv[1:]
    try:
        rscript = _resolve_rscript(arguments.r_executable)
        helper = _fixed_helper(arguments.helper)
        completed = subprocess.run(
            [os.fspath(rscript), "--vanilla", os.fspath(helper), *helper_argv],
            check=False,
            shell=False,
        )
    except (OSError, ValueError) as error:
        print(f"R helper configuration error: {error}", file=sys.stderr)
        return 2
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
