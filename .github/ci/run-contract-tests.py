#!/usr/bin/env python3
"""Run the CI contract suite without leaving pytest or bytecode caches."""

from __future__ import annotations

import os
import subprocess
import sys


def main() -> int:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            "-m",
            "pytest",
            "-q",
            "-p",
            "no:cacheprovider",
            ".github/ci/tests",
        ],
        check=False,
        env=environment,
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
