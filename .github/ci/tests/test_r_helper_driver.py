import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / ".github" / "ci" / "run_r_helper.py"


def make_executable(path, body="raise SystemExit(0)\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_launcher_uses_absolute_rscript_sibling_and_fixed_helper(tmp_path):
    received = tmp_path / "received.json"
    r_binary = make_executable(tmp_path / "R home" / "bin" / "R")
    make_executable(
        r_binary.with_name("Rscript"),
        "import json, os, pathlib, sys\n"
        "pathlib.Path(os.environ['RECEIVED']).write_text(json.dumps(sys.argv[1:]))\n"
        "raise SystemExit(17)\n",
    )
    environment = os.environ.copy()
    environment["RECEIVED"] = os.fspath(received)

    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(SCRIPT),
            "--r-executable",
            os.fspath(r_binary),
            "--helper",
            "verify-mkl",
            "--",
            "--policy=path with spaces/policy.yml",
            "--evidence=literal;not-shell.json",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert completed.returncode == 17, completed.stderr
    arguments = json.loads(received.read_text(encoding="utf-8"))
    assert arguments[0] == "--vanilla"
    assert arguments[1] == os.fspath((ROOT / ".github/ci/verify-mkl.R").resolve())
    assert arguments[2:] == [
        "--policy=path with spaces/policy.yml",
        "--evidence=literal;not-shell.json",
    ]


@pytest.mark.parametrize("helper", ["../arbitrary.R", "unknown", "prepare-rhub-dependencies.R"])
def test_launcher_rejects_non_allowlisted_helper(tmp_path, helper):
    r_binary = make_executable(tmp_path / "bin" / "R")
    make_executable(r_binary.with_name("Rscript"))

    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(SCRIPT),
            "--r-executable",
            os.fspath(r_binary),
            "--helper",
            helper,
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
