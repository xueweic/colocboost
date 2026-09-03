import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / ".github" / "ci" / "run_unit_driver.py"


def make_executable(path, body="raise SystemExit(0)\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    path.chmod(0o755)
    return path


def invoke(r_executable, tmp_path, filter_value=None):
    received = tmp_path / "received args.json"
    environment = os.environ.copy()
    environment["UNIT_DRIVER_RECEIVED"] = os.fspath(received)
    arguments = [
            sys.executable,
            os.fspath(SCRIPT),
            "--r-executable",
            os.fspath(r_executable),
            "--package",
            "package path's value",
            "--load-package",
            "installed",
            "--context",
            "r-cmd-check-installed",
            "--policy",
            "policy path's value.yml",
            "--output",
            "output path's value.json",
        ]
    if filter_value is not None:
        arguments.extend(["--filter", filter_value])
    completed = subprocess.run(
        arguments,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    return completed, received


def test_unit_driver_uses_only_absolute_rscript_sibling_and_preserves_arguments(
    tmp_path,
):
    r_binary = make_executable(tmp_path / "R home with spaces" / "bin" / "R")
    recorder = (
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "Path(os.environ['UNIT_DRIVER_RECEIVED']).write_text(json.dumps(sys.argv[1:]))\n"
        "raise SystemExit(23)\n"
    )
    rscript = make_executable(r_binary.with_name("Rscript"), recorder)

    completed, received = invoke(r_binary, tmp_path)

    assert completed.returncode == 23, completed.stderr
    arguments = json.loads(received.read_text(encoding="utf-8"))
    assert arguments == [
        "--vanilla",
        os.fspath((ROOT / ".github" / "ci" / "run-unit-tests.R").resolve()),
        "--package=package path's value",
        "--load-package=installed",
        "--context=r-cmd-check-installed",
        "--policy=policy path's value.yml",
        "--output=output path's value.json",
    ]
    assert rscript.is_file()


def test_unit_driver_forwards_an_optional_filter_as_one_literal_argument(tmp_path):
    r_binary = make_executable(tmp_path / "R home" / "bin" / "R")
    recorder = (
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "Path(os.environ['UNIT_DRIVER_RECEIVED']).write_text(json.dumps(sys.argv[1:]))\n"
    )
    make_executable(r_binary.with_name("Rscript"), recorder)

    completed, received = invoke(r_binary, tmp_path, "^utils$; literal $(no-shell)")

    assert completed.returncode == 0, completed.stderr
    arguments = json.loads(received.read_text(encoding="utf-8"))
    assert arguments[-1] == "--filter=^utils$; literal $(no-shell)"


@pytest.mark.parametrize("kind", ["relative", "symlink", "nonexec", "wrong-basename"])
def test_unit_driver_rejects_unbound_r_executable(tmp_path, kind):
    if kind == "relative":
        r_binary = Path("R")
    elif kind == "symlink":
        target = make_executable(tmp_path / "real" / "R")
        r_binary = tmp_path / "linked-R"
        r_binary.symlink_to(target)
    elif kind == "nonexec":
        r_binary = tmp_path / "bin" / "R"
        r_binary.parent.mkdir()
        r_binary.write_text("not executable", encoding="utf-8")
    else:
        r_binary = make_executable(tmp_path / "bin" / "R-custom")

    completed, received = invoke(r_binary, tmp_path)

    assert completed.returncode == 2
    assert "R executable" in completed.stderr
    assert not received.exists()


@pytest.mark.parametrize("kind", ["missing", "symlink", "nonexec"])
def test_unit_driver_rejects_invalid_rscript_sibling(tmp_path, kind):
    r_binary = make_executable(tmp_path / "bin" / "R")
    rscript = r_binary.with_name("Rscript")
    if kind == "symlink":
        target = make_executable(tmp_path / "real-Rscript")
        rscript.symlink_to(target)
    elif kind == "nonexec":
        rscript.write_text("not executable", encoding="utf-8")

    completed, received = invoke(r_binary, tmp_path)

    assert completed.returncode == 2
    assert "Rscript sibling" in completed.stderr
    assert not received.exists()
