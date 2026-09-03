import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "verify_special_check.py"
MANIFEST = CI_DIR / "check-matrix.yml"
SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
DIGEST = "c" * 64


def invoke(tmp_path, *, environment_id="r-release-linux-x86-64", profile="full-documentation", log=None, exit_code=0):
    check_root = tmp_path / "input"
    check = check_root / "colocboost.Rcheck"
    check.mkdir(parents=True)
    if log is None:
        log = (CI_DIR / "tests/fixtures/release-full-doc-00check.log").read_text(
            encoding="utf-8"
        )
    (check / "00check.log").write_text(log, encoding="utf-8")
    native = tmp_path / "native.json"
    native.write_text(json.dumps({"wrapper_exit_code": exit_code}) + "\n", encoding="utf-8")
    runtime = tmp_path / "runtime.json"
    runtime.write_text(json.dumps({
        "environment_id": environment_id,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": DIGEST,
    }) + "\n", encoding="utf-8")
    output = tmp_path / "special.json"
    github_output = tmp_path / "github-output"
    completed = subprocess.run([
        sys.executable, "-B", os.fspath(SCRIPT), f"--manifest={MANIFEST}",
        f"--environment-id={environment_id}", f"--profile={profile}",
        f"--check-root={check_root}", f"--native-evidence={native}",
        f"--runtime-evidence={runtime}", f"--source-sha={SOURCE_SHA}",
        f"--event-sha={EVENT_SHA}", f"--tarball-sha256={DIGEST}",
        f"--output={output}", f"--github-output={github_output}",
    ], check=False, capture_output=True, text=True)
    return completed, output, github_output, check


def test_full_documentation_requires_executed_manual_and_vignettes(tmp_path):
    completed, output, github_output, _ = invoke(tmp_path)
    assert completed.returncode == 0, completed.stderr
    proof = json.loads(output.read_text(encoding="utf-8"))
    assert proof["profile"] == "full-documentation"
    assert proof["manual_executed"] is True
    assert proof["vignettes_executed"] is True
    assert proof["raw_exit_code"] == 0
    assert proof["effective_exit_code"] == 0
    assert "effective_exit_code=0\n" in github_output.read_text(encoding="utf-8")


@pytest.mark.parametrize("missing", ["manual", "package", "rebuilding"])
def test_full_documentation_rejects_missing_or_skipped_stage(tmp_path, missing):
    lines = {
        "manual": "* checking PDF version of manual ... OK",
        "package": "* checking package vignettes ... OK",
        "rebuilding": "* checking re-building of vignette outputs ... OK",
    }
    lines[missing] = lines[missing].replace("OK", "SKIPPED")
    completed, *_ = invoke(tmp_path, log="\n".join(lines.values()) + "\n* DONE\nStatus: OK\n")
    assert completed.returncode != 0


def test_standard_no_documentation_rejects_documentation_execution(tmp_path):
    completed, *_ = invoke(
        tmp_path,
        environment_id="r-devel-linux-x86-64-debian-clang",
        profile="standard-no-documentation",
    )
    assert completed.returncode != 0


def test_standard_no_documentation_accepts_clean_log_without_doc_stages(tmp_path):
    completed, output, _, _ = invoke(
        tmp_path,
        environment_id="r-devel-linux-x86-64-debian-clang",
        profile="standard-no-documentation",
        log="* checking examples ... OK\n* DONE\nStatus: OK\n",
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(output.read_text())["effective_exit_code"] == 0


@pytest.mark.parametrize("kind", ["nonzero", "duplicate", "symlink", "stale-runtime"])
def test_fails_closed_on_exit_tree_or_identity_ambiguity(tmp_path, kind):
    completed, output, github_output, check = invoke(
        tmp_path, exit_code=1 if kind == "nonzero" else 0
    )
    if kind == "duplicate":
        other = check.parent / "other.Rcheck"
        other.mkdir()
        (other / "00check.log").write_text("* DONE\nStatus: OK\n")
    elif kind == "symlink":
        target = check / "target"
        target.write_text("x")
        try:
            (check / "linked").symlink_to(target)
        except OSError:
            pytest.skip("symlinks unavailable")
    elif kind == "stale-runtime":
        runtime = tmp_path / "runtime.json"
        document = json.loads(runtime.read_text())
        document["source_sha"] = "d" * 40
        runtime.write_text(json.dumps(document))
    if kind != "nonzero":
        completed = subprocess.run(completed.args, check=False, capture_output=True, text=True)
    assert completed.returncode != 0
