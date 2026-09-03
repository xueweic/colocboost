import json
import os
import stat
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402
from run_native_check import run_native_check  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def make_executable(path, body="raise SystemExit(0)\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def make_contract(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    tarball = source / "colocboost_1.0.9.tar.gz"
    tarball.write_bytes(b"verified immutable bytes")
    metadata = source / "source-metadata.json"
    create_metadata(tarball, metadata, source_sha=SOURCE_SHA, event_sha=EVENT_SHA)
    return tarball, metadata


def invoke(tmp_path, *, environment_id="mkl", wrapper_body=None):
    tarball, metadata = make_contract(tmp_path)
    if wrapper_body is None:
        wrapper_body = (
            "import pathlib, sys\n"
            "parent = pathlib.Path(sys.argv[1])\n"
            "assert len(list(parent.glob('*.tar.gz'))) == 1\n"
            "check = parent / 'colocboost.Rcheck'\n"
            "check.mkdir()\n"
            "(check / '00check.log').write_text('* DONE\\nStatus: OK\\n')\n"
        )
        if environment_id == "nosuggests":
            wrapper_body += (
                "(check / 'tests').mkdir()\n"
                "(check / 'tests' / 'testthat.Rout').write_text('1 test passed\\n')\n"
            )
    wrapper = make_executable(tmp_path / "bin" / "r-check", wrapper_body)
    system_r = make_executable(tmp_path / "bin" / "R")
    check_library = tmp_path / "check-library"
    work_dir = tmp_path / "check-work"
    output_log = tmp_path / "published" / "00check.log"
    evidence = tmp_path / "evidence.json"
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)
    environment["R_LIBS_USER"] = os.fspath(check_library)
    row = {
        "id": environment_id,
        "driver": "native-wrapper",
        "wrapper_path": os.fspath(wrapper),
        "system_r": os.fspath(system_r),
        "wrapper_input": "tarball-parent",
    }
    result = run_native_check(
        row,
        wrapper=wrapper,
        required_r_executable=system_r,
        tarball=tarball,
        metadata=metadata,
        expected_source_sha=SOURCE_SHA,
        expected_event_sha=EVENT_SHA,
        work_dir=work_dir,
        check_library=check_library,
        evidence_path=evidence,
        log_output=output_log,
        environment=environment,
    )
    return result, work_dir, check_library, output_log, evidence


def test_runs_wrapper_in_fresh_tree_and_publishes_exact_log_and_evidence(tmp_path):
    code, work_dir, check_library, output_log, evidence = invoke(tmp_path)

    assert code == 0
    assert check_library.is_dir()
    assert output_log.read_text() == "* DONE\nStatus: OK\n"
    copied = list((work_dir / "input").glob("*.tar.gz"))
    assert len(copied) == 1
    document = json.loads(evidence.read_text())
    assert document["wrapper_exit_code"] == 0
    assert document["check_log"].endswith("colocboost.Rcheck/00check.log")
    assert document["environment"]["R_PROFILE_USER"] == os.devnull
    assert document["environment"]["R_ENVIRON_USER"] == os.devnull


def test_preserves_nonzero_wrapper_exit_when_a_unique_log_exists(tmp_path):
    body = (
        "import pathlib, sys\n"
        "parent = pathlib.Path(sys.argv[1])\n"
        "check = parent / 'colocboost.Rcheck'\n"
        "check.mkdir()\n"
        "(check / '00check.log').write_text('* DONE\\nStatus: 1 ERROR\\n')\n"
        "raise SystemExit(23)\n"
    )
    code, *_ = invoke(tmp_path, wrapper_body=body)
    assert code == 23


@pytest.mark.parametrize("kind", ["missing", "duplicate", "symlink", "empty-rout"])
def test_native_check_fails_closed_on_ambiguous_outputs(tmp_path, kind):
    body = "import pathlib, sys\nparent = pathlib.Path(sys.argv[1])\n"
    if kind == "missing":
        body += "raise SystemExit(0)\n"
        environment_id = "mkl"
    elif kind == "duplicate":
        body += (
            "\nfor name in ('one.Rcheck', 'two.Rcheck'):\n"
            "    check = parent / name\n"
            "    check.mkdir()\n"
            "    (check / '00check.log').write_text('log')\n"
        )
        environment_id = "mkl"
    elif kind == "symlink":
        body += (
            "target = parent / 'target'\n"
            "target.write_text('log')\n"
            "check = parent / 'colocboost.Rcheck'\n"
            "check.mkdir()\n"
            "(check / '00check.log').symlink_to(target)\n"
        )
        environment_id = "mkl"
    else:
        body += (
            "check = parent / 'colocboost.Rcheck'\n"
            "(check / 'tests').mkdir(parents=True)\n"
            "(check / '00check.log').write_text('log')\n"
            "(check / 'tests' / 'testthat.Rout').write_text('')\n"
        )
        environment_id = "nosuggests"

    with pytest.raises(ValueError, match="00check.log|testthat.Rout"):
        invoke(tmp_path, environment_id=environment_id, wrapper_body=body)


def test_rejects_reused_work_or_check_library(tmp_path):
    tarball, metadata = make_contract(tmp_path)
    wrapper = make_executable(tmp_path / "bin" / "r-check")
    system_r = make_executable(tmp_path / "bin" / "R")
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    (work_dir / "stale").write_text("old")
    check_library = tmp_path / "check-lib"
    environment["R_LIBS_USER"] = os.fspath(check_library)

    with pytest.raises(ValueError, match="initially empty"):
        run_native_check(
            {
                "id": "mkl",
                "driver": "native-wrapper",
                "wrapper_path": os.fspath(wrapper),
                "system_r": os.fspath(system_r),
                "wrapper_input": "tarball-parent",
            },
            wrapper=wrapper,
            required_r_executable=system_r,
            tarball=tarball,
            metadata=metadata,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
            work_dir=work_dir,
            check_library=check_library,
            evidence_path=tmp_path / "evidence.json",
            log_output=tmp_path / "log",
            environment=environment,
        )


def test_native_check_disables_only_user_startup_files(tmp_path):
    body = (
        "import json, os, pathlib, sys\n"
        "parent = pathlib.Path(sys.argv[1])\n"
        "check = parent / 'colocboost.Rcheck'\n"
        "check.mkdir()\n"
        "(check / '00check.log').write_text('* DONE\\nStatus: OK\\n')\n"
        "(parent / 'startup.json').write_text(json.dumps({k: os.environ.get(k) for k in ('R_PROFILE_USER', 'R_ENVIRON_USER', 'R_PROFILE', 'R_ENVIRON')}))\n"
    )
    code, work, *_ = invoke(tmp_path, wrapper_body=body)
    assert code == 0
    captured = json.loads((work / "input" / "startup.json").read_text())
    assert captured["R_PROFILE_USER"] == os.devnull
    assert captured["R_ENVIRON_USER"] == os.devnull
    assert captured["R_PROFILE"] is None
    assert captured["R_ENVIRON"] is None
