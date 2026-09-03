import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402
from run_r_binary_check import _same_portable_path, run_r_binary_check  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def make_executable(tmp_path):
    path = tmp_path / "selected R with spaces" / (
        "R.exe" if os.name == "nt" else "R"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"fake R executable\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def make_contract(tmp_path):
    source = tmp_path / "source"
    source.mkdir(exist_ok=True)
    tarball = source / "colocboost_1.0.9.tar.gz"
    tarball.write_bytes(b"verified immutable source bytes")
    metadata = source / "source-metadata.json"
    document = create_metadata(
        tarball, metadata, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
    )
    return tarball, metadata, document


def invoke(tmp_path, monkeypatch, *, exit_code=0, output_mode="one", **overrides):
    tarball, metadata, metadata_document = make_contract(tmp_path)
    r_binary = make_executable(tmp_path)
    check_library = tmp_path / "check library"
    check_library.mkdir(exist_ok=True)
    (check_library / "prepared-dependencies").write_text("bound", encoding="utf-8")
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        result_dir = Path(kwargs["cwd"])
        if output_mode != "missing":
            names = ("colocboost.Rcheck",) if output_mode != "duplicate" else (
                "one.Rcheck",
                "two.Rcheck",
            )
            for name in names:
                check = result_dir / name
                check.mkdir()
                log = check / "00check.log"
                if output_mode == "symlink":
                    target = result_dir / "target.log"
                    target.write_text("target\n", encoding="utf-8")
                    try:
                        log.symlink_to(target)
                    except OSError:
                        pytest.skip("symlinks unavailable")
                else:
                    log.write_text(
                        "* DONE\nStatus: OK\n"
                        if exit_code == 0
                        else "Status: 1 ERROR\n",
                        encoding="utf-8",
                    )
        return subprocess.CompletedProcess(command, exit_code)

    monkeypatch.setattr(subprocess, "run", fake_run)
    arguments = {
        "manifest_row": {"id": "windows-release", "driver": "r-binary"},
        "environment_id": "windows-release",
        "r_executable": r_binary,
        "tarball": tarball,
        "metadata": metadata,
        "expected_source_sha": SOURCE_SHA,
        "expected_event_sha": EVENT_SHA,
        "work_dir": tmp_path / "check work",
        "check_library": check_library,
        "evidence_path": tmp_path / "published" / "check-proof.json",
        "log_output": tmp_path / "published" / "00check.log",
        "environment": {"PATH": os.environ.get("PATH", "")},
    }
    arguments.update(overrides)
    code = run_r_binary_check(**arguments)
    return code, calls, arguments, metadata_document


def test_runs_exact_r_cmd_check_and_atomically_publishes_bound_proof(
    tmp_path, monkeypatch
):
    code, calls, arguments, metadata = invoke(tmp_path, monkeypatch)

    assert code == 0
    assert len(calls) == 1
    command, options = calls[0]
    assert command[:6] == [
        os.fspath(arguments["r_executable"]),
        "CMD",
        "check",
        "--as-cran",
        "--no-manual",
        "--no-build-vignettes",
    ]
    assert len(command) == 7 and command[-1].endswith(".tar.gz")
    assert options["shell"] is False
    assert options["check"] is False
    assert options["env"]["R_LIBS_USER"] == os.fspath(arguments["check_library"])
    assert options["env"]["R_PROFILE_USER"] == os.devnull
    assert options["env"]["R_ENVIRON_USER"] == os.devnull
    proof = json.loads(Path(arguments["evidence_path"]).read_text(encoding="utf-8"))
    assert proof["kind"] == "r-binary-check-proof"
    assert proof["status"] == "passed"
    assert proof["environment_id"] == "windows-release"
    assert proof["source_sha"] == SOURCE_SHA
    assert proof["event_sha"] == EVENT_SHA
    assert proof["tarball"]["sha256"] == metadata["sha256"]
    assert proof["r_executable"]["path"] == os.fspath(arguments["r_executable"])
    assert proof["exit_code"] == 0
    assert proof["startup_environment"] == {
        "R_LIBS_USER": os.fspath(arguments["check_library"]),
        "R_PROFILE_USER": os.devnull,
        "R_ENVIRON_USER": os.devnull,
    }
    assert Path(arguments["log_output"]).read_text(encoding="utf-8") == (
        "* DONE\nStatus: OK\n"
    )


def test_preserves_nonzero_r_check_exit_when_contract_output_exists(
    tmp_path, monkeypatch
):
    code, _, arguments, _ = invoke(tmp_path, monkeypatch, exit_code=19)
    proof = json.loads(Path(arguments["evidence_path"]).read_text(encoding="utf-8"))
    assert code == 19
    assert proof["exit_code"] == 19
    assert proof["status"] == "failed"
    assert "1 ERROR" in Path(arguments["log_output"]).read_text(encoding="utf-8")


def test_full_documentation_requires_an_explicit_opt_in(tmp_path, monkeypatch):
    _, calls, _, _ = invoke(tmp_path, monkeypatch, allow_full_documentation=True)
    command = calls[0][0]
    assert "--as-cran" not in command
    assert "--no-manual" not in command
    assert "--no-build-vignettes" not in command


@pytest.mark.parametrize("output_mode", ["missing", "duplicate", "symlink"])
def test_rejects_missing_ambiguous_or_symlinked_check_log(
    tmp_path, monkeypatch, output_mode
):
    with pytest.raises(ValueError, match="exactly one|regular non-symlink"):
        invoke(tmp_path, monkeypatch, output_mode=output_mode)


def test_rejects_reused_work_directory(tmp_path, monkeypatch):
    path = tmp_path / "check work"
    path.mkdir()
    (path / "stale").write_text("old", encoding="utf-8")
    with pytest.raises(ValueError, match="initially empty"):
        invoke(tmp_path, monkeypatch, work_dir=path)


@pytest.mark.parametrize("mode", ["missing", "empty", "symlink"])
def test_requires_a_nonempty_prepared_check_library(tmp_path, monkeypatch, mode):
    library = tmp_path / "prepared library override"
    if mode == "empty":
        library.mkdir()
    elif mode == "symlink":
        target = tmp_path / "library target"
        target.mkdir()
        (target / "dependency").write_text("installed", encoding="utf-8")
        try:
            library.symlink_to(target, target_is_directory=True)
        except OSError:
            pytest.skip("symlinks unavailable")
    with pytest.raises(ValueError, match="prepared check library"):
        invoke(tmp_path, monkeypatch, check_library=library)


def test_rejects_environment_manifest_and_metadata_identity_mismatches(
    tmp_path, monkeypatch
):
    with pytest.raises(ValueError, match="environment ID"):
        invoke(tmp_path, monkeypatch, environment_id="windows-devel")
    with pytest.raises(ValueError, match="environment ID"):
        invoke(tmp_path, monkeypatch, environment_id="windows-release\ninjected")
    with pytest.raises(ValueError, match="source_sha"):
        invoke(tmp_path, monkeypatch, expected_source_sha="c" * 40)


def test_rejects_symlinked_r_and_overlapping_directories(tmp_path, monkeypatch):
    real_r = make_executable(tmp_path)
    linked_r = tmp_path / "linked-R"
    try:
        linked_r.symlink_to(real_r)
    except OSError:
        pytest.skip("symlinks unavailable")
    with pytest.raises(ValueError, match="R executable"):
        invoke(tmp_path, monkeypatch, r_executable=linked_r)

    library = tmp_path / "shared"
    library.mkdir()
    (library / "dependency").write_text("installed", encoding="utf-8")
    work = library / "work"
    with pytest.raises(ValueError, match="independent"):
        invoke(tmp_path, monkeypatch, work_dir=work, check_library=library)


def test_windows_path_spelling_is_compared_portably():
    assert _same_portable_path(r"C:\R\bin\R.exe", "C:/R/bin/R.exe")
    assert not _same_portable_path(r"D:\R\bin\R.exe", "C:/R/bin/R.exe")
