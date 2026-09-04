import json
import hashlib
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402
import run_driver as driver  # noqa: E402
from run_driver import capture_environment_evidence, run_driver  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
TARBALL_TOKEN = "{tarball}"
TARBALL_PARENT_TOKEN = "{tarball-parent}"


def make_executable(path, body="raise SystemExit(0)\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def make_artifact(tmp_path):
    tarball = tmp_path / "colocboost_1.0.0.tar.gz"
    tarball.write_bytes(b"one immutable source package")
    metadata = tmp_path / "source-metadata.json"
    create_metadata(tarball, metadata, source_sha=SOURCE_SHA, event_sha=EVENT_SHA)
    return tarball, metadata


def invoke(
    row,
    requested_driver,
    executable,
    argv,
    tarball,
    metadata,
    **kwargs,
):
    if requested_driver == "native-wrapper" and "required_r_executable" not in kwargs:
        system_r = make_executable(Path(tarball).parent / "system-r" / "R")
        selected_environment = os.environ.copy()
        selected_environment["PATH"] = os.fspath(system_r.parent)
        kwargs["required_r_executable"] = system_r
        kwargs.setdefault("environment", selected_environment)
    return run_driver(
        row,
        requested_driver=requested_driver,
        executable=executable,
        argv=argv,
        tarball=tarball,
        metadata=metadata,
        expected_source_sha=SOURCE_SHA,
        expected_event_sha=EVENT_SHA,
        **kwargs,
    )


def test_native_wrapper_preserves_absolute_path_and_argv_with_spaces(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    output = tmp_path / "captured argv.json"
    wrapper = make_executable(
        tmp_path / "directory with spaces" / "native wrapper",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]))\n",
    )
    argv = [str(output), "argument with spaces", TARBALL_TOKEN, "literal;not-shell"]

    exit_code = invoke(
        {"id": "mkl", "driver": "native-wrapper"},
        "native-wrapper",
        wrapper,
        argv,
        tarball,
        metadata,
    )

    assert exit_code == 0
    assert json.loads(output.read_text()) == [
        "argument with spaces",
        str(tarball.absolute()),
        "literal;not-shell",
    ]


def test_native_shell_wrapper_may_be_nonexecutable_when_hash_bound(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    marker = tmp_path / "called"
    wrapper = tmp_path / "bin" / "r-check"
    wrapper.parent.mkdir(parents=True)
    wrapper.write_text(f"#!/bin/sh\nprintf called > '{marker}'\n")
    system_r = make_executable(wrapper.parent / "R")
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)
    row = {
        "id": "gcc16",
        "driver": "native-wrapper",
        "wrapper_path": os.fspath(wrapper),
        "wrapper_sha256": hashlib.sha256(wrapper.read_bytes()).hexdigest(),
        "system_r": os.fspath(system_r),
        "wrapper_input": "tarball-parent",
    }
    assert run_driver(
        row, requested_driver="native-wrapper", executable=wrapper,
        argv=[TARBALL_PARENT_TOKEN], tarball=tarball, metadata=metadata,
        expected_source_sha=SOURCE_SHA, expected_event_sha=EVENT_SHA,
        required_r_executable=system_r, environment=environment,
    ) == 0
    assert marker.read_text() == "called"


def test_native_shell_wrapper_resolves_trusted_system_shell_symlink(tmp_path, monkeypatch):
    target = make_executable(tmp_path / "bin" / "dash")
    shell = tmp_path / "bin" / "sh"
    shell.symlink_to(target)
    wrapper = tmp_path / "bin" / "wrapper"
    wrapper.write_text("#!/bin/sh\n")
    monkeypatch.setattr(driver, "SYSTEM_SHELL_PATH", shell)

    command, invoked = driver._build_command("native-wrapper", wrapper, [], None)
    assert command == [os.fspath(target), os.fspath(wrapper)]
    assert invoked == target


def test_native_wrapper_binds_verified_single_tarball_parent_and_system_r(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    output = tmp_path / "captured-parent.json"
    wrapper = make_executable(
        tmp_path / "bin" / "wrapper",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]))\n",
    )
    system_r = make_executable(tmp_path / "bin" / "R")
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)
    evidence = tmp_path / "evidence.json"

    exit_code = run_driver(
        {
            "id": "mkl",
            "driver": "native-wrapper",
            "wrapper_path": os.fspath(wrapper),
            "system_r": os.fspath(system_r),
            "wrapper_input": "tarball-parent",
        },
        requested_driver="native-wrapper",
        executable=wrapper,
        argv=[os.fspath(output), TARBALL_PARENT_TOKEN],
        tarball=tarball,
        metadata=metadata,
        expected_source_sha=SOURCE_SHA,
        expected_event_sha=EVENT_SHA,
        required_r_executable=system_r,
        evidence_path=evidence,
        environment=environment,
    )

    assert exit_code == 0
    assert json.loads(output.read_text()) == [os.fspath(tarball.parent)]
    proof = json.loads(evidence.read_text())
    assert proof["required_r_executable"]["path"] == os.fspath(system_r)
    assert proof["path_r_resolution"] == os.fspath(system_r)


def test_native_wrapper_binds_declared_hash_and_check_args(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "bin" / "r-check")
    system_r = make_executable(tmp_path / "bin" / "R")
    wrapper_sha256 = hashlib.sha256(wrapper.read_bytes()).hexdigest()
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)
    environment["CHECK_ARGS"] = "--no-manual --no-build-vignettes"
    row = {
        "id": "r-devel-linux-x86-64-debian-clang",
        "driver": "native-wrapper",
        "wrapper_path": os.fspath(wrapper),
        "wrapper_sha256": wrapper_sha256,
        "system_r": os.fspath(system_r),
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
    }

    assert run_driver(
        row,
        requested_driver="native-wrapper",
        executable=wrapper,
        argv=[TARBALL_PARENT_TOKEN],
        tarball=tarball,
        metadata=metadata,
        expected_source_sha=SOURCE_SHA,
        expected_event_sha=EVENT_SHA,
        required_r_executable=system_r,
        environment=environment,
    ) == 0

    bad_hash = dict(row, wrapper_sha256="0" * 64)
    with pytest.raises(ValueError, match="wrapper.*sha256"):
        run_driver(
            bad_hash,
            requested_driver="native-wrapper",
            executable=wrapper,
            argv=[TARBALL_PARENT_TOKEN],
            tarball=tarball,
            metadata=metadata,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
            required_r_executable=system_r,
            environment=environment,
        )

    wrong_args_environment = dict(environment, CHECK_ARGS="--no-manual")
    with pytest.raises(ValueError, match="CHECK_ARGS"):
        run_driver(
            row,
            requested_driver="native-wrapper",
            executable=wrapper,
            argv=[TARBALL_PARENT_TOKEN],
            tarball=tarball,
            metadata=metadata,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
            required_r_executable=system_r,
            environment=wrong_args_environment,
        )


@pytest.mark.parametrize("extra_kind", ["regular", "symlink", "fifo"])
def test_tarball_parent_rejects_every_ambiguous_or_nonregular_tarball(tmp_path, extra_kind):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "bin" / "wrapper")
    system_r = make_executable(tmp_path / "bin" / "R")
    extra = tmp_path / "extra.tar.gz"
    if extra_kind == "regular":
        extra.write_bytes(b"extra")
    elif extra_kind == "symlink":
        extra.symlink_to(tarball)
    else:
        os.mkfifo(extra)
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(system_r.parent)

    with pytest.raises(ValueError, match="exactly one regular non-symlink"):
        run_driver(
            {
                "id": "mkl",
                "driver": "native-wrapper",
                "wrapper_path": os.fspath(wrapper),
                "system_r": os.fspath(system_r),
                "wrapper_input": "tarball-parent",
            },
            requested_driver="native-wrapper",
            executable=wrapper,
            argv=[TARBALL_PARENT_TOKEN],
            tarball=tarball,
            metadata=metadata,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
            required_r_executable=system_r,
            environment=environment,
        )


def test_native_wrapper_rejects_path_r_that_is_not_the_required_system_r(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "wrapper-bin" / "wrapper")
    required_r = make_executable(tmp_path / "system" / "R")
    shadow_r = make_executable(tmp_path / "pixi" / "R")
    environment = os.environ.copy()
    environment["PATH"] = os.fspath(shadow_r.parent)

    with pytest.raises(ValueError, match="PATH.*required system R"):
        run_driver(
            {
                "id": "mkl",
                "driver": "native-wrapper",
                "wrapper_path": os.fspath(wrapper),
                "system_r": os.fspath(required_r),
                "wrapper_input": "tarball-parent",
            },
            requested_driver="native-wrapper",
            executable=wrapper,
            argv=[TARBALL_PARENT_TOKEN],
            tarball=tarball,
            metadata=metadata,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
            required_r_executable=required_r,
            environment=environment,
        )


@pytest.mark.parametrize("kind", ["missing", "relative", "non-executable", "directory"])
def test_driver_rejects_invalid_executable_paths_before_invocation(tmp_path, kind):
    tarball, metadata = make_artifact(tmp_path)
    marker = tmp_path / "invoked"
    if kind == "missing":
        executable = tmp_path / "missing-wrapper"
    elif kind == "relative":
        executable = Path("relative-wrapper")
    elif kind == "directory":
        executable = tmp_path / "wrapper-directory"
        executable.mkdir()
    else:
        executable = tmp_path / "non-executable-wrapper"
        executable.write_text("do not run")

    with pytest.raises(ValueError, match="absolute|regular|executable"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            executable,
            [str(marker), TARBALL_TOKEN],
            tarball,
            metadata,
        )
    assert not marker.exists()


def test_driver_rejects_symlink_executable(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    target = make_executable(tmp_path / "real-wrapper", "raise SystemExit(0)\n")
    symlink = tmp_path / "linked-wrapper"
    symlink.symlink_to(target)

    with pytest.raises(ValueError, match="regular.*non-symlink"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            symlink,
            [TARBALL_TOKEN],
            tarball,
            metadata,
        )


def test_driver_rejects_manifest_and_requested_driver_mismatch(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "wrapper", "raise SystemExit(0)\n")

    with pytest.raises(ValueError, match="driver mismatch"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "r-binary",
            wrapper,
            [TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="r-cmd",
        )


def test_driver_rejects_unknown_manifest_driver(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "wrapper", "raise SystemExit(0)\n")

    with pytest.raises(ValueError, match="unsupported driver"):
        invoke(
            {"id": "bad", "driver": "shell-command"},
            "shell-command",
            wrapper,
            [TARBALL_TOKEN],
            tarball,
            metadata,
        )


def test_native_wrapper_exit_status_is_propagated(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "wrapper", "raise SystemExit(23)\n")

    assert (
        invoke(
            {"id": "valgrind", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            [TARBALL_TOKEN],
            tarball,
            metadata,
        )
        == 23
    )


def test_invalid_r_binary_never_falls_back_to_path_r(tmp_path, monkeypatch):
    tarball, metadata = make_artifact(tmp_path)
    path_bin = tmp_path / "path-bin"
    marker = tmp_path / "path-r-was-invoked"
    make_executable(
        path_bin / "R",
        f"import pathlib\npathlib.Path({str(marker)!r}).write_text('bad')\n",
    )
    monkeypatch.setenv("PATH", str(path_bin))

    with pytest.raises(ValueError, match="regular"):
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            tmp_path / "missing" / "R",
            ["script.R", TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="rscript",
        )
    assert not marker.exists()


def test_rscript_uses_verified_sibling_and_preserves_argv(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    output = tmp_path / "rscript-argv.json"
    r_binary = make_executable(
        tmp_path / "R installation with spaces" / "bin" / "R",
        "raise SystemExit(99)\n",
    )
    make_executable(
        r_binary.with_name("Rscript"),
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]))\n",
    )
    argv = [str(output), "script with spaces.R", TARBALL_TOKEN, "arg with spaces"]

    assert (
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            argv,
            tarball,
            metadata,
            r_entrypoint="rscript",
        )
        == 0
    )
    assert json.loads(output.read_text()) == [
        "script with spaces.R",
        str(tarball.absolute()),
        "arg with spaces",
    ]


def test_rscript_evidence_keeps_selected_r_and_invoked_sibling_identities(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    evidence_path = tmp_path / "rscript-evidence.json"
    r_binary = make_executable(
        tmp_path / "R-home" / "bin" / "R", "raise SystemExit(99)\n"
    )
    rscript = make_executable(
        r_binary.with_name("Rscript"), "raise SystemExit(0)\n"
    )

    assert (
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            ["script.R", TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="rscript",
            evidence_path=evidence_path,
        )
        == 0
    )
    evidence = json.loads(evidence_path.read_text())
    assert evidence["executable"]["path"] == str(r_binary)
    assert evidence["invoked_executable"]["path"] == str(rscript)


def test_r_cmd_uses_selected_r_executable(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    output = tmp_path / "r-cmd-argv.json"
    r_binary = make_executable(
        tmp_path / "R-home" / "bin" / "R",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[2]).write_text(json.dumps(sys.argv[1:]))\n",
    )

    assert (
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            [str(output), "check", "--as-cran", TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="r-cmd",
        )
        == 0
    )
    assert json.loads(output.read_text()) == [
        "CMD",
        str(output),
        "check",
        "--as-cran",
        str(tarball),
    ]


@pytest.mark.parametrize("sibling_state", ["missing", "non-executable"])
def test_rscript_requires_executable_sibling(tmp_path, sibling_state):
    tarball, metadata = make_artifact(tmp_path)
    r_binary = make_executable(tmp_path / "R-home" / "bin" / "R", "raise SystemExit(0)\n")
    if sibling_state == "non-executable":
        r_binary.with_name("Rscript").write_text("do not run")

    with pytest.raises(ValueError, match="Rscript.*regular|Rscript.*executable"):
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            ["script.R", TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="rscript",
        )


def test_r_binary_rejects_wrong_entrypoint_and_non_r_basename(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    not_r = make_executable(tmp_path / "python", "raise SystemExit(0)\n")

    with pytest.raises(ValueError, match="R or R.exe"):
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            not_r,
            [TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="r-cmd",
        )

    r_binary = make_executable(tmp_path / "R", "raise SystemExit(0)\n")
    with pytest.raises(ValueError, match="r_entrypoint"):
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            [TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="shell",
        )


def test_native_wrapper_rejects_r_entrypoint_switch(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    wrapper = make_executable(tmp_path / "wrapper", "raise SystemExit(0)\n")

    with pytest.raises(ValueError, match="r_entrypoint"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            [TARBALL_TOKEN],
            tarball,
            metadata,
            r_entrypoint="r-cmd",
        )


def test_corrupted_tarball_blocks_driver_execution(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    marker = tmp_path / "invoked"
    wrapper = make_executable(
        tmp_path / "wrapper",
        f"import pathlib\npathlib.Path({str(marker)!r}).write_text('bad')\n",
    )
    tarball.write_bytes(b"changed after prepare")

    with pytest.raises(ValueError, match="size|sha256"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            [TARBALL_TOKEN],
            tarball,
            metadata,
        )
    assert not marker.exists()


def test_environment_evidence_is_whitelisted_and_records_executable_identity(
    tmp_path, monkeypatch
):
    executable = make_executable(tmp_path / "wrapper", "raise SystemExit(0)\n")
    monkeypatch.setenv("PATH", "/approved/bin")
    monkeypatch.setenv("R_HOME", "/approved/R")
    monkeypatch.setenv("R_LIBS", "/approved/R/library")
    monkeypatch.setenv("R_LIBS_USER", "/approved/user/library")
    monkeypatch.setenv("R_LIBS_SITE", "/approved/site/library")
    monkeypatch.setenv("LD_LIBRARY_PATH", "/approved/lib")
    monkeypatch.setenv("DYLD_LIBRARY_PATH", "/approved/dyld")
    monkeypatch.setenv("CC", "clang")
    monkeypatch.setenv("CXX", "clang++")
    monkeypatch.setenv("FC", "flang")
    monkeypatch.setenv("LANG", "C.UTF-8")
    monkeypatch.setenv("LC_ALL", "C")
    monkeypatch.delenv("LC_CTYPE", raising=False)
    monkeypatch.delenv("F77", raising=False)
    monkeypatch.delenv("LIBPATH", raising=False)
    monkeypatch.delenv("SHLIB_PATH", raising=False)
    monkeypatch.setenv("SUPER_SECRET_TOKEN", "must-not-leak")

    evidence = capture_environment_evidence(executable)

    assert evidence["environment"] == {
        "PATH": "/approved/bin",
        "R_HOME": "/approved/R",
        "R_LIBS": "/approved/R/library",
        "R_LIBS_USER": "/approved/user/library",
        "R_LIBS_SITE": "/approved/site/library",
        "R_PROFILE_USER": None,
        "R_ENVIRON_USER": None,
        "LD_LIBRARY_PATH": "/approved/lib",
        "DYLD_LIBRARY_PATH": "/approved/dyld",
        "LIBPATH": None,
        "SHLIB_PATH": None,
        "CC": "clang",
        "CXX": "clang++",
        "FC": "flang",
        "F77": None,
        "LANG": "C.UTF-8",
        "LC_ALL": "C",
        "LC_CTYPE": None,
    }
    assert evidence["executable"]["path"] == str(executable)
    assert evidence["executable"]["resolved_path"] == str(executable.resolve())
    assert len(evidence["executable"]["sha256"]) == 64
    assert "SUPER_SECRET_TOKEN" not in json.dumps(evidence)
    assert evidence["dynamic_libraries"]["method"] in {
        "proc-self-maps",
        "unsupported",
    }


def test_evidence_is_written_before_driver_runs(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    evidence_path = tmp_path / "evidence.json"
    wrapper = make_executable(
        tmp_path / "wrapper",
        "import json, pathlib, sys\n"
        "evidence = pathlib.Path(sys.argv[1])\n"
        "raise SystemExit(0 if json.loads(evidence.read_text()) else 31)\n",
    )

    assert (
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            [str(evidence_path), TARBALL_TOKEN],
            tarball,
            metadata,
            evidence_path=evidence_path,
        )
        == 0
    )


def test_driver_cli_enforces_manifest_wrapper_identity(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    output = tmp_path / "cli argv.json"
    evidence = tmp_path / "cli evidence.json"
    wrapper = make_executable(
        tmp_path / "cli wrapper with spaces",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]))\n"
        "raise SystemExit(19)\n",
    )
    system_r = make_executable(tmp_path / "system R" / "bin" / "R")
    environment = os.environ.copy()
    environment["PATH"] = os.pathsep.join(
        [os.fspath(system_r.parent), environment.get("PATH", "")]
    )
    command_argv = [str(output), "arg with spaces", TARBALL_TOKEN]

    completed = subprocess.run(
        [
            sys.executable,
            str(CI_DIR / "run_driver.py"),
            "--manifest",
            str(CI_DIR / "check-matrix.yml"),
            "--environment-id",
            "noomp",
            "--driver",
            "native-wrapper",
            "--executable",
            str(wrapper),
            "--tarball",
            str(tarball),
            "--metadata",
            str(metadata),
            "--source-sha",
            SOURCE_SHA,
            "--event-sha",
            EVENT_SHA,
            "--evidence",
            str(evidence),
            "--required-r-executable",
            str(system_r),
            "--",
            *command_argv,
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert completed.returncode == 2
    assert "does not match manifest wrapper_path" in completed.stderr
    assert not output.exists()
    assert not evidence.exists()


def test_driver_cli_rejects_unknown_environment_without_invocation(tmp_path):
    tarball, metadata = make_artifact(tmp_path)
    marker = tmp_path / "invoked"
    wrapper = make_executable(
        tmp_path / "wrapper",
        f"import pathlib\npathlib.Path({str(marker)!r}).write_text('bad')\n",
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(CI_DIR / "run_driver.py"),
            "--manifest",
            str(CI_DIR / "check-matrix.yml"),
            "--environment-id",
            "not-in-manifest",
            "--driver",
            "native-wrapper",
            "--executable",
            str(wrapper),
            "--tarball",
            str(tarball),
            "--metadata",
            str(metadata),
            "--source-sha",
            SOURCE_SHA,
            "--event-sha",
            EVENT_SHA,
            "--evidence",
            str(tmp_path / "unused-evidence.json"),
            "--",
            TARBALL_TOKEN,
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert "environment_id" in completed.stderr
    assert not marker.exists()


@pytest.mark.parametrize(
    "argv_factory",
    [
        lambda tarball_b: [str(tarball_b)],
        lambda tarball_b: [TARBALL_TOKEN, TARBALL_TOKEN],
        lambda tarball_b: [f"--tarball={TARBALL_TOKEN}"],
    ],
    ids=(
        "different-tarball-no-placeholder",
        "multiple-placeholders",
        "embedded-placeholder",
    ),
)
def test_native_wrapper_rejects_unbound_or_ambiguous_tarball_argv(
    tmp_path, argv_factory
):
    tarball_a, metadata_a = make_artifact(tmp_path)
    tarball_b = tmp_path / "different-package.tar.gz"
    tarball_b.write_bytes(b"unverified source package")
    marker = tmp_path / "invoked"
    wrapper = make_executable(
        tmp_path / "wrapper",
        f"import pathlib\npathlib.Path({str(marker)!r}).write_text('bad')\n",
    )

    with pytest.raises(ValueError, match="exactly one literal.*\\{tarball\\}"):
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            argv_factory(tarball_b),
            tarball_a,
            metadata_a,
        )
    assert not marker.exists()


def test_native_wrapper_substitutes_verified_a_and_never_substitutes_b(tmp_path):
    tarball_a, metadata_a = make_artifact(tmp_path)
    tarball_b = tmp_path / "different-package.tar.gz"
    tarball_b.write_bytes(b"unverified source package")
    output = tmp_path / "received.json"
    wrapper = make_executable(
        tmp_path / "wrapper",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]))\n",
    )

    assert (
        invoke(
            {"id": "mkl", "driver": "native-wrapper"},
            "native-wrapper",
            wrapper,
            [str(output), TARBALL_TOKEN],
            tarball_a,
            metadata_a,
        )
        == 0
    )
    received = json.loads(output.read_text())
    assert received == [str(tarball_a.absolute())]
    assert str(tarball_b) not in received


def test_r_binary_substitutes_absolute_lexical_verified_tarball_path(
    tmp_path, monkeypatch
):
    tarball_a, metadata_a = make_artifact(tmp_path)
    output = tmp_path / "r-received.json"
    r_binary = make_executable(
        tmp_path / "R-home" / "bin" / "R",
        "import json, pathlib, sys\n"
        "pathlib.Path(sys.argv[2]).write_text(json.dumps(sys.argv[3:]))\n",
    )
    monkeypatch.chdir(tmp_path)

    assert (
        invoke(
            {"id": "r-release", "driver": "r-binary"},
            "r-binary",
            r_binary,
            [str(output), TARBALL_TOKEN],
            Path(tarball_a.name),
            Path(metadata_a.name),
            r_entrypoint="r-cmd",
        )
        == 0
    )
    assert json.loads(output.read_text()) == [str(tarball_a)]
