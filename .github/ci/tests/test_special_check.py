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
sys.path.insert(0, os.fspath(Path(__file__).resolve().parent))
from test_runtime_evidence import valid_document  # noqa: E402


ACTIVE = {
    "clang-asan": ("clang-asan", "/usr/local/bin/r-check", "9732ef12761fd6ecd4631dd6ce0861fbbbd68fb33f6e7643745ec36de456b4f9", "/opt/R/devel-asan/bin/R", "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes"),
    "clang-ubsan": ("clang-ubsan", "/usr/local/bin/r-check", "9732ef12761fd6ecd4631dd6ce0861fbbbd68fb33f6e7643745ec36de456b4f9", "/opt/R/devel-asan/bin/R", "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes"),
    "donttest": ("donttest", "/usr/local/bin/r-check", "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090", "/opt/R/devel/bin/R", "--no-manual --no-build-vignettes"),
    "gcc-asan": ("gcc-asan", "/usr/local/bin/r-check", "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090", "/opt/R/devel/bin/R", "--no-manual --no-build-vignettes"),
    "gcc-ubsan": ("gcc-ubsan", "/usr/local/bin/r-check", "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090", "/opt/R/devel/bin/R", "--no-manual --no-build-vignettes"),
    "nold": ("nold", "/usr/local/bin/r-check", "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090", "/opt/R/devel-nold/bin/R", "--no-manual --no-build-vignettes"),
    "valgrind": ("valgrind", "/usr/local/bin/r-check", "79261a338b0a381a157cf2025ef40f389b3906fb4711d14b7a72314f5316cbc8", "/opt/R/devel-valgrind/bin/R", "--use-valgrind --extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes"),
    "vnu": ("vnu", "/usr/local/bin/r-check", "0de8ba373ec3bbe81b84122857d071dda4eee72e035d3180464607837c0bb089", "/opt/R/release/bin/R", "--no-manual --no-build-vignettes"),
}


def clean_log():
    return (
        "* checking whether package 'colocboost' can be installed ... [11s/11s] OK\n"
        "* checking tests ... OK\n"
        "* checking examples ... OK\n"
        "* DONE\nStatus: OK\n"
    )


def invoke(
    tmp_path, *, environment_id="r-release-linux-x86-64",
    profile="full-documentation", log=None, exit_code=0, files=None,
    native_mutation=None, runtime_environment=None,
):
    check_root = tmp_path / "input"
    check = check_root / "colocboost.Rcheck"
    check.mkdir(parents=True)
    if log is None:
        log = (CI_DIR / "tests/fixtures/release-full-doc-00check.log").read_text(
            encoding="utf-8"
        )
    (check / "00check.log").write_text(log, encoding="utf-8")
    if files:
        for relative, content in files.items():
            path = check / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content if isinstance(content, bytes) else content.encode())
    native = tmp_path / "native.json"
    if environment_id in ACTIVE:
        _, wrapper, wrapper_sha, system_r, check_args = ACTIVE[environment_id]
    else:
        wrapper = "/usr/local/bin/r-check"
        wrapper_sha = "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090"
        system_r = "/opt/R/release/bin/R" if environment_id == "r-release-linux-x86-64" else "/opt/R/devel/bin/R"
        check_args = "" if environment_id == "r-release-linux-x86-64" else "--no-manual --no-build-vignettes"
    native_document = {
        "wrapper_exit_code": exit_code,
        "wrapper_sha256": wrapper_sha,
        "check_args": check_args,
        "check_log": os.fspath(check / "00check.log"),
        "executable": {"path": wrapper, "sha256": wrapper_sha},
        "required_r_executable": {
            "path": system_r,
            "resolved_path": "/opt/R/4.6.1/bin/R" if environment_id == "vnu" else system_r,
        },
        "path_r_resolution": system_r,
    }
    if native_mutation:
        field, value = native_mutation
        native_document[field] = value
    native.write_text(json.dumps(native_document) + "\n", encoding="utf-8")
    runtime = tmp_path / "runtime.json"
    if environment_id in ACTIVE:
        runtime_document = valid_document(environment_id)
        if runtime_environment is not None:
            runtime_document["environment"].update(runtime_environment)
    else:
        runtime_document = {
            "environment_id": environment_id,
            "source_sha": SOURCE_SHA,
            "event_sha": EVENT_SHA,
            "tarball_sha256": DIGEST,
            "environment": runtime_environment or {},
        }
    runtime.write_text(json.dumps(runtime_document) + "\n", encoding="utf-8")
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


@pytest.mark.parametrize("environment_id", ["clang-asan", "clang-ubsan", "valgrind"])
def test_known_wrapper_no_rout_fail_defect_is_reconciled_only_after_clean_full_scan(tmp_path, environment_id):
    files = {"tests/testthat.Rout": "all tests passed\n"}
    if environment_id == "valgrind":
        files["tests/testthat.Rout"] = "==42== ERROR SUMMARY: 0 errors from 0 contexts\n"
    completed, output, github_output, _ = invoke(
        tmp_path, environment_id=environment_id, profile=environment_id,
        log=clean_log(), exit_code=1, files=files,
    )
    assert completed.returncode == 0, completed.stderr
    proof = json.loads(output.read_text())
    assert proof["raw_exit_code"] == 1
    assert proof["effective_exit_code"] == 0
    assert proof["reconciliation"] == "rhub-empty-rout-fail-grep"
    assert "effective_exit_code=0\n" in github_output.read_text()


@pytest.mark.parametrize(
    ("environment_id", "bad_text"),
    [
        ("clang-asan", "ERROR: AddressSanitizer: heap-use-after-free\n"),
        ("clang-ubsan", "UndefinedBehaviorSanitizer: runtime error\n"),
        ("gcc-asan", "SUMMARY: AddressSanitizer: buffer-overflow\n"),
        ("gcc-ubsan", "runtime error: signed integer overflow\n"),
        ("valgrind", "==42== Invalid read of size 8\n==42== ERROR SUMMARY: 1 errors from 1 contexts\n"),
    ],
)
def test_special_scan_rejects_diagnostics_anywhere_in_complete_tree(tmp_path, environment_id, bad_text):
    completed, *_ = invoke(
        tmp_path, environment_id=environment_id, profile=environment_id,
        log=clean_log(), files={"tests/nested/testthat.Rout": bad_text},
    )
    assert completed.returncode != 0


@pytest.mark.parametrize(
    ("environment_id", "hidden_diagnostic", "clean_transcript"),
    [
        ("clang-asan", "runtime error: hidden overflow\n", "clean\n"),
        (
            "valgrind",
            "==42== Invalid write of size 8\n",
            "==42== ERROR SUMMARY: 0 errors from 0 contexts\n",
        ),
    ],
)
def test_special_scan_covers_non_rout_files_in_the_entire_check_tree(
    tmp_path, environment_id, hidden_diagnostic, clean_transcript
):
    completed, *_ = invoke(
        tmp_path,
        environment_id=environment_id,
        profile=environment_id,
        log=clean_log(),
        files={
            "tests/testthat.Rout": clean_transcript,
            "nested/worker-diagnostic.log": hidden_diagnostic,
        },
    )

    assert completed.returncode != 0


@pytest.mark.parametrize("kind", ["raw-two", "rout-fail", "no-rout", "incomplete", "wrong-wrapper", "wrong-log-path"])
def test_reconciliation_fails_closed_without_exact_known_defect_signature(tmp_path, kind):
    files = {"tests/testthat.Rout": "clean\n"}
    exit_code = 2 if kind == "raw-two" else 1
    log = clean_log()
    mutation = None
    if kind == "rout-fail":
        files["tests/testthat.Rout.fail"] = "clean but failed\n"
    elif kind == "no-rout":
        files = {}
    elif kind == "incomplete":
        log = "* checking tests ... OK\n* DONE\nStatus: OK\n"
    elif kind == "wrong-wrapper":
        mutation = ("wrapper_sha256", "0" * 64)
    elif kind == "wrong-log-path":
        mutation = ("check_log", "/tmp/elsewhere/00check.log")
    completed, *_ = invoke(
        tmp_path, environment_id="clang-asan", profile="clang-asan",
        log=log, exit_code=exit_code, files=files, native_mutation=mutation,
    )
    assert completed.returncode != 0


def test_donttest_requires_policy_and_nonempty_examples_outputs(tmp_path):
    files = {
        "colocboost-Ex.R": "example(\"colocboost\")\n",
        "colocboost-Ex.Rout": "examples completed\n",
        "00_pkg_src/colocboost/man/colocboost.Rd": "\\name{colocboost}\n",
    }
    completed, output, *_ = invoke(
        tmp_path, environment_id="donttest", profile="donttest", log=clean_log(),
        files=files, runtime_environment={"_R_CHECK_DONTTEST_EXAMPLES_": "true"},
    )
    assert completed.returncode == 0, completed.stderr
    proof = json.loads(output.read_text())
    assert proof["donttest_policy_enabled"] is True
    assert proof["examples_executed"] is True
    assert proof["donttest_block_count"] == 0


@pytest.mark.parametrize("kind", ["missing-policy", "false-policy", "missing-source", "empty-rout", "skipped"])
def test_donttest_rejects_unproven_expanded_examples(tmp_path, kind):
    files = {
        "colocboost-Ex.R": "code\n",
        "colocboost-Ex.Rout": "output\n",
        "00_pkg_src/colocboost/man/colocboost.Rd": "\\name{colocboost}\n",
    }
    runtime_environment = {"_R_CHECK_DONTTEST_EXAMPLES_": "true"}
    log = clean_log()
    if kind == "missing-policy":
        runtime_environment = {"_R_CHECK_DONTTEST_EXAMPLES_": None}
    elif kind == "false-policy":
        runtime_environment = {"_R_CHECK_DONTTEST_EXAMPLES_": "false"}
    elif kind == "missing-source":
        del files["colocboost-Ex.R"]
    elif kind == "empty-rout":
        files["colocboost-Ex.Rout"] = ""
    else:
        log = clean_log().replace("* checking examples ... OK", "* checking examples ... SKIPPED")
    completed, *_ = invoke(
        tmp_path, environment_id="donttest", profile="donttest", log=log,
        files=files, runtime_environment=runtime_environment,
    )
    assert completed.returncode != 0


def test_donttest_rejects_a_source_package_with_donttest_blocks(tmp_path):
    completed, *_ = invoke(
        tmp_path,
        environment_id="donttest",
        profile="donttest",
        log=clean_log(),
        files={
            "colocboost-Ex.R": "code\n",
            "colocboost-Ex.Rout": "output\n",
            "00_pkg_src/colocboost/man/colocboost.Rd": (
                "\\name{colocboost}\n\\examples{\\donttest{stop('must run')}}\n"
            ),
        },
        runtime_environment={"_R_CHECK_DONTTEST_EXAMPLES_": "true"},
    )

    assert completed.returncode != 0


def test_valgrind_requires_a_clean_instrumented_transcript(tmp_path):
    completed, *_ = invoke(
        tmp_path, environment_id="valgrind", profile="valgrind", log=clean_log(),
        files={"tests/testthat.Rout": "ordinary R output only\n"},
    )
    assert completed.returncode != 0


@pytest.mark.parametrize(
    "executed_stage",
    [
        "* checking PDF version of manual ... OK\n",
        "* checking re-building of vignette outputs ... OK\n",
    ],
)
def test_active_no_documentation_profiles_reject_full_doc_stage(
    tmp_path, executed_stage
):
    completed, *_ = invoke(
        tmp_path, environment_id="gcc-asan", profile="gcc-asan",
        log=executed_stage + clean_log(),
        files={"tests/testthat.Rout": "clean\n"},
    )
    assert completed.returncode != 0


def test_no_build_vignettes_allows_incidental_vignette_checks_but_not_rebuild(tmp_path):
    log = clean_log().replace(
        "* DONE\n",
        "* checking package vignettes ... OK\n"
        "* checking running R code from vignettes ...\n"
        "  'Input_Data_Format.Rmd' using 'UTF-8'... OK\n"
        "* checking re-building of vignette outputs ... SKIPPED\n"
        "* DONE\n",
    )

    completed, output, *_ = invoke(
        tmp_path,
        environment_id="gcc-asan",
        profile="gcc-asan",
        log=log,
        files={"tests/testthat.Rout": "clean\n"},
    )

    assert completed.returncode == 0, completed.stderr
    proof = json.loads(output.read_text())
    assert proof["manual_executed"] is False
    assert proof["vignettes_executed"] is False


def test_real_task6_install_tests_and_examples_heading_shapes_are_accepted(tmp_path):
    headings = [
        "* checking whether package 'colocboost' can be installed ... [11s/11s] OK",
        "* checking examples ... OK",
        "* checking tests ... [28s/28s] OK",
    ]
    completed, *_ = invoke(
        tmp_path, environment_id="gcc-asan", profile="gcc-asan",
        log="\n".join(headings) + "\n* DONE\nStatus: OK\n",
        files={"tests/testthat.Rout": "clean\n"},
    )
    assert completed.returncode == 0, completed.stderr


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
