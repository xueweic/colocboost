import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "finalize_result.py"
sys.path.insert(0, str(CI_DIR))

from result_contract import validate_result  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
DIGEST = "c" * 64
ENVIRONMENT_ID = "r-release-linux-x86-64"


def unit_result(status="pass", **updates):
    document = {
        "schema_version": 1,
        "result_kind": "unit",
        "environment_id": ENVIRONMENT_ID,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": DIGEST,
        "status": status,
        "tests": {"pass": 5, "fail": 0, "error": 0, "warning": 0, "skip": 0},
    }
    document.update(updates)
    return document


def applicability_result():
    return {
        "schema_version": 1,
        "result_kind": "applicability",
        "environment_id": ENVIRONMENT_ID,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": DIGEST,
        "status": "not-applicable",
        "predicate": "no-native-source",
        "evidence": "The source tarball contains no package-owned native code.",
    }


def run_finalizer(
    result,
    *,
    expected_kind="unit",
    environment_id=ENVIRONMENT_ID,
    source_sha=SOURCE_SHA,
    event_sha=EVENT_SHA,
    digest=DIGEST,
    job_status="success",
    producer_outcome="success",
    timeout=10,
):
    return subprocess.run(
        [
            sys.executable,
            os.fspath(SCRIPT),
            "--result",
            os.fspath(result),
            "--expected-kind",
            expected_kind,
            "--environment-id",
            environment_id,
            "--source-sha",
            source_sha,
            "--event-sha",
            event_sha,
            "--tarball-sha256",
            digest,
            "--job-status",
            job_status,
            "--producer-outcome",
            producer_outcome,
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def write_json(path, document):
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return path


def assert_infrastructure_failure(path):
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_result(document)
    assert document["result_kind"] == "infrastructure"
    assert document["status"] == "fail"
    assert document["environment_id"] == ENVIRONMENT_ID
    assert document["source_sha"] == SOURCE_SHA
    assert document["event_sha"] == EVENT_SHA
    assert document["tarball_sha256"] == DIGEST
    assert document["classification"] == "result-finalization"
    assert document["message"]
    return document


def test_valid_green_result_is_preserved_only_when_both_outcomes_succeed(tmp_path):
    result = write_json(tmp_path / "result.json", unit_result())
    before = result.read_bytes()

    completed = run_finalizer(result)

    assert completed.returncode == 0, completed.stderr
    assert result.read_bytes() == before


def test_not_applicable_is_the_green_applicability_terminal_state(tmp_path):
    result = write_json(tmp_path / "result.json", applicability_result())
    before = result.read_bytes()

    completed = run_finalizer(result, expected_kind="applicability")

    assert completed.returncode == 0, completed.stderr
    assert result.read_bytes() == before


@pytest.mark.parametrize("status", ["fail", "missing", "cancelled", "uncovered"])
def test_valid_identity_matching_terminal_non_green_result_is_preserved(
    tmp_path, status
):
    result = write_json(tmp_path / "result.json", unit_result(status=status))
    before = result.read_bytes()

    completed = run_finalizer(
        result, job_status="failure", producer_outcome="skipped"
    )

    assert completed.returncode == 1
    assert result.read_bytes() == before


@pytest.mark.parametrize(
    ("job_status", "producer_outcome"),
    [
        ("failure", "success"),
        ("cancelled", "success"),
        ("success", "failure"),
        ("success", "cancelled"),
        ("success", "skipped"),
    ],
)
def test_green_result_with_non_success_outcome_is_replaced_by_infrastructure_failure(
    tmp_path, job_status, producer_outcome
):
    result = write_json(tmp_path / "result.json", unit_result())

    completed = run_finalizer(
        result, job_status=job_status, producer_outcome=producer_outcome
    )

    assert completed.returncode == 1
    document = assert_infrastructure_failure(result)
    assert job_status in document["message"] or producer_outcome in document["message"]


@pytest.mark.parametrize("case", ["missing", "provisional", "malformed", "duplicate"])
def test_missing_provisional_or_malformed_result_is_finalized_as_infrastructure(
    tmp_path, case
):
    result = tmp_path / "result.json"
    if case == "provisional":
        document = unit_result(status="provisional")
        document.pop("tests")
        write_json(result, document)
    elif case == "malformed":
        result.write_text("{not json", encoding="utf-8")
    elif case == "duplicate":
        result.write_text(
            '{"schema_version":1,"schema_version":1,"result_kind":"unit"}\n',
            encoding="utf-8",
        )

    completed = run_finalizer(result)

    assert completed.returncode == 1
    assert_infrastructure_failure(result)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("result_kind", "package-check"),
        ("environment_id", "wrong-environment"),
        ("source_sha", "d" * 40),
        ("event_sha", "e" * 40),
        ("tarball_sha256", "f" * 64),
    ],
)
def test_wrong_result_identity_is_replaced_by_infrastructure_failure(
    tmp_path, field, value
):
    document = unit_result()
    document[field] = value
    if field == "result_kind":
        document.pop("tests")
        document["check"] = {
            "errors": 0,
            "warnings": 0,
            "notes": 0,
            "log_path": "00check.log",
        }
    result = write_json(tmp_path / "result.json", document)

    completed = run_finalizer(result)

    assert completed.returncode == 1
    assert_infrastructure_failure(result)


def test_valid_infrastructure_failure_from_the_producer_is_preserved(tmp_path):
    document = {
        "schema_version": 1,
        "result_kind": "infrastructure",
        "environment_id": ENVIRONMENT_ID,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": DIGEST,
        "status": "fail",
        "classification": "runner",
        "message": "R was unavailable.",
    }
    result = write_json(tmp_path / "result.json", document)
    before = result.read_bytes()

    completed = run_finalizer(result, job_status="failure", producer_outcome="failure")

    assert completed.returncode == 1
    assert result.read_bytes() == before


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFO is unavailable")
@pytest.mark.parametrize("kind", ["symlink", "fifo"])
def test_non_regular_result_path_is_replaced_without_following_or_blocking(
    tmp_path, kind
):
    result = tmp_path / "result.json"
    target = tmp_path / "target.json"
    target.write_text("do not modify\n", encoding="utf-8")
    if kind == "symlink":
        result.symlink_to(target)
    else:
        os.mkfifo(result)

    completed = run_finalizer(result, timeout=5)

    assert completed.returncode == 1
    assert result.is_file() and not result.is_symlink()
    assert target.read_text(encoding="utf-8") == "do not modify\n"
    assert_infrastructure_failure(result)


def test_invalid_caller_identity_never_creates_a_result(tmp_path):
    result = tmp_path / "result.json"

    completed = run_finalizer(result, source_sha="A" * 40)

    assert completed.returncode == 2
    assert not result.exists()
