from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


CI_DIR = Path(__file__).resolve().parents[1]
ROOT = CI_DIR.parents[1]
MANIFEST_PATH = CI_DIR / "check-matrix.yml"
POLICY_PATH = CI_DIR / "check-policy.yml"
AGGREGATOR = CI_DIR / "aggregate_results.py"
UNIT_ADAPTER = CI_DIR / "adapt_unit_result.py"
sys.path.insert(0, str(CI_DIR))

from result_contract import validate_result  # noqa: E402
from validate_manifest import (  # noqa: E402
    expected_result_keys,
    load_manifest,
    validate_manifest,
    validate_policy,
)


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
TARBALL_SHA256 = "c" * 64


@pytest.fixture(scope="module")
def manifest():
    return validate_manifest(load_manifest(MANIFEST_PATH))


@pytest.fixture(scope="module")
def policy():
    return validate_policy(load_manifest(POLICY_PATH))


def package_result(environment_id: str, *, status="pass", **counts):
    document = {
        "schema_version": 1,
        "result_kind": "package-check",
        "environment_id": environment_id,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": TARBALL_SHA256,
        "status": status,
        "check": {
            "errors": counts.get("errors", 0),
            "warnings": counts.get("warnings", 0),
            "notes": counts.get("notes", 0),
            "log_path": f"logs/{environment_id}/00check.log",
        },
    }
    return document


def applicability_result(environment_id: str, *, status="not-applicable"):
    return {
        "schema_version": 1,
        "result_kind": "applicability",
        "environment_id": environment_id,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": TARBALL_SHA256,
        "status": status,
        "predicate": "no-native-source | exact",
        "evidence": "first line\nsecond line",
    }


def unit_result(environment_id: str, *, status="pass", skip=None, **counts):
    if skip is None:
        skip = 5 if environment_id == "nosuggests" else 0
    return {
        "schema_version": 1,
        "result_kind": "unit",
        "environment_id": environment_id,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": TARBALL_SHA256,
        "status": status,
        "tests": {
            "pass": counts.get("pass_count", 3),
            "fail": counts.get("fail", 0),
            "error": counts.get("error", 0),
            "warning": counts.get("warning", 0),
            "skip": skip,
        },
    }


def all_green_documents(manifest):
    rows = {row["id"]: row for row in manifest["coverage"]}
    documents = []
    for result_kind, environment_id in expected_result_keys(manifest):
        if result_kind == "package-check":
            documents.append(package_result(environment_id))
        elif result_kind == "applicability":
            assert rows[environment_id]["state"] == "not-applicable"
            documents.append(applicability_result(environment_id))
        else:
            documents.append(unit_result(environment_id))
    assert len(documents) == 60
    return documents


def write_documents(root: Path, documents):
    paths = []
    for index, document in enumerate(documents):
        path = root / f"artifact-{index:02d}" / "result.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(document) + "\n", encoding="utf-8")
        paths.append(path)
    return paths


def run_aggregator(
    tmp_path: Path,
    manifest,
    *,
    documents=None,
    results_root=None,
    expected_source_sha=SOURCE_SHA,
    expected_event_sha=EVENT_SHA,
    expected_tarball_sha256=TARBALL_SHA256,
):
    if results_root is None:
        results_root = tmp_path / "results"
    if documents is not None:
        paths = write_documents(results_root, documents)
    else:
        paths = []
    summary = tmp_path / "summary.md"
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(AGGREGATOR),
            f"--manifest={MANIFEST_PATH}",
            f"--results-root={results_root}",
            f"--expected-source-sha={expected_source_sha}",
            f"--expected-event-sha={expected_event_sha}",
            f"--expected-tarball-sha256={expected_tarball_sha256}",
            f"--summary={summary}",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return completed, summary, paths


def assert_failed_with_complete_table(completed, summary):
    assert completed.returncode != 0, completed.stdout + completed.stderr
    text = summary.read_text(encoding="utf-8")
    table_rows = [line for line in text.splitlines() if line.startswith("| ")]
    assert len(table_rows) == 62  # header, separator, and all 60 stable rows
    return text


def test_all_60_results_pass_and_markdown_is_escaped(tmp_path, manifest):
    documents = all_green_documents(manifest)
    completed, summary, _ = run_aggregator(tmp_path, manifest, documents=documents)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    text = summary.read_text(encoding="utf-8")
    assert "Overall: PASS" in text
    assert "60/60" in text
    assert "no-native-source \\| exact" in text
    assert "first line<br>second line" in text
    assert len([line for line in text.splitlines() if line.startswith("| ")]) == 62


def test_non_result_sidecars_are_ignored(tmp_path, manifest):
    documents = all_green_documents(manifest)
    results_root = tmp_path / "results"
    write_documents(results_root, documents)
    (results_root / "diagnostic.json").write_text("not json", encoding="utf-8")

    completed, summary, _ = run_aggregator(
        tmp_path, manifest, documents=None, results_root=results_root
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Overall: PASS" in summary.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "mutation, expected_text",
    [
        ("missing", "missing"),
        ("duplicate", "duplicate"),
        ("schema-invalid", "schema"),
        ("wrong-source", "source_sha"),
        ("wrong-event", "event_sha"),
        ("wrong-tarball", "tarball_sha256"),
        ("cancelled", "cancelled"),
        ("skipped", "skipped"),
        ("fail", "fail"),
        ("uncovered", "uncovered"),
        ("provisional", "provisional"),
        ("unknown", "unexpected"),
        ("infrastructure", "infrastructure"),
    ],
)
def test_every_artifact_failure_mode_is_fail_closed(
    tmp_path, manifest, mutation, expected_text
):
    documents = all_green_documents(manifest)
    target = next(
        document
        for document in documents
        if document["result_kind"] == "package-check"
    )
    extras = []

    if mutation == "missing":
        documents.remove(target)
    elif mutation == "duplicate":
        extras.append(copy.deepcopy(target))
    elif mutation == "schema-invalid":
        target["unexpected"] = True
    elif mutation == "wrong-source":
        target["source_sha"] = "d" * 40
    elif mutation == "wrong-event":
        target["event_sha"] = "d" * 40
    elif mutation == "wrong-tarball":
        target["tarball_sha256"] = "d" * 64
    elif mutation in {"cancelled", "fail", "uncovered"}:
        target["status"] = mutation
    elif mutation == "skipped":
        target["status"] = "skipped"
    elif mutation == "provisional":
        target["status"] = "provisional"
        target["tarball_sha256"] = None
        target.pop("check")
    elif mutation == "unknown":
        extra = package_result("unknown-environment")
        extras.append(extra)
    elif mutation == "infrastructure":
        extras.append(
            {
                "schema_version": 1,
                "result_kind": "infrastructure",
                "environment_id": target["environment_id"],
                "source_sha": SOURCE_SHA,
                "event_sha": EVENT_SHA,
                "tarball_sha256": TARBALL_SHA256,
                "status": "fail",
                "classification": "runner",
                "message": "runner failed",
            }
        )

    completed, summary, paths = run_aggregator(
        tmp_path, manifest, documents=documents + extras
    )
    if mutation == "duplicate":
        # The ordinary writer already puts the duplicate under a separate directory.
        assert len(paths) == 61
    text = assert_failed_with_complete_table(completed, summary)
    assert expected_text.lower() in text.lower()


def test_missing_results_root_and_missing_prepare_digest_render_all_rows(
    tmp_path, manifest
):
    missing_root = tmp_path / "does-not-exist"
    completed, summary, _ = run_aggregator(
        tmp_path,
        manifest,
        results_root=missing_root,
        expected_tarball_sha256="",
    )

    text = assert_failed_with_complete_table(completed, summary)
    assert "results root" in text.lower()
    assert "tarball" in text.lower()
    assert "0/60" in text


@pytest.mark.parametrize(
    "kind, mutation, expected_text",
    [
        ("package-check", {"errors": 1}, "errors"),
        ("package-check", {"warnings": 1}, "warnings"),
        ("unit", {"pass_count": 0}, "pass"),
        ("unit", {"fail": 1}, "fail"),
        ("unit", {"error": 1}, "error"),
        ("unit", {"warning": 1}, "warning"),
        ("source-skip", {"skip": 1}, "skip"),
        ("installed-skip", {"skip": 4}, "skip"),
    ],
)
def test_green_status_cannot_hide_blocking_counts(
    tmp_path, manifest, kind, mutation, expected_text
):
    documents = all_green_documents(manifest)
    if kind == "package-check":
        target = next(x for x in documents if x["result_kind"] == kind)
        target["check"].update(mutation)
    elif kind == "unit":
        target = next(
            x
            for x in documents
            if x["result_kind"] == "unit" and x["environment_id"] != "nosuggests"
        )
        target["tests"].update(
            {"pass" if key == "pass_count" else key: value for key, value in mutation.items()}
        )
    elif kind == "source-skip":
        target = next(
            x
            for x in documents
            if x["result_kind"] == "unit" and x["environment_id"] != "nosuggests"
        )
        target["tests"].update(mutation)
    else:
        target = next(
            x
            for x in documents
            if x["result_kind"] == "unit" and x["environment_id"] == "nosuggests"
        )
        target["tests"].update(mutation)

    completed, summary, _ = run_aggregator(tmp_path, manifest, documents=documents)
    text = assert_failed_with_complete_table(completed, summary)
    assert expected_text.lower() in text.lower()


def test_current_empty_note_policy_requires_zero_notes(tmp_path, manifest):
    documents = all_green_documents(manifest)
    target = next(x for x in documents if x["result_kind"] == "package-check")
    target["check"]["notes"] = 1

    completed, summary, _ = run_aggregator(tmp_path, manifest, documents=documents)

    text = assert_failed_with_complete_table(completed, summary)
    assert "notes" in text.lower()


def detailed_case(case_class="pass", *, file="test-pass.R", title="pass case", reason=None):
    names = [
        "error",
        "failure",
        "warning",
        "empty",
        "skip_like_success",
        "skip",
        "pass",
    ]
    observation_class = case_class.replace("-", "_")
    observations = {name: 0 for name in names}
    observations[observation_class] = 1
    observations["total"] = 1
    return {
        "class": case_class,
        "message": (
            "success"
            if case_class == "pass"
            else f"Reason: {reason}" if case_class == "skip" else "synthetic failure"
        ),
        "call": None,
        "file": file,
        "line": 1,
        "test_title": title,
        "test_context": None,
        "reason": reason,
        "observations": observations,
    }


def detailed_report(
    context,
    mode,
    cases,
    *,
    status="pass",
    matched=None,
    unused=None,
    violations=None,
    filter_value=None,
):
    classes = [case["class"] for case in cases]
    summary = {
        "total": len(cases),
        "pass": classes.count("pass"),
        "failure": classes.count("failure"),
        "error": classes.count("error"),
        "warning": classes.count("warning"),
        "empty": classes.count("empty"),
        "skip_like_success": classes.count("skip-like-success"),
        "skip": classes.count("skip"),
    }
    observation_names = [
        "total",
        "error",
        "failure",
        "warning",
        "empty",
        "skip_like_success",
        "skip",
        "pass",
    ]
    summary["observations"] = {
        name: sum(case["observations"][name] for case in cases)
        for name in observation_names
    }
    return {
        "schema_version": 1,
        "kind": "unit-tests",
        "context": context,
        "status": status,
        "environment": {
            "id": context,
            "load_package": mode,
            "package": "/tmp/colocboost",
            "r_version": "4.5.1",
            "suite": "full" if filter_value is None else "filtered",
            "filter": filter_value,
        },
        "summary": summary,
        "cases": cases,
        "matched_allowances": [] if matched is None else matched,
        "unused_allowances": [] if unused is None else unused,
        "violations": [] if violations is None else violations,
    }


def public_waiver(waiver):
    return {**waiver, "observed_count": waiver["expected_count"]}


def valid_nosuggests_report(policy):
    waivers = policy["unit_tests"]["allowed_skips"]
    cases = [detailed_case()]
    for waiver in waivers:
        cases.append(
            detailed_case(
                "skip",
                file=waiver["file"],
                title=waiver["test_title"],
                reason=waiver["reason"],
            )
        )
    return detailed_report(
        "r-cmd-check-installed",
        "installed",
        cases,
        matched=[public_waiver(waiver) for waiver in waivers],
    )


_AUTO_SIDECARS = object()


def valid_mkl_sidecars():
    return [
        detailed_report(
            "mkl",
            "source",
            [detailed_case(file=filename, title=f"{filename} passes")],
            filter_value=filter_value,
        )
        for filename, filter_value in (
            ("test_utils.R", "^utils$"),
            ("test_Xref.R", "^Xref$"),
        )
    ]


def run_adapter(tmp_path, document, environment_id, sidecars=_AUTO_SIDECARS):
    diagnostic = tmp_path / "diagnostic.json"
    output = tmp_path / "result.json"
    diagnostic.write_text(json.dumps(document) + "\n", encoding="utf-8")
    arguments = [
            sys.executable,
            "-B",
            os.fspath(UNIT_ADAPTER),
            f"--diagnostic={diagnostic}",
            f"--manifest={MANIFEST_PATH}",
            f"--policy={POLICY_PATH}",
            f"--output={output}",
            f"--environment-id={environment_id}",
            f"--source-sha={SOURCE_SHA}",
            f"--event-sha={EVENT_SHA}",
            f"--tarball-sha256={TARBALL_SHA256}",
        ]
    if sidecars is _AUTO_SIDECARS:
        sidecars = valid_mkl_sidecars() if environment_id == "mkl" else []
    for index, sidecar in enumerate(sidecars):
        sidecar_path = tmp_path / f"sidecar-{index}.json"
        sidecar_path.write_text(json.dumps(sidecar) + "\n", encoding="utf-8")
        arguments.append(f"--sidecar={sidecar_path}")
    completed = subprocess.run(
        arguments,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))
    validate_result(result)
    return completed, result


def test_unit_adapter_accepts_valid_source_report(tmp_path):
    report = detailed_report("mkl", "source", [detailed_case()])
    completed, result = run_adapter(tmp_path, report, "mkl")

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert result["result_kind"] == "unit"
    assert result["environment_id"] == "mkl"
    assert result["status"] == "pass"
    assert result["tests"] == {
        "pass": 1,
        "fail": 0,
        "error": 0,
        "warning": 0,
        "skip": 0,
    }


@pytest.mark.parametrize(
    "mutation",
    [
        "missing",
        "extra",
        "duplicate",
        "wrong-filter",
        "wrong-file",
        "wrong-mode",
        "skip",
        "zero-pass",
        "violation",
    ],
)
def test_mkl_unit_adapter_requires_exact_green_targeted_sidecars(tmp_path, mutation):
    report = detailed_report("mkl", "source", [detailed_case()])
    sidecars = valid_mkl_sidecars()
    if mutation == "missing":
        sidecars.pop()
    elif mutation == "extra":
        sidecars.append(
            detailed_report(
                "mkl",
                "source",
                [detailed_case(file="test_other.R")],
                filter_value="^other$",
            )
        )
    elif mutation == "duplicate":
        sidecars[1] = copy.deepcopy(sidecars[0])
    elif mutation == "wrong-filter":
        sidecars[0]["environment"]["filter"] = "utils"
    elif mutation == "wrong-file":
        sidecars[0]["cases"][0]["file"] = "test_other.R"
    elif mutation == "wrong-mode":
        sidecars[0]["environment"]["load_package"] = "installed"
    elif mutation == "skip":
        sidecars[0] = detailed_report(
            "mkl",
            "source",
            [detailed_case("skip", file="test_utils.R", reason="bad")],
            filter_value="^utils$",
        )
    elif mutation == "zero-pass":
        sidecars[0] = detailed_report(
            "mkl", "source", [], filter_value="^utils$"
        )
    else:
        sidecars[0]["violations"] = [
            {
                "type": "infrastructure-error",
                "message": "sidecar runner failed",
            }
        ]

    completed, result = run_adapter(tmp_path, report, "mkl", sidecars=sidecars)

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"
    assert result["status"] == "fail"


def test_non_mkl_lane_rejects_undeclared_sidecars(tmp_path):
    report = detailed_report("atlas", "source", [detailed_case()])
    sidecar = detailed_report(
        "atlas",
        "source",
        [detailed_case(file="test_utils.R")],
        filter_value="^utils$",
    )

    completed, result = run_adapter(tmp_path, report, "atlas", sidecars=[sidecar])

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"


def test_unit_adapter_accepts_exact_installed_nosuggests_waivers(tmp_path, policy):
    report = valid_nosuggests_report(policy)
    completed, result = run_adapter(tmp_path, report, "nosuggests")

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert result["environment_id"] == "nosuggests"
    assert report["context"] == "r-cmd-check-installed"
    assert result["status"] == "pass"
    assert result["tests"]["pass"] == 1
    assert result["tests"]["skip"] == 5


@pytest.mark.parametrize(
    "mutation",
    [
        "extra-root-field",
        "wrong-context",
        "wrong-environment-context",
        "wrong-mode",
        "summary-count-mismatch",
        "observation-count-mismatch",
        "case-count-mismatch",
        "violations",
        "report-status-fail",
        "source-skip",
        "source-allowance",
        "filtered-suite",
    ],
)
def test_unit_adapter_never_passes_malformed_or_non_green_source(
    tmp_path, mutation
):
    report = detailed_report("mkl", "source", [detailed_case()])
    if mutation == "extra-root-field":
        report["unexpected"] = True
    elif mutation == "wrong-context":
        report["context"] = "other"
        report["environment"]["id"] = "other"
    elif mutation == "wrong-environment-context":
        report["environment"]["id"] = "other"
    elif mutation == "wrong-mode":
        report["environment"]["load_package"] = "installed"
    elif mutation == "summary-count-mismatch":
        report["summary"]["pass"] = 2
    elif mutation == "observation-count-mismatch":
        report["summary"]["observations"]["pass"] = 2
    elif mutation == "case-count-mismatch":
        report["summary"]["total"] = 2
    elif mutation == "violations":
        report["violations"] = [{"type": "failure"}]
    elif mutation == "report-status-fail":
        report = detailed_report(
            "mkl", "source", [detailed_case("failure")], status="fail"
        )
    elif mutation == "source-skip":
        report = detailed_report(
            "mkl", "source", [detailed_case(), detailed_case("skip", reason="x")]
        )
    elif mutation == "source-allowance":
        report["matched_allowances"] = [{"id": "impossible"}]
    elif mutation == "filtered-suite":
        report["environment"]["suite"] = "filtered"
        report["environment"]["filter"] = "test_utils"

    completed, result = run_adapter(tmp_path, report, "mkl")

    assert completed.returncode != 0
    assert result["status"] != "pass"
    if mutation == "report-status-fail":
        assert result["result_kind"] == "unit"
        assert result["tests"]["fail"] == 1
    else:
        assert result["result_kind"] == "infrastructure"


@pytest.mark.parametrize(
    "mutation",
    ["missing-allowance", "wrong-allowance", "unused-allowance", "wrong-skip-case"],
)
def test_unit_adapter_requires_exact_installed_allowances(
    tmp_path, policy, mutation
):
    report = valid_nosuggests_report(policy)
    if mutation == "missing-allowance":
        report["matched_allowances"].pop()
    elif mutation == "wrong-allowance":
        report["matched_allowances"][0]["reason"] = "changed reason"
    elif mutation == "unused-allowance":
        report["unused_allowances"] = [report["matched_allowances"].pop()]
    elif mutation == "wrong-skip-case":
        report["cases"][1]["reason"] = "changed reason"

    completed, result = run_adapter(tmp_path, report, "nosuggests")

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"
    assert result["status"] == "fail"


@pytest.mark.parametrize(
    "mutation",
    [
        "bool-count",
        "float-count",
        "negative-count",
        "extra-summary-field",
        "extra-case-field",
        "case-observation-total",
        "global-observation-total",
        "empty-suite",
        "status-error",
    ],
)
def test_unit_adapter_rejects_invalid_nested_types_and_totals(tmp_path, mutation):
    report = detailed_report("mkl", "source", [detailed_case()])
    if mutation == "bool-count":
        report["summary"]["pass"] = True
    elif mutation == "float-count":
        report["summary"]["pass"] = 1.0
    elif mutation == "negative-count":
        report["summary"]["pass"] = -1
    elif mutation == "extra-summary-field":
        report["summary"]["unexpected"] = 0
    elif mutation == "extra-case-field":
        report["cases"][0]["unexpected"] = 0
    elif mutation == "case-observation-total":
        report["cases"][0]["observations"]["total"] = 2
    elif mutation == "global-observation-total":
        report["summary"]["observations"]["total"] = 2
    elif mutation == "empty-suite":
        report = detailed_report("mkl", "source", [])
    elif mutation == "status-error":
        report["status"] = "error"

    completed, result = run_adapter(tmp_path, report, "mkl")

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"
    assert result["status"] == "fail"


def test_unit_adapter_rejects_malformed_violation_objects_even_when_non_green(
    tmp_path,
):
    report = detailed_report(
        "mkl",
        "source",
        [detailed_case("failure")],
        status="fail",
        violations=[{"type": "failure", "unexpected": "not runner output"}],
    )

    completed, result = run_adapter(tmp_path, report, "mkl")

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"


def test_unit_adapter_rejects_skip_message_that_disagrees_with_reason(
    tmp_path, policy
):
    report = valid_nosuggests_report(policy)
    report["cases"][1]["message"] = "Reason: changed reason"

    completed, result = run_adapter(tmp_path, report, "nosuggests")

    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"


def test_unit_adapter_rejects_duplicate_json_keys(tmp_path):
    diagnostic = tmp_path / "diagnostic.json"
    output = tmp_path / "result.json"
    diagnostic.write_text(
        '{"schema_version":1,"schema_version":1,"kind":"unit-tests"}\n',
        encoding="utf-8",
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(UNIT_ADAPTER),
            f"--diagnostic={diagnostic}",
            f"--manifest={MANIFEST_PATH}",
            f"--policy={POLICY_PATH}",
            f"--output={output}",
            "--environment-id=mkl",
            f"--source-sha={SOURCE_SHA}",
            f"--event-sha={EVENT_SHA}",
            f"--tarball-sha256={TARBALL_SHA256}",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))

    validate_result(result)
    assert completed.returncode != 0
    assert result["result_kind"] == "infrastructure"
