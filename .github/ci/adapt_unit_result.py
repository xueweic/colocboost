#!/usr/bin/env python3
"""Strictly adapt a detailed unit-test diagnostic to the aggregate contract."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from result_contract import finalize, validate_result, write_provisional
from validate_manifest import load_manifest, validate_manifest, validate_policy


ROOT_FIELDS = {
    "schema_version",
    "kind",
    "context",
    "status",
    "environment",
    "summary",
    "cases",
    "matched_allowances",
    "unused_allowances",
    "violations",
}
ENVIRONMENT_FIELDS = {
    "id",
    "load_package",
    "package",
    "r_version",
    "suite",
    "filter",
}
SUMMARY_COUNT_FIELDS = {
    "total",
    "pass",
    "failure",
    "error",
    "warning",
    "empty",
    "skip_like_success",
    "skip",
}
OBSERVATION_FIELDS = {
    "total",
    "error",
    "failure",
    "warning",
    "empty",
    "skip_like_success",
    "skip",
    "pass",
}
CASE_FIELDS = {
    "class",
    "message",
    "call",
    "file",
    "line",
    "test_title",
    "test_context",
    "reason",
    "observations",
}
CASE_CLASSES = {
    "error",
    "failure",
    "warning",
    "empty",
    "skip-like-success",
    "skip",
    "pass",
}
WAIVER_FIELDS = {
    "id",
    "context",
    "file",
    "test_title",
    "reason",
    "expected_count",
    "rationale",
    "expires",
    "observed_count",
}
STRICT_VIOLATION_FIELDS = {"type", "file", "line", "test_title", "message"}
VIOLATION_FIELDS = {
    "unexpected-skip": {
        "type",
        "file",
        "line",
        "test_title",
        "reason",
        "message",
    },
    "skip-count-mismatch": {
        "type",
        "id",
        "expected_count",
        "observed_count",
        "message",
    },
    "expired-allowance": {
        "type",
        "id",
        "expires",
        "observed_count",
        "message",
    },
    "unused-allowance": {"type", "id", "message"},
    "cli-error": {"type", "message"},
    "filesystem-error": {"type", "message"},
    "policy-error": {"type", "message"},
    "infrastructure-error": {"type", "message"},
    "execution-error": {"type", "message"},
}
STRICT_SUMMARY_FIELDS = {
    "failure",
    "error",
    "warning",
    "empty",
    "skip_like_success",
}


def _reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _load_json(path: Path) -> Any:
    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream, object_pairs_hook=_reject_duplicate_keys)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"could not read diagnostic JSON {path}: {error}") from error


def _require_exact_fields(value: Any, expected: set[str], label: str) -> Mapping:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"{label} fields differ: missing={missing}, extra={extra}")
    return value


def _require_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty string")
    return value


def _require_nullable_string(value: Any, label: str) -> None:
    if value is not None and not isinstance(value, str):
        raise ValueError(f"{label} must be a string or null")


def _require_count(value: Any, label: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _validate_count_object(
    value: Any, expected: set[str], label: str
) -> dict[str, int]:
    count_object = _require_exact_fields(value, expected, label)
    return {
        field: _require_count(count_object[field], f"{label}.{field}")
        for field in expected
    }


def _class_from_observations(observations: Mapping[str, int]) -> str:
    for field, case_class in (
        ("error", "error"),
        ("failure", "failure"),
        ("warning", "warning"),
        ("empty", "empty"),
        ("skip_like_success", "skip-like-success"),
        ("skip", "skip"),
    ):
        if observations[field] > 0:
            return case_class
    return "pass"


def _validate_public_allowance(value: Any, label: str) -> None:
    allowance = _require_exact_fields(value, WAIVER_FIELDS, label)
    for field in (
        "id",
        "context",
        "file",
        "test_title",
        "reason",
        "rationale",
        "expires",
    ):
        _require_string(allowance[field], f"{label}.{field}")
    _require_count(allowance["expected_count"], f"{label}.expected_count")
    _require_count(allowance["observed_count"], f"{label}.observed_count")
    if allowance["expected_count"] < 1:
        raise ValueError(f"{label}.expected_count must be positive")


def _validate_violation(value: Any, label: str) -> None:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    violation_type = value.get("type")
    if violation_type in {"error", "failure", "warning", "empty", "skip-like-success"}:
        expected = STRICT_VIOLATION_FIELDS
    else:
        expected = VIOLATION_FIELDS.get(violation_type)
    if expected is None:
        raise ValueError(f"{label}.type is invalid")
    violation = _require_exact_fields(value, expected, label)
    _require_string(violation["type"], f"{label}.type")
    _require_string(violation["message"], f"{label}.message")
    if "file" in violation:
        _require_nullable_string(violation["file"], f"{label}.file")
    if "line" in violation and violation["line"] is not None and (
        type(violation["line"]) is not int or violation["line"] < 1
    ):
        raise ValueError(f"{label}.line is invalid")
    for field in ("test_title", "reason", "id", "expires"):
        if field in violation:
            _require_string(violation[field], f"{label}.{field}")
    for field in ("expected_count", "observed_count"):
        if field in violation:
            _require_count(violation[field], f"{label}.{field}")


def validate_diagnostic(document: Any) -> Mapping[str, Any]:
    """Validate the closed Task 4 diagnostic and its redundant counts."""

    report = _require_exact_fields(document, ROOT_FIELDS, "diagnostic")
    if type(report["schema_version"]) is not int or report["schema_version"] != 1:
        raise ValueError("diagnostic.schema_version must be integer 1")
    if report["kind"] != "unit-tests":
        raise ValueError("diagnostic.kind must be 'unit-tests'")
    _require_string(report["context"], "diagnostic.context")
    if report["status"] not in {"pass", "fail", "error"}:
        raise ValueError("diagnostic.status is invalid")

    environment = _require_exact_fields(
        report["environment"], ENVIRONMENT_FIELDS, "diagnostic.environment"
    )
    for field in ("id", "load_package", "package", "r_version", "suite"):
        _require_string(environment[field], f"diagnostic.environment.{field}")
    if environment["load_package"] not in {"source", "installed"}:
        raise ValueError("diagnostic.environment.load_package is invalid")
    if environment["suite"] not in {"full", "filtered"}:
        raise ValueError("diagnostic.environment.suite is invalid")
    _require_nullable_string(environment["filter"], "diagnostic.environment.filter")
    if environment["suite"] == "full" and environment["filter"] is not None:
        raise ValueError("a full-suite diagnostic must have a null filter")
    if environment["suite"] == "filtered" and not environment["filter"]:
        raise ValueError("a filtered diagnostic must name its filter")
    if environment["id"] != report["context"]:
        raise ValueError("diagnostic environment.id does not match context")

    summary = _require_exact_fields(
        report["summary"], SUMMARY_COUNT_FIELDS | {"observations"}, "diagnostic.summary"
    )
    summary_counts = {
        field: _require_count(summary[field], f"diagnostic.summary.{field}")
        for field in SUMMARY_COUNT_FIELDS
    }
    if summary_counts["total"] != sum(
        summary_counts[field] for field in SUMMARY_COUNT_FIELDS - {"total"}
    ):
        raise ValueError("diagnostic summary class counts do not sum to total")
    summary_observations = _validate_count_object(
        summary["observations"], OBSERVATION_FIELDS, "diagnostic.summary.observations"
    )
    if summary_observations["total"] != sum(
        summary_observations[field] for field in OBSERVATION_FIELDS - {"total"}
    ):
        raise ValueError("diagnostic summary observations do not sum to total")

    cases = report["cases"]
    if not isinstance(cases, list):
        raise ValueError("diagnostic.cases must be a list")
    if len(cases) != summary_counts["total"]:
        raise ValueError("diagnostic case count does not match summary.total")
    class_counts = Counter()
    accumulated_observations = Counter({field: 0 for field in OBSERVATION_FIELDS})
    for index, raw_case in enumerate(cases):
        case = _require_exact_fields(raw_case, CASE_FIELDS, f"diagnostic.cases[{index}]")
        case_class = case["class"]
        if case_class not in CASE_CLASSES:
            raise ValueError(f"diagnostic.cases[{index}].class is invalid")
        _require_string(case["message"], f"diagnostic.cases[{index}].message")
        _require_nullable_string(case["call"], f"diagnostic.cases[{index}].call")
        _require_nullable_string(case["file"], f"diagnostic.cases[{index}].file")
        _require_nullable_string(
            case["test_context"], f"diagnostic.cases[{index}].test_context"
        )
        _require_nullable_string(case["reason"], f"diagnostic.cases[{index}].reason")
        _require_string(case["test_title"], f"diagnostic.cases[{index}].test_title")
        if case["line"] is not None and (
            type(case["line"]) is not int or case["line"] < 1
        ):
            raise ValueError(f"diagnostic.cases[{index}].line is invalid")
        observations = _validate_count_object(
            case["observations"],
            OBSERVATION_FIELDS,
            f"diagnostic.cases[{index}].observations",
        )
        if observations["total"] < 1 or observations["total"] != sum(
            observations[field] for field in OBSERVATION_FIELDS - {"total"}
        ):
            raise ValueError(
                f"diagnostic.cases[{index}] observations do not sum to a positive total"
            )
        if _class_from_observations(observations) != case_class:
            raise ValueError(
                f"diagnostic.cases[{index}] class contradicts its observations"
            )
        if case_class == "skip":
            reason = _require_string(case["reason"], f"diagnostic.cases[{index}].reason")
            _require_string(case["file"], f"diagnostic.cases[{index}].file")
            if case["message"] != f"Reason: {reason}":
                raise ValueError(
                    f"diagnostic.cases[{index}] skip message contradicts reason"
                )
        class_counts[case_class] += 1
        accumulated_observations.update(observations)

    class_to_summary = {
        "pass": "pass",
        "failure": "failure",
        "error": "error",
        "warning": "warning",
        "empty": "empty",
        "skip-like-success": "skip_like_success",
        "skip": "skip",
    }
    for case_class, summary_field in class_to_summary.items():
        if class_counts[case_class] != summary_counts[summary_field]:
            raise ValueError(
                f"diagnostic summary.{summary_field} contradicts case classes"
            )
    if dict(accumulated_observations) != summary_observations:
        raise ValueError("diagnostic summary observations contradict case observations")

    for field in ("matched_allowances", "unused_allowances"):
        value = report[field]
        if not isinstance(value, list):
            raise ValueError(f"diagnostic.{field} must be a list")
        for index, allowance in enumerate(value):
            _validate_public_allowance(
                allowance, f"diagnostic.{field}[{index}]"
            )
    if not isinstance(report["violations"], list):
        raise ValueError("diagnostic.violations must be a list")
    for index, violation in enumerate(report["violations"]):
        _validate_violation(violation, f"diagnostic.violations[{index}]")

    return report


def _expected_public_allowances(policy: Mapping, context: str) -> list[dict]:
    return [
        {**waiver, "observed_count": waiver["expected_count"]}
        for waiver in policy["unit_tests"]["allowed_skips"]
        if waiver["context"] == context
    ]


def validate_lane_report(
    report: Mapping[str, Any], lane: Mapping[str, Any], policy: Mapping[str, Any]
) -> None:
    """Bind a valid detailed report to one exact manifest unit lane."""

    environment = report["environment"]
    summary = report["summary"]
    if report["context"] != lane["runner_context"]:
        raise ValueError("diagnostic context does not match unit-lane runner_context")
    if environment["load_package"] != lane["mode"]:
        raise ValueError("diagnostic load mode does not match unit lane")
    if lane["suite"] != "full" or environment["suite"] != "full":
        raise ValueError("aggregate unit lanes require an unfiltered full-suite diagnostic")
    if environment["filter"] is not None:
        raise ValueError("aggregate unit lanes cannot adapt a filtered diagnostic")

    if report["status"] == "error":
        raise ValueError("unit runner reported an infrastructure error")
    if report["status"] == "pass":
        if report["violations"] or report["unused_allowances"]:
            raise ValueError("passing diagnostic contains violations or unused allowances")
        if any(summary[field] != 0 for field in STRICT_SUMMARY_FIELDS):
            raise ValueError("passing diagnostic contains a strict non-green test class")
        if summary["pass"] < 1:
            raise ValueError("passing diagnostic contains no passing test cases")

    if lane["mode"] == "source":
        if summary["skip"] != 0:
            raise ValueError("source-mode unit lane contains skips")
        if report["matched_allowances"] or report["unused_allowances"]:
            raise ValueError("source-mode unit lane contains skip allowances")
        return

    expected_allowances = _expected_public_allowances(policy, lane["runner_context"])
    if report["unused_allowances"]:
        raise ValueError("installed unit lane contains unused allowances")
    if report["matched_allowances"] != expected_allowances:
        raise ValueError("installed unit lane allowances do not match current policy")
    expected_skips = sum(item["observed_count"] for item in expected_allowances)
    if summary["skip"] != expected_skips:
        raise ValueError("installed unit lane skip count does not match current policy")

    observed_skip_keys = Counter(
        (case["file"], case["test_title"], case["reason"])
        for case in report["cases"]
        if case["class"] == "skip"
    )
    expected_skip_keys = Counter()
    for waiver in policy["unit_tests"]["allowed_skips"]:
        if waiver["context"] == lane["runner_context"]:
            expected_skip_keys[
                (waiver["file"], waiver["test_title"], waiver["reason"])
            ] += waiver["expected_count"]
    if observed_skip_keys != expected_skip_keys:
        raise ValueError("installed unit lane skip cases do not match current policy")


def _identity(arguments, result_kind: str) -> dict[str, Any]:
    return {
        "result_kind": result_kind,
        "environment_id": arguments.environment_id,
        "source_sha": arguments.source_sha,
        "event_sha": arguments.event_sha,
        "tarball_sha256": arguments.tarball_sha256,
    }


def _emit(arguments, result_kind: str, status: str, details: Mapping[str, Any]):
    output = Path(arguments.output)
    identity = _identity(arguments, result_kind)
    provisional = {"schema_version": 1, "status": "provisional", **identity}
    validate_result(provisional)
    write_provisional(output, identity)
    return finalize(output, status, details)


def _emit_infrastructure(arguments, message: str) -> None:
    _emit(
        arguments,
        "infrastructure",
        "fail",
        {"classification": "unit-adapter", "message": message},
    )


def _parse_args(argv: list[str]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--diagnostic", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    arguments = _parse_args(argv)
    try:
        # Validate caller-owned identity before it is used for any result kind.
        validate_result(
            {
                "schema_version": 1,
                "status": "provisional",
                **_identity(arguments, "unit"),
            }
        )
    except (TypeError, ValueError) as error:
        print(f"unit adapter identity invalid: {error}", file=sys.stderr)
        return 2

    try:
        manifest = validate_manifest(load_manifest(arguments.manifest))
        policy = validate_policy(load_manifest(arguments.policy))
        lanes = [
            lane
            for lane in manifest["unit_lanes"]
            if lane["environment_id"] == arguments.environment_id
        ]
        if len(lanes) != 1:
            raise ValueError("environment-id does not identify one declared unit lane")
        report = validate_diagnostic(_load_json(Path(arguments.diagnostic)))
        validate_lane_report(report, lanes[0], policy)
    except (OSError, TypeError, ValueError) as error:
        try:
            _emit_infrastructure(arguments, str(error))
        except (OSError, TypeError, ValueError) as write_error:
            print(f"unit adapter could not write failure result: {write_error}", file=sys.stderr)
            return 2
        print(f"unit adapter rejected diagnostic: {error}", file=sys.stderr)
        return 2

    summary = report["summary"]
    counts = {
        "pass": summary["pass"],
        "fail": summary["failure"]
        + summary["empty"]
        + summary["skip_like_success"],
        "error": summary["error"],
        "warning": summary["warning"],
        "skip": summary["skip"],
    }
    status = "pass" if report["status"] == "pass" else "fail"
    try:
        _emit(arguments, "unit", status, {"tests": counts})
    except (OSError, TypeError, ValueError) as error:
        print(f"unit adapter could not write result: {error}", file=sys.stderr)
        return 2
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
