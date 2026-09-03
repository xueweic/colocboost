#!/usr/bin/env python3
"""Fail-closed aggregation for the complete manifest-declared result inventory."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from result_contract import validate_result
from validate_manifest import (
    expected_result_keys,
    load_manifest,
    validate_manifest,
    validate_policy,
)


SHA = re.compile(r"^[0-9a-f]{40}$")
DIGEST = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class ExpectedResult:
    result_kind: str
    environment_id: str
    declared: str
    coverage_state: str | None
    unit_mode: str | None
    runner_context: str | None

    @property
    def key(self) -> tuple[str, str]:
        return (self.result_kind, self.environment_id)


@dataclass
class Artifact:
    path: Path
    document: Mapping[str, Any] | None
    error: str | None


def _reject_duplicate_keys(pairs):
    document = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate JSON key {key!r}")
        document[key] = value
    return document


def _read_result(path: Path) -> Any:
    if path.is_symlink() or not path.is_file():
        raise ValueError("result artifact must be a regular non-symlink file")
    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream, object_pairs_hook=_reject_duplicate_keys)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"invalid JSON: {error}") from error


def _build_expected(manifest: Mapping[str, Any]) -> list[ExpectedResult]:
    coverage_by_id = {row["id"]: row for row in manifest["coverage"]}
    lanes_by_id = {lane["environment_id"]: lane for lane in manifest["unit_lanes"]}
    expected = []
    for result_kind, environment_id in expected_result_keys(manifest):
        if result_kind == "unit":
            lane = lanes_by_id[environment_id]
            expected.append(
                ExpectedResult(
                    result_kind=result_kind,
                    environment_id=environment_id,
                    declared=f"unit/{lane['mode']}/full",
                    coverage_state=None,
                    unit_mode=lane["mode"],
                    runner_context=lane["runner_context"],
                )
            )
        else:
            row = coverage_by_id[environment_id]
            expected.append(
                ExpectedResult(
                    result_kind=result_kind,
                    environment_id=environment_id,
                    declared=f"{row['group']}/{row['state']}",
                    coverage_state=row["state"],
                    unit_mode=None,
                    runner_context=None,
                )
            )
    if len(expected) != 60:
        raise ValueError(f"manifest yielded {len(expected)} results rather than 60")
    return expected


def _recover_key(document: Any) -> tuple[str, str] | None:
    if not isinstance(document, Mapping):
        return None
    result_kind = document.get("result_kind")
    environment_id = document.get("environment_id")
    if not isinstance(result_kind, str) or not isinstance(environment_id, str):
        return None
    return result_kind, environment_id


def _scan_artifacts(
    results_root: Path, expected_keys: set[tuple[str, str]]
) -> tuple[dict[tuple[str, str], list[Artifact]], list[str]]:
    by_key: dict[tuple[str, str], list[Artifact]] = defaultdict(list)
    issues = []
    if not results_root.exists():
        return by_key, [f"results root does not exist: {results_root}"]
    if results_root.is_symlink() or not results_root.is_dir():
        return by_key, [f"results root is not a regular directory: {results_root}"]

    try:
        candidates = sorted(results_root.rglob("result.json"))
    except OSError as error:
        return by_key, [f"could not scan results root {results_root}: {error}"]

    for path in candidates:
        try:
            document = _read_result(path)
        except ValueError as error:
            issues.append(f"{path}: {error}")
            continue
        key = _recover_key(document)
        if key is None:
            issues.append(f"{path}: schema invalid: result key is not recoverable")
            continue
        if key not in expected_keys:
            label = "infrastructure result" if key[0] == "infrastructure" else "unexpected result"
            issues.append(f"{path}: {label} {key[0]}/{key[1]}")
            continue
        try:
            validate_result(document)
        except (TypeError, ValueError) as error:
            by_key[key].append(Artifact(path, document, f"schema invalid: {error}"))
        else:
            by_key[key].append(Artifact(path, document, None))
    return by_key, issues


def _validate_expected_identity(arguments) -> list[str]:
    issues = []
    for field in ("expected_source_sha", "expected_event_sha"):
        value = getattr(arguments, field)
        if SHA.fullmatch(value) is None:
            issues.append(f"{field} must be a 40-character lowercase Git SHA")
    if DIGEST.fullmatch(arguments.expected_tarball_sha256) is None:
        issues.append(
            "expected_tarball_sha256 must be exactly 64 lowercase hexadecimal characters"
        )
    return issues


def _installed_skip_count(
    policy: Mapping[str, Any] | None, runner_context: str | None
) -> int | None:
    if policy is None or runner_context is None:
        return None
    return sum(
        waiver["expected_count"]
        for waiver in policy["unit_tests"]["allowed_skips"]
        if waiver["context"] == runner_context
    )


def _semantic_errors(
    expected: ExpectedResult,
    document: Mapping[str, Any],
    arguments,
    policy: Mapping[str, Any] | None,
) -> list[str]:
    errors = []
    for document_field, argument_field in (
        ("source_sha", "expected_source_sha"),
        ("event_sha", "expected_event_sha"),
        ("tarball_sha256", "expected_tarball_sha256"),
    ):
        if document[document_field] != getattr(arguments, argument_field):
            errors.append(f"{document_field} does not match expected identity")

    if expected.result_kind == "package-check":
        if document["status"] != "pass":
            errors.append(f"status is {document['status']}, expected pass")
        check = document.get("check")
        if isinstance(check, Mapping):
            if check["errors"] != 0:
                errors.append(f"errors={check['errors']}, expected 0")
            if check["warnings"] != 0:
                errors.append(f"warnings={check['warnings']}, expected 0")
            # The committed policy is deliberately empty. A future reviewed
            # allowance must update this aggregation contract as well.
            if check["notes"] != 0:
                errors.append(f"notes={check['notes']}, committed policy expects 0")
    elif expected.result_kind == "applicability":
        if document["status"] != "not-applicable":
            errors.append(
                f"status is {document['status']}, expected not-applicable"
            )
    else:
        if document["status"] != "pass":
            errors.append(f"status is {document['status']}, expected pass")
        tests = document.get("tests")
        if isinstance(tests, Mapping):
            if tests["pass"] < 1:
                errors.append("unit pass count must be positive")
            for field in ("fail", "error", "warning"):
                if tests[field] != 0:
                    errors.append(f"unit {field}={tests[field]}, expected 0")
            if expected.unit_mode == "source" and tests["skip"] != 0:
                errors.append(f"source unit skip={tests['skip']}, expected 0")
            elif expected.unit_mode == "installed":
                allowed = _installed_skip_count(policy, expected.runner_context)
                if allowed is None:
                    errors.append("installed unit skip policy is unavailable")
                elif tests["skip"] != allowed:
                    errors.append(
                        f"installed unit skip={tests['skip']}, expected {allowed}"
                    )
    return errors


def _details(document: Mapping[str, Any] | None) -> str:
    if document is None:
        return "missing"
    kind = document.get("result_kind")
    if kind == "package-check" and isinstance(document.get("check"), Mapping):
        check = document["check"]
        return (
            f"errors={check.get('errors')}; warnings={check.get('warnings')}; "
            f"notes={check.get('notes')}; log={check.get('log_path')}"
        )
    if kind == "unit" and isinstance(document.get("tests"), Mapping):
        tests = document["tests"]
        return "; ".join(f"{field}={tests.get(field)}" for field in (
            "pass", "fail", "error", "warning", "skip"
        ))
    if kind == "applicability":
        return f"predicate={document.get('predicate')}; evidence={document.get('evidence')}"
    return f"status={document.get('status')}"


def _escape_markdown(value: Any) -> str:
    return str(value).replace("\r\n", "\n").replace("\r", "\n").replace("|", "\\|").replace("\n", "<br>")


def _write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _render_summary(
    expected_rows: list[ExpectedResult],
    row_results: list[tuple[ExpectedResult, Mapping[str, Any] | None, list[str]]],
    global_issues: list[str],
) -> tuple[str, bool]:
    passed = sum(not errors for _, _, errors in row_results)
    overall_pass = passed == len(expected_rows) and not global_issues
    lines = [
        "# CRAN preflight summary",
        "",
        f"**Overall: {'PASS' if overall_pass else 'FAIL'}**",
        "",
        f"Coverage: **{passed}/{len(expected_rows)}** required results passed.",
        "",
        "| Result kind | Environment | Declared | Observed | Details |",
        "| --- | --- | --- | --- | --- |",
    ]
    for expected, document, errors in row_results:
        observed = "missing" if document is None else str(document.get("status", "invalid"))
        details = _details(document)
        if errors:
            details = f"FAIL: {'; '.join(errors)}; {details}"
        lines.append(
            "| "
            + " | ".join(
                _escape_markdown(value)
                for value in (
                    expected.result_kind,
                    expected.environment_id,
                    expected.declared,
                    observed,
                    details,
                )
            )
            + " |"
        )
    if global_issues:
        lines.extend(["", "## Artifact and configuration issues", ""])
        lines.extend(f"- {_escape_markdown(issue)}" for issue in global_issues)
    lines.append("")
    return "\n".join(lines), overall_pass


def _parse_args(argv: list[str]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--expected-source-sha", required=True)
    parser.add_argument("--expected-event-sha", required=True)
    parser.add_argument("--expected-tarball-sha256", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    arguments = _parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(arguments.manifest))
        expected_rows = _build_expected(manifest)
    except (KeyError, OSError, TypeError, ValueError) as error:
        print(f"aggregate manifest invalid: {error}", file=sys.stderr)
        return 2

    global_issues = _validate_expected_identity(arguments)
    policy = None
    policy_path = Path(arguments.manifest).with_name("check-policy.yml")
    try:
        policy = validate_policy(load_manifest(policy_path))
    except (KeyError, OSError, TypeError, ValueError) as error:
        global_issues.append(f"committed check policy invalid: {error}")

    expected_keys = {row.key for row in expected_rows}
    artifacts, scan_issues = _scan_artifacts(
        Path(arguments.results_root), expected_keys
    )
    global_issues.extend(scan_issues)

    identity_invalid = bool(_validate_expected_identity(arguments))
    row_results = []
    for expected in expected_rows:
        candidates = artifacts.get(expected.key, [])
        document = candidates[0].document if candidates else None
        errors = []
        if identity_invalid:
            errors.append("expected identity is invalid")
        if not candidates:
            errors.append("missing required result")
        elif len(candidates) > 1:
            errors.append(f"duplicate result artifacts ({len(candidates)})")
        if len(candidates) == 1:
            artifact = candidates[0]
            if artifact.error:
                errors.append(artifact.error)
            elif artifact.document is not None:
                errors.extend(
                    _semantic_errors(expected, artifact.document, arguments, policy)
                )
        row_results.append((expected, document, errors))

    summary, overall_pass = _render_summary(
        expected_rows, row_results, global_issues
    )
    try:
        _write_atomic(Path(arguments.summary), summary)
    except OSError as error:
        print(f"could not write aggregate summary: {error}", file=sys.stderr)
        return 2
    print(
        f"aggregate {'passed' if overall_pass else 'failed'}: "
        f"{sum(not errors for _, _, errors in row_results)}/{len(expected_rows)}"
    )
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
