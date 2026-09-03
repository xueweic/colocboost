#!/usr/bin/env python3
"""Finalize one producer result without allowing missing work to look green."""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from result_contract import finalize, validate_result, write_provisional


EXPECTED_KINDS = ("unit", "package-check", "applicability")
JOB_STATUSES = ("success", "failure", "cancelled")
PRODUCER_OUTCOMES = ("success", "failure", "cancelled", "skipped")


def _reject_duplicate_keys(pairs):
    document = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate JSON key {key!r}")
        document[key] = value
    return document


def _expected_identity(arguments, result_kind: str) -> dict[str, Any]:
    return {
        "result_kind": result_kind,
        "environment_id": arguments.environment_id,
        "source_sha": arguments.source_sha,
        "event_sha": arguments.event_sha,
        "tarball_sha256": arguments.tarball_sha256,
    }


def _validate_caller_identity(arguments) -> None:
    validate_result(
        {
            "schema_version": 1,
            "status": "provisional",
            **_expected_identity(arguments, arguments.expected_kind),
        }
    )


def _read_regular_result(path: Path) -> Mapping[str, Any]:
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError("result is missing or inaccessible") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("result must be a regular non-symlink file")
    try:
        with path.open(encoding="utf-8") as stream:
            document = json.load(stream, object_pairs_hook=_reject_duplicate_keys)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"result JSON is malformed: {error}") from error
    validate_result(document)
    return document


def _validate_result_identity(document: Mapping[str, Any], arguments) -> None:
    result_kind = document["result_kind"]
    if result_kind not in {arguments.expected_kind, "infrastructure"}:
        raise ValueError(
            f"result kind {result_kind!r} does not match {arguments.expected_kind!r}"
        )
    for field, expected in (
        ("environment_id", arguments.environment_id),
        ("source_sha", arguments.source_sha),
        ("event_sha", arguments.event_sha),
        ("tarball_sha256", arguments.tarball_sha256),
    ):
        if document[field] != expected:
            raise ValueError(f"result {field} does not match expected identity")
    if document["status"] == "provisional":
        raise ValueError("result remained provisional")
    if result_kind == "infrastructure" and document["status"] != "fail":
        raise ValueError("infrastructure result must have fail status")


def _is_green(document: Mapping[str, Any], expected_kind: str) -> bool:
    if document["result_kind"] != expected_kind:
        return False
    expected_status = (
        "not-applicable" if expected_kind == "applicability" else "pass"
    )
    return document["status"] == expected_status


def _write_infrastructure_failure(path: Path, arguments, message: str) -> None:
    identity = _expected_identity(arguments, "infrastructure")
    write_provisional(path, identity)
    finalize(
        path,
        "fail",
        {"classification": "result-finalization", "message": message},
    )


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Finalize a producer result and fail closed on missing evidence."
    )
    parser.add_argument("--result", required=True)
    parser.add_argument("--expected-kind", required=True, choices=EXPECTED_KINDS)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--job-status", required=True, choices=JOB_STATUSES)
    parser.add_argument(
        "--producer-outcome", required=True, choices=PRODUCER_OUTCOMES
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _argument_parser().parse_args(argv)
    try:
        _validate_caller_identity(arguments)
    except (TypeError, ValueError) as error:
        print(f"result finalizer identity error: {error}", file=sys.stderr)
        return 2

    result_path = Path(arguments.result)
    try:
        document = _read_regular_result(result_path)
        _validate_result_identity(document, arguments)
    except (TypeError, ValueError) as error:
        try:
            _write_infrastructure_failure(
                result_path, arguments, f"Producer result was unusable: {error}"
            )
        except (OSError, TypeError, ValueError) as write_error:
            print(
                f"result finalizer could not write infrastructure failure: {write_error}",
                file=sys.stderr,
            )
            return 2
        print(f"result finalizer failed closed: {error}", file=sys.stderr)
        return 1

    if not _is_green(document, arguments.expected_kind):
        print(
            f"result finalizer preserved non-green {document['result_kind']}/"
            f"{document['status']}",
            file=sys.stderr,
        )
        return 1

    if arguments.job_status != "success" or arguments.producer_outcome != "success":
        message = (
            "Producer wrote a green result but execution was not successful: "
            f"job-status={arguments.job_status}, "
            f"producer-outcome={arguments.producer_outcome}"
        )
        try:
            _write_infrastructure_failure(result_path, arguments, message)
        except (OSError, TypeError, ValueError) as error:
            print(
                f"result finalizer could not write infrastructure failure: {error}",
                file=sys.stderr,
            )
            return 2
        print(message, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
