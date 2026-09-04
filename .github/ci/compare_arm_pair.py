#!/usr/bin/env python3
"""Make the sole linux-arm64 result from two raw native check proofs."""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from result_contract import finalize, write_provisional

ARM_INDEX = "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:74f3781764e90b6c4bf8e300c788b0f853e867289585f020dd3b842180c08bf6"
CHILD_IMAGES = {
    "amd64": "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:ef57b4fbbcc76a173e53622143d59f8d4a97978481237e81153af23b03c66e9c",
    "arm64": "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:bcfe37df999b716bd41ffc139d242c251c06b79d82f3632b4e8c65b6eb13f050",
}
CHECK_ARGS = ["--as-cran", "--no-manual", "--no-build-vignettes"]
CONFIGURATION = {"args": CHECK_ARGS, "r": "/usr/bin/R"}
CHILD_FIELDS = {
    "environment_id", "architecture", "image", "index", "system_r",
    "source_sha", "event_sha", "tarball_sha256", "status", "check_status",
    "check_exit_code", "check_args", "configuration", "tarball", "native",
}
NATIVE_FIELDS = {"clean", "uname_m", "r_platform", "r_home", "image", "index"}


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    document: dict[str, Any] = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate JSON key: {key}")
        document[key] = value
    return document


def _load(path: str | os.PathLike[str]) -> Mapping[str, Any]:
    value = Path(path)
    status = value.lstat()
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"ARM evidence is not a regular non-symlink file: {value}")
    with value.open(encoding="utf-8") as stream:
        document = json.load(stream, object_pairs_hook=_reject_duplicate_keys)
    if not isinstance(document, Mapping):
        raise ValueError(f"ARM evidence is not an object: {value}")
    return document


def compare_pair(
    amd64_path: str | os.PathLike[str],
    arm64_path: str | os.PathLike[str],
    *,
    result_path: str | os.PathLike[str],
    evidence_path: str | os.PathLike[str],
    source_sha: str,
    event_sha: str,
    tarball_sha256: str,
    expected_index: str = ARM_INDEX,
) -> Mapping[str, Any]:
    amd64 = _load(amd64_path)
    arm64 = _load(arm64_path)
    errors: list[str] = []
    documents = {"amd64": amd64, "arm64": arm64}
    for label, document in documents.items():
        missing = sorted(CHILD_FIELDS - set(document))
        extra = sorted(set(document) - CHILD_FIELDS)
        if missing or extra:
            errors.append(f"{label} evidence fields differ: missing={missing}, extra={extra}")
        if document.get("environment_id") != "linux-arm64":
            errors.append(f"{label} environment identity is not linux-arm64")
        for key, expected_value in (("source_sha", source_sha), ("event_sha", event_sha), ("tarball_sha256", tarball_sha256)):
            if document.get(key) != expected_value:
                errors.append(f"{label} {key} does not match shared identity")
        if document.get("architecture") != label:
            errors.append(f"{label} native architecture identity is missing or wrong")
        if document.get("index") != expected_index:
            errors.append(f"{label} image index does not match pinned manifest")
        if document.get("system_r") != "/usr/bin/R":
            errors.append(f"{label} system R is not /usr/bin/R")
        if document.get("check_args") != CHECK_ARGS:
            errors.append(f"{label} check arguments are not exact")
        if document.get("configuration") != CONFIGURATION:
            errors.append(f"{label} check configuration is not exact")
        if document.get("tarball") != "shared-source":
            errors.append(f"{label} did not consume the shared tarball")
        if document.get("status") != "pass" or document.get("check_status") != "pass":
            errors.append(f"{label} package check is not clean")
        if type(document.get("check_exit_code")) is not int or document.get("check_exit_code") != 0:
            errors.append(f"{label} package check exit code is not zero")
        if document.get("image") != CHILD_IMAGES[label]:
            errors.append(f"{label} image is not a pinned rcheckserver child")
        if not isinstance(document.get("native"), Mapping):
            errors.append(f"{label} native runtime proof is missing")
        else:
            native = document["native"]
            native_missing = sorted(NATIVE_FIELDS - set(native))
            native_extra = sorted(set(native) - NATIVE_FIELDS)
            if native_missing or native_extra:
                errors.append(f"{label} native fields differ: missing={native_missing}, extra={native_extra}")
            if native.get("clean") is not True:
                errors.append(f"{label} native runtime proof is not clean")
            if native.get("image") != CHILD_IMAGES[label] or native.get("index") != expected_index:
                errors.append(f"{label} native image identity does not match its child")
            expected_machine = "x86_64" if label == "amd64" else "aarch64"
            if native.get("uname_m") not in {expected_machine, "amd64" if label == "amd64" else "arm64"}:
                errors.append(f"{label} uname architecture does not match child")
            if not native.get("r_platform") or expected_machine not in str(native.get("r_platform")):
                errors.append(f"{label} R platform evidence does not match child")
            if not str(native.get("r_home", "")).startswith("/"):
                errors.append(f"{label} R home is not an absolute native path")
    comparison = {
        "schema_version": 1,
        "environment_id": "linux-arm64",
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "amd64": amd64,
        "arm64": arm64,
        "errors": errors,
        "architecture_only": not errors,
    }
    output = Path(evidence_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, output)
    write_provisional(result_path, {"result_kind": "package-check", "environment_id": "linux-arm64", "source_sha": source_sha, "event_sha": event_sha, "tarball_sha256": tarball_sha256})
    finalize(result_path, "pass" if not errors else "fail", {"check": {"errors": len(errors), "warnings": 0, "notes": 0, "log_path": os.fspath(evidence_path)}})
    return comparison


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amd64", required=True)
    parser.add_argument("--arm64", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--expected-index", default=ARM_INDEX)
    args = parser.parse_args(argv)
    try:
        comparison = compare_pair(**vars(args))
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
        print(f"ARM pair comparison failed closed: {error}", file=sys.stderr)
        return 2
    return 0 if not comparison["errors"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
