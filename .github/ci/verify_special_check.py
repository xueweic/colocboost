#!/usr/bin/env python3
"""Fail-closed semantic validation for CRAN check trees and special lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

from artifact_contract import append_github_outputs
from validate_manifest import load_manifest, validate_manifest


FULL_DOCUMENTATION = "full-documentation"
NO_DOCUMENTATION = "standard-no-documentation"
PROFILES = {FULL_DOCUMENTATION, NO_DOCUMENTATION}
DOC_HEADINGS = {
    "manual": re.compile(r"^\* checking PDF version of manual \.\.\.(?: \[[^ ]+\])? OK$"),
    "package": re.compile(r"^\* checking package vignettes \.\.\.(?: \[[^ ]+\])? OK$"),
    "rebuilding": re.compile(r"^\* checking re-building of vignette outputs \.\.\.(?: \[[^ ]+\])? OK$"),
}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _read_json(path_value, label):
    path = Path(path_value)
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is missing") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file")
    try:
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is malformed") from error
    if not isinstance(document, Mapping):
        raise ValueError(f"{label} must be an object")
    return document


def _scan_tree(root_value):
    root = Path(root_value)
    if not root.is_absolute():
        raise ValueError("check root must be absolute")
    try:
        root_status = root.lstat()
    except OSError as error:
        raise ValueError("check root is missing") from error
    if stat.S_ISLNK(root_status.st_mode) or not stat.S_ISDIR(root_status.st_mode):
        raise ValueError("check root must be a regular non-symlink directory")
    for current, directories, files in os.walk(root, followlinks=False):
        for name in [*directories, *files]:
            path = Path(current) / name
            status = path.lstat()
            if stat.S_ISLNK(status.st_mode):
                raise ValueError("check tree must not contain symlinks")
            if name in directories and not stat.S_ISDIR(status.st_mode):
                raise ValueError("check tree contains a non-directory entry")
            if name in files and not stat.S_ISREG(status.st_mode):
                raise ValueError("check tree contains a non-regular file")
    candidates = [entry for entry in root.iterdir() if entry.name.endswith(".Rcheck")]
    if len(candidates) != 1:
        raise ValueError("check root must contain exactly one *.Rcheck directory")
    check = candidates[0]
    if not check.is_dir():
        raise ValueError("*.Rcheck must be a directory")
    log = check / "00check.log"
    try:
        status = log.lstat()
    except OSError as error:
        raise ValueError("check tree lacks 00check.log") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("00check.log must be a regular non-symlink file")
    return check, log


def _atomic_write(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(document, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def verify_special_check(
    manifest_row,
    *,
    environment_id,
    profile,
    check_root,
    native_evidence,
    runtime_evidence,
    source_sha,
    event_sha,
    tarball_sha256,
    output,
    github_output=None,
):
    if profile not in PROFILES:
        raise ValueError("unsupported special check profile")
    if manifest_row.get("id") != environment_id:
        raise ValueError("environment ID does not match manifest")
    expected_profile = (
        FULL_DOCUMENTATION if manifest_row.get("check_args") == [] else NO_DOCUMENTATION
    )
    if profile != expected_profile:
        raise ValueError("check profile does not match manifest check arguments")
    native = _read_json(native_evidence, "native evidence")
    raw_exit = native.get("wrapper_exit_code")
    if type(raw_exit) is not int or raw_exit != 0:
        raise ValueError("native wrapper exit is not a clean zero")
    runtime = _read_json(runtime_evidence, "runtime evidence")
    for field, expected in {
        "environment_id": environment_id,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
    }.items():
        if runtime.get(field) != expected:
            raise ValueError(f"runtime evidence {field} does not match")

    check, log = _scan_tree(check_root)
    log_bytes = log.read_bytes()
    try:
        text = log_bytes.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("00check.log must be UTF-8") from error
    if not text.endswith("* DONE\nStatus: OK\n"):
        raise ValueError("00check.log does not end in a strict successful footer")
    lines = text.splitlines()
    matched = {
        name: sum(bool(pattern.fullmatch(line)) for line in lines)
        for name, pattern in DOC_HEADINGS.items()
    }
    if profile == FULL_DOCUMENTATION:
        if any(count != 1 for count in matched.values()):
            raise ValueError("full documentation stages were not each executed exactly once")
    elif any(count for count in matched.values()) or any(
        line.startswith((
            "* checking PDF version of manual ...",
            "* checking package vignettes ...",
            "* checking running R code from vignettes ...",
            "* checking re-building of vignette outputs ...",
        ))
        for line in lines
    ):
        raise ValueError("no-documentation profile unexpectedly executed a doc stage")

    proof = {
        "schema_version": 1,
        "kind": "special-check-proof",
        "status": "pass",
        "environment_id": environment_id,
        "profile": profile,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "check_directory": os.fspath(check),
        "check_log_sha256": hashlib.sha256(log_bytes).hexdigest(),
        "manual_executed": matched["manual"] == 1,
        "vignettes_executed": matched["package"] == matched["rebuilding"] == 1,
        "raw_exit_code": raw_exit,
        "effective_exit_code": raw_exit,
        "reconciliation": None,
    }
    _atomic_write(output, proof)
    if github_output is not None:
        append_github_outputs(github_output, {"effective_exit_code": str(raw_exit)})
    return proof


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--profile", required=True, choices=sorted(PROFILES))
    parser.add_argument("--check-root", required=True)
    parser.add_argument("--native-evidence", required=True)
    parser.add_argument("--runtime-evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--github-output")
    args = parser.parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(args.manifest))
        rows = [row for row in manifest["coverage"] if row["id"] == args.environment_id]
        if len(rows) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        verify_special_check(
            rows[0], environment_id=args.environment_id, profile=args.profile,
            check_root=args.check_root, native_evidence=args.native_evidence,
            runtime_evidence=args.runtime_evidence, source_sha=args.source_sha,
            event_sha=args.event_sha, tarball_sha256=args.tarball_sha256,
            output=args.output, github_output=args.github_output,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"special check evidence error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
