#!/usr/bin/env python3
"""Verify immutable active R-hub helper files and emit identity-bound proof."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
from collections.abc import Mapping
from pathlib import Path

from run_driver import _atomic_write_json
from validate_manifest import load_manifest, validate_manifest


GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
PURPOSES = {"valgrind-suppression", "vnu-dispatcher"}


def _hash(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _regular(path_value, label: str, *, executable: bool = False) -> Path:
    path = Path(path_value)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is missing") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file")
    if executable and not os.access(path, os.X_OK):
        raise ValueError(f"{label} must be executable")
    return path


def _identity(value: str, pattern: re.Pattern[str], label: str) -> None:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ValueError(f"{label} has invalid identity syntax")


def verify_bound_file(
    manifest_row: Mapping,
    *,
    environment_id: str,
    purpose: str,
    path,
    source_sha: str,
    event_sha: str,
    tarball_sha256: str,
    output,
    reference=None,
):
    if purpose not in PURPOSES:
        raise ValueError("unsupported file contract purpose")
    if manifest_row.get("id") != environment_id:
        raise ValueError("file contract environment does not match manifest")
    _identity(source_sha, GIT_SHA, "source_sha")
    _identity(event_sha, GIT_SHA, "event_sha")
    _identity(tarball_sha256, SHA256, "tarball_sha256")

    if purpose == "vnu-dispatcher":
        expected_path = manifest_row.get("vnu_path")
        expected_hash = manifest_row.get("vnu_sha256")
        candidate = _regular(path, "vnu dispatcher", executable=True)
        if os.fspath(candidate) != expected_path:
            raise ValueError("vnu dispatcher path does not match manifest")
        data = candidate.read_bytes()
        digest = _hash(data)
        if digest != expected_hash:
            raise ValueError("vnu dispatcher sha256 does not match manifest")
        suffix_bytes = None
        suffix_hash = None
        is_executable = True
    else:
        expected_path = manifest_row.get("suppression_path")
        suffix_bytes = manifest_row.get("suppression_suffix_bytes")
        expected_hash = manifest_row.get("suppression_suffix_sha256")
        if type(suffix_bytes) is not int or suffix_bytes <= 0:
            raise ValueError("suppression suffix byte count is invalid")
        candidate = _regular(path, "Valgrind suppression")
        if os.fspath(candidate) != expected_path:
            raise ValueError("Valgrind suppression path does not match manifest")
        reference_path = _regular(
            reference
            if reference is not None
            else Path(__file__).with_name("fixtures") / "rhub-valgrind.supp",
            "committed Valgrind suppression reference",
        )
        reference_data = reference_path.read_bytes()
        if len(reference_data) != suffix_bytes or _hash(reference_data) != expected_hash:
            raise ValueError("committed Valgrind suppression reference does not match manifest")
        data = candidate.read_bytes()
        if len(data) < suffix_bytes or data[-suffix_bytes:] != reference_data:
            raise ValueError("Valgrind default suppression lacks the exact committed suffix")
        digest = _hash(data)
        suffix_hash = _hash(data[-suffix_bytes:])
        if suffix_hash != expected_hash:
            raise ValueError("Valgrind suppression suffix sha256 does not match manifest")
        is_executable = bool(os.access(candidate, os.X_OK))

    proof = {
        "schema_version": 1,
        "kind": "bound-file-proof",
        "status": "pass",
        "environment_id": environment_id,
        "purpose": purpose,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "path": os.fspath(candidate),
        "size": len(data),
        "sha256": digest,
        "executable": is_executable,
        "suffix_bytes": suffix_bytes,
        "suffix_sha256": suffix_hash,
    }
    _atomic_write_json(Path(output), proof)
    return proof


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--purpose", required=True, choices=sorted(PURPOSES))
    parser.add_argument("--path", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(args.manifest))
        rows = [row for row in manifest["coverage"] if row["id"] == args.environment_id]
        if len(rows) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        verify_bound_file(
            rows[0], environment_id=args.environment_id, purpose=args.purpose,
            path=args.path, source_sha=args.source_sha, event_sha=args.event_sha,
            tarball_sha256=args.tarball_sha256, output=args.output,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"file contract error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
