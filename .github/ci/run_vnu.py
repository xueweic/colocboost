#!/usr/bin/env python3
"""Run the checksum-bound R-hub VNU dispatcher on freshly extracted source."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

from artifact_contract import verify_tarball
from extract_source import extract_verified_source
from run_driver import _atomic_write_json, _require_executable, _sha256
from validate_manifest import load_manifest, validate_manifest


FILE_PROOF_FIELDS = {
    "schema_version", "kind", "status", "environment_id", "purpose",
    "source_sha", "event_sha", "tarball_sha256", "path", "size", "sha256",
    "executable", "suffix_bytes", "suffix_sha256",
}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _directory(path_value, label, *, nonempty=False):
    path = Path(path_value)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is missing") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink directory")
    if nonempty and not any(path.iterdir()):
        raise ValueError(f"{label} must be nonempty")
    return path


def _read_file_proof(path_value, expected):
    path = Path(path_value)
    try:
        status = path.lstat()
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("VNU dispatcher file proof is unreadable") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("VNU dispatcher file proof must be regular and non-symlink")
    if not isinstance(document, Mapping) or set(document) != FILE_PROOF_FIELDS:
        raise ValueError("VNU dispatcher file proof fields are not closed")
    if (
        type(document["schema_version"]) is not int
        or document["schema_version"] != 1
        or document["kind"] != "bound-file-proof"
        or document["status"] != "pass"
        or document["purpose"] != "vnu-dispatcher"
        or document["executable"] is not True
        or document["suffix_bytes"] is not None
        or document["suffix_sha256"] is not None
    ):
        raise ValueError("VNU dispatcher file proof is not terminal green")
    for field, value in expected.items():
        if document.get(field) != value:
            raise ValueError(f"VNU dispatcher file proof {field} does not match")
    return document


def run_vnu(
    manifest_row: Mapping,
    *,
    environment_id: str,
    dispatcher,
    r_executable,
    tarball,
    metadata,
    source_sha: str,
    event_sha: str,
    tarball_sha256: str,
    work_dir,
    library,
    file_evidence,
    output,
    environment: Mapping[str, str] | None = None,
):
    if manifest_row.get("id") != environment_id or environment_id != "vnu":
        raise ValueError("VNU environment does not match manifest")
    metadata_document = verify_tarball(
        tarball, metadata, expected_source_sha=source_sha,
        expected_event_sha=event_sha,
    )
    if metadata_document["sha256"] != tarball_sha256:
        raise ValueError("VNU tarball digest does not match verified metadata")
    selected_r = _require_executable(r_executable, label="selected VNU R")
    selected_dispatcher = _require_executable(dispatcher, label="VNU dispatcher")
    if os.fspath(selected_r) != manifest_row.get("system_r"):
        raise ValueError("selected VNU R does not match manifest")
    if os.fspath(selected_dispatcher) != manifest_row.get("vnu_path"):
        raise ValueError("VNU dispatcher path does not match manifest")
    dispatcher_hash = _sha256(selected_dispatcher)
    if dispatcher_hash != manifest_row.get("vnu_sha256"):
        raise ValueError("VNU dispatcher sha256 does not match manifest")
    _read_file_proof(file_evidence, {
        "environment_id": environment_id,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "path": os.fspath(selected_dispatcher),
        "sha256": dispatcher_hash,
    })

    selected_library = _directory(library, "VNU check library", nonempty=True)
    target = Path(work_dir)
    if not target.is_absolute():
        raise ValueError("VNU source work directory must be absolute")
    target.parent.mkdir(parents=True, exist_ok=True)
    output_file = target.parent / f".{target.name}-extract-output"
    output_file.unlink(missing_ok=True)
    package_root = extract_verified_source(
        os.fspath(tarball), os.fspath(metadata), source_sha=source_sha,
        event_sha=event_sha, destination=os.fspath(target),
        github_output=os.fspath(output_file),
    )
    output_file.unlink(missing_ok=True)

    selected_environment = dict(os.environ if environment is None else environment)
    if selected_environment.get("R_LIBS_USER") != os.fspath(selected_library):
        raise ValueError("R_LIBS_USER must name the VNU check library")
    path_r = shutil.which("R", path=selected_environment.get("PATH"))
    if path_r is None or Path(os.path.abspath(path_r)) != selected_r:
        raise ValueError("PATH must resolve the selected VNU R")
    selected_environment["R_PROFILE_USER"] = os.devnull
    selected_environment["R_ENVIRON_USER"] = os.devnull
    completed = subprocess.run(
        [os.fspath(selected_dispatcher)], cwd=package_root,
        env=selected_environment, check=False, shell=False,
        capture_output=True, text=True, encoding="utf-8", errors="strict",
    )
    if completed.returncode != 0:
        raise ValueError(f"VNU dispatcher exited {completed.returncode}")
    if not completed.stdout.strip() and not completed.stderr.strip():
        raise ValueError("VNU dispatcher produced no validator transcript")
    html = package_root / "pkg.html"
    try:
        html_status = html.lstat()
    except OSError as error:
        raise ValueError("VNU dispatcher did not create pkg.html") from error
    if (
        stat.S_ISLNK(html_status.st_mode)
        or not stat.S_ISREG(html_status.st_mode)
        or html_status.st_size <= 0
    ):
        raise ValueError("VNU pkg.html must be a nonempty regular non-symlink file")
    html_data = html.read_bytes()
    proof = {
        "schema_version": 1,
        "kind": "vnu-proof",
        "status": "pass",
        "environment_id": environment_id,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "r_executable": os.fspath(selected_r),
        "r_resolved": os.fspath(selected_r.resolve(strict=True)),
        "dispatcher": os.fspath(selected_dispatcher),
        "dispatcher_sha256": dispatcher_hash,
        "dispatcher_exit_code": completed.returncode,
        "package_root": os.fspath(package_root),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "pkg_html": os.fspath(html),
        "pkg_html_size": len(html_data),
        "pkg_html_sha256": hashlib.sha256(html_data).hexdigest(),
        "zero_bad_entries": True,
    }
    _atomic_write_json(Path(output), proof)
    return proof


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--dispatcher", required=True)
    parser.add_argument("--r-executable", required=True)
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--file-evidence", required=True)
    parser.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(args.manifest))
        rows = [row for row in manifest["coverage"] if row["id"] == args.environment_id]
        if len(rows) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        run_vnu(
            rows[0], environment_id=args.environment_id,
            dispatcher=args.dispatcher, r_executable=args.r_executable,
            tarball=args.tarball, metadata=args.metadata, source_sha=args.source_sha,
            event_sha=args.event_sha, tarball_sha256=args.tarball_sha256,
            work_dir=args.work_dir, library=args.library,
            file_evidence=args.file_evidence, output=args.output,
        )
    except (OSError, TypeError, UnicodeError, ValueError) as error:
        print(f"VNU proof error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
