#!/usr/bin/env python3
"""Run one verified source tarball through one exact cross-platform R binary."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from artifact_contract import _atomic_write_json, verify_tarball
from validate_manifest import load_manifest, validate_manifest


_ENVIRONMENT_ID = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
_NO_DOCUMENTATION_ARGS = ("--as-cran", "--no-manual", "--no-build-vignettes")


def _single_line(value: str, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty string")
    if any(character in value for character in ("\r", "\n", "\x00")):
        raise ValueError(f"{label} must be an unambiguous single-line value")
    return value


def _same_portable_path(left: str | os.PathLike[str], right: str | os.PathLike[str]) -> bool:
    """Compare lexical paths across Windows slash/case spellings."""

    left_text = os.fspath(left).replace("\\", "/")
    right_text = os.fspath(right).replace("\\", "/")
    windows = bool(re.match(r"^[A-Za-z]:/", left_text)) or bool(
        re.match(r"^[A-Za-z]:/", right_text)
    )
    if windows:
        left_text = left_text.casefold()
        right_text = right_text.casefold()
    return left_text == right_text


def _require_r_executable(value: str | os.PathLike[str]) -> tuple[Path, Path]:
    raw = os.fspath(value)
    _single_line(raw, "R executable path")
    lexical = Path(raw)
    if not lexical.is_absolute():
        raise ValueError("R executable path must be absolute")
    try:
        status = lexical.lstat()
    except OSError as error:
        raise ValueError("R executable must be an existing regular file") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("R executable must be a regular non-symlink file")
    if not os.access(lexical, os.X_OK):
        raise ValueError("R executable must be executable")
    if lexical.name.lower() not in {"r", "r.exe"}:
        raise ValueError("R executable basename must be R or R.exe")
    try:
        resolved = lexical.resolve(strict=True)
        resolved_status = resolved.stat()
    except OSError as error:
        raise ValueError("resolved R executable is not accessible") from error
    if not stat.S_ISREG(resolved_status.st_mode) or not os.access(resolved, os.X_OK):
        raise ValueError("resolved R executable must be a regular executable file")
    return lexical, resolved


def _fresh_directory(value: str | os.PathLike[str], label: str) -> Path:
    raw = os.fspath(value)
    _single_line(raw, label)
    path = Path(raw)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    if path.exists():
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
            raise ValueError(f"{label} must be a regular non-symlink directory")
        if any(path.iterdir()):
            raise ValueError(f"{label} must be initially empty")
    else:
        path.mkdir(parents=True)
    return path


def _prepared_library(value: str | os.PathLike[str]) -> Path:
    raw = os.fspath(value)
    _single_line(raw, "prepared check library")
    path = Path(raw)
    if not path.is_absolute():
        raise ValueError("prepared check library must be absolute")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError("prepared check library must already exist") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError("prepared check library must be a regular non-symlink directory")
    if not any(path.iterdir()):
        raise ValueError("prepared check library must be nonempty")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as output:
            temporary = Path(output.name)
            with source.open("rb") as incoming:
                shutil.copyfileobj(incoming, output)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _unique_check_log(result_directory: Path) -> tuple[Path, Path]:
    candidates = list(result_directory.glob("*.Rcheck"))
    if len(candidates) != 1:
        raise ValueError("R CMD check must produce exactly one *.Rcheck directory")
    check_directory = candidates[0]
    try:
        check_status = check_directory.lstat()
    except OSError as error:
        raise ValueError("R CMD check output is not accessible") from error
    if stat.S_ISLNK(check_status.st_mode) or not stat.S_ISDIR(check_status.st_mode):
        raise ValueError("*.Rcheck must be a regular non-symlink directory")
    check_log = check_directory / "00check.log"
    try:
        log_status = check_log.lstat()
    except OSError as error:
        raise ValueError(
            "R CMD check must produce exactly one *.Rcheck/00check.log"
        ) from error
    if stat.S_ISLNK(log_status.st_mode) or not stat.S_ISREG(log_status.st_mode):
        raise ValueError("00check.log must be a regular non-symlink file")
    return check_directory, check_log


def run_r_binary_check(
    manifest_row: Mapping[str, Any],
    *,
    environment_id: str,
    r_executable: str | os.PathLike[str],
    tarball: str | os.PathLike[str],
    metadata: str | os.PathLike[str],
    expected_source_sha: str,
    expected_event_sha: str,
    work_dir: str | os.PathLike[str],
    check_library: str | os.PathLike[str],
    evidence_path: str | os.PathLike[str],
    log_output: str | os.PathLike[str],
    environment: Mapping[str, str] | None = None,
    allow_full_documentation: bool = False,
) -> int:
    """Execute R CMD check and publish its identity-bound evidence and exact log."""

    if not isinstance(manifest_row, Mapping):
        raise ValueError("manifest row must be an object")
    selected_id = _single_line(environment_id, "environment ID")
    if _ENVIRONMENT_ID.fullmatch(selected_id) is None:
        raise ValueError("environment ID contains unsupported characters")
    if manifest_row.get("id") != selected_id:
        raise ValueError("environment ID does not match the manifest row")
    if manifest_row.get("driver") != "r-binary":
        raise ValueError("environment ID must select an r-binary manifest row")
    if not isinstance(allow_full_documentation, bool):
        raise ValueError("allow_full_documentation must be a boolean")

    lexical_r, resolved_r = _require_r_executable(r_executable)
    declared_r = manifest_row.get("system_r", manifest_row.get("r_executable"))
    if declared_r is not None and not _same_portable_path(declared_r, lexical_r):
        raise ValueError("R executable does not match the manifest row")
    metadata_document = verify_tarball(
        tarball,
        metadata,
        expected_source_sha=expected_source_sha,
        expected_event_sha=expected_event_sha,
    )

    work = _fresh_directory(work_dir, "check work directory")
    library = _prepared_library(check_library)
    if work == library or work in library.parents or library in work.parents:
        raise ValueError("check work directory and check library must be independent")
    input_directory = work / "input"
    result_directory = work / "results"
    input_directory.mkdir()
    result_directory.mkdir()
    copied_tarball = input_directory / metadata_document["filename"]
    _copy_atomic(Path(tarball), copied_tarball)
    if _sha256(copied_tarball) != metadata_document["sha256"]:
        raise ValueError("copied tarball digest does not match verified metadata")

    check_args = [] if allow_full_documentation else list(_NO_DOCUMENTATION_ARGS)
    declared_args = manifest_row.get("check_args")
    if declared_args is not None and declared_args != check_args:
        raise ValueError("R CMD check arguments do not match the manifest row")
    command = [
        os.fspath(lexical_r),
        "CMD",
        "check",
        *check_args,
        os.fspath(copied_tarball),
    ]
    process_environment = dict(os.environ if environment is None else environment)
    process_environment["R_LIBS_USER"] = os.fspath(library)
    process_environment["R_PROFILE_USER"] = os.devnull
    process_environment["R_ENVIRON_USER"] = os.devnull
    completed = subprocess.run(
        command,
        cwd=result_directory,
        env=process_environment,
        check=False,
        shell=False,
    )
    check_directory, check_log = _unique_check_log(result_directory)
    _copy_atomic(check_log, Path(log_output))

    r_status = lexical_r.stat()
    log_status = check_log.stat()
    proof: dict[str, Any] = {
        "schema_version": 1,
        "kind": "r-binary-check-proof",
        "status": "passed" if completed.returncode == 0 else "failed",
        "environment_id": selected_id,
        "source_sha": metadata_document["source_sha"],
        "event_sha": metadata_document["event_sha"],
        "tarball": {
            "filename": metadata_document["filename"],
            "size": metadata_document["size"],
            "sha256": metadata_document["sha256"],
            "verified_input_path": os.fspath(Path(tarball).absolute()),
            "copied_path": os.fspath(copied_tarball),
        },
        "r_executable": {
            "path": os.fspath(lexical_r),
            "resolved_path": os.fspath(resolved_r),
            "size": r_status.st_size,
            "sha256": _sha256(lexical_r),
        },
        "check_arguments": check_args,
        "full_documentation": allow_full_documentation,
        "check_library": os.fspath(library),
        "startup_environment": {
            "R_LIBS_USER": os.fspath(library),
            "R_PROFILE_USER": os.devnull,
            "R_ENVIRON_USER": os.devnull,
        },
        "check_directory": os.fspath(check_directory),
        "check_log": {
            "path": os.fspath(check_log),
            "size": log_status.st_size,
            "sha256": _sha256(check_log),
        },
        "exit_code": completed.returncode,
    }
    _atomic_write_json(Path(evidence_path), proof)
    return completed.returncode


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--r-executable", required=True)
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--check-library", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--allow-full-documentation", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(arguments.manifest))
        matches = [
            row
            for row in manifest["coverage"]
            if row["id"] == arguments.environment_id
        ]
        if len(matches) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        return run_r_binary_check(
            matches[0],
            environment_id=arguments.environment_id,
            r_executable=arguments.r_executable,
            tarball=arguments.tarball,
            metadata=arguments.metadata,
            expected_source_sha=arguments.source_sha,
            expected_event_sha=arguments.event_sha,
            work_dir=arguments.work_dir,
            check_library=arguments.check_library,
            evidence_path=arguments.evidence,
            log_output=arguments.log,
            allow_full_documentation=arguments.allow_full_documentation,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"r-binary check error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
