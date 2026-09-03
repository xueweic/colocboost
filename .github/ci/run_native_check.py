#!/usr/bin/env python3
"""Run one R-hub native wrapper in a fresh, identity-bound check tree."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import stat
import sys
from collections.abc import Mapping
from pathlib import Path

from artifact_contract import verify_tarball
from run_driver import _atomic_write_json, run_driver
from validate_manifest import load_manifest, validate_manifest


def _fresh_directory(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    path = path.absolute()
    if path.exists():
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
            raise ValueError(f"{label} must be a regular non-symlink directory")
        if any(path.iterdir()):
            raise ValueError(f"{label} must be initially empty")
    else:
        path.mkdir(parents=True)
    return path


def _regular_unique(paths: list[Path], label: str, *, nonempty: bool = False) -> Path:
    if len(paths) != 1:
        raise ValueError(f"native wrapper must produce exactly one {label}")
    path = paths[0]
    try:
        status = path.lstat()
        parent_status = path.parent.lstat()
    except OSError as error:
        raise ValueError(f"native wrapper {label} is not accessible") from error
    if (
        stat.S_ISLNK(status.st_mode)
        or not stat.S_ISREG(status.st_mode)
        or stat.S_ISLNK(parent_status.st_mode)
        or not stat.S_ISDIR(parent_status.st_mode)
    ):
        raise ValueError(f"native wrapper {label} must be a regular non-symlink file")
    if nonempty and status.st_size == 0:
        raise ValueError(f"native wrapper {label} must be nonzero")
    return path


def _copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        with source.open("rb") as incoming, temporary.open("xb") as outgoing:
            shutil.copyfileobj(incoming, outgoing)
            outgoing.flush()
            os.fsync(outgoing.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def run_native_check(
    manifest_row: Mapping,
    *,
    wrapper: str | os.PathLike[str],
    required_r_executable: str | os.PathLike[str],
    tarball: str | os.PathLike[str],
    metadata: str | os.PathLike[str],
    expected_source_sha: str,
    expected_event_sha: str,
    work_dir: str | os.PathLike[str],
    check_library: str | os.PathLike[str],
    evidence_path: str | os.PathLike[str],
    log_output: str | os.PathLike[str],
    environment: Mapping[str, str] | None = None,
) -> int:
    verify_tarball(
        tarball,
        metadata,
        expected_source_sha=expected_source_sha,
        expected_event_sha=expected_event_sha,
    )
    work = _fresh_directory(work_dir, "check work directory")
    library = _fresh_directory(check_library, "check library")
    if work == library or work in library.parents or library in work.parents:
        raise ValueError("check work directory and check library must be independent")
    selected_environment = dict(os.environ if environment is None else environment)
    if selected_environment.get("R_LIBS_USER") != os.fspath(library):
        raise ValueError("R_LIBS_USER must name the fresh check library")
    selected_environment["R_PROFILE_USER"] = os.devnull
    selected_environment["R_ENVIRON_USER"] = os.devnull

    input_dir = work / "input"
    input_dir.mkdir()
    copied_tarball = input_dir / Path(tarball).name
    shutil.copyfile(tarball, copied_tarball, follow_symlinks=False)
    exit_code = run_driver(
        manifest_row,
        requested_driver="native-wrapper",
        executable=wrapper,
        argv=["{tarball-parent}"],
        tarball=copied_tarball,
        metadata=metadata,
        expected_source_sha=expected_source_sha,
        expected_event_sha=expected_event_sha,
        required_r_executable=required_r_executable,
        evidence_path=evidence_path,
        environment=selected_environment,
        cwd=work,
    )
    evidence_file = Path(evidence_path)
    try:
        evidence = json.loads(evidence_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("native wrapper evidence is unreadable") from error
    evidence["wrapper_exit_code"] = exit_code
    _atomic_write_json(evidence_file, evidence)

    check_log = _regular_unique(
        list(input_dir.glob("*.Rcheck/00check.log")), "*.Rcheck/00check.log"
    )
    evidence["check_log"] = os.fspath(check_log)
    if manifest_row.get("id") == "nosuggests":
        rout = _regular_unique(
            list(input_dir.glob("*.Rcheck/tests/testthat.Rout")),
            "*.Rcheck/tests/testthat.Rout",
            nonempty=True,
        )
        evidence["testthat_rout"] = {
            "path": os.fspath(rout),
            "size": rout.stat().st_size,
        }
    _atomic_write_json(evidence_file, evidence)
    _copy_atomic(check_log, Path(log_output))
    return exit_code


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--wrapper", required=True)
    parser.add_argument("--required-r-executable", required=True)
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--check-library", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--log", required=True)
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
        if len(matches) != 1 or matches[0]["driver"] != "native-wrapper":
            raise ValueError("environment-id must select one native-wrapper row")
        return run_native_check(
            matches[0],
            wrapper=arguments.wrapper,
            required_r_executable=arguments.required_r_executable,
            tarball=arguments.tarball,
            metadata=arguments.metadata,
            expected_source_sha=arguments.source_sha,
            expected_event_sha=arguments.event_sha,
            work_dir=arguments.work_dir,
            check_library=arguments.check_library,
            evidence_path=arguments.evidence,
            log_output=arguments.log,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"native check error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
