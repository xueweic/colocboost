#!/usr/bin/env python3
"""Discover and verify the exact downloaded preflight source contract."""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from pathlib import Path

from artifact_contract import (
    append_github_outputs,
    github_output_values,
    verify_tarball,
)


_METADATA_NAME = "source-metadata.json"


def _require_unambiguous_directory(value: str) -> Path:
    if not value or any(character in value for character in ("\r", "\n", "\x00")):
        raise ValueError("source directory must be a nonempty unambiguous path")
    directory = Path(value).absolute()
    try:
        status = directory.lstat()
    except OSError as error:
        raise ValueError("source directory must be an existing directory") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError("source directory must be a regular non-symlink directory")
    return directory


def _require_regular_non_symlink(path: Path, label: str) -> None:
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} must be an existing regular non-symlink file") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file")


def verify_source_directory(
    source_directory: str,
    *,
    source_sha: str,
    event_sha: str,
    github_output: str,
) -> dict:
    """Verify exactly one metadata file and one raw source tarball."""

    directory = _require_unambiguous_directory(source_directory)
    try:
        entries = list(directory.iterdir())
    except OSError as error:
        raise ValueError("source directory is not readable") from error

    tarballs = [entry for entry in entries if entry.name.endswith(".tar.gz")]
    if len(tarballs) != 1:
        raise ValueError("source directory must contain exactly one .tar.gz file")
    tarball = tarballs[0]
    _require_regular_non_symlink(tarball, "source tarball")

    metadata = directory / _METADATA_NAME
    _require_regular_non_symlink(metadata, "source metadata")
    expected_names = {tarball.name, _METADATA_NAME}
    if {entry.name for entry in entries} != expected_names or len(entries) != 2:
        raise ValueError(
            "source directory must contain exactly one tarball and its metadata"
        )

    document = verify_tarball(
        tarball,
        metadata,
        expected_source_sha=source_sha,
        expected_event_sha=event_sha,
    )
    append_github_outputs(
        github_output,
        github_output_values(document, tarball, metadata),
    )
    return document


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Discover and verify one downloaded source artifact contract."
    )
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--github-output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _argument_parser().parse_args(argv)
    try:
        document = verify_source_directory(
            arguments.source_dir,
            source_sha=arguments.source_sha,
            event_sha=arguments.event_sha,
            github_output=arguments.github_output,
        )
    except (OSError, ValueError) as error:
        print(f"source verification error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
