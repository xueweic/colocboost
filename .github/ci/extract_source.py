#!/usr/bin/env python3
"""Safely extract the one verified source package used by source-mode tests."""

from __future__ import annotations

import argparse
import os
import shutil
import stat
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

from artifact_contract import append_github_outputs, verify_tarball


def _destination(value: str) -> Path:
    if not value or any(character in value for character in ("\r", "\n", "\x00")):
        raise ValueError("destination must be a nonempty unambiguous path")
    path = Path(value)
    if not path.is_absolute():
        raise ValueError("destination must be absolute")
    path = path.absolute()
    try:
        parent = path.parent.resolve(strict=True)
    except OSError as error:
        raise ValueError("destination parent must exist") from error
    if path != parent / path.name:
        raise ValueError("destination contains ambiguous parent components")
    if path.exists():
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
            raise ValueError("destination must be a regular non-symlink directory")
        if any(path.iterdir()):
            raise ValueError("destination must be initially empty")
    return path


def _validated_members(archive: tarfile.TarFile) -> tuple[list[tarfile.TarInfo], str]:
    members = archive.getmembers()
    if not members:
        raise ValueError("source tarball is empty")
    seen = set()
    roots = set()
    description_found = False
    for member in members:
        name = member.name
        if not name or "\\" in name or "\x00" in name:
            raise ValueError("source tarball contains an invalid member name")
        path = PurePosixPath(name)
        if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
            raise ValueError("source tarball contains path traversal")
        normalized = path.as_posix()
        if normalized in seen:
            raise ValueError("source tarball contains duplicate members")
        seen.add(normalized)
        roots.add(path.parts[0])
        if not (member.isdir() or member.isfile()):
            raise ValueError("source tarball contains a link or special file")
        if len(path.parts) == 2 and path.parts[1] == "DESCRIPTION":
            if not member.isfile():
                raise ValueError("package DESCRIPTION must be a regular file")
            description_found = True
    if len(roots) != 1:
        raise ValueError("source tarball must contain exactly one top-level directory")
    if not description_found:
        raise ValueError("source tarball lacks a top-level package DESCRIPTION")
    return members, next(iter(roots))


def extract_verified_source(
    tarball: str,
    metadata: str,
    *,
    source_sha: str,
    event_sha: str,
    destination: str,
    github_output: str,
) -> Path:
    verify_tarball(
        tarball,
        metadata,
        expected_source_sha=source_sha,
        expected_event_sha=event_sha,
    )
    target = _destination(destination)
    temporary = Path(tempfile.mkdtemp(prefix=f".{target.name}.", dir=target.parent))
    committed = False
    try:
        with tarfile.open(tarball, mode="r:gz") as archive:
            members, root = _validated_members(archive)
            for member in members:
                output = temporary.joinpath(*PurePosixPath(member.name).parts)
                if member.isdir():
                    output.mkdir(parents=True, exist_ok=True)
                    continue
                output.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    raise ValueError(f"could not read source member {member.name}")
                with source, output.open("xb") as stream:
                    shutil.copyfileobj(source, stream)
                output.chmod(member.mode & 0o777)
        package = temporary / root
        if target.exists():
            target.rmdir()
        os.replace(temporary, target)
        committed = True
        package = target / root
        append_github_outputs(github_output, {"package_path": os.fspath(package)})
        return package
    except (OSError, tarfile.TarError) as error:
        raise ValueError(f"could not safely extract source tarball: {error}") from error
    finally:
        if not committed:
            shutil.rmtree(temporary, ignore_errors=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--github-output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        package = extract_verified_source(
            arguments.tarball,
            arguments.metadata,
            source_sha=arguments.source_sha,
            event_sha=arguments.event_sha,
            destination=arguments.destination,
            github_output=arguments.github_output,
        )
    except (OSError, ValueError) as error:
        print(f"source extraction error: {error}", file=sys.stderr)
        return 2
    print(os.fspath(package))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
