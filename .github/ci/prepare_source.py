#!/usr/bin/env python3
"""Build the one preflight source tarball from an exact clean Git commit."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

from artifact_contract import (
    append_github_outputs,
    create_metadata,
    github_output_values,
)


_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def _require_no_controls(value: str, label: str) -> None:
    if not value or any(character in value for character in ("\r", "\n", "\x00")):
        raise ValueError(f"{label} must be a nonempty unambiguous path")


def _require_repository(value: str) -> Path:
    _require_no_controls(value, "repository")
    lexical = Path(value).absolute()
    try:
        status = lexical.lstat()
    except OSError as error:
        raise ValueError("repository must be an existing directory") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError("repository must be a regular non-symlink directory")
    repository = lexical.resolve(strict=True)
    completed = subprocess.run(
        ["git", "-C", os.fspath(repository), "rev-parse", "--show-toplevel"],
        check=False,
        capture_output=True,
        text=True,
        shell=False,
    )
    if completed.returncode != 0:
        raise ValueError("repository is not a Git working tree")
    reported_lines = completed.stdout.splitlines()
    if len(reported_lines) != 1:
        raise ValueError("Git reported an ambiguous repository root")
    try:
        reported_root = Path(reported_lines[0]).resolve(strict=True)
    except OSError as error:
        raise ValueError("Git repository root is not accessible") from error
    if reported_root != repository:
        raise ValueError("repository must name the exact Git working-tree root")
    return repository


def _require_output_directory(value: str, repository: Path) -> Path:
    _require_no_controls(value, "output directory")
    lexical = Path(value).absolute()
    try:
        status = lexical.lstat()
    except OSError as error:
        raise ValueError("output directory must already exist") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError("output directory must be a regular non-symlink directory")
    output = lexical.resolve(strict=True)
    if output == repository or output.is_relative_to(repository):
        raise ValueError("output directory must be outside the repository")
    try:
        if any(output.iterdir()):
            raise ValueError("output directory must initially be empty")
    except OSError as error:
        raise ValueError("output directory is not readable") from error
    return output


def _git_output(repository: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", os.fspath(repository), *arguments],
        check=False,
        capture_output=True,
        text=True,
        shell=False,
    )
    if completed.returncode != 0:
        raise ValueError(f"git {' '.join(arguments)} failed")
    return completed.stdout


def _require_exact_clean_head(repository: Path, event_sha: str) -> str:
    if _GIT_SHA.fullmatch(event_sha) is None:
        raise ValueError("event SHA must be a 40-character lowercase Git SHA")
    head_lines = _git_output(repository, "rev-parse", "HEAD").splitlines()
    if len(head_lines) != 1 or _GIT_SHA.fullmatch(head_lines[0]) is None:
        raise ValueError("Git HEAD is not a single full commit SHA")
    head = head_lines[0]
    if head != event_sha:
        raise ValueError("Git HEAD does not match the workflow event SHA")
    dirty = _git_output(
        repository, "status", "--porcelain=v1", "--untracked-files=all"
    )
    if dirty:
        raise ValueError("repository must be completely clean before source export")
    return head


def _require_pixi_r() -> tuple[Path, Path]:
    prefix_value = os.environ.get("CONDA_PREFIX")
    if not prefix_value:
        raise ValueError("CONDA_PREFIX must identify the active Pixi environment")
    _require_no_controls(prefix_value, "Pixi prefix")
    try:
        prefix = Path(prefix_value).resolve(strict=True)
    except OSError as error:
        raise ValueError("Pixi prefix does not exist") from error
    r_from_path = shutil.which("R")
    if r_from_path is None:
        raise ValueError("R is unavailable in the active Pixi environment")
    try:
        r_binary = Path(r_from_path).resolve(strict=True)
        status = r_binary.stat()
    except OSError as error:
        raise ValueError("Pixi R executable is not accessible") from error
    if not stat.S_ISREG(status.st_mode) or not os.access(r_binary, os.X_OK):
        raise ValueError("Pixi R must be a regular executable file")
    if not r_binary.is_relative_to(prefix):
        raise ValueError("R executable is outside the active Pixi prefix")

    completed = subprocess.run(
        [os.fspath(r_binary), "RHOME"],
        check=False,
        capture_output=True,
        text=True,
        shell=False,
    )
    r_home_lines = completed.stdout.splitlines()
    if completed.returncode != 0 or len(r_home_lines) != 1 or not r_home_lines[0]:
        raise ValueError("Pixi R did not report one unambiguous RHOME")
    try:
        r_home = Path(r_home_lines[0]).resolve(strict=True)
    except OSError as error:
        raise ValueError("Pixi RHOME does not exist") from error
    if not r_home.is_dir() or not r_home.is_relative_to(prefix):
        raise ValueError("RHOME is outside the active Pixi prefix")
    return r_binary, prefix


def _require_pixi_pandoc(prefix: Path) -> Path:
    pandoc_from_path = shutil.which("pandoc")
    if pandoc_from_path is None:
        raise ValueError("Pandoc is unavailable in the active Pixi environment")
    try:
        pandoc_binary = Path(pandoc_from_path).resolve(strict=True)
        status = pandoc_binary.stat()
    except OSError as error:
        raise ValueError("Pixi Pandoc executable is not accessible") from error
    if not stat.S_ISREG(status.st_mode) or not os.access(pandoc_binary, os.X_OK):
        raise ValueError("Pixi Pandoc must be a regular executable file")
    if not pandoc_binary.is_relative_to(prefix):
        raise ValueError("Pandoc executable is outside the active Pixi prefix")
    return pandoc_binary.parent


def _export_head(repository: Path, head: str, destination: Path) -> None:
    archive_path = destination.parent / "source.tar"
    with archive_path.open("wb") as archive:
        completed = subprocess.run(
            [
                "git",
                "-C",
                os.fspath(repository),
                "archive",
                "--format=tar",
                head,
            ],
            check=False,
            stdout=archive,
            stderr=subprocess.PIPE,
            shell=False,
        )
    if completed.returncode != 0:
        raise ValueError("git archive HEAD failed")
    destination.mkdir()
    with tarfile.open(archive_path, mode="r:") as source_archive:
        for member in source_archive.getmembers():
            member_path = PurePosixPath(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise ValueError("git archive contains an unsafe member path")
        source_archive.extractall(destination, filter="data")


def _build_tarball(
    r_binary: Path,
    pandoc_directory: Path,
    source: Path,
    build_directory: Path,
) -> Path:
    environment = os.environ.copy()
    environment["RSTUDIO_PANDOC"] = os.fspath(pandoc_directory)
    completed = subprocess.run(
        [
            os.fspath(r_binary),
            "--vanilla",
            "CMD",
            "build",
            os.fspath(source),
        ],
        cwd=build_directory,
        check=False,
        env=environment,
        shell=False,
    )
    if completed.returncode != 0:
        raise ValueError(f"R CMD build failed with exit code {completed.returncode}")
    candidates = list(build_directory.glob("*.tar.gz"))
    if len(candidates) != 1:
        raise ValueError(
            "R CMD build must produce exactly one regular non-symlink .tar.gz"
        )
    candidate = candidates[0]
    try:
        status = candidate.lstat()
    except OSError as error:
        raise ValueError(
            "R CMD build must produce exactly one regular non-symlink .tar.gz"
        ) from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(
            "R CMD build must produce exactly one regular non-symlink .tar.gz"
        )
    return candidate


def prepare_source(
    repository_value: str,
    output_value: str,
    event_sha: str,
    github_output: str,
) -> dict:
    repository = _require_repository(repository_value)
    output = _require_output_directory(output_value, repository)
    head = _require_exact_clean_head(repository, event_sha)
    r_binary, pixi_prefix = _require_pixi_r()
    pandoc_directory = _require_pixi_pandoc(pixi_prefix)

    output_tarball: Path | None = None
    metadata_path = output / "source-metadata.json"
    try:
        with tempfile.TemporaryDirectory(
            prefix=".prepare-source-", dir=output.parent
        ) as temporary_value:
            temporary = Path(temporary_value)
            source = temporary / "source"
            build_directory = temporary / "build"
            build_directory.mkdir()
            _export_head(repository, head, source)
            built_tarball = _build_tarball(
                r_binary, pandoc_directory, source, build_directory
            )
            output_tarball = output / built_tarball.name
            os.replace(built_tarball, output_tarball)

        document = create_metadata(
            output_tarball,
            metadata_path,
            source_sha=head,
            event_sha=event_sha,
        )
        append_github_outputs(
            github_output,
            github_output_values(document, output_tarball, metadata_path),
        )
        return document
    except Exception:
        metadata_path.unlink(missing_ok=True)
        if output_tarball is not None:
            output_tarball.unlink(missing_ok=True)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build and describe one clean Git-archive source tarball."
    )
    parser.add_argument("--repository", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--github-output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _argument_parser().parse_args(argv)
    try:
        prepare_source(
            arguments.repository,
            arguments.output_dir,
            arguments.event_sha,
            arguments.github_output,
        )
    except (OSError, subprocess.SubprocessError, tarfile.TarError, ValueError) as error:
        print(f"source preparation error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
