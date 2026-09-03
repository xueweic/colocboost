"""Create and verify metadata for the single CRAN preflight source tarball."""

from __future__ import annotations

import hashlib
import argparse
import json
import os
import re
import stat
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any


_METADATA_FIELDS = frozenset(
    {"source_sha", "event_sha", "filename", "size", "sha256"}
)
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _validate_git_sha(value: Any, field: str) -> None:
    if not isinstance(value, str) or _GIT_SHA.fullmatch(value) is None:
        raise ValueError(f"{field} must be a 40-character lowercase Git SHA")


def _validate_filename(value: Any) -> None:
    if not isinstance(value, str) or not value:
        raise ValueError("filename must be a non-empty basename")
    if value in {".", ".."}:
        raise ValueError("filename must not contain traversal")
    if PurePosixPath(value).name != value or PureWindowsPath(value).name != value:
        raise ValueError("filename must be a basename without traversal")


def validate_metadata(document: Mapping[str, Any]) -> Mapping[str, Any]:
    """Validate a closed source-tarball metadata document."""

    if not isinstance(document, Mapping):
        raise ValueError("metadata must be an object")
    missing = _METADATA_FIELDS - set(document)
    if missing:
        raise ValueError(f"metadata is missing {sorted(missing)[0]}")
    extra = set(document) - _METADATA_FIELDS
    if extra:
        raise ValueError(f"metadata contains unexpected property {sorted(extra)[0]}")

    _validate_git_sha(document["source_sha"], "source_sha")
    _validate_git_sha(document["event_sha"], "event_sha")
    _validate_filename(document["filename"])

    size = document["size"]
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise ValueError("size must be a nonnegative integer")
    digest = document["sha256"]
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise ValueError(
            "sha256 must be exactly 64 lowercase hexadecimal characters"
        )
    return document


def _require_regular_tarball(path: Path) -> os.stat_result:
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"tarball must be an existing regular file: {path}") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"tarball must be a regular non-symlink file: {path}")
    return status


def _require_regular_metadata(path: Path) -> os.stat_result:
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(
            f"metadata must be an existing regular non-symlink file: {path}"
        ) from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(
            f"metadata must be a regular non-symlink file: {path}"
        )
    return status


def _hash_file(path: Path) -> tuple[int, str]:
    _require_regular_tarball(path)
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(document, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def create_metadata(
    tarball: str | os.PathLike[str],
    metadata_path: str | os.PathLike[str],
    *,
    source_sha: str,
    event_sha: str,
) -> dict[str, Any]:
    """Hash raw tarball bytes and atomically write their source metadata."""

    _validate_git_sha(source_sha, "source_sha")
    _validate_git_sha(event_sha, "event_sha")
    tarball_path = Path(tarball)
    size, digest = _hash_file(tarball_path)
    document = {
        "source_sha": source_sha,
        "event_sha": event_sha,
        "filename": tarball_path.name,
        "size": size,
        "sha256": digest,
    }
    validate_metadata(document)
    _atomic_write_json(Path(metadata_path), document)
    return document


def verify_tarball(
    tarball: str | os.PathLike[str],
    metadata_path: str | os.PathLike[str],
    *,
    expected_source_sha: str,
    expected_event_sha: str,
) -> dict[str, Any]:
    """Recompute raw-byte integrity and match the expected workflow identity."""

    _validate_git_sha(expected_source_sha, "expected_source_sha")
    _validate_git_sha(expected_event_sha, "expected_event_sha")
    tarball_path = Path(tarball)
    _require_regular_tarball(tarball_path)
    metadata_file = Path(metadata_path)
    _require_regular_metadata(metadata_file)
    try:
        document = json.loads(metadata_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"metadata is not readable valid JSON: {metadata_path}") from error
    validate_metadata(document)

    if document["source_sha"] != expected_source_sha:
        raise ValueError("metadata source_sha does not match expected source_sha")
    if document["event_sha"] != expected_event_sha:
        raise ValueError("metadata event_sha does not match expected event_sha")
    if document["filename"] != tarball_path.name:
        raise ValueError("metadata filename does not match the tarball filename")

    size, digest = _hash_file(tarball_path)
    if document["size"] != size:
        raise ValueError("metadata size does not match the tarball bytes")
    if document["sha256"] != digest:
        raise ValueError("metadata sha256 does not match the tarball bytes")
    return dict(document)


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create or verify raw-byte source-tarball metadata."
    )
    parser.add_argument("operation", choices=("create", "verify"))
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _argument_parser().parse_args(argv)
    try:
        if args.operation == "create":
            document = create_metadata(
                args.tarball,
                args.metadata,
                source_sha=args.source_sha,
                event_sha=args.event_sha,
            )
        else:
            document = verify_tarball(
                args.tarball,
                args.metadata,
                expected_source_sha=args.source_sha,
                expected_event_sha=args.event_sha,
            )
    except ValueError as error:
        print(f"artifact contract error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
