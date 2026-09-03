"""Validation and atomic writes for CRAN preflight result documents."""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA_PATH = Path(__file__).with_name("result.schema.json")
_SCHEMA = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
_COMMON_REQUIRED = tuple(_SCHEMA["required"])
_ALLOWED_PROPERTIES = frozenset(_SCHEMA["properties"])
_RESULT_KINDS = frozenset(_SCHEMA["properties"]["result_kind"]["enum"])
_STATUSES = frozenset(_SCHEMA["properties"]["status"]["enum"])
_SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
_DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_DETAIL_FIELDS = {
    "unit": ("tests",),
    "package-check": ("check",),
    "applicability": ("predicate", "evidence"),
    "infrastructure": ("classification", "message"),
}


def _require_string(document: Mapping[str, Any], field: str) -> None:
    value = document.get(field)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{field} must be a non-empty string")


def _validate_count_object(
    document: Mapping[str, Any], field: str, required: tuple[str, ...]
) -> None:
    value = document.get(field)
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")

    missing = set(required) - set(value)
    if missing:
        raise ValueError(f"{field} is missing {sorted(missing)[0]}")
    extra = set(value) - set(required)
    if extra:
        raise ValueError(f"{field} contains unexpected property {sorted(extra)[0]}")

    for count_name in required:
        if count_name == "log_path":
            if not isinstance(value[count_name], str) or not value[count_name]:
                raise ValueError(f"{field}.{count_name} must be a non-empty string")
            continue
        count = value[count_name]
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"{field}.{count_name} must be a nonnegative integer")


def validate_result(document: Mapping[str, Any]) -> Mapping[str, Any]:
    """Validate *document* against the committed result contract.

    This standard-library validator mirrors ``result.schema.json`` so CI callers
    do not need to install a third-party JSON Schema implementation merely to
    emit or aggregate result artifacts.
    """

    if not isinstance(document, Mapping):
        raise ValueError("result must be an object")

    missing = set(_COMMON_REQUIRED) - set(document)
    if missing:
        raise ValueError(f"result is missing {sorted(missing)[0]}")
    extra = set(document) - _ALLOWED_PROPERTIES
    if extra:
        raise ValueError(f"result contains unexpected property {sorted(extra)[0]}")

    schema_version = document["schema_version"]
    if isinstance(schema_version, bool) or schema_version != 1:
        raise ValueError("schema_version must be 1")

    result_kind = document["result_kind"]
    if result_kind not in _RESULT_KINDS:
        raise ValueError(f"result_kind is invalid: {result_kind!r}")

    status = document["status"]
    if status not in _STATUSES:
        raise ValueError(f"status is invalid: {status!r}")

    _require_string(document, "environment_id")
    for field in ("source_sha", "event_sha"):
        value = document[field]
        if not isinstance(value, str) or _SHA_PATTERN.fullmatch(value) is None:
            raise ValueError(f"{field} must be a 40-character lowercase Git SHA")

    digest = document["tarball_sha256"]
    if digest is None:
        if status != "provisional":
            raise ValueError("tarball_sha256 may be null only for a provisional result")
    elif not isinstance(digest, str) or _DIGEST_PATTERN.fullmatch(digest) is None:
        raise ValueError("tarball_sha256 must be exactly 64 lowercase hexadecimal characters")

    if "tests" in document:
        _validate_count_object(
            document, "tests", ("pass", "fail", "error", "warning", "skip")
        )
    if "check" in document:
        _validate_count_object(
            document, "check", ("errors", "warnings", "notes", "log_path")
        )
    for field in ("predicate", "evidence", "classification", "message"):
        if field in document:
            _require_string(document, field)

    if status != "provisional":
        required_details = _DETAIL_FIELDS[result_kind]
        missing_details = [field for field in required_details if field not in document]
        if missing_details:
            raise ValueError(
                f"{result_kind} result is missing {missing_details[0]}"
            )

    return document


def _atomic_write(path: Path, document: Mapping[str, Any]) -> None:
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


def write_provisional(path: str | os.PathLike[str], identity: Mapping[str, Any]):
    """Validate and atomically write a provisional identity document."""

    document = dict(identity)
    document["schema_version"] = 1
    document["status"] = "provisional"
    validate_result(document)
    _atomic_write(Path(path), document)
    return document


def finalize(
    path: str | os.PathLike[str], status: str, details: Mapping[str, Any]
):
    """Validate and atomically replace a provisional result with its final form."""

    result_path = Path(path)
    current = json.loads(result_path.read_text(encoding="utf-8"))
    validate_result(current)
    if current["status"] != "provisional":
        raise ValueError("only a provisional result can be finalized")
    if status == "provisional":
        raise ValueError("final status cannot be provisional")
    if not isinstance(details, Mapping):
        raise ValueError("details must be an object")

    document = dict(current)
    document.update(details)
    document["status"] = status
    validate_result(document)
    _atomic_write(result_path, document)
    return document
