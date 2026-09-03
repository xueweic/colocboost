#!/usr/bin/env python3
"""Validate runtime evidence captured inside the selected R process."""

from __future__ import annotations

import argparse
import json
import math
import re
import stat
import sys
from collections.abc import Mapping
from pathlib import Path

from validate_manifest import load_manifest, validate_manifest


FIELDS = {
    "schema_version", "kind", "status", "environment_id", "runtime_profile",
    "source_sha", "event_sha", "tarball_sha256", "r_executable", "r_resolved",
    "r_version", "r_status", "r_platform", "os", "os_release",
    "distribution_id", "distribution_version", "architecture", "compiler",
    "cc", "cc_path", "cc_version", "cxx", "cxx_path", "cxx_version",
    "fc", "fc_path", "fc_version", "locale", "session_info", "ext_soft_version",
    "la_library", "la_version", "blas_libs", "matrix_dimension",
    "matrix_checksum", "maps_method", "loaded_libraries", "long_double",
    "environment",
}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _single_line(value, label, *, allow_empty=False):
    if (
        not isinstance(value, str)
        or (not allow_empty and not value)
        or any(character in value for character in "\r\n\x00")
    ):
        raise ValueError(f"{label} must be a single-line string")


def _classify_r(document):
    status = document["r_status"].lower()
    version = document["r_version"].lower()
    if "under development" in status or "under development" in version:
        return "devel"
    if "patched" in status or "patched" in version:
        return "patched"
    return "release"


def validate_runtime_evidence(
    document,
    *,
    manifest,
    environment_id,
    source_sha,
    event_sha,
    tarball_sha256,
):
    if not isinstance(document, Mapping) or set(document) != FIELDS:
        raise ValueError("runtime evidence fields are not closed")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("schema_version must be integer 1")
    if document["kind"] != "r-runtime-proof" or document["status"] != "pass":
        raise ValueError("runtime evidence is not a passing R runtime proof")
    rows = [row for row in manifest["coverage"] if row["id"] == environment_id]
    if len(rows) != 1:
        raise ValueError("runtime environment must select exactly one manifest row")
    row = rows[0]
    expected = {
        "environment_id": environment_id,
        "runtime_profile": row.get("runtime_profile"),
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "r_executable": row.get("system_r"),
        "r_resolved": row.get("system_r"),
        "os": row.get("expected_os"),
        "distribution_id": row.get("expected_distribution"),
        "distribution_version": row.get("expected_distribution_version"),
        "architecture": row.get("expected_architecture"),
    }
    if any(value is None for value in expected.values()):
        raise ValueError("manifest row lacks a closed runtime binding")
    for field, value in expected.items():
        if document[field] != value:
            raise ValueError(f"runtime {field} does not match manifest/caller")
    for field in (
        "r_executable", "r_resolved", "r_version", "r_platform", "os",
        "os_release", "distribution_id", "distribution_version", "architecture",
        "compiler", "cc", "cc_path", "cc_version", "cxx", "cxx_path",
        "cxx_version", "fc", "fc_path", "fc_version", "locale",
        "session_info", "la_library", "la_version", "blas_libs",
    ):
        _single_line(document[field], field)
    _single_line(document["r_status"], "r_status", allow_empty=True)
    expected_kind = row["expected_r_kind"]
    status_lower = document["r_status"].lower()
    if (
        _classify_r(document) != expected_kind
        or expected_kind == "devel" and "under development" not in status_lower
        or expected_kind == "patched" and "patched" not in status_lower
        or expected_kind == "release" and status_lower
    ):
        raise ValueError("runtime R release kind does not match manifest")

    profile = document["runtime_profile"]
    for field in ("cc_path", "cxx_path", "fc_path"):
        if not document[field].startswith("/"):
            raise ValueError(f"runtime {field} is not an absolute compiler path")
    if profile == "clang22":
        if not re.search(r"clang(?:\+\+)?-?22(?:\b|$)", document["cc"], re.I):
            raise ValueError("clang22 runtime does not prove clang 22")
        if not re.search(r"clang(?:\+\+)?-?22(?:\b|$)", document["cxx"], re.I):
            raise ValueError("clang22 runtime does not prove clang++ 22")
        if not re.search(r"flang(?:-new)?-?22(?:\b|$)", document["fc"], re.I):
            raise ValueError("clang22 runtime does not prove flang 22")
        for field, pattern in (
            ("cc_path", r"clang-?22$"),
            ("cxx_path", r"clang\+\+-?22$"),
            ("fc_path", r"flang(?:-new)?-?22$"),
        ):
            if re.search(pattern, document[field], re.I) is None:
                raise ValueError("clang22 runtime compiler executable/version is wrong")
        for field, family in (
            ("cc_version", "clang"), ("cxx_version", "clang"),
            ("fc_version", "flang"),
        ):
            value = document[field]
            if family not in value.lower() or re.search(r"\b22(?:\.|\b)", value) is None:
                raise ValueError("clang22 runtime compiler executable/version is wrong")
    elif profile == "ubuntu-gcc16":
        if not re.search(r"gcc-?16(?:\b|$)", document["cc"], re.I):
            raise ValueError("ubuntu-gcc16 runtime does not prove GCC 16")
        if re.search(r"/gcc-?16$", document["cc_path"], re.I) is None:
            raise ValueError("ubuntu-gcc16 runtime compiler executable is wrong")
        if "gcc" not in document["cc_version"].lower() or re.search(
            r"\b16(?:\.|\b)", document["cc_version"]
        ) is None:
            raise ValueError("ubuntu-gcc16 runtime compiler version is wrong")
    elif profile == "gcc16":
        if "gcc" not in document["cc"].lower() or "gcc" not in document["cc_version"].lower():
            raise ValueError("gcc16 runtime does not prove GCC")
        if re.search(r"\b16(?:\.|\b)", document["cc_version"]) is None:
            raise ValueError("gcc16 runtime compiler version is wrong")

    dimension = document["matrix_dimension"]
    checksum = document["matrix_checksum"]
    if type(dimension) is not int or dimension < 2:
        raise ValueError("matrix operation dimension is invalid")
    if type(checksum) not in {int, float} or isinstance(checksum, bool) or not math.isfinite(checksum) or checksum == 0:
        raise ValueError("matrix operation checksum is invalid")
    if document["maps_method"] != "proc-self-maps":
        raise ValueError("runtime libraries were not captured from this R process")
    libraries = document["loaded_libraries"]
    if (
        not isinstance(libraries, list) or not libraries
        or any(not isinstance(item, str) or not item.startswith("/") or "\x00" in item for item in libraries)
        or len(libraries) != len(set(libraries))
    ):
        raise ValueError("loaded library evidence is invalid")
    if type(document["long_double"]) is not bool:
        raise ValueError("long_double must be boolean")
    if not isinstance(document["ext_soft_version"], Mapping):
        raise ValueError("ext_soft_version must be an object")
    environment = document["environment"]
    if not isinstance(environment, Mapping) or any(
        not isinstance(key, str)
        or value is not None and (
            not isinstance(value, str) or any(character in value for character in "\r\n\x00")
        )
        for key, value in environment.items()
    ):
        raise ValueError("runtime environment evidence is invalid")

    lowered = "\n".join(libraries).lower()
    if profile == "atlas":
        if re.search(r"lib(?:mkl|openblas|blis)", lowered):
            raise ValueError("ATLAS runtime loaded a forbidden BLAS")
        if re.search(r"lib(?:s?atlas)(?:\.so(?:\.\d+)*)?(?:$|\s|/)", lowered) is None:
            raise ValueError("ATLAS runtime did not load libsatlas")
        if "atlas" not in document["la_library"].lower():
            raise ValueError("ATLAS La_library identity is missing")
    return document


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    args = parser.parse_args(argv)
    try:
        path = Path(args.evidence)
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
            raise ValueError("evidence must be a regular non-symlink file")
        manifest = load_manifest(args.manifest)
        validate_manifest(manifest)
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
        validate_runtime_evidence(
            document,
            manifest=manifest,
            environment_id=args.environment_id,
            source_sha=args.source_sha,
            event_sha=args.event_sha,
            tarball_sha256=args.tarball_sha256,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        print(f"runtime evidence error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
