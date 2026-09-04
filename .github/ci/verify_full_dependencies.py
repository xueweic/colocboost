#!/usr/bin/env python3
"""Validate identity-bound full dependency preparation evidence."""

from __future__ import annotations

import argparse
import json
import os
import posixpath
import re
import stat
import sys
from collections.abc import Mapping
from pathlib import Path

from validate_manifest import load_manifest, validate_manifest


FIELDS = {
    "schema_version", "kind", "status", "environment_id", "purpose",
    "dependency_policy", "source_sha", "event_sha", "tarball_sha256",
    "r_executable", "r_resolved", "r_version", "r_platform", "library",
    "tarball", "hard_dependencies", "suggested_dependencies",
    "vignette_builders", "tooling", "installed_packages", "availability",
    "package_origins", "plan_only",
}
HARD = ["R", "Rfast", "matrixStats"]
SUGGESTS = ["testthat", "knitr", "rmarkdown", "ashr", "MASS", "susieR"]
VIGNETTE_BUILDERS = ["knitr"]
TOOLING = {
    "unit": ["jsonlite", "yaml"],
    "check": ["jsonlite", "yaml", "V8"],
}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _string_list(value, label, *, allow_empty=False):
    if (
        not isinstance(value, list)
        or (not allow_empty and not value)
        or any(not isinstance(item, str) or not item or "\x00" in item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be a unique string list")


def _single_line(value, label):
    if not isinstance(value, str) or not value or any(c in value for c in "\r\n\x00"):
        raise ValueError(f"{label} must be a nonempty single-line string")


def _absolute(value, label):
    _single_line(value, label)
    if not (Path(value).is_absolute() or (len(value) >= 3 and value[1:3] in {":\\", ":/"})):
        raise ValueError(f"{label} must be absolute")


def _normal_portable_path(value, label):
    _absolute(value, label)
    portable = value.replace("\\", "/")
    pieces = portable.split("/")
    if ".." in pieces:
        raise ValueError(f"{label} must not contain parent traversal")
    normal = posixpath.normpath(portable)
    if re_drive := (len(normal) >= 3 and normal[1:3] == ":/"):
        normal = normal.casefold()
    return normal, re_drive


def _same_portable_path(left, right, label):
    left_normal, left_windows = _normal_portable_path(left, label)
    right_normal, right_windows = _normal_portable_path(right, label)
    if left_windows != right_windows:
        return False
    return left_normal == right_normal


def _strict_portable_descendant(child, parent, label):
    child_normal, child_windows = _normal_portable_path(child, label)
    parent_normal, parent_windows = _normal_portable_path(parent, "library")
    if child_windows != parent_windows:
        return False
    return child_normal.startswith(parent_normal.rstrip("/") + "/")


def _approved_native_r_resolution(declared, resolved):
    if _same_portable_path(declared, resolved, "r_resolved"):
        return True
    declared_normal, declared_windows = _normal_portable_path(
        declared, "r_executable"
    )
    resolved_normal, resolved_windows = _normal_portable_path(
        resolved, "r_resolved"
    )
    if declared_windows or resolved_windows:
        return False
    return (
        re.fullmatch(r"/opt/R/[^/]+/bin/R", declared_normal) is not None
        and re.fullmatch(
            r"/opt/R/[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9._-]+)?/bin/R",
            resolved_normal,
        )
        is not None
    )


def validate_full_dependency_evidence(
    document,
    *,
    manifest,
    environment_id,
    purpose,
    source_sha,
    event_sha,
    tarball_sha256,
    r_executable,
):
    if not isinstance(document, Mapping) or set(document) != FIELDS:
        raise ValueError("full dependency evidence fields are not closed")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("schema_version must be integer 1")
    if document["kind"] != "full-dependency-proof":
        raise ValueError("full dependency evidence kind is invalid")
    if document["status"] != "pass" or document["plan_only"] is not False:
        raise ValueError("only an actual passing dependency run is acceptable")
    if purpose not in TOOLING or document["purpose"] != purpose:
        raise ValueError("dependency purpose does not match")
    if document["environment_id"] != environment_id:
        raise ValueError("dependency environment does not match")

    coverage = [row for row in manifest["coverage"] if row["id"] == environment_id]
    if len(coverage) != 1 or coverage[0]["state"] not in {"direct", "proxy"}:
        raise ValueError("dependency environment is not an executable manifest row")
    if purpose == "unit":
        lanes = [
            lane for lane in manifest["unit_lanes"]
            if lane["environment_id"] == environment_id
        ]
        if len(lanes) != 1 or lanes[0]["mode"] != "source" or lanes[0]["suite"] != "full":
            raise ValueError("dependency environment is not a full source unit lane")

    expected_identity = {
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
    }
    for field, expected in expected_identity.items():
        _single_line(expected, f"expected {field}")
        if document[field] != expected:
            raise ValueError(f"dependency evidence {field} does not match")
    if not _same_portable_path(document["r_executable"], r_executable, "r_executable"):
        raise ValueError("dependency evidence r_executable does not match")
    declared_r = coverage[0].get("system_r")
    if declared_r is not None and not _same_portable_path(
        declared_r, r_executable, "manifest system_r"
    ):
        raise ValueError("selected R does not match manifest system_r")

    for field in ("r_executable", "r_resolved", "library", "tarball"):
        _absolute(document[field], field)
    if (
        coverage[0].get("system_r") is not None
        and coverage[0].get("expected_os") == "linux"
        and not _approved_native_r_resolution(
            document["r_executable"], document["r_resolved"]
        )
    ):
        raise ValueError("resolved R does not match native manifest system_r")
    for field in ("r_version", "r_platform"):
        _single_line(document[field], field)
    if document["dependency_policy"] != "all":
        raise ValueError("dependency policy must be all")
    if type(document["plan_only"]) is not bool:
        raise ValueError("plan_only must be boolean")

    for field in (
        "hard_dependencies", "suggested_dependencies", "vignette_builders",
        "tooling", "installed_packages",
    ):
        _string_list(document[field], field, allow_empty=False)
    if document["hard_dependencies"] != HARD:
        raise ValueError("hard dependency set does not match DESCRIPTION")
    if document["suggested_dependencies"] != SUGGESTS:
        raise ValueError("Suggests set does not match DESCRIPTION")
    if document["vignette_builders"] != VIGNETTE_BUILDERS:
        raise ValueError("VignetteBuilder set does not match DESCRIPTION")
    if document["tooling"] != TOOLING[purpose]:
        raise ValueError("dependency tooling does not match lane purpose")

    availability = document["availability"]
    expected_available = set(HARD[1:] + SUGGESTS + TOOLING[purpose])
    if (
        not isinstance(availability, Mapping)
        or set(availability) != expected_available
        or any(type(value) is not bool for value in availability.values())
        or not all(availability.values())
    ):
        raise ValueError("effective dependency availability is incomplete")
    installed = set(document["installed_packages"])
    required_installed = set(
        HARD[1:] + [name for name in SUGGESTS if name != "MASS"] + TOOLING[purpose]
    )
    if not required_installed <= installed:
        raise ValueError("fresh lane library lacks required installed packages")
    origins = document["package_origins"]
    if (
        not isinstance(origins, Mapping)
        or set(origins) != expected_available
        or any(not isinstance(value, str) or not value for value in origins.values())
    ):
        raise ValueError("dependency package origins are incomplete")
    for name in required_installed:
        if name not in origins or not _strict_portable_descendant(
            origins[name], document["library"], f"{name} package origin"
        ):
            raise ValueError(f"{name} was not loaded from the fresh lane library")
    return document


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--purpose", required=True, choices=sorted(TOOLING))
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--r-executable", required=True)
    args = parser.parse_args(argv)
    try:
        path = Path(args.evidence)
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
            raise ValueError("evidence must be a regular non-symlink file")
        manifest = load_manifest(args.manifest)
        validate_manifest(manifest)
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
        validate_full_dependency_evidence(
            document,
            manifest=manifest,
            environment_id=args.environment_id,
            purpose=args.purpose,
            source_sha=args.source_sha,
            event_sha=args.event_sha,
            tarball_sha256=args.tarball_sha256,
            r_executable=args.r_executable,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        print(f"full dependency evidence error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
