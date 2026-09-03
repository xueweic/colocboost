#!/usr/bin/env python3
"""Validate actual R-hub dependency preparation evidence, never a plan."""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from collections.abc import Mapping
from pathlib import Path


FIELDS = {
    "schema_version",
    "kind",
    "status",
    "environment_id",
    "dependency_policy",
    "ci_tooling",
    "hard_dependencies",
    "suggested_dependencies",
    "recognized_testing_frameworks",
    "vignette_builders",
    "selected_soft_dependencies",
    "selected_refs",
    "excluded_suggests",
    "install_package",
    "plan_only",
    "library",
    "tarball",
    "installed_packages",
    "availability",
    "target_installed",
    "installed_tests",
}
HARD = ["R", "Rfast", "matrixStats"]
SUGGESTS = ["testthat", "knitr", "rmarkdown", "ashr", "MASS", "susieR"]
TESTING = ["testthat", "RUnit", "tinytest"]


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
        or (not value and not allow_empty)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be a unique string list")


def validate_dependency_evidence(document, environment_id: str):
    if not isinstance(document, Mapping) or set(document) != FIELDS:
        raise ValueError("dependency evidence fields are not closed")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("schema_version must be integer 1")
    if document["kind"] != "rhub-dependency-proof":
        raise ValueError("dependency evidence kind is invalid")
    if document["status"] != "pass" or document["plan_only"] is not False:
        raise ValueError("only an actual passing dependency run is acceptable")
    if environment_id not in {"mkl", "nosuggests"} or document["environment_id"] != environment_id:
        raise ValueError("dependency evidence environment does not match")
    for field in (
        "ci_tooling",
        "hard_dependencies",
        "suggested_dependencies",
        "recognized_testing_frameworks",
        "vignette_builders",
        "selected_soft_dependencies",
        "selected_refs",
        "excluded_suggests",
        "installed_packages",
    ):
        _string_list(
            document[field],
            field,
            allow_empty=field in {"ci_tooling", "excluded_suggests"},
        )
    if document["hard_dependencies"] != HARD:
        raise ValueError("hard dependency set does not match DESCRIPTION")
    if document["suggested_dependencies"] != SUGGESTS:
        raise ValueError("Suggests set does not match DESCRIPTION")
    if document["recognized_testing_frameworks"] != TESTING:
        raise ValueError("testing framework policy does not match R-hub")
    if document["vignette_builders"] != ["knitr"]:
        raise ValueError("VignetteBuilder set does not match DESCRIPTION")
    if not all(
        isinstance(document[field], str)
        and document[field].startswith("/")
        and not any(character in document[field] for character in "\r\n\x00")
        for field in ("library", "tarball")
    ):
        raise ValueError("dependency evidence paths must be absolute and single-line")
    for field in ("install_package", "plan_only", "target_installed", "installed_tests"):
        if type(document[field]) is not bool:
            raise ValueError(f"{field} must be boolean")
    availability = document["availability"]
    if not isinstance(availability, Mapping) or any(
        not isinstance(key, str) or type(value) is not bool
        for key, value in availability.items()
    ):
        raise ValueError("availability must be a boolean mapping")

    installed = set(document["installed_packages"])
    if environment_id == "nosuggests":
        selected = ["testthat", "knitr"]
        forbidden = ["devtools", "ashr", "susieR"]
        if document["dependency_policy"] != "hard-plus-testing-and-vignette-builder":
            raise ValueError("noSuggests dependency policy is invalid")
        if document["ci_tooling"]:
            raise ValueError("noSuggests must not add CI tooling packages")
        if document["selected_soft_dependencies"] != selected:
            raise ValueError("noSuggests selected soft dependencies are invalid")
        if document["excluded_suggests"] != ["rmarkdown", "ashr", "MASS", "susieR"]:
            raise ValueError("noSuggests excluded Suggests are invalid")
        if len(document["selected_refs"]) != len(selected):
            raise ValueError("noSuggests selected refs are incomplete")
        if not {"Rfast", "matrixStats", "testthat", "knitr", "colocboost"} <= installed:
            raise ValueError("noSuggests installed package set is incomplete")
        if installed & set(document["excluded_suggests"]):
            raise ValueError("noSuggests installed an excluded Suggests package")
        if installed & set(forbidden):
            raise ValueError("noSuggests installed a forbidden package")
        expected_availability = {"Rfast", "matrixStats", *selected, *forbidden}
        if set(availability) != expected_availability:
            raise ValueError("noSuggests availability proof is incomplete")
        if not all(availability[name] for name in ["Rfast", "matrixStats", *selected]):
            raise ValueError("a hard or allowed noSuggests package is unavailable")
        if any(availability[name] for name in forbidden):
            raise ValueError("a forbidden noSuggests package is available")
        if not (
            document["install_package"]
            and document["target_installed"]
            and document["installed_tests"]
        ):
            raise ValueError("noSuggests target or installed tests are unproven")
    else:
        if document["dependency_policy"] != "all":
            raise ValueError("MKL dependency policy is invalid")
        if document["ci_tooling"] != ["devtools", "jsonlite", "yaml"]:
            raise ValueError("MKL CI tooling is incomplete")
        if document["selected_soft_dependencies"] != SUGGESTS:
            raise ValueError("MKL must select all Suggests")
        if document["excluded_suggests"]:
            raise ValueError("MKL must not exclude Suggests")
        expected_available = {"Rfast", "matrixStats", *SUGGESTS, "devtools", "jsonlite", "yaml"}
        if set(availability) != expected_available or not all(availability.values()):
            raise ValueError("MKL effective dependency availability is incomplete")
        if not {"Rfast", "matrixStats"} <= installed:
            raise ValueError("MKL lane library lacks package hard dependencies")
        if document["install_package"] or document["target_installed"] or document["installed_tests"]:
            raise ValueError("MKL source-mode prep must not install the target")
    return document


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--environment-id", required=True, choices=("mkl", "nosuggests"))
    arguments = parser.parse_args(argv)
    try:
        path = Path(arguments.evidence)
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
            raise ValueError("evidence must be a regular non-symlink file")
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
        validate_dependency_evidence(document, arguments.environment_id)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        print(f"dependency evidence error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
