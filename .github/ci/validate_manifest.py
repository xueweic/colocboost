#!/usr/bin/env python3
"""Validate the committed CRAN preflight coverage manifest."""

from __future__ import annotations

import re
import sys
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml


PRIMARY_STATES = {
    "r-devel-linux-x86_64-debian-clang": "proxy",
    "r-devel-linux-x86_64-debian-gcc": "proxy",
    "r-devel-linux-x86_64-fedora-clang": "proxy",
    "r-devel-linux-x86_64-fedora-gcc": "direct",
    "r-devel-windows-x86_64": "proxy",
    "r-patched-linux-x86_64": "proxy",
    "r-release-linux-x86_64": "proxy",
    "r-release-macos-arm64": "proxy",
    "r-release-macos-x86_64": "proxy",
    "r-release-windows-x86_64": "proxy",
    "r-oldrel-macos-arm64": "proxy",
    "r-oldrel-macos-x86_64": "proxy",
    "r-oldrel-windows-x86_64": "proxy",
}
ADDITIONAL_STATES = {
    "ATLAS": "direct",
    "BLIS": "proxy",
    "BLAS": "not-applicable",
    "C23": "not-applicable",
    "Intel": "not-applicable",
    "LTO": "not-applicable",
    "M1mac": "proxy",
    "MKL": "direct",
    "OpenBLAS": "proxy",
    "Strict": "not-applicable",
    "clang-ASAN": "direct",
    "clang-UBSAN": "direct",
    "donttest": "direct",
    "gcc-ASAN": "direct",
    "gcc-UBSAN": "direct",
    "gcc": "not-applicable",
    "gcc15": "not-applicable",
    "noLD": "direct",
    "noOMP": "proxy",
    "noRemap": "not-applicable",
    "noSuggests": "direct",
    "valgrind": "direct",
    "0len": "not-applicable",
    "rchk": "not-applicable",
    "rcnst": "proxy",
    "rlibro": "proxy",
    "musl": "direct",
    "linux-arm64": "direct",
    "vnu": "direct",
}
EXPECTED_STATES = {"primary": PRIMARY_STATES, "additional": ADDITIONAL_STATES}
ARTIFACT_ID = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
RHUB_IMAGE = re.compile(
    r"^ghcr\.io/r-hub/containers/(?P<name>[a-z0-9-]+)@sha256:[0-9a-f]{64}$"
)
DEPRECATED_RHUB_IMAGES = {"c23", "gcc15", "intel", "noremap", "ubuntu-gcc12"}
PRIMARY_PROOF = {
    "source-tarball-sha256",
    "r-executable",
    "r-version",
    "operating-system",
    "compiler",
    "architecture",
    "locale",
}
ACCEPTANCE_SPIKE_PROOFS = {
    "BLIS": {"loaded-blis", "blis-version", "single-thread"},
    "noOMP": {
        "no-openmp-compile-flags",
        "no-openmp-link-flags",
        "no-loaded-openmp-runtime",
        "dependency-openmp-audit",
    },
}
SPECIAL_IMAGE_PROOFS = {
    "MKL": {
        "container-system-r",
        "session-info",
        "blas-identity",
        "lapack-identity",
        "la-library",
        "la-version",
        "process-mappings",
        "loaded-mkl",
        "blas-operation",
        "mkl-verbose",
        "single-thread",
    },
    "noSuggests": {
        "container-system-r",
        "depends-only-policy",
        "allowed-test-frameworks",
        "allowed-vignette-builders",
        "nonexempt-suggests-absent",
        "nonzero-installed-test-count",
    },
    "valgrind": {
        "container-system-r",
        "opt-r-devel-valgrind",
        "use-valgrind",
        "valgrind-runtime",
        "suppression-and-error-scan",
    },
    "vnu": {
        "container-system-r",
        "vnu-special-dispatch",
        "nu-validator-executed",
        "zero-bad-entries",
        "validator-output",
    },
}
APPROVED_RHUB_IMAGE = {
    "r-devel-linux-x86_64-debian-clang": "clang22",
    "r-devel-linux-x86_64-debian-gcc": "ubuntu-gcc16",
    "r-devel-linux-x86_64-fedora-clang": "clang22",
    "r-devel-linux-x86_64-fedora-gcc": "gcc16",
    "r-patched-linux-x86_64": "ubuntu-next",
    "r-release-linux-x86_64": "ubuntu-release",
    "ATLAS": "atlas",
    "MKL": "mkl",
    "OpenBLAS": "gcc16",
    "clang-ASAN": "clang-asan",
    "clang-UBSAN": "clang-ubsan",
    "donttest": "donttest",
    "gcc-ASAN": "gcc-asan",
    "gcc-UBSAN": "gcc-asan",
    "noLD": "nold",
    "noSuggests": "nosuggests",
    "valgrind": "valgrind",
    "rcnst": "ubuntu-clang",
    "vnu": "vnu",
}
DRIVER_ENDPOINT = {"native-wrapper": "image", "r-binary": "runner"}


def _require_nonempty_string(value: Any, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")


def _require_string_list(value: Any, label: str) -> None:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be a non-empty list of unique strings")


def _valid_endpoint(value: Any) -> bool:
    if isinstance(value, str):
        return bool(value.strip())
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(item, str) and bool(item.strip()) for item in value)
        and len(value) == len(set(value))
    )


def _expected_id(group: str, cran_name: str) -> str:
    del group
    return cran_name.lower().replace("_", "-")


def _validate_predicate(row: Mapping[str, Any]) -> None:
    predicate = row.get("predicate")
    if not isinstance(predicate, Mapping):
        raise ValueError(f"{row['cran_name']} predicate must be an object")
    if set(predicate) != {"input", "command"}:
        raise ValueError(
            f"{row['cran_name']} predicate must contain only input and command"
        )
    if predicate["input"] != "built-source-tarball":
        raise ValueError(f"{row['cran_name']} predicate must inspect built-source-tarball")
    command = predicate["command"]
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(part, str) or not part for part in command)
        or "{tarball}" not in command
    ):
        raise ValueError(
            f"{row['cran_name']} predicate command must be an executable argv list using {{tarball}}"
        )


def _validate_row(row: Any, index: int) -> None:
    if not isinstance(row, Mapping):
        raise ValueError(f"coverage row {index} must be an object")
    required = {"id", "group", "cran_name", "state", "driver", "proof"}
    missing = required - set(row)
    if missing:
        raise ValueError(f"coverage row {index} is missing {sorted(missing)[0]}")

    for field in ("id", "group", "cran_name", "state", "driver"):
        _require_nonempty_string(row[field], f"coverage row {index}.{field}")
    if ARTIFACT_ID.fullmatch(row["id"]) is None:
        raise ValueError(f"{row['id']!r} is not a lowercase artifact-safe id")
    if row["group"] not in EXPECTED_STATES:
        raise ValueError(f"{row['id']} has invalid group {row['group']!r}")
    if row["state"] not in {"direct", "proxy", "not-applicable", "uncovered"}:
        raise ValueError(f"{row['id']} has invalid state {row['state']!r}")
    if row["state"] == "uncovered":
        raise ValueError(f"{row['id']} is uncovered")

    expected_id = _expected_id(row["group"], row["cran_name"])
    if row["id"] != expected_id:
        raise ValueError(
            f"{row['cran_name']} id must remain the stable artifact-safe id {expected_id}"
        )

    _require_string_list(row["proof"], f"{row['id']}.proof")
    proof = set(row["proof"])
    if "source-tarball-sha256" not in proof:
        raise ValueError(f"{row['id']}.proof must include source-tarball-sha256")

    driver = row["driver"]
    if driver not in DRIVER_ENDPOINT:
        raise ValueError(f"{row['id']} has invalid driver {driver!r}")
    endpoint = DRIVER_ENDPOINT[driver]
    other_endpoint = "runner" if endpoint == "image" else "image"
    if endpoint not in row or not _valid_endpoint(row[endpoint]):
        raise ValueError(f"{driver} row {row['id']} requires a valid {endpoint}")
    if other_endpoint in row:
        raise ValueError(f"{driver} row {row['id']} must not declare {other_endpoint}")

    if endpoint == "image" and isinstance(row["image"], str):
        image = row["image"]
        if image.startswith("ghcr.io/r-hub/containers/"):
            image_name = image.split("/", maxsplit=4)[-1].split("@", maxsplit=1)[0]
            if image_name in DEPRECATED_RHUB_IMAGES:
                raise ValueError(f"{row['id']} uses deprecated R-hub image {image_name}")
            match = RHUB_IMAGE.fullmatch(image)
            if match is None:
                raise ValueError(f"{row['id']} R-hub image must use an immutable digest")
            approved = APPROVED_RHUB_IMAGE.get(row["cran_name"])
            if approved != match.group("name"):
                raise ValueError(
                    f"{row['id']} must use approved R-hub image {approved}, "
                    f"not {match.group('name')}"
                )
            if "container-system-r" not in proof:
                raise ValueError(
                    f"{row['id']}.proof must include container-system-r"
                )

    if row["group"] == "primary":
        missing_proof = PRIMARY_PROOF - proof
        if missing_proof:
            raise ValueError(
                f"{row['id']} primary proof is missing {sorted(missing_proof)[0]}"
            )
    if row["state"] == "direct" and "r-executable" not in proof:
        raise ValueError(f"{row['id']} direct row lacks an identity assertion")
    if row["state"] == "proxy":
        if not proof:
            raise ValueError(f"{row['id']}.proof is required for proxy coverage")
        _require_nonempty_string(row.get("limitation"), f"{row['id']}.limitation")
    if row["state"] == "not-applicable":
        _validate_predicate(row)

    required_spike_proof = ACCEPTANCE_SPIKE_PROOFS.get(row["cran_name"])
    if required_spike_proof and not required_spike_proof <= proof:
        missing_spike = sorted(required_spike_proof - proof)[0]
        raise ValueError(
            f"{row['cran_name']} acceptance spike proof is missing {missing_spike}"
        )
    required_special_proof = SPECIAL_IMAGE_PROOFS.get(row["cran_name"])
    if required_special_proof and not required_special_proof <= proof:
        missing_special = sorted(required_special_proof - proof)[0]
        raise ValueError(
            f"{row['cran_name']} special image proof is missing {missing_special}"
        )


def load_manifest(path: str | Path) -> Mapping[str, Any]:
    """Load a YAML manifest without accepting an empty document."""

    manifest_path = Path(path)
    try:
        data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise ValueError(f"could not load manifest {manifest_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ValueError("manifest must be an object")
    return data


def validate_manifest(data: Mapping[str, Any]) -> Mapping[str, Any]:
    """Fail closed unless *data* declares the exact approved inventory."""

    if not isinstance(data, Mapping):
        raise ValueError("manifest must be an object")
    if set(data) != {"version", "coverage"}:
        raise ValueError("manifest must contain only version and coverage")
    if isinstance(data["version"], bool) or data["version"] != 1:
        raise ValueError("manifest version must be 1")
    rows = data["coverage"]
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        raise ValueError("coverage must be a list")

    ids = [row.get("id") for row in rows if isinstance(row, Mapping)]
    duplicates = sorted(
        name
        for name, count in Counter(ids).items()
        if isinstance(name, str) and count > 1
    )
    if duplicates:
        raise ValueError(f"duplicate id {duplicates[0]}")

    for index, row in enumerate(rows):
        _validate_row(row, index)

    names = [(row["group"], row["cran_name"]) for row in rows]
    duplicate_names = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicate_names:
        raise ValueError(f"duplicate inventory entry {duplicate_names[0]}")

    actual_by_group = {
        group: {row["cran_name"]: row["state"] for row in rows if row["group"] == group}
        for group in EXPECTED_STATES
    }
    for group, expected in EXPECTED_STATES.items():
        actual = actual_by_group[group]
        if actual != expected:
            missing = sorted(set(expected) - set(actual))
            extra = sorted(set(actual) - set(expected))
            wrong = sorted(
                name
                for name in set(actual) & set(expected)
                if actual[name] != expected[name]
            )
            detail = (
                f"missing={missing}, extra={extra}, wrong_state={wrong}"
            )
            raise ValueError(f"{group} inventory does not match approved mapping: {detail}")

    return data


def expected_result_ids(data: Mapping[str, Any]) -> list[str]:
    """Return all stable result IDs after validating the complete inventory."""

    validate_manifest(data)
    return [row["id"] for row in data["coverage"]]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {Path(argv[0]).name} MANIFEST", file=sys.stderr)
        return 2
    try:
        data = load_manifest(argv[1])
        validate_manifest(data)
    except ValueError as error:
        print(f"manifest invalid: {error}", file=sys.stderr)
        return 1

    rows = data["coverage"]
    primary = sum(row["group"] == "primary" for row in rows)
    additional = sum(row["group"] == "additional" for row in rows)
    uncovered = sum(row["state"] == "uncovered" for row in rows)
    print(f"manifest valid: primary={primary}, additional={additional}, uncovered={uncovered}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
