import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "verify_full_dependencies.py"
MANIFEST = CI_DIR / "check-matrix.yml"
SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
TARBALL_SHA256 = "c" * 64


def test_dependency_prep_binds_active_r_by_r_home_not_launcher_path():
    source = (CI_DIR / "prepare-full-dependencies.R").read_text(encoding="utf-8")
    assert "selected_r_home <- system2(" in source
    assert 'selected_r, "RHOME"' in source
    assert "active_r_home <- absolute(R.home()" in source
    assert 'file.path(R.home("bin")' not in source


def test_full_dependency_prep_preserves_single_vignette_builder_as_json_array():
    source = (CI_DIR / "prepare-full-dependencies.R").read_text(encoding="utf-8")
    assert "vignette_builders = I(builders)" in source


def valid_document(environment_id="r-devel-linux-x86-64-debian-clang", purpose="unit"):
    selected_r = {
        "r-release-linux-x86-64": "/opt/R/release/bin/R",
        "r-release-macos-arm64": "/Library/Frameworks/R.framework/Resources/bin/R",
        "r-release-windows-x86-64": "C:/R/bin/R.exe",
    }.get(environment_id, "/opt/R/devel/bin/R")
    resolved_r = (
        "/Library/Frameworks/R.framework/Versions/4.6-arm64/Resources/bin/R"
        if environment_id == "r-release-macos-arm64"
        else selected_r
    )
    return {
        "schema_version": 1,
        "kind": "full-dependency-proof",
        "status": "pass",
        "environment_id": environment_id,
        "purpose": purpose,
        "dependency_policy": "all",
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": TARBALL_SHA256,
        "r_executable": selected_r,
        "r_resolved": resolved_r,
        "r_version": "R Under development (unstable) (2026-09-01 r99999)",
        "r_platform": "x86_64-pc-linux-gnu",
        "library": f"/tmp/{environment_id}-{purpose}",
        "tarball": "/tmp/source/colocboost_1.0.9.tar.gz",
        "hard_dependencies": ["R", "Rfast", "matrixStats"],
        "suggested_dependencies": [
            "testthat", "knitr", "rmarkdown", "ashr", "MASS", "susieR"
        ],
        "vignette_builders": ["knitr"],
        "tooling": ["jsonlite", "yaml"],
        "installed_packages": [
            "Rfast", "matrixStats", "testthat", "knitr", "rmarkdown", "ashr",
            "susieR", "jsonlite", "yaml"
        ],
        "availability": {
            name: True
            for name in [
                "Rfast", "matrixStats", "testthat", "knitr", "rmarkdown",
                "ashr", "MASS", "susieR", "jsonlite", "yaml"
            ]
        },
        "package_origins": {
            name: (
                "/Library/Frameworks/R.framework/Resources/library/MASS"
                if name == "MASS"
                else f"/tmp/{environment_id}-{purpose}/{name}"
            )
            for name in [
                "Rfast", "matrixStats", "testthat", "knitr", "rmarkdown",
                "ashr", "MASS", "susieR", "jsonlite", "yaml"
            ]
        },
        "plan_only": False,
    }


def invoke(tmp_path, document, *, environment_id=None, purpose=None, r_executable=None):
    evidence = tmp_path / "evidence.json"
    evidence.write_text(json.dumps(document) + "\n", encoding="utf-8")
    return subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(SCRIPT),
            f"--manifest={MANIFEST}",
            f"--environment-id={environment_id or document['environment_id']}",
            f"--purpose={purpose or document['purpose']}",
            f"--evidence={evidence}",
            f"--source-sha={SOURCE_SHA}",
            f"--event-sha={EVENT_SHA}",
            f"--tarball-sha256={TARBALL_SHA256}",
            f"--r-executable={r_executable or document['r_executable']}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    "environment_id",
    [
        "r-devel-linux-x86-64-debian-clang",
        "r-release-linux-x86-64",
        "r-release-macos-arm64",
        "r-release-windows-x86-64",
        "atlas",
    ],
)
@pytest.mark.parametrize("purpose", ["unit", "check"])
def test_accepts_closed_identity_bound_full_dependency_evidence(
    tmp_path, environment_id, purpose
):
    document = valid_document(environment_id, purpose)
    completed = invoke(tmp_path, document)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("status", "plan"),
        ("plan_only", True),
        ("environment_id", "atlas"),
        ("purpose", "check"),
        ("source_sha", "d" * 40),
        ("event_sha", "d" * 40),
        ("tarball_sha256", "d" * 64),
        ("r_executable", "/tmp/R"),
        ("r_resolved", "/tmp/R"),
        ("dependency_policy", "hard"),
    ],
)
def test_rejects_stale_cross_lane_or_nonactual_evidence(tmp_path, field, value):
    document = valid_document()
    document[field] = value
    completed = invoke(
        tmp_path,
        document,
        environment_id="r-devel-linux-x86-64-debian-clang",
        purpose="unit",
        r_executable="/opt/R/devel/bin/R",
    )
    assert completed.returncode != 0


@pytest.mark.parametrize(
    "mutation",
    ["missing-available", "false-available", "missing-installed", "wrong-origin", "extra-field", "bad-type"],
)
def test_rejects_incomplete_or_open_dependency_proof(tmp_path, mutation):
    document = valid_document()
    if mutation == "missing-available":
        document["availability"].pop("MASS")
    elif mutation == "false-available":
        document["availability"]["susieR"] = False
    elif mutation == "missing-installed":
        document["installed_packages"].remove("Rfast")
    elif mutation == "wrong-origin":
        document["package_origins"]["ashr"] = "/ambient/library/ashr"
    elif mutation == "extra-field":
        document["untrusted"] = True
    else:
        document["plan_only"] = 0
    completed = invoke(tmp_path, document)
    assert completed.returncode != 0


def test_rejects_duplicate_json_keys(tmp_path):
    document = valid_document()
    evidence = tmp_path / "evidence.json"
    text = json.dumps(document)
    evidence.write_text(text[:-1] + ', "status": "pass"}\n', encoding="utf-8")
    completed = subprocess.run(
        [
            sys.executable, "-B", os.fspath(SCRIPT), f"--manifest={MANIFEST}",
            f"--environment-id={document['environment_id']}", "--purpose=unit",
            f"--evidence={evidence}", f"--source-sha={SOURCE_SHA}",
            f"--event-sha={EVENT_SHA}", f"--tarball-sha256={TARBALL_SHA256}",
            f"--r-executable={document['r_executable']}",
        ], check=False, capture_output=True, text=True,
    )
    assert completed.returncode != 0


def test_accepts_canonical_windows_paths_with_native_separator_evidence(tmp_path):
    document = valid_document("r-release-windows-x86-64", "check")
    document["r_executable"] = r"C:\R\bin\R.exe"
    document["r_resolved"] = r"C:\R\bin\R.exe"
    document["library"] = r"C:\runner\fresh-check-library"
    document["tarball"] = r"C:\runner\source\colocboost_1.0.9.tar.gz"
    for name in document["package_origins"]:
        if name != "MASS":
            document["package_origins"][name] = (
                rf"C:\runner\fresh-check-library\{name}"
            )
    completed = invoke(
        tmp_path,
        document,
        environment_id="r-release-windows-x86-64",
        purpose="check",
        r_executable="C:/R/bin/R.exe",
    )
    assert completed.returncode == 0, completed.stderr


def test_accepts_versioned_resolution_of_native_opt_r_alias(tmp_path):
    document = valid_document("atlas", "unit")
    assert document["r_executable"] == "/opt/R/devel/bin/R"
    document["r_resolved"] = "/opt/R/4.6.1/bin/R"

    completed = invoke(tmp_path, document)

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "resolved_r",
    [
        "/tmp/R",
        "/opt/R/not-a-version/bin/R",
        "/opt/R/4.6.1/bin/Rscript",
        "/opt/R/4.6.1/../other/bin/R",
    ],
)
def test_rejects_unbound_resolution_of_native_opt_r_alias(tmp_path, resolved_r):
    document = valid_document("atlas", "unit")
    document["r_resolved"] = resolved_r

    completed = invoke(tmp_path, document)

    assert completed.returncode != 0


@pytest.mark.parametrize(
    "origin",
    [
        "/tmp/r-devel-linux-x86-64-debian-clang-unit",
        "/tmp/r-devel-linux-x86-64-debian-clang-unit/../ambient/ashr",
    ],
)
def test_rejects_non_descendant_or_traversing_package_origins(tmp_path, origin):
    document = valid_document()
    document["package_origins"]["ashr"] = origin
    completed = invoke(tmp_path, document)
    assert completed.returncode != 0
    assert (
        "fresh lane library" in completed.stderr
        or "parent traversal" in completed.stderr
    )
