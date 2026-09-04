import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "verify_dependency_evidence.py"


def test_rhub_dependency_prep_preserves_single_vignette_builder_as_json_array():
    source = (CI_DIR / "prepare-rhub-dependencies.R").read_text(encoding="utf-8")
    assert "vignette_builders = I(vignette_builders)" in source


def valid_document(environment_id="nosuggests"):
    suggestions = ["testthat", "knitr", "rmarkdown", "ashr", "MASS", "susieR"]
    if environment_id == "nosuggests":
        selected = ["testthat", "knitr"]
        excluded = ["rmarkdown", "ashr", "MASS", "susieR"]
        installed = ["Rfast", "matrixStats", "testthat", "knitr", "colocboost"]
        availability = {
            "Rfast": True,
            "matrixStats": True,
            "testthat": True,
            "knitr": True,
            "devtools": False,
            "ashr": False,
            "susieR": False,
        }
        return {
            "schema_version": 1,
            "kind": "rhub-dependency-proof",
            "status": "pass",
            "environment_id": environment_id,
            "dependency_policy": "hard-plus-testing-and-vignette-builder",
            "ci_tooling": [],
            "hard_dependencies": ["R", "Rfast", "matrixStats"],
            "suggested_dependencies": suggestions,
            "recognized_testing_frameworks": ["testthat", "RUnit", "tinytest"],
            "vignette_builders": ["knitr"],
            "selected_soft_dependencies": selected,
            "selected_refs": ["testthat@>= 3.0.0", "knitr"],
            "excluded_suggests": excluded,
            "install_package": True,
            "plan_only": False,
            "library": "/tmp/check/unit-library",
            "tarball": "/tmp/source/colocboost_1.0.9.tar.gz",
            "installed_packages": installed,
            "availability": availability,
            "target_installed": True,
            "installed_tests": True,
        }
    return {
        "schema_version": 1,
        "kind": "rhub-dependency-proof",
        "status": "pass",
        "environment_id": "mkl",
        "dependency_policy": "all",
        "ci_tooling": ["devtools", "jsonlite", "yaml"],
        "hard_dependencies": ["R", "Rfast", "matrixStats"],
        "suggested_dependencies": suggestions,
        "recognized_testing_frameworks": ["testthat", "RUnit", "tinytest"],
        "vignette_builders": ["knitr"],
        "selected_soft_dependencies": suggestions,
        "selected_refs": suggestions,
        "excluded_suggests": [],
        "install_package": False,
        "plan_only": False,
        "library": "/tmp/check/unit-library",
        "tarball": "/tmp/source/colocboost_1.0.9.tar.gz",
        "installed_packages": ["Rfast", "matrixStats", *suggestions],
        "availability": {
            name: True
            for name in [
                "Rfast", "matrixStats", *suggestions, "devtools", "jsonlite", "yaml"
            ]
        },
        "target_installed": False,
        "installed_tests": False,
    }


def invoke(tmp_path, document, environment_id="nosuggests"):
    evidence = tmp_path / "evidence.json"
    evidence.write_text(json.dumps(document) + "\n", encoding="utf-8")
    return subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(SCRIPT),
            f"--evidence={evidence}",
            f"--environment-id={environment_id}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("environment_id", ["mkl", "nosuggests"])
def test_accepts_closed_actual_dependency_proof(tmp_path, environment_id):
    completed = invoke(tmp_path, valid_document(environment_id), environment_id)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "plan",
        "wrong-policy",
        "missing-installed",
        "target-not-installed",
        "tests-missing",
        "forbidden-available",
        "excluded-installed",
        "allowed-unavailable",
        "extra-field",
    ],
)
def test_nosuggests_rejects_plans_and_unproven_runtime_state(tmp_path, mutation):
    document = valid_document()
    if mutation == "plan":
        document["plan_only"] = True
        document["status"] = "plan"
    elif mutation == "wrong-policy":
        document["dependency_policy"] = "all"
    elif mutation == "missing-installed":
        document["installed_packages"].remove("testthat")
    elif mutation == "target-not-installed":
        document["target_installed"] = False
    elif mutation == "tests-missing":
        document["installed_tests"] = False
    elif mutation == "forbidden-available":
        document["availability"]["ashr"] = True
    elif mutation == "excluded-installed":
        document["installed_packages"].append("rmarkdown")
    elif mutation == "allowed-unavailable":
        document["availability"]["knitr"] = False
    else:
        document["untrusted"] = True

    completed = invoke(tmp_path, document)
    assert completed.returncode != 0
