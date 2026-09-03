import copy
import re
import sys
from pathlib import Path

import pytest
import yaml


CI_DIR = Path(__file__).resolve().parents[1]
MANIFEST_PATH = CI_DIR / "check-matrix.yml"
POLICY_PATH = CI_DIR / "check-policy.yml"
sys.path.insert(0, str(CI_DIR))

from validate_manifest import (  # noqa: E402
    expected_result_ids,
    load_manifest,
    validate_manifest,
)


PRIMARY_NAMES = {
    "r-devel-linux-x86_64-debian-clang",
    "r-devel-linux-x86_64-debian-gcc",
    "r-devel-linux-x86_64-fedora-clang",
    "r-devel-linux-x86_64-fedora-gcc",
    "r-devel-windows-x86_64",
    "r-patched-linux-x86_64",
    "r-release-linux-x86_64",
    "r-release-macos-arm64",
    "r-release-macos-x86_64",
    "r-release-windows-x86_64",
    "r-oldrel-macos-arm64",
    "r-oldrel-macos-x86_64",
    "r-oldrel-windows-x86_64",
}
ADDITIONAL_NAMES = {
    "ATLAS",
    "BLIS",
    "BLAS",
    "C23",
    "Intel",
    "LTO",
    "M1mac",
    "MKL",
    "OpenBLAS",
    "Strict",
    "clang-ASAN",
    "clang-UBSAN",
    "donttest",
    "gcc-ASAN",
    "gcc-UBSAN",
    "gcc",
    "gcc15",
    "noLD",
    "noOMP",
    "noRemap",
    "noSuggests",
    "valgrind",
    "0len",
    "rchk",
    "rcnst",
    "rlibro",
    "musl",
    "linux-arm64",
    "vnu",
}
@pytest.fixture
def manifest():
    return load_manifest(MANIFEST_PATH)


def row_for(data, cran_name):
    return next(row for row in data["coverage"] if row["cran_name"] == cran_name)


def test_committed_manifest_has_exact_inventory_and_totals(manifest):
    assert validate_manifest(manifest) == manifest
    primary = [row for row in manifest["coverage"] if row["group"] == "primary"]
    additional = [
        row for row in manifest["coverage"] if row["group"] == "additional"
    ]

    assert {row["cran_name"] for row in primary} == PRIMARY_NAMES
    assert {row["cran_name"] for row in additional} == ADDITIONAL_NAMES
    assert {state: sum(row["state"] == state for row in primary) for state in {"direct", "proxy"}} == {
        "direct": 1,
        "proxy": 12,
    }
    assert {
        state: sum(row["state"] == state for row in additional)
        for state in {"direct", "proxy", "not-applicable", "uncovered"}
    } == {"direct": 13, "proxy": 6, "not-applicable": 10, "uncovered": 0}


@pytest.mark.parametrize("cran_name", sorted(PRIMARY_NAMES | ADDITIONAL_NAMES))
def test_removing_any_exact_inventory_name_is_rejected(manifest, cran_name):
    mutated = copy.deepcopy(manifest)
    mutated["coverage"] = [
        row for row in mutated["coverage"] if row["cran_name"] != cran_name
    ]

    with pytest.raises(ValueError, match="inventory"):
        validate_manifest(mutated)


def test_result_ids_are_complete_unique_and_artifact_safe(manifest):
    result_ids = expected_result_ids(manifest)

    assert result_ids == [row["id"] for row in manifest["coverage"]]
    assert len(result_ids) == len(set(result_ids)) == 42
    assert all(re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", value) for value in result_ids)


@pytest.mark.parametrize("bad_id", ["MKL", "mkl_result", "mkl result", "-mkl", "mkl-"])
def test_unsafe_artifact_ids_are_rejected(manifest, bad_id):
    mutated = copy.deepcopy(manifest)
    mutated["coverage"][0]["id"] = bad_id

    with pytest.raises(ValueError, match="artifact-safe"):
        validate_manifest(mutated)


def test_duplicate_ids_are_rejected(manifest):
    mutated = copy.deepcopy(manifest)
    mutated["coverage"][1]["id"] = mutated["coverage"][0]["id"]

    with pytest.raises(ValueError, match="duplicate id"):
        validate_manifest(mutated)


@pytest.mark.parametrize("field", ["proof", "limitation"])
def test_proxy_rows_require_proof_and_visible_limitation(manifest, field):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "BLIS").pop(field)

    with pytest.raises(ValueError, match=field):
        validate_manifest(mutated)


@pytest.mark.parametrize("mutation", ["missing", "no-command", "working-tree"])
def test_not_applicable_rows_require_executable_built_tarball_predicates(
    manifest, mutation
):
    mutated = copy.deepcopy(manifest)
    row = row_for(mutated, "C23")
    if mutation == "missing":
        row.pop("predicate")
    elif mutation == "no-command":
        row["predicate"].pop("command")
    else:
        row["predicate"]["input"] = "working-tree"

    with pytest.raises(ValueError, match="predicate"):
        validate_manifest(mutated)


def test_direct_rows_require_identity_assertions(manifest):
    mutated = copy.deepcopy(manifest)
    row = row_for(mutated, "MKL")
    row["proof"].remove("r-executable")

    with pytest.raises(ValueError, match="identity assertion"):
        validate_manifest(mutated)


def test_primary_rows_require_complete_identity_and_tarball_proof(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "r-devel-linux-x86_64-fedora-gcc")["proof"].remove("locale")

    with pytest.raises(ValueError, match="primary proof"):
        validate_manifest(mutated)


def test_uncovered_rows_are_rejected(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "BLIS")["state"] = "uncovered"

    with pytest.raises(ValueError, match="uncovered"):
        validate_manifest(mutated)


@pytest.mark.parametrize(
    ("cran_name", "missing_proof"),
    [("BLIS", "loaded-blis"), ("BLIS", "single-thread"), ("noOMP", "no-loaded-openmp-runtime")],
)
def test_acceptance_spikes_require_hard_runtime_proofs(
    manifest, cran_name, missing_proof
):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, cran_name)["proof"].remove(missing_proof)

    with pytest.raises(ValueError, match="acceptance spike"):
        validate_manifest(mutated)


@pytest.mark.parametrize(
    ("driver", "endpoint"),
    [("native-wrapper", "runner"), ("r-binary", "image")],
)
def test_driver_endpoint_types_cannot_be_substituted(manifest, driver, endpoint):
    mutated = copy.deepcopy(manifest)
    row = next(row for row in mutated["coverage"] if row["driver"] == driver)
    row[endpoint] = row.pop("image" if endpoint == "runner" else "runner")

    with pytest.raises(ValueError, match=driver):
        validate_manifest(mutated)


def test_active_rhub_images_must_be_immutable(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "MKL")["image"] = "ghcr.io/r-hub/containers/mkl:latest"

    with pytest.raises(ValueError, match="immutable"):
        validate_manifest(mutated)


def test_deprecated_rhub_images_cannot_be_declared(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "MKL")["image"] = (
        "ghcr.io/r-hub/containers/intel@sha256:" + "a" * 64
    )

    with pytest.raises(ValueError, match="deprecated"):
        validate_manifest(mutated)


def test_approved_rhub_image_mapping_cannot_be_substituted(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "MKL")["image"] = row_for(mutated, "ATLAS")["image"]

    with pytest.raises(ValueError, match="approved R-hub image"):
        validate_manifest(mutated)


def test_rhub_images_require_container_system_r_proof(manifest):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, "MKL")["proof"].remove("container-system-r")

    with pytest.raises(ValueError, match="container-system-r"):
        validate_manifest(mutated)


@pytest.mark.parametrize(
    ("cran_name", "missing_proof"),
    [
        ("MKL", "mkl-verbose"),
        ("noSuggests", "depends-only-policy"),
        ("valgrind", "use-valgrind"),
        ("vnu", "vnu-special-dispatch"),
    ],
)
def test_special_rhub_images_retain_native_setup_and_postprocessing(
    manifest, cran_name, missing_proof
):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, cran_name)["proof"].remove(missing_proof)

    with pytest.raises(ValueError, match="special image"):
        validate_manifest(mutated)


def test_policy_contains_only_the_five_exact_installed_skip_waivers():
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    waivers = policy["unit_tests"]["allowed_skips"]

    assert len(waivers) == 5
    assert all(set(waiver) == {
        "id",
        "context",
        "file",
        "test_title",
        "reason",
        "expected_count",
        "rationale",
        "expires",
    } for waiver in waivers)
    assert {waiver["context"] for waiver in waivers} == {"r-cmd-check-installed"}
    assert {waiver["expected_count"] for waiver in waivers} == {1}
    assert all(re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", waiver["id"]) for waiver in waivers)
    assert all("regex" not in key.lower() for waiver in waivers for key in waiver)
    assert {(waiver["file"], waiver["test_title"], waiver["reason"]) for waiver in waivers} == {
        (
            "test_model.R",
            "colocboost_init_data correctly initializes data",
            "colocboost_init_data not directly accessible",
        ),
        (
            "test_model.R",
            "colocboost correctly maps focal outcome to keep_variables with dict_keep_variables",
            "colocboost_init_data not directly accessible for integration test",
        ),
        (
            "test_model.R",
            "colocboost_assemble processes model results",
            "colocboost_assemble not directly accessible",
        ),
        (
            "test_model.R",
            "colocboost_workhorse performs boosting iterations",
            "colocboost_workhorse not directly accessible",
        ),
        (
            "test_utils.R",
            "colocboost_init_data handles complex dictionary mappings",
            "colocboost_init_data not directly accessible",
        ),
    }


def test_policy_fails_all_notes_and_uses_strict_mkl_library_patterns():
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    mkl = policy["numerical_backends"]["mkl"]

    assert policy["r_cmd_check"]["allowed_notes"] == []
    assert set(mkl["required_library_patterns"]) == {
        r"(?i)(?:^|/)libmkl_intel_lp64(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libmkl_core(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libmkl_sequential(?:\.so(?:\.\d+)*)?(?:$|\s)",
    }
    assert set(mkl["forbidden_library_patterns"]) == {
        r"(?i)(?:^|/)libmkl_(?:intel|gnu|tbb)_thread(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libopenblas(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)lib(?:s?atlas)(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libblis(?:\.so(?:\.\d+)*)?(?:$|\s)",
    }
