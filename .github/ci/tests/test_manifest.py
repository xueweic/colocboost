import copy
import re
import sys
from datetime import date
from pathlib import Path

import pytest
import yaml


CI_DIR = Path(__file__).resolve().parents[1]
MANIFEST_PATH = CI_DIR / "check-matrix.yml"
POLICY_PATH = CI_DIR / "check-policy.yml"
sys.path.insert(0, str(CI_DIR))

import validate_manifest as validator_module  # noqa: E402
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

UNIT_COVERAGE_IDS = [
    "r-devel-linux-x86-64-debian-clang",
    "r-devel-linux-x86-64-debian-gcc",
    "r-devel-linux-x86-64-fedora-clang",
    "r-devel-linux-x86-64-fedora-gcc",
    "r-devel-windows-x86-64",
    "r-patched-linux-x86-64",
    "r-release-linux-x86-64",
    "r-release-macos-arm64",
    "r-release-macos-x86-64",
    "r-release-windows-x86-64",
    "r-oldrel-macos-arm64",
    "r-oldrel-macos-x86-64",
    "r-oldrel-windows-x86-64",
    "atlas",
    "blis",
    "mkl",
    "openblas",
    "nosuggests",
]

# Independent golden contract: source contexts are deliberately distinct from
# the sole installed-waiver context. Targeted MKL runs are diagnostic sidecars,
# not additional aggregate results.
GOLDEN_UNIT_LANES = [
    {
        "environment_id": coverage_id,
        "coverage_id": coverage_id,
        "mode": "installed" if coverage_id == "nosuggests" else "source",
        "suite": "full",
        "runner_context": (
            "r-cmd-check-installed" if coverage_id == "nosuggests" else coverage_id
        ),
        "required_sidecar_tests": (
            ["test_utils.R", "test_Xref.R"] if coverage_id == "mkl" else []
        ),
    }
    for coverage_id in UNIT_COVERAGE_IDS
]

# Independent golden contract: no value here is derived from the YAML under test.
GOLDEN_ROWS = {
    "r-devel-linux-x86_64-debian-clang": ("primary", "r-devel-linux-x86-64-debian-clang", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e"),
    "r-devel-linux-x86_64-debian-gcc": ("primary", "r-devel-linux-x86-64-debian-gcc", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-gcc16@sha256:2e9576e51ad17a706887b7e06fc4057388765226ec795f3c7af0e7348fb8fbf1"),
    "r-devel-linux-x86_64-fedora-clang": ("primary", "r-devel-linux-x86-64-fedora-clang", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e"),
    "r-devel-linux-x86_64-fedora-gcc": ("primary", "r-devel-linux-x86-64-fedora-gcc", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2"),
    "r-devel-windows-x86_64": ("primary", "r-devel-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "r-patched-linux-x86_64": ("primary", "r-patched-linux-x86-64", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-next@sha256:1c29eed93b0aa05fe147e476e6b3510a09fb0476460eeefbbe99a4a1eb2c3bd9"),
    "r-release-linux-x86_64": ("primary", "r-release-linux-x86-64", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-release@sha256:714722b7ecb4307fbf88a707f83fb045181482c0c8ce6086a0e26860416ca9c2"),
    "r-release-macos-arm64": ("primary", "r-release-macos-arm64", "proxy", "r-binary", "runner", "macos-15"),
    "r-release-macos-x86_64": ("primary", "r-release-macos-x86-64", "proxy", "r-binary", "runner", "macos-15-intel"),
    "r-release-windows-x86_64": ("primary", "r-release-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "r-oldrel-macos-arm64": ("primary", "r-oldrel-macos-arm64", "proxy", "r-binary", "runner", "macos-15"),
    "r-oldrel-macos-x86_64": ("primary", "r-oldrel-macos-x86-64", "proxy", "r-binary", "runner", "macos-15-intel"),
    "r-oldrel-windows-x86_64": ("primary", "r-oldrel-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "ATLAS": ("additional", "atlas", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/atlas@sha256:7597f7d0b6b2f009ae7bb425391523d8f4388223238db50d7dfb1572add63a88"),
    "BLIS": ("additional", "blis", "proxy", "native-wrapper", "image", "local/colocboost-blis-proxy"),
    "BLAS": ("additional", "blas", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "C23": ("additional", "c23", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "Intel": ("additional", "intel", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "LTO": ("additional", "lto", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "M1mac": ("additional", "m1mac", "proxy", "r-binary", "runner", "macos-15"),
    "MKL": ("additional", "mkl", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/mkl@sha256:d84847130b3ae0b9b0402208ee17360ab527df1c448c44ba045b44524b17620a"),
    "OpenBLAS": ("additional", "openblas", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2"),
    "Strict": ("additional", "strict", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "clang-ASAN": ("additional", "clang-asan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang-asan@sha256:dfab3d2274151577eb705d2be9acd3391798ba4b56d384ca9df84544e5b6be96"),
    "clang-UBSAN": ("additional", "clang-ubsan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang-ubsan@sha256:a58b00b52e9c4b210c3474eac4c28dc18bf70c912d75bb45e466410452a694e6"),
    "donttest": ("additional", "donttest", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/donttest@sha256:fd1942c8b8627d7e1d80582a023b2792acfa158a6edfeeeba3edd27a52c57967"),
    "gcc-ASAN": ("additional", "gcc-asan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f"),
    "gcc-UBSAN": ("additional", "gcc-ubsan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f"),
    "gcc": ("additional", "gcc", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "gcc15": ("additional", "gcc15", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "noLD": ("additional", "nold", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/nold@sha256:9dac23acd0e610b5fb5f82957c8136ffc6acf216a328a067795241542b3dc091"),
    "noOMP": ("additional", "noomp", "proxy", "native-wrapper", "image", "local/colocboost-noomp-proxy"),
    "noRemap": ("additional", "noremap", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "noSuggests": ("additional", "nosuggests", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/nosuggests@sha256:588bd73470d657d0d2560a90e5fe036d191a89fde394b2c67c6423b71d7a12df"),
    "valgrind": ("additional", "valgrind", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/valgrind@sha256:98dfda5016513c269e33c8155ea57745d977dad4b1a703b28db82fcb1237ef9c"),
    "0len": ("additional", "0len", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "rchk": ("additional", "rchk", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "rcnst": ("additional", "rcnst", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-clang@sha256:b66e5f86ce6f8fa3e6afd497fabfdb3aeac79c74ee3525e8a8814de79aa2bc82"),
    "rlibro": ("additional", "rlibro", "proxy", "native-wrapper", "image", "local/colocboost-rlibro-proxy"),
    "musl": ("additional", "musl", "direct", "native-wrapper", "image", "cran-linked-public-musl-reproduction"),
    "linux-arm64": ("additional", "linux-arm64", "direct", "r-binary", "runner", ("ubuntu-24.04", "ubuntu-24.04-arm")),
    "vnu": ("additional", "vnu", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/vnu@sha256:03f45d5944fc092627cae2e944f16f037fed9795f8768c90d9511799098a7884"),
}

GOLDEN_COVERAGE_RESULT_KEYS = [
    (
        "applicability" if row[2] == "not-applicable" else "package-check",
        row[1],
    )
    for row in GOLDEN_ROWS.values()
]
GOLDEN_UNIT_RESULT_KEYS = [("unit", environment_id) for environment_id in UNIT_COVERAGE_IDS]
GOLDEN_RESULT_KEYS = GOLDEN_COVERAGE_RESULT_KEYS + GOLDEN_UNIT_RESULT_KEYS


def proof_set(tokens):
    return frozenset(tokens.split())


GOLDEN_PROOFS = {
    "r-devel-linux-x86_64-debian-clang": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler clang-22 flang-22 libcxx architecture locale"),
    "r-devel-linux-x86_64-debian-gcc": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler gcc-16 architecture locale"),
    "r-devel-linux-x86_64-fedora-clang": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler clang-22 flang-22 libcxx architecture locale"),
    "r-devel-linux-x86_64-fedora-gcc": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system fedora-44 compiler gcc-version architecture locale"),
    "r-devel-windows-x86_64": proof_set("source-tarball-sha256 r-executable r-version r-devel-revision operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "r-patched-linux-x86_64": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-revision r-patched operating-system compiler architecture locale"),
    "r-release-linux-x86_64": proof_set("source-tarball-sha256 container-system-r r-executable r-version r-release operating-system compiler architecture locale"),
    "r-release-macos-arm64": proof_set("source-tarball-sha256 r-executable r-version r-release operating-system macos-version compiler apple-and-gnu-compilers architecture arm64 locale"),
    "r-release-macos-x86_64": proof_set("source-tarball-sha256 r-executable r-version r-release operating-system macos-version compiler apple-and-gnu-compilers architecture x86-64 locale"),
    "r-release-windows-x86_64": proof_set("source-tarball-sha256 r-executable r-version r-release operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "r-oldrel-macos-arm64": proof_set("source-tarball-sha256 r-executable r-version r-oldrel operating-system macos-version compiler apple-and-gnu-compilers architecture arm64 locale"),
    "r-oldrel-macos-x86_64": proof_set("source-tarball-sha256 r-executable r-version r-oldrel operating-system macos-version compiler apple-and-gnu-compilers architecture x86-64 locale"),
    "r-oldrel-windows-x86_64": proof_set("source-tarball-sha256 r-executable r-version r-oldrel operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "ATLAS": proof_set("source-tarball-sha256 container-system-r r-executable extsoftversion-blas la-library process-mappings loaded-atlas"),
    "BLIS": proof_set("source-tarball-sha256 r-executable r-devel-revision operating-system compiler architecture loaded-blis blis-version single-thread"),
    "BLAS": proof_set("source-tarball-sha256 no-package-owned-native-source no-linkingto no-compilation-marker no-direct-native-call"),
    "C23": proof_set("source-tarball-sha256 no-package-c-source obsolete-special-area"),
    "Intel": proof_set("source-tarball-sha256 no-package-owned-native-source obsolete-special-area"),
    "LTO": proof_set("source-tarball-sha256 no-package-owned-native-object no-compilation-target"),
    "M1mac": proof_set("source-tarball-sha256 r-executable r-version r-devel-revision operating-system macos-version compiler architecture arm64"),
    "MKL": proof_set("source-tarball-sha256 container-system-r r-executable r-version operating-system compiler architecture session-info blas-identity lapack-identity la-library la-version process-mappings loaded-mkl blas-operation mkl-verbose single-thread"),
    "OpenBLAS": proof_set("source-tarball-sha256 container-system-r r-executable r-devel-revision operating-system compiler loaded-openblas single-thread"),
    "Strict": proof_set("source-tarball-sha256 no-package-owned-compilation-target strict-r-headers-not-applicable"),
    "clang-ASAN": proof_set("source-tarball-sha256 container-system-r r-executable sanitizer-compiler-flags sanitizer-linker-flags asan-runtime sanitizer-output-scan"),
    "clang-UBSAN": proof_set("source-tarball-sha256 container-system-r r-executable sanitizer-compiler-flags sanitizer-linker-flags ubsan-runtime undefined-behavior-output-scan"),
    "donttest": proof_set("source-tarball-sha256 container-system-r r-executable check-donttest-examples-true expanded-examples-executed"),
    "gcc-ASAN": proof_set("source-tarball-sha256 container-system-r r-executable gcc-asan-flags asan-preload asan-runtime sanitizer-output-scan"),
    "gcc-UBSAN": proof_set("source-tarball-sha256 container-system-r r-executable gcc-ubsan-flags ubsan-preload ubsan-runtime sanitizer-output-scan"),
    "gcc": proof_set("source-tarball-sha256 no-package-owned-compiled-source"),
    "gcc15": proof_set("source-tarball-sha256 no-package-owned-native-source obsolete-special-area"),
    "noLD": proof_set("source-tarball-sha256 container-system-r r-executable opt-r-devel-nold long-double-disabled"),
    "noOMP": proof_set("source-tarball-sha256 r-executable r-devel-revision no-openmp-compile-flags no-openmp-link-flags no-loaded-openmp-runtime dependency-openmp-audit"),
    "noRemap": proof_set("source-tarball-sha256 no-package-cpp-source obsolete-special-area"),
    "noSuggests": proof_set("source-tarball-sha256 container-system-r r-executable depends-only-policy allowed-test-frameworks allowed-vignette-builders ashr-absent susier-absent nonzero-installed-test-count"),
    "valgrind": proof_set("source-tarball-sha256 container-system-r r-executable opt-r-devel-valgrind use-valgrind valgrind-runtime suppression-and-error-scan"),
    "0len": proof_set("source-tarball-sha256 no-direct-native-interface removed-official-experimental-support"),
    "rchk": proof_set("source-tarball-sha256 no-package-c-cpp-object compiled-dependencies-visible-limitation"),
    "rcnst": proof_set("source-tarball-sha256 container-system-r r-executable r-devel-revision r-compile-pkgs-1 r-jit-strategy-4 r-check-constants-5 constant-corruption-diagnostics"),
    "rlibro": proof_set("source-tarball-sha256 r-executable nonroot-uid readonly-bind-mount write-failure installed-library-path package-check-executed"),
    "musl": proof_set("source-tarball-sha256 r-executable r-version musl-libc alpine-version locale architecture"),
    "linux-arm64": proof_set("source-tarball-sha256 r-executable native-amd64 native-arm64 identical-rcheckserver-image identical-tarball-and-configuration architecture-only-comparison"),
    "vnu": proof_set("source-tarball-sha256 container-system-r r-executable vnu-special-dispatch nu-validator-executed zero-bad-entries validator-output"),
}

PREDICATE_PREFIX = ("python", ".github/ci/evaluate_native_features.py", "--tarball", "{tarball}", "--rule")
GOLDEN_PREDICATES = {
    "BLAS": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-source-linkingto-compilation-or-direct-call",)),
    "C23": ("built-source-tarball", PREDICATE_PREFIX + ("no-c-source",)),
    "Intel": ("built-source-tarball", PREDICATE_PREFIX + ("no-package-owned-native-source",)),
    "LTO": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-object",)),
    "Strict": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-compilation-target",)),
    "gcc": ("built-source-tarball", PREDICATE_PREFIX + ("no-compiled-source",)),
    "gcc15": ("built-source-tarball", PREDICATE_PREFIX + ("no-package-owned-native-source",)),
    "noRemap": ("built-source-tarball", PREDICATE_PREFIX + ("no-cpp-source",)),
    "0len": ("built-source-tarball", PREDICATE_PREFIX + ("no-direct-native-call",)),
    "rchk": ("built-source-tarball", PREDICATE_PREFIX + ("no-c-cpp-object",)),
}

GOLDEN_WAIVERS = {
    "installed-model-init-data": ("r-cmd-check-installed", "test_model.R", "colocboost_init_data correctly initializes data", "colocboost_init_data not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-dictionary-mapping": ("r-cmd-check-installed", "test_model.R", "colocboost correctly maps focal outcome to keep_variables with dict_keep_variables", "colocboost_init_data not directly accessible for integration test", 1, "The legacy installed-package integration branch cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-assemble": ("r-cmd-check-installed", "test_model.R", "colocboost_assemble processes model results", "colocboost_assemble not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-workhorse": ("r-cmd-check-installed", "test_model.R", "colocboost_workhorse performs boosting iterations", "colocboost_workhorse not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-utils-dictionary-mapping": ("r-cmd-check-installed", "test_utils.R", "colocboost_init_data handles complex dictionary mappings", "colocboost_init_data not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
}


def predicate_contract(data):
    return {
        row["cran_name"]: (
            row["predicate"]["input"],
            tuple(row["predicate"]["command"]),
        )
        for row in data["coverage"]
        if row["state"] == "not-applicable"
    }


def test_all_ten_committed_predicates_match_independent_golden_contract(manifest):
    assert len(GOLDEN_PREDICATES) == 10
    assert predicate_contract(manifest) == GOLDEN_PREDICATES


@pytest.mark.parametrize("cran_name", sorted(GOLDEN_PREDICATES))
def test_each_not_applicable_predicate_rejects_coordinated_drift(
    manifest, cran_name
):
    mutated = copy.deepcopy(manifest)
    row_for(mutated, cran_name)["predicate"] = {
        "input": "built-source-tarball",
        "command": ["true", "{tarball}"],
    }

    with pytest.raises(ValueError, match="predicate"):
        validate_manifest(mutated)
@pytest.fixture
def manifest():
    return load_manifest(MANIFEST_PATH)


def row_for(data, cran_name):
    return next(row for row in data["coverage"] if row["cran_name"] == cran_name)


def unit_lane_for(data, environment_id):
    return next(
        lane
        for lane in data["unit_lanes"]
        if lane["environment_id"] == environment_id
    )


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


def test_manifest_top_level_is_closed_over_unit_inventory(manifest):
    assert set(manifest) == {"version", "coverage", "unit_lanes"}

    missing = copy.deepcopy(manifest)
    missing.pop("unit_lanes")
    with pytest.raises(ValueError, match="top-level|only version"):
        validate_manifest(missing)

    extra = copy.deepcopy(manifest)
    extra["undeclared"] = []
    with pytest.raises(ValueError, match="top-level|only version"):
        validate_manifest(extra)


def test_committed_unit_lanes_match_independent_golden_contract(manifest):
    assert manifest["unit_lanes"] == GOLDEN_UNIT_LANES
    assert len(manifest["unit_lanes"]) == 18
    assert all(
        set(lane)
        == {
            "environment_id",
            "coverage_id",
            "mode",
            "suite",
            "runner_context",
            "required_sidecar_tests",
        }
        for lane in manifest["unit_lanes"]
    )
    assert {
        lane["runner_context"] for lane in manifest["unit_lanes"][:-1]
    }.isdisjoint({"r-cmd-check-installed"})


def test_result_inventory_apis_preserve_42_ids_and_declare_60_composite_keys(
    manifest,
):
    assert expected_result_ids(manifest) == [
        row["id"] for row in manifest["coverage"]
    ]
    assert len(expected_result_ids(manifest)) == 42

    coverage_keys = validator_module.expected_coverage_result_keys(manifest)
    unit_keys = validator_module.expected_unit_result_keys(manifest)
    result_keys = validator_module.expected_result_keys(manifest)

    assert coverage_keys == GOLDEN_COVERAGE_RESULT_KEYS
    assert unit_keys == GOLDEN_UNIT_RESULT_KEYS
    assert result_keys == GOLDEN_RESULT_KEYS
    assert len(coverage_keys) == 42
    assert sum(kind == "package-check" for kind, _ in coverage_keys) == 32
    assert sum(kind == "applicability" for kind, _ in coverage_keys) == 10
    assert len(unit_keys) == 18
    assert len(result_keys) == len(set(result_keys)) == 60


def test_coverage_order_is_part_of_the_aggregate_contract(manifest):
    mutated = copy.deepcopy(manifest)
    mutated["coverage"][0], mutated["coverage"][1] = (
        mutated["coverage"][1],
        mutated["coverage"][0],
    )

    with pytest.raises(ValueError, match="coverage order"):
        validate_manifest(mutated)


def test_manifest_cli_reports_coverage_unit_and_composite_totals(capsys):
    assert validator_module.main(["validate_manifest.py", str(MANIFEST_PATH)]) == 0

    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == (
        "manifest valid: primary=13, additional=29, uncovered=0, "
        "unit=18, results=60\n"
    )


@pytest.mark.parametrize("mutation", ["remove", "extra"])
def test_unit_lane_inventory_rejects_removals_and_extras(manifest, mutation):
    mutated = copy.deepcopy(manifest)
    if mutation == "remove":
        mutated["unit_lanes"].pop()
    else:
        lane = copy.deepcopy(mutated["unit_lanes"][0])
        lane["environment_id"] = "extra-unit-lane"
        mutated["unit_lanes"].append(lane)

    with pytest.raises(ValueError, match="unit lane|unit_lanes"):
        validate_manifest(mutated)


def test_duplicate_unit_environment_ids_are_rejected(manifest):
    mutated = copy.deepcopy(manifest)
    mutated["unit_lanes"][1]["environment_id"] = mutated["unit_lanes"][0][
        "environment_id"
    ]

    with pytest.raises(ValueError, match="duplicate unit environment_id"):
        validate_manifest(mutated)


@pytest.mark.parametrize("field", ["runner_context", "required_sidecar_tests"])
def test_unit_lane_objects_are_closed_and_exact(manifest, field):
    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, "mkl").pop(field)

    with pytest.raises(ValueError, match="unit lane.*fields"):
        validate_manifest(mutated)

    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, "mkl")["undeclared"] = True
    with pytest.raises(ValueError, match="unit lane.*fields"):
        validate_manifest(mutated)


@pytest.mark.parametrize("bad_id", ["MKL", "mkl_unit", "mkl unit", "-mkl", "mkl-"])
def test_unsafe_unit_environment_ids_are_rejected(manifest, bad_id):
    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, "mkl")["environment_id"] = bad_id

    with pytest.raises(ValueError, match="unit environment_id.*artifact-safe"):
        validate_manifest(mutated)


def test_non_string_unit_environment_id_fails_closed(manifest):
    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, "mkl")["environment_id"] = ["mkl"]

    with pytest.raises(ValueError, match="unit lane.*environment_id"):
        validate_manifest(mutated)


@pytest.mark.parametrize("coverage_id", ["missing-coverage", "blas"])
def test_unit_lanes_reference_only_direct_or_proxy_coverage(
    manifest, coverage_id
):
    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, "mkl")["coverage_id"] = coverage_id

    with pytest.raises(ValueError, match="direct or proxy coverage"):
        validate_manifest(mutated)


@pytest.mark.parametrize(
    ("environment_id", "field", "value"),
    [
        ("mkl", "mode", "installed"),
        ("mkl", "mode", "not-applicable"),
        ("mkl", "suite", "targeted"),
        ("mkl", "runner_context", "r-cmd-check-installed"),
        ("mkl", "required_sidecar_tests", []),
        ("mkl", "required_sidecar_tests", ["test_utils.R"]),
        ("mkl", "required_sidecar_tests", ["test_Xref.R", "test_utils.R"]),
        ("atlas", "required_sidecar_tests", ["test_utils.R"]),
        ("nosuggests", "mode", "source"),
        ("nosuggests", "runner_context", "nosuggests"),
    ],
)
def test_unit_lane_mode_context_suite_and_sidecars_cannot_drift(
    manifest, environment_id, field, value
):
    mutated = copy.deepcopy(manifest)
    unit_lane_for(mutated, environment_id)[field] = value

    with pytest.raises(ValueError, match="unit lane.*approved contract"):
        validate_manifest(mutated)


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
        r"(?i)libopenblas",
        r"(?i)(?:^|/)lib(?:s?atlas)(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libblis(?:\.so(?:\.\d+)*)?(?:$|\s)",
    }


def manifest_core(data):
    result = {}
    for row in data["coverage"]:
        endpoint_name = "image" if "image" in row else "runner"
        endpoint = row[endpoint_name]
        if isinstance(endpoint, list):
            endpoint = tuple(endpoint)
        result[row["cran_name"]] = (
            row["group"],
            row["id"],
            row["state"],
            row["driver"],
            endpoint_name,
            endpoint,
        )
    return result


def test_every_committed_row_matches_independent_golden_contract(manifest):
    assert manifest_core(manifest) == GOLDEN_ROWS
    assert set(GOLDEN_PROOFS) == set(GOLDEN_ROWS)
    for row in manifest["coverage"]:
        assert GOLDEN_PROOFS[row["cran_name"]] <= set(row["proof"])


@pytest.mark.parametrize(
    "mutation",
    [
        "arbitrary-same-name-digest",
        "local-mkl-image",
        "mkl-driver-endpoint",
        "openblas-proof",
        "c23-constant-predicate",
        "state-swap",
        "id-remap",
    ],
)
def test_coordinated_contract_mutations_are_rejected(manifest, mutation):
    mutated = copy.deepcopy(manifest)
    if mutation == "arbitrary-same-name-digest":
        row = row_for(mutated, "MKL")
        row["image"] = row["image"][:-1] + ("b" if row["image"][-1] != "b" else "c")
    elif mutation == "local-mkl-image":
        row_for(mutated, "MKL")["image"] = "local/mkl"
    elif mutation == "mkl-driver-endpoint":
        row = row_for(mutated, "MKL")
        row["driver"] = "r-binary"
        row["runner"] = "ubuntu-24.04"
        row.pop("image")
    elif mutation == "openblas-proof":
        row_for(mutated, "OpenBLAS")["proof"].remove("loaded-openblas")
    elif mutation == "c23-constant-predicate":
        row_for(mutated, "C23")["predicate"]["command"] = ["true", "{tarball}"]
    elif mutation == "state-swap":
        row_for(mutated, "ATLAS")["state"] = "proxy"
        row_for(mutated, "ATLAS")["limitation"] = "mutated"
        row_for(mutated, "BLIS")["state"] = "direct"
    else:
        row_for(mutated, "MKL")["id"] = "mkl-remapped"

    with pytest.raises(ValueError):
        validate_manifest(mutated)


@pytest.mark.parametrize("cran_name", sorted(GOLDEN_PROOFS))
def test_every_row_rejects_removal_of_a_required_proof(manifest, cran_name):
    mutated = copy.deepcopy(manifest)
    removed = sorted(GOLDEN_PROOFS[cran_name])[0]
    row_for(mutated, cran_name)["proof"].remove(removed)

    with pytest.raises(ValueError):
        validate_manifest(mutated)


def waiver_contract(policy):
    return {
        waiver["id"]: (
            waiver["context"],
            waiver["file"],
            waiver["test_title"],
            waiver["reason"],
            waiver["expected_count"],
            waiver["rationale"],
            waiver["expires"],
        )
        for waiver in policy["unit_tests"]["allowed_skips"]
    }


def test_current_policy_matches_independent_golden_and_is_valid():
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))

    assert waiver_contract(policy) == GOLDEN_WAIVERS
    assert validator_module.validate_policy(policy, today=date(2026, 9, 3)) == policy


def test_expired_waiver_is_rejected():
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    policy["unit_tests"]["allowed_skips"][0]["expires"] = "2026-09-02"

    with pytest.raises(ValueError, match="expired"):
        validator_module.validate_policy(policy, today=date(2026, 9, 3))


@pytest.mark.parametrize("expires", ["2027-3-3", "not-a-date", 20270303])
def test_malformed_waiver_expiry_is_rejected(expires):
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    policy["unit_tests"]["allowed_skips"][0]["expires"] = expires

    with pytest.raises(ValueError, match="ISO expiry"):
        validator_module.validate_policy(policy, today=date(2026, 9, 3))


@pytest.mark.parametrize(
    "mapped_name",
    ["libopenblas.so.0", "libopenblasp-r0.3.26.so", "libopenblas64_.so.0"],
)
def test_openblas_forbidden_pattern_catches_common_mapped_names(mapped_name):
    policy = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    patterns = policy["numerical_backends"]["mkl"]["forbidden_library_patterns"]

    assert any(re.search(pattern, f"/usr/lib/{mapped_name}") for pattern in patterns)
