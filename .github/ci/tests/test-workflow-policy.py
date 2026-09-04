import re
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
WORKFLOW_PATH = ROOT / ".github" / "workflows" / "cran-preflight.yml"
FORK_GUARD = "github.repository == 'xueweic/colocboost'"
CHECKOUT = "actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1"
SETUP_PIXI = "prefix-dev/setup-pixi@d3f436a425481402e6a95a1d1fc10331c708cd9e"
UPLOAD = "actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a"
DOWNLOAD = "actions/download-artifact@37930b1c2abaa49bbe596cd826c3c89aef350131"
SETUP_R = "r-lib/actions/setup-r@465b7d8e732ca3921382b1674c59bada9cbf3399"
ALLOWED_ACTIONS = {CHECKOUT, SETUP_PIXI, SETUP_R, UPLOAD, DOWNLOAD}
MKL_IMAGE = "ghcr.io/r-hub/containers/mkl@sha256:d84847130b3ae0b9b0402208ee17360ab527df1c448c44ba045b44524b17620a"
NOSUGGESTS_IMAGE = "ghcr.io/r-hub/containers/nosuggests@sha256:588bd73470d657d0d2560a90e5fe036d191a89fde394b2c67c6423b71d7a12df"
WRAPPER_SHA256 = {
    "mkl": "680c40e9c70355c4dd47c4404b4107fc980748efd08da62a1784aeb5a787d87e",
    "nosuggests": "a8559fbb931019e59aa683aa631deeaa0d14fe96719d058d9af4e9c95ba991c0",
}


def load_workflow():
    assert WORKFLOW_PATH.is_file(), "fork-scoped preflight workflow is missing"
    text = WORKFLOW_PATH.read_text(encoding="utf-8")
    document = yaml.safe_load(text)
    assert isinstance(document, dict)
    assert "on" in document, 'workflow must quote the top-level "on" key'
    return text, document


def step_with_id(job, step_id):
    matches = [step for step in job["steps"] if step.get("id") == step_id]
    assert len(matches) == 1, step_id
    return matches[0]


def steps_using(job, action):
    return [step for step in job["steps"] if step.get("uses") == action]


def test_triggers_permissions_concurrency_and_job_inventory_are_exact():
    _, workflow = load_workflow()

    assert workflow["on"] == {
        "push": {"branches": ["**"]},
        "pull_request": None,
    }
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["env"] == {"PYTHONDONTWRITEBYTECODE": "1"}
    assert workflow["concurrency"] == {
        "group": "cran-preflight-${{ github.event_name }}-${{ github.ref }}",
        "cancel-in-progress": True,
    }
    assert set(workflow["jobs"]) == {
        "prepare", "primary-linux", "primary-platform", "mkl", "nosuggests",
        "active-rhub", "remaining-docker", "linux-arm64", "linux-arm64-compare",
        "applicability-native", "summary-gate",
    }
    assert workflow["jobs"]["primary-linux"]["strategy"]["fail-fast"] is False
    assert workflow["jobs"]["primary-platform"]["strategy"]["fail-fast"] is False
    assert workflow["jobs"]["active-rhub"]["strategy"]["fail-fast"] is False


def test_every_job_is_inert_outside_the_fork_and_dependents_stop_when_cancelled():
    _, workflow = load_workflow()

    assert workflow["jobs"]["prepare"]["if"] == "${{ " + FORK_GUARD + " }}"
    for name in ("primary-linux", "primary-platform", "mkl", "nosuggests", "active-rhub", "remaining-docker", "linux-arm64", "linux-arm64-compare", "applicability-native"):
        assert workflow["jobs"][name]["if"] == (
            "${{ always() && !cancelled() && " + FORK_GUARD + " }}"
        )
        expected_needs = ["prepare", "linux-arm64"] if name == "linux-arm64-compare" else "prepare"
        assert workflow["jobs"][name]["needs"] == expected_needs
    summary = workflow["jobs"]["summary-gate"]
    assert summary["if"] == (
        "${{ always() && !cancelled() && " + FORK_GUARD + " }}"
    )
    assert summary["needs"] == [
        "prepare", "primary-linux", "primary-platform", "mkl", "nosuggests",
        "active-rhub", "applicability-native", "remaining-docker", "linux-arm64", "linux-arm64-compare",
    ]


def test_proxy_dockerfiles_are_passed_to_validator_as_absolute_paths():
    _, workflow = load_workflow()
    build = step_with_id(workflow["jobs"]["remaining-docker"], "build-image")

    assert build["env"]["DOCKERFILE"] == (
        "${{ github.workspace }}/${{ matrix.dockerfile }}"
    )


def test_as_cran_platform_checks_disable_only_remote_incoming_version_lookup():
    _, workflow = load_workflow()
    platform = workflow["jobs"]["primary-platform"]

    assert platform["env"]["_R_CHECK_CRAN_INCOMING_REMOTE_"] == "false"


def test_primary_linux_matrix_is_exact_and_binds_native_runtime_contracts():
    _, workflow = load_workflow()
    job = workflow["jobs"]["primary-linux"]
    rows = job["strategy"]["matrix"]["include"]
    expected = {
        "r-devel-linux-x86-64-debian-clang": ("ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e", "/opt/R/devel/bin/R", "clang22", "standard-no-documentation", "--no-manual --no-build-vignettes"),
        "r-devel-linux-x86-64-debian-gcc": ("ghcr.io/r-hub/containers/ubuntu-gcc16@sha256:2e9576e51ad17a706887b7e06fc4057388765226ec795f3c7af0e7348fb8fbf1", "/opt/R/devel/bin/R", "ubuntu-gcc16", "standard-no-documentation", "--no-manual --no-build-vignettes"),
        "r-devel-linux-x86-64-fedora-clang": ("ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e", "/opt/R/devel/bin/R", "clang22", "standard-no-documentation", "--no-manual --no-build-vignettes"),
        "r-devel-linux-x86-64-fedora-gcc": ("ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2", "/opt/R/devel-gcc16/bin/R", "gcc16", "standard-no-documentation", "--no-manual --no-build-vignettes"),
        "r-patched-linux-x86-64": ("ghcr.io/r-hub/containers/ubuntu-next@sha256:1c29eed93b0aa05fe147e476e6b3510a09fb0476460eeefbbe99a4a1eb2c3bd9", "/opt/R/next/bin/R", "ubuntu-next", "standard-no-documentation", "--no-manual --no-build-vignettes"),
        "r-release-linux-x86-64": ("ghcr.io/r-hub/containers/ubuntu-release@sha256:714722b7ecb4307fbf88a707f83fb045181482c0c8ce6086a0e26860416ca9c2", "/opt/R/release/bin/R", "ubuntu-release", "full-documentation", ""),
    }
    assert len(rows) == len(expected)
    assert len({row["environment_id"] for row in rows}) == len(rows)
    assert {row["environment_id"] for row in rows} == set(expected)
    for row in rows:
        image, system_r, profile, check_profile, check_args = expected[row["environment_id"]]
        assert row["image"] == image
        assert row["system_r"] == system_r
        assert row["runtime_profile"] == profile
        assert row["check_profile"] == check_profile
        assert row["check_args"] == check_args
        assert row["wrapper"] == "/usr/local/bin/r-check"
        assert row["wrapper_sha256"] == "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090"
    assert job["container"] == {"image": "${{ matrix.image }}", "options": "--user 0"}
    assert job["defaults"] == {"run": {"shell": "bash"}}


def test_primary_linux_children_have_two_independent_results_and_full_diagnostics():
    _, workflow = load_workflow()
    job = workflow["jobs"]["primary-linux"]
    ids = [step.get("id") for step in job["steps"]]
    for step_id in (
        "verify-source", "extract-source", "prepare-unit-dependencies",
        "verify-unit-dependencies", "probe-runtime", "verify-runtime", "unit-full",
        "adapt-unit", "unit-result-gate", "native-check", "package-result-gate",
        "finalize-unit", "finalize-check", "upload-unit", "upload-check",
        "upload-diagnostics", "producer-gate",
    ):
        assert ids.count(step_id) == 1
    package_gate = step_with_id(job, "package-result-gate")
    assert "UNIT_GATE_OUTCOME" not in package_gate.get("env", {})
    assert "steps.unit-full.outcome" not in package_gate.get("env", {}).values()
    assert package_gate["env"]["SEMANTIC_PROOF_OUTCOME"] == (
        "${{ steps.verify-check-semantics.outcome }}"
    )
    unit_gate = step_with_id(job, "unit-result-gate")
    assert "NATIVE_CHECK_OUTCOME" not in unit_gate.get("env", {})
    for step_id in ("finalize-unit", "finalize-check", "upload-unit", "upload-check", "upload-diagnostics", "producer-gate"):
        assert step_with_id(job, step_id)["if"] == "${{ always() }}"
    diagnostics = step_with_id(job, "upload-diagnostics")
    assert diagnostics["with"]["include-hidden-files"] is True
    assert diagnostics["with"]["if-no-files-found"] == "error"
    assert "*.Rcheck" in diagnostics["with"]["path"]
    assert "ci-verify-runtime" in step_with_id(job, "verify-runtime")["run"]
    release_proof = step_with_id(job, "verify-check-semantics")
    assert "ci-special-check" in release_proof["run"]
    assert '"$CHECK_PROFILE"' in release_proof["run"]
    assert "matrix.environment_id" not in release_proof["run"]
    ids = [step.get("id") for step in job["steps"]]
    assert ids.index("verify-check-semantics") < ids.index("package-result-gate")
    assert ids.index("finalize-unit") < ids.index("upload-unit")
    assert ids.index("finalize-check") < ids.index("upload-check")
    assert max(ids.index("upload-unit"), ids.index("upload-check")) < ids.index(
        "producer-gate"
    )
    assert step_with_id(job, "upload-unit")["with"] == {
        "name": "cran-preflight-result-unit-${{ matrix.environment_id }}",
        "path": "${{ runner.temp }}/primary-linux-results/unit/result.json",
        "if-no-files-found": "error",
    }
    assert step_with_id(job, "upload-check")["with"] == {
        "name": "cran-preflight-result-package-check-${{ matrix.environment_id }}",
        "path": "${{ runner.temp }}/primary-linux-results/package-check/result.json",
        "if-no-files-found": "error",
    }


def test_primary_platform_matrix_is_exact_and_binds_setup_r_identity():
    _, workflow = load_workflow()
    job = workflow["jobs"]["primary-platform"]
    rows = job["strategy"]["matrix"]["include"]
    expected = {
        "r-devel-windows-x86-64": (
            "windows-2022", "devel", "C:/R/bin/R.exe", "devel", "windows",
            "x86_64", True,
        ),
        "r-release-windows-x86-64": (
            "windows-2022", "release", "C:/R/bin/R.exe", "release", "windows",
            "x86_64", True,
        ),
        "r-oldrel-windows-x86-64": (
            "windows-2022", "oldrel-1", "C:/R/bin/R.exe", "release", "windows",
            "x86_64", True,
        ),
        "r-release-macos-arm64": (
            "macos-15", "release",
            "/Library/Frameworks/R.framework/Resources/bin/R", "release",
            "macos", "aarch64", True,
        ),
        "r-oldrel-macos-arm64": (
            "macos-15", "oldrel-1",
            "/Library/Frameworks/R.framework/Resources/bin/R", "release",
            "macos", "aarch64", True,
        ),
        "r-release-macos-x86-64": (
            "macos-15-intel", "release",
            "/Library/Frameworks/R.framework/Resources/bin/R", "release",
            "macos", "x86_64", True,
        ),
        "r-oldrel-macos-x86-64": (
            "macos-15-intel", "oldrel-1",
            "/Library/Frameworks/R.framework/Resources/bin/R", "release",
            "macos", "x86_64", True,
        ),
        "m1mac": (
            "macos-15", "devel",
            "/Library/Frameworks/R.framework/Resources/bin/R", "devel",
            "macos", "aarch64", False,
        ),
    }
    assert len(rows) == len(expected) == 8
    for row in rows:
        assert (
            row["runner"], row["setup_r_selector"], row["system_r"],
            row["expected_r_kind"], row["expected_os"],
            row["expected_architecture"], row["unit_enabled"],
        ) == expected[row["environment_id"]]
        assert row["check_args"] == "--as-cran --no-manual --no-build-vignettes"
    assert job["runs-on"] == "${{ matrix.runner }}"
    assert job["defaults"] == {"run": {"shell": "pwsh"}}
    setup_r = step_with_id(job, "setup-r")
    assert setup_r["uses"] == SETUP_R
    assert setup_r["continue-on-error"] is True
    assert setup_r["with"]["r-version"] == "${{ matrix.setup_r_selector }}"
    assert setup_r["with"]["use-public-rspm"] is True
    identity = step_with_id(job, "verify-platform-r")
    assert identity["env"]["SETUP_R_VERSION"] == (
        "${{ steps.setup-r.outputs.installed-r-version }}"
    )
    assert '"-B", ".github/ci/verify_platform_r.py"' in identity["run"]
    assert "pixi run --locked python @arguments" in identity["run"]
    assert "ci-platform-r" not in identity["run"]


def test_primary_platform_children_are_cross_platform_and_fail_closed():
    _, workflow = load_workflow()
    job = workflow["jobs"]["primary-platform"]
    ids = [step.get("id") for step in job["steps"]]
    for step_id in (
        "initialize", "download-source", "verify-source", "extract-source",
        "verify-platform-r", "prepare-unit-dependencies",
        "verify-unit-dependencies", "unit-full", "adapt-unit",
        "unit-result-gate", "prepare-check-dependencies",
        "verify-check-dependencies", "r-binary-check", "package-result-gate",
        "finalize-unit", "finalize-check", "upload-unit", "upload-check",
        "upload-diagnostics", "producer-gate",
    ):
        assert ids.count(step_id) == 1
    for step in job["steps"]:
        run = step.get("run", "")
        assert "set -" not in run
        assert "export " not in run
        assert " \\\n" not in run
    assert '"-B", ".github/ci/verify_source.py"' in step_with_id(
        job, "verify-source"
    )["run"]
    assert '"-B", ".github/ci/extract_source.py"' in step_with_id(
        job, "extract-source"
    )["run"]
    assert '"-B", ".github/ci/run_r_binary_check.py"' in step_with_id(
        job, "r-binary-check"
    )["run"]
    unit_gate = step_with_id(job, "unit-result-gate")
    package_gate = step_with_id(job, "package-result-gate")
    assert unit_gate["env"]["SETUP_R_OUTCOME"] == "${{ steps.setup-r.outcome }}"
    assert package_gate["env"]["SETUP_R_OUTCOME"] == "${{ steps.setup-r.outcome }}"
    assert "BINARY_CHECK_OUTCOME" not in unit_gate["env"]
    assert "UNIT_OUTCOME" not in package_gate["env"]
    assert "ADAPTER_OUTCOME" not in package_gate["env"]
    for step_id in (
        "prepare-unit-dependencies", "verify-unit-dependencies", "unit-full",
        "adapt-unit", "unit-result-gate",
    ):
        assert step_with_id(job, step_id)["if"] == "${{ matrix.unit_enabled }}"
    assert step_with_id(job, "finalize-unit")["if"] == (
        "${{ always() && matrix.unit_enabled }}"
    )
    assert step_with_id(job, "upload-unit")["if"] == (
        "${{ always() && matrix.unit_enabled }}"
    )
    for step_id in ("finalize-check", "upload-check", "upload-diagnostics", "producer-gate"):
        assert step_with_id(job, step_id)["if"] == "${{ always() }}"
    for step_id in (
        "prepare-check-dependencies", "verify-check-dependencies",
        "r-binary-check", "package-result-gate",
    ):
        assert "if" not in step_with_id(job, step_id)
    diagnostics = step_with_id(job, "upload-diagnostics")
    assert diagnostics["with"]["include-hidden-files"] is True
    assert diagnostics["with"]["if-no-files-found"] == "error"
    assert "*.Rcheck" in diagnostics["with"]["path"]
    producer_gate = step_with_id(job, "producer-gate")
    assert producer_gate["env"]["UNIT_ENABLED"] == "${{ matrix.unit_enabled }}"
    assert 'if ($env:UNIT_ENABLED -eq "true")' in producer_gate["run"]
    assert 'elseif ($env:UNIT_ENABLED -eq "false")' in producer_gate["run"]
    assert "else {\n  exit 1" in producer_gate["run"]


def test_task11_final_six_producers_produce_exactly_60_terminal_keys():
    _, workflow = load_workflow()
    keys = set()
    for job_name in ("primary-linux", "primary-platform"):
        rows = workflow["jobs"][job_name]["strategy"]["matrix"]["include"]
        for row in rows:
            keys.add(("package-check", row["environment_id"]))
            if row.get("unit_enabled", True):
                keys.add(("unit", row["environment_id"]))
    for environment_id in ("mkl", "nosuggests"):
        keys.add(("package-check", environment_id))
        keys.add(("unit", environment_id))
    active_rows = workflow["jobs"]["active-rhub"]["strategy"]["matrix"]["include"]
    for row in active_rows:
        keys.add(("package-check", row["environment_id"]))
        if row["unit_enabled"]:
            keys.add(("unit", row["environment_id"]))
    for row in workflow["jobs"]["applicability-native"]["strategy"]["matrix"]["include"]:
        keys.add(("applicability", row["environment_id"]))
    for row in workflow["jobs"]["remaining-docker"]["strategy"]["matrix"]["include"]:
        keys.add(("package-check", row["environment_id"]))
        if row["unit_enabled"]:
            keys.add(("unit", row["environment_id"]))
    keys.add(("package-check", "linux-arm64"))
    assert len(keys) == 60
    assert sum(kind == "package-check" for kind, _ in keys) == 32
    assert sum(kind == "unit" for kind, _ in keys) == 18
    assert sum(kind == "applicability" for kind, _ in keys) == 10
    assert ("package-check", "m1mac") in keys
    assert ("unit", "m1mac") not in keys


def test_all_actions_are_immutable_allowlisted_and_checkout_is_safe():
    _, workflow = load_workflow()

    for job in workflow["jobs"].values():
        action_steps = [step for step in job["steps"] if "uses" in step]
        for step in action_steps:
            assert step["uses"] in ALLOWED_ACTIONS
            assert re.search(r"@[0-9a-f]{40}$", step["uses"])

        checkout_steps = steps_using(job, CHECKOUT)
        assert len(checkout_steps) == 1
        assert checkout_steps[0]["with"] == {
            "ref": "${{ github.sha }}",
            "fetch-depth": 1,
            "persist-credentials": False,
        }

        setup_steps = steps_using(job, SETUP_PIXI)
        assert len(setup_steps) == 1
        assert setup_steps[0]["with"] == {
            "pixi-version": "v0.67.2",
            "locked": True,
            "cache": True,
            "cache-write": "${{ github.event_name == 'push' }}",
            "activate-environment": False,
        }


def test_special_producers_use_pinned_job_containers_and_system_r_bindings():
    _, workflow = load_workflow()
    expected = {
        "mkl": (MKL_IMAGE, "/opt/R/devel-mkl/bin/R"),
        "nosuggests": (NOSUGGESTS_IMAGE, "/opt/R/devel/bin/R"),
    }
    for name, (image, system_r) in expected.items():
        job = workflow["jobs"][name]
        assert job["runs-on"] == "ubuntu-24.04"
        assert job["container"] == {"image": image, "options": "--user 0"}
        assert job["defaults"] == {"run": {"shell": "bash"}}
        assert job["env"]["SYSTEM_R"] == system_r
        assert job["env"]["NATIVE_WRAPPER"] == "/usr/local/bin/r-check"
        assert job["env"]["WRAPPER_SHA256"] == WRAPPER_SHA256[name]
        assert all(
            "runner.temp" not in str(value) for value in job.get("env", {}).values()
        )
        initialize = step_with_id(job, "initialize")["run"]
        assert "$GITHUB_ENV" in initialize
        assert "$RUNNER_TEMP" in initialize
        for variable in (
            "UNIT_LIBRARY",
            "CHECK_LIBRARY",
            "SOURCE_DOWNLOAD_DIR",
            "SOURCE_TREE",
            "DIAGNOSTIC_DIR",
            "RESULT_DIR",
            "CHECK_WORK_DIR",
        ):
            assert f"{variable}=" in initialize


def test_active_rhub_matrix_is_exact_and_immutable():
    _, workflow = load_workflow()
    job = workflow["jobs"]["active-rhub"]
    rows = job["strategy"]["matrix"]["include"]
    generic = "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090"
    clang = "9732ef12761fd6ecd4631dd6ce0861fbbbd68fb33f6e7643745ec36de456b4f9"
    expected = {
        "atlas": ("ghcr.io/r-hub/containers/atlas@sha256:7597f7d0b6b2f009ae7bb425391523d8f4388223238db50d7dfb1572add63a88", "/opt/R/devel/bin/R", generic, "atlas", "--no-manual --no-build-vignettes", True, ""),
        "openblas": ("ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2", "/opt/R/devel-gcc16/bin/R", generic, "openblas", "--no-manual --no-build-vignettes", True, ""),
        "rcnst": ("ghcr.io/r-hub/containers/ubuntu-clang@sha256:b66e5f86ce6f8fa3e6afd497fabfdb3aeac79c74ee3525e8a8814de79aa2bc82", "/opt/R/devel/bin/R", generic, "rcnst", "--no-manual --no-build-vignettes", False, ""),
        "clang-asan": ("ghcr.io/r-hub/containers/clang-asan@sha256:dfab3d2274151577eb705d2be9acd3391798ba4b56d384ca9df84544e5b6be96", "/opt/R/devel-asan/bin/R", clang, "clang-asan", "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes", False, ""),
        "clang-ubsan": ("ghcr.io/r-hub/containers/clang-ubsan@sha256:a58b00b52e9c4b210c3474eac4c28dc18bf70c912d75bb45e466410452a694e6", "/opt/R/devel-asan/bin/R", clang, "clang-ubsan", "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes", False, ""),
        "donttest": ("ghcr.io/r-hub/containers/donttest@sha256:fd1942c8b8627d7e1d80582a023b2792acfa158a6edfeeeba3edd27a52c57967", "/opt/R/devel/bin/R", generic, "donttest", "--no-manual --no-build-vignettes", False, ""),
        "gcc-asan": ("ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f", "/opt/R/devel/bin/R", generic, "gcc-asan", "--no-manual --no-build-vignettes", False, ""),
        "gcc-ubsan": ("ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f", "/opt/R/devel/bin/R", generic, "gcc-ubsan", "--no-manual --no-build-vignettes", False, ""),
        "nold": ("ghcr.io/r-hub/containers/nold@sha256:9dac23acd0e610b5fb5f82957c8136ffc6acf216a328a067795241542b3dc091", "/opt/R/devel-nold/bin/R", generic, "nold", "--no-manual --no-build-vignettes", False, ""),
        "valgrind": ("ghcr.io/r-hub/containers/valgrind@sha256:98dfda5016513c269e33c8155ea57745d977dad4b1a703b28db82fcb1237ef9c", "/opt/R/devel-valgrind/bin/R", "79261a338b0a381a157cf2025ef40f389b3906fb4711d14b7a72314f5316cbc8", "valgrind", "--use-valgrind --extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes", False, "valgrind-suppression"),
        "vnu": ("ghcr.io/r-hub/containers/vnu@sha256:03f45d5944fc092627cae2e944f16f037fed9795f8768c90d9511799098a7884", "/opt/R/release/bin/R", "0de8ba373ec3bbe81b84122857d071dda4eee72e035d3180464607837c0bb089", "vnu", "--no-manual --no-build-vignettes", False, "vnu-dispatcher"),
    }
    assert len(rows) == 11
    assert len({row["environment_id"] for row in rows}) == 11
    assert {row["environment_id"] for row in rows} == set(expected)
    for row in rows:
        assert (
            row["image"], row["system_r"], row["wrapper_sha256"],
            row["runtime_profile"], row["check_args"], row["unit_enabled"],
            row["file_purpose"],
        ) == expected[row["environment_id"]]
        assert row["wrapper"] == "/usr/local/bin/r-check"
    assert job["runs-on"] == "ubuntu-24.04"
    assert job["container"] == {"image": "${{ matrix.image }}", "options": "--user 0"}
    assert job["defaults"] == {"run": {"shell": "bash"}}


def test_active_rhub_orders_identity_raw_check_semantics_parse_and_special_dispatch():
    _, workflow = load_workflow()
    job = workflow["jobs"]["active-rhub"]
    ids = [step.get("id") for step in job["steps"]]
    for step_id in (
        "initialize", "download-source", "verify-source", "extract-source",
        "verify-special-file", "probe-runtime", "verify-runtime",
        "native-check", "verify-check-semantics", "run-vnu", "parse-check",
        "package-result-gate", "finalize-check", "upload-check",
        "upload-diagnostics", "producer-gate",
    ):
        assert ids.count(step_id) == 1
    assert ids.index("verify-special-file") < ids.index("native-check")
    assert ids.index("native-check") < ids.index("verify-check-semantics")
    assert ids.index("verify-check-semantics") < ids.index("run-vnu")
    assert ids.index("verify-check-semantics") < ids.index("parse-check")
    native = step_with_id(job, "native-check")
    assert "ci-native-check" in native["run"]
    assert "raw_exit_code=" in native["run"]
    semantics = step_with_id(job, "verify-check-semantics")
    assert "ci-special-check" in semantics["run"]
    parser = step_with_id(job, "parse-check")
    assert parser["env"]["EFFECTIVE_EXIT_CODE"] == (
        "${{ steps.verify-check-semantics.outputs.effective_exit_code }}"
    )
    assert '"$EFFECTIVE_EXIT_CODE"' in parser["run"]
    special_file = step_with_id(job, "verify-special-file")
    assert "ci-file-contract" in special_file["run"]
    assert special_file["if"] == "${{ matrix.file_purpose != '' }}"
    vnu = step_with_id(job, "run-vnu")
    assert vnu["if"] == "${{ matrix.environment_id == 'vnu' }}"
    assert "ci-run-vnu" in vnu["run"]
    assert "/usr/local/bin/vnu.sh" in vnu["run"]


def test_active_atlas_alone_emits_independent_full_source_unit_result():
    _, workflow = load_workflow()
    job = workflow["jobs"]["active-rhub"]
    for step_id in (
        "prepare-unit-dependencies", "verify-unit-dependencies", "unit-full",
        "adapt-unit", "unit-result-gate", "finalize-unit", "upload-unit",
    ):
        step = step_with_id(job, step_id)
        assert "matrix.unit_enabled" in step["if"]
    unit = step_with_id(job, "unit-full")
    assert unit["env"]["PACKAGE_PATH"] == "${{ steps.extract-source.outputs.package_path }}"
    assert "ci-unit source" in unit["run"]
    assert "--filter" not in unit["run"]
    package_gate = step_with_id(job, "package-result-gate")
    assert "UNIT_OUTCOME" not in package_gate["env"]
    assert "ADAPTER_OUTCOME" not in package_gate["env"]
    assert "NATIVE_CHECK_OUTCOME" not in package_gate["env"]
    assert "steps.native-check.outcome" not in package_gate["env"].values()
    assert package_gate["env"]["SEMANTIC_PROOF_OUTCOME"] == (
        "${{ steps.verify-check-semantics.outcome }}"
    )
    unit_gate = step_with_id(job, "unit-result-gate")
    assert "NATIVE_CHECK_OUTCOME" not in unit_gate["env"]


def test_active_results_and_complete_diagnostics_are_fail_closed():
    _, workflow = load_workflow()
    job = workflow["jobs"]["active-rhub"]
    for step_id in ("finalize-check", "upload-check", "upload-diagnostics", "producer-gate"):
        assert step_with_id(job, step_id)["if"] == "${{ always() }}"
    assert step_with_id(job, "finalize-unit")["if"] == "${{ always() && matrix.unit_enabled }}"
    assert step_with_id(job, "upload-unit")["if"] == "${{ always() && matrix.unit_enabled }}"
    diagnostics = step_with_id(job, "upload-diagnostics")
    assert diagnostics["with"]["include-hidden-files"] is True
    assert diagnostics["with"]["if-no-files-found"] == "error"
    assert "*.Rcheck" in diagnostics["with"]["path"]
    assert step_with_id(job, "upload-check")["with"]["name"] == (
        "cran-preflight-result-package-check-${{ matrix.environment_id }}"
    )
    assert step_with_id(job, "upload-unit")["with"]["name"] == (
        "cran-preflight-result-unit-${{ matrix.environment_id }}"
    )
    ids = [step.get("id") for step in job["steps"]]
    assert ids.index("finalize-check") < ids.index("upload-check")
    assert ids.index("finalize-unit") < ids.index("upload-unit")
    assert max(ids.index("upload-check"), ids.index("upload-unit"), ids.index("upload-diagnostics")) < ids.index("producer-gate")


def test_special_producers_consume_verified_tarball_and_not_checkout_package():
    _, workflow = load_workflow()
    for name in ("mkl", "nosuggests"):
        job = workflow["jobs"][name]
        download = step_with_id(job, "download-source")
        assert download["uses"] == DOWNLOAD
        assert download["with"]["name"] == "cran-preflight-source"
        verify = step_with_id(job, "verify-source")
        extract = step_with_id(job, "extract-source")
        assert "ci-verify-source" in verify["run"]
        assert "ci-extract-source" in extract["run"]
        assert extract["env"]["TARBALL"] == "${{ steps.verify-source.outputs.tarball_path }}"
        assert extract["env"]["METADATA"] == "${{ steps.verify-source.outputs.metadata_path }}"
        unit_steps = [
            step
            for step in job["steps"]
            if step.get("id") in {"unit-utils", "unit-xref", "unit-full"}
        ]
        assert unit_steps
        for step in unit_steps:
            expected_package = (
                "${{ steps.extract-source.outputs.package_path }}"
                if name == "mkl"
                else "${{ runner.temp }}/nosuggests-unit-library/colocboost"
            )
            assert step["env"]["PACKAGE_PATH"] == expected_package
            assert '"$PACKAGE_PATH"' in step["run"]
            assert " --package=." not in step["run"]


def test_mkl_order_and_targeted_sidecars_are_strictly_declared():
    _, workflow = load_workflow()
    mkl = workflow["jobs"]["mkl"]
    ids = [step.get("id") for step in mkl["steps"]]
    ordered = [
        "prepare-dependencies",
        "verify-dependencies",
        "verify-mkl",
        "unit-utils",
        "unit-xref",
        "unit-full",
        "adapt-unit",
        "native-check",
    ]
    assert [ids.index(step_id) for step_id in ordered] == sorted(
        ids.index(step_id) for step_id in ordered
    )
    assert '"^utils$"' in step_with_id(mkl, "unit-utils")["run"]
    assert '"^Xref$"' in step_with_id(mkl, "unit-xref")["run"]
    adapter_run = step_with_id(mkl, "adapt-unit")["run"]
    assert '"$DIAGNOSTIC_DIR/unit-utils.json"' in adapter_run
    assert '"$DIAGNOSTIC_DIR/unit-xref.json"' in adapter_run
    assert "verify-mkl" in step_with_id(mkl, "verify-mkl")["run"]
    verify_mkl = step_with_id(mkl, "verify-mkl")["run"]
    assert "source /opt/intel/oneapi/setvars.sh" in verify_mkl
    normalized_verify_mkl = " ".join(verify_mkl.replace("\\\n", " ").split())
    assert (
        'pixi run --locked ci-verify-mkl "$SYSTEM_R" '
        '"$GITHUB_WORKSPACE/.github/ci/check-policy.yml" '
        '"$DIAGNOSTIC_DIR/mkl-runtime.json"'
    ) in normalized_verify_mkl


def test_nosuggests_proves_runtime_policy_installed_tests_and_native_log():
    text, workflow = load_workflow()
    job = workflow["jobs"]["nosuggests"]
    ids = [step.get("id") for step in job["steps"]]
    ordered = [
        "prepare-dependencies",
        "verify-dependencies",
        "unit-full",
        "adapt-unit",
        "native-check",
    ]
    assert [ids.index(step_id) for step_id in ordered] == sorted(
        ids.index(step_id) for step_id in ordered
    )
    assert "verify-dependencies" in step_with_id(job, "verify-dependencies")["run"]
    native = step_with_id(job, "native-check")
    assert "ci-native-check" in native["run"]
    assert "ci-package-result" in native["run"]
    assert "testthat.Rout" in text
    assert "devtools" not in job["env"]


def test_each_special_producer_finalizes_and_uploads_two_independent_results():
    _, workflow = load_workflow()
    for name in ("mkl", "nosuggests"):
        job = workflow["jobs"][name]
        for step_id in ("finalize-unit", "finalize-check"):
            step = step_with_id(job, step_id)
            assert step["if"] == "${{ always() }}"
            assert "ci-finalize-result" in step["run"]
        uploads = [
            step for step in job["steps"]
            if step.get("uses") == UPLOAD and step["with"]["name"].startswith("cran-preflight-result-")
        ]
        assert {step["with"]["name"] for step in uploads} == {
            f"cran-preflight-result-package-check-{name}",
            f"cran-preflight-result-unit-{name}",
        }
        assert all(step["if"] == "${{ always() }}" for step in uploads)
        assert all(step["with"]["if-no-files-found"] == "error" for step in uploads)
        assert step_with_id(job, "producer-gate")["if"] == "${{ always() }}"
        ids = [step.get("id") for step in job["steps"]]
        assert ids.index("finalize-unit") < ids.index("upload-unit")
        assert ids.index("finalize-check") < ids.index("upload-check")
        assert max(
            ids.index("upload-unit"),
            ids.index("upload-check"),
            ids.index("upload-diagnostics"),
        ) < ids.index("producer-gate")


def test_special_diagnostics_uploads_include_complete_native_check_tree():
    _, workflow = load_workflow()
    expected = {
        "mkl": "${{ runner.temp }}/mkl-native-check/input/*.Rcheck",
        "nosuggests": (
            "${{ runner.temp }}/nosuggests-native-check/input/*.Rcheck"
        ),
    }
    for name, check_tree in expected.items():
        upload = step_with_id(workflow["jobs"][name], "upload-diagnostics")
        assert upload["if"] == "${{ always() }}"
        assert upload["continue-on-error"] is True
        assert upload["uses"] == UPLOAD
        assert upload["with"]["path"].splitlines() == [
            f"${{{{ runner.temp }}}}/{name}-diagnostics",
            check_tree,
        ]
        assert upload["with"]["if-no-files-found"] == "error"
        assert upload["with"]["include-hidden-files"] is True


def test_package_and_unit_gates_are_independent():
    _, workflow = load_workflow()
    for name in ("mkl", "nosuggests"):
        package_gate = step_with_id(workflow["jobs"][name], "package-result-gate")
        assert "UNIT_GATE_OUTCOME" not in package_gate["env"]
        assert "steps.unit-full.outcome" not in package_gate["env"].values()
        assert "steps.adapt-unit.outcome" not in package_gate["env"].values()
        assert package_gate["env"]["NATIVE_CHECK_OUTCOME"] == (
            "${{ steps.native-check.outcome }}"
        )


def test_prepare_builds_one_source_with_local_r45_and_uploads_only_its_contract():
    _, workflow = load_workflow()
    prepare = workflow["jobs"]["prepare"]

    assert prepare["runs-on"] == "ubuntu-24.04"
    assert prepare["outputs"] == {
        "source_sha": "${{ steps.prepare-source.outputs.source_sha }}",
        "event_sha": "${{ steps.prepare-source.outputs.event_sha }}",
        "tarball_filename": "${{ steps.prepare-source.outputs.tarball_filename }}",
        "tarball_sha256": "${{ steps.prepare-source.outputs.tarball_sha256 }}",
    }
    assert all("path" not in output for output in prepare["outputs"])
    run_steps = [step for step in prepare["steps"] if "run" in step]
    assert any("pixi run --locked ci-validate" in step["run"] for step in run_steps)
    assert any(
        "pixi run --locked ci-contract-tests" in step["run"]
        for step in run_steps
    )
    source = step_with_id(prepare, "prepare-source")
    assert source["env"] == {"EVENT_SHA": "${{ github.sha }}"}
    assert "pixi run --locked --environment local-r45 ci-prepare-source" in source["run"]
    assert '"$EVENT_SHA"' in source["run"]
    assert '"$GITHUB_OUTPUT"' in source["run"]

    uploads = steps_using(prepare, UPLOAD)
    assert len(uploads) == 1
    upload = uploads[0]
    assert upload["with"]["name"] == "cran-preflight-source"
    assert upload["with"]["if-no-files-found"] == "error"
    assert "steps.prepare-source.outputs.tarball_path" in upload["with"]["path"]
    assert "steps.prepare-source.outputs.metadata_path" in upload["with"]["path"]


def test_summary_downloads_without_flattening_and_verifies_raw_tarball_bytes():
    text, workflow = load_workflow()
    summary = workflow["jobs"]["summary-gate"]
    assert summary["runs-on"] == "ubuntu-24.04"

    source_download = step_with_id(summary, "download-source")
    assert source_download["uses"] == DOWNLOAD
    assert source_download["continue-on-error"] is True
    assert source_download["with"]["name"] == "cran-preflight-source"

    result_download = step_with_id(summary, "download-results")
    assert result_download["uses"] == DOWNLOAD
    assert result_download["continue-on-error"] is True
    assert result_download["with"]["pattern"] == "cran-preflight-result-*"
    assert result_download["with"]["merge-multiple"] is False

    verify = step_with_id(summary, "verify-source")
    assert verify["continue-on-error"] is True
    assert verify["env"] == {"EVENT_SHA": "${{ github.sha }}"}
    assert "ci-verify-source" in verify["run"]
    assert '"$RUNNER_TEMP/cran-preflight-source-download"' in verify["run"]
    assert '"$GITHUB_OUTPUT"' in verify["run"]
    assert "find " not in verify["run"]
    assert "mapfile" not in verify["run"]

    aggregate = step_with_id(summary, "aggregate")
    assert "continue-on-error" not in aggregate
    assert "ci-summary" in aggregate["run"]
    assert aggregate["env"]["EXPECTED_TARBALL_SHA256"] == (
        "${{ steps.verify-source.outputs.tarball_sha256 }}"
    )
    assert "artifact-digest" not in text
    assert "subset" not in aggregate["run"].lower()
    assert "smoke" not in aggregate["run"].lower()


def test_summary_publication_and_explicit_outcome_gate_run_after_failures():
    _, workflow = load_workflow()
    summary = workflow["jobs"]["summary-gate"]
    publish = step_with_id(summary, "publish-summary")
    final_gate = step_with_id(summary, "final-gate")

    assert publish["if"] == "${{ always() }}"
    assert '"$GITHUB_STEP_SUMMARY"' in publish["run"]
    assert final_gate["if"] == "${{ always() }}"
    assert final_gate["env"] == {
        "PREPARE_RESULT": "${{ needs.prepare.result }}",
        "PRIMARY_LINUX_RESULT": "${{ needs.primary-linux.result }}",
        "PRIMARY_PLATFORM_RESULT": "${{ needs.primary-platform.result }}",
        "MKL_RESULT": "${{ needs.mkl.result }}",
        "NOSUGGESTS_RESULT": "${{ needs.nosuggests.result }}",
        "ACTIVE_RHUB_RESULT": "${{ needs.active-rhub.result }}",
        "APPLICABILITY_RESULT": "${{ needs.applicability-native.result }}",
        "REMAINING_DOCKER_RESULT": "${{ needs.remaining-docker.result }}",
        "LINUX_ARM64_RESULT": "${{ needs.linux-arm64.result }}",
        "LINUX_ARM64_COMPARE_RESULT": "${{ needs.linux-arm64-compare.result }}",
        "SOURCE_DOWNLOAD_OUTCOME": "${{ steps.download-source.outcome }}",
        "RESULT_DOWNLOAD_OUTCOME": "${{ steps.download-results.outcome }}",
        "SOURCE_VERIFY_OUTCOME": "${{ steps.verify-source.outcome }}",
        "AGGREGATE_OUTCOME": "${{ steps.aggregate.outcome }}",
    }
    for variable in final_gate["env"]:
        assert f'"${variable}"' in final_gate["run"]


def test_workflow_has_no_write_credentials_or_untrusted_context_in_shell():
    text, workflow = load_workflow()
    lowered = text.lower()

    for forbidden in (
        "pull_request_target",
        "secrets.",
        "github.token",
        "permissions: write",
    ):
        assert forbidden not in lowered

    for job in workflow["jobs"].values():
        for step in job["steps"]:
            if "run" in step:
                assert not re.search(
                    r"\$\{\{\s*(?:github|steps|needs)\.", step["run"]
                )
