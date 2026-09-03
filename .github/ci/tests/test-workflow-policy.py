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
ALLOWED_ACTIONS = {CHECKOUT, SETUP_PIXI, UPLOAD, DOWNLOAD}
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
    assert set(workflow["jobs"]) == {"prepare", "mkl", "nosuggests", "summary-gate"}
    assert all("strategy" not in job for job in workflow["jobs"].values())


def test_every_job_is_inert_outside_the_fork_and_summary_always_runs():
    _, workflow = load_workflow()

    assert workflow["jobs"]["prepare"]["if"] == "${{ " + FORK_GUARD + " }}"
    for name in ("mkl", "nosuggests"):
        assert workflow["jobs"][name]["if"] == (
            "${{ always() && " + FORK_GUARD + " }}"
        )
        assert workflow["jobs"][name]["needs"] == "prepare"
    summary = workflow["jobs"]["summary-gate"]
    assert summary["if"] == "${{ always() && " + FORK_GUARD + " }}"
    assert summary["needs"] == ["prepare", "mkl", "nosuggests"]


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
        "MKL_RESULT": "${{ needs.mkl.result }}",
        "NOSUGGESTS_RESULT": "${{ needs.nosuggests.result }}",
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
