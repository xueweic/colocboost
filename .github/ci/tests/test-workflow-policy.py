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
    assert set(workflow["jobs"]) == {"prepare", "summary-gate"}
    assert all("strategy" not in job for job in workflow["jobs"].values())


def test_every_job_is_inert_outside_the_fork_and_summary_always_runs():
    _, workflow = load_workflow()

    assert workflow["jobs"]["prepare"]["if"] == "${{ " + FORK_GUARD + " }}"
    summary = workflow["jobs"]["summary-gate"]
    assert summary["if"] == "${{ always() && " + FORK_GUARD + " }}"
    assert summary["needs"] == "prepare"


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
        }


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
