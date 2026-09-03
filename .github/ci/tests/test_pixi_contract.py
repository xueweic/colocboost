import hashlib
import os
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[3]
PIXI_MANIFEST = ROOT / "pixi.toml"
PIXI_LOCK = ROOT / "pixi.lock"
LEGACY_MANIFEST = ROOT / ".github" / "environment" / "pixi.toml"

HELPER_PLATFORMS = {
    "linux-64",
    "linux-aarch64",
    "osx-64",
    "osx-arm64",
    "win-64",
}
UNIX_PLATFORMS = HELPER_PLATFORMS - {"win-64"}
HELPER_DEPENDENCIES = {
    "python",
    "pyyaml",
    "jsonschema",
    "pytest",
}
R_DEPENDENCIES = {
    "r-base",
    "r-testthat",
    "r-devtools",
    "r-rcmdcheck",
    "r-jsonlite",
    "r-rfast",
    "r-matrixstats",
    "r-knitr",
    "r-rmarkdown",
    "r-ashr",
    "r-mass",
    "r-susier",
    "r-yaml",
}
LEGACY_BODY_SHA256 = (
    "4137647ef8016144e2b98e4c94b529ae9c20133a92cdfec9254c803568e7ab61"
)


@pytest.fixture(scope="module")
def manifest():
    assert PIXI_MANIFEST.is_file(), "the canonical root pixi.toml is missing"
    return tomllib.loads(PIXI_MANIFEST.read_text(encoding="utf-8"))


def task_command(task):
    return task if isinstance(task, str) else task["cmd"]


def task_args(task):
    return task.get("args", []) if isinstance(task, dict) else []


def strip_leading_comments(text):
    lines = text.splitlines(keepends=True)
    while lines and (not lines[0].strip() or lines[0].startswith("#")):
        lines.pop(0)
    return "".join(lines)


def test_root_manifest_declares_five_platform_helper(manifest):
    assert "project" not in manifest
    assert manifest["workspace"]["name"] == "colocboost-ci"
    assert set(manifest["workspace"]["platforms"]) == HELPER_PLATFORMS
    assert set(manifest["dependencies"]) == HELPER_DEPENDENCIES
    assert manifest["dependencies"]["python"] == "3.12.*"


def test_unix_helpers_are_not_exposed_on_windows(manifest):
    assert "jq" not in manifest["dependencies"]
    assert "coreutils" not in manifest["dependencies"]
    targets = manifest["target"]
    assert "win-64" not in targets
    for platform in UNIX_PLATFORMS:
        assert set(targets[platform]["dependencies"]) == {"jq", "coreutils"}


@pytest.mark.parametrize(
    ("feature", "r_version"),
    [("local-r44", "4.4.*"), ("local-r45", "4.5.*")],
)
def test_local_r_features_are_complete_and_unix_only(manifest, feature, r_version):
    feature_config = manifest["feature"][feature]
    assert set(feature_config["platforms"]) == UNIX_PLATFORMS
    assert set(feature_config["dependencies"]) == R_DEPENDENCIES
    assert feature_config["dependencies"]["r-base"] == r_version
    environment = manifest["environments"][feature]
    assert environment["features"] == [feature]
    assert not environment.get("no-default-feature", False)


def test_top_level_task_interfaces_are_explicit(manifest):
    tasks = manifest["tasks"]
    assert set(tasks) == {
        "ci-validate",
        "ci-contract-tests",
        "ci-prepare",
        "ci-unit",
        "ci-check",
        "ci-summary",
    }

    assert task_args(tasks["ci-unit"]) == ["mode", "output", "environment_id"]
    assert task_args(tasks["ci-prepare"]) == [
        "tarball",
        "metadata",
        "source_sha",
        "event_sha",
    ]
    assert task_args(tasks["ci-check"]) == [
        "environment_id",
        "driver",
        "executable",
        "tarball",
        "metadata",
        "source_sha",
        "event_sha",
        "evidence",
        {"arg": "r_entrypoint", "default": ""},
    ]

    commands = {name: task_command(task) for name, task in tasks.items()}
    assert commands["ci-validate"] == (
        "python .github/ci/validate_manifest.py .github/ci/check-matrix.yml"
    )
    assert commands["ci-contract-tests"] == "pytest -q .github/ci/tests"
    assert ".github/ci/artifact_contract.py create" in commands["ci-prepare"]
    assert ".github/ci/run-unit-tests.R" in commands["ci-unit"]
    assert ".github/ci/run_driver.py" in commands["ci-check"]
    assert commands["ci-check"].rstrip().endswith("--")
    assert commands["ci-summary"] == "python .github/ci/aggregate_results.py"
    assert all(task.get("cwd") == "." for task in tasks.values())


@pytest.mark.parametrize("feature", ["local-r44", "local-r45"])
def test_local_r_tasks_keep_strict_policy_and_full_check(manifest, feature):
    tasks = manifest["feature"][feature]["tasks"]
    assert set(tasks) == {"test", "check", "ci-r-contract-tests"}

    test_command = task_command(tasks["test"])
    assert ".github/ci/run-unit-tests.R" in test_command
    assert f"--context={feature}" in test_command
    assert "--policy=.github/ci/check-policy.yml" in test_command
    assert f".pixi/{feature}/unit-results.json" in test_command

    assert task_command(tasks["ci-r-contract-tests"]) == (
        "Rscript .github/ci/tests/test-unit-policy.R"
    )
    check_command = task_command(tasks["check"])
    assert tasks["check"].get("env") == {"LC_ALL": "C"}
    assert check_command.startswith("Rscript --vanilla ")
    assert "devtools::check" in check_command
    assert "manual = TRUE" in check_command
    assert "build_args = character()" in check_command
    assert "--no-manual" not in check_command
    assert "--no-build-vignettes" not in check_command
    assert "--ignore-vignettes" not in check_command


def test_lock_covers_every_declared_environment_platform(manifest):
    assert PIXI_LOCK.is_file(), "pixi.lock is missing"
    lock = yaml.safe_load(PIXI_LOCK.read_text(encoding="utf-8"))
    locked_environments = lock["environments"]
    assert set(locked_environments["default"]["packages"]) == HELPER_PLATFORMS
    for environment in ("local-r44", "local-r45"):
        assert set(locked_environments[environment]["packages"]) == UNIX_PLATFORMS


def test_legacy_manifest_body_is_byte_for_byte_unchanged():
    text = LEGACY_MANIFEST.read_text(encoding="utf-8")
    leading_comments = text[: -len(strip_leading_comments(text))]
    assert "legacy" in leading_comments.lower()
    assert "release/pkgdown" in leading_comments.lower()
    assert "migrated and verified" in leading_comments.lower()
    body = strip_leading_comments(text).encode()
    assert hashlib.sha256(body).hexdigest() == LEGACY_BODY_SHA256


def run_pixi_dry(task, *arguments):
    pixi = shutil.which("pixi") or "/Users/xueweic/.pixi/bin/pixi"
    if not Path(pixi).is_file():
        pytest.skip("Pixi is not installed")
    return subprocess.run(
        [
            pixi,
            "run",
            "--locked",
            "--dry-run",
            "--manifest-path",
            os.fspath(PIXI_MANIFEST),
            task,
            *arguments,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_pixi_renders_typed_task_arguments_and_passthrough():
    unit = run_pixi_dry("ci-unit", "source", "unit.json", "local-r44")
    assert unit.returncode == 0, unit.stderr
    unit_output = unit.stdout + unit.stderr
    assert "--load-package=source" in unit_output
    assert "--output=unit.json" in unit_output
    assert "--context=local-r44" in unit_output

    prepare = run_pixi_dry(
        "ci-prepare",
        "package.tar.gz",
        "metadata.json",
        "a" * 40,
        "b" * 40,
    )
    assert prepare.returncode == 0, prepare.stderr
    prepare_output = prepare.stdout + prepare.stderr
    assert "--tarball=package.tar.gz" in prepare_output
    assert "--metadata=metadata.json" in prepare_output

    check = run_pixi_dry(
        "ci-check",
        "r-release-linux-x86-64",
        "r-binary",
        "/opt/R/bin/R",
        "package.tar.gz",
        "metadata.json",
        "a" * 40,
        "b" * 40,
        "evidence.json",
        "r-cmd",
        "--",
        "check",
        "--as-cran",
        "{tarball}",
    )
    assert check.returncode == 0, check.stderr
    check_output = check.stdout + check.stderr
    assert "--r-entrypoint=r-cmd" in check_output
    assert "-- check --as-cran {tarball}" in check_output

    missing_argument = run_pixi_dry("ci-unit", "source", "unit.json")
    assert missing_argument.returncode != 0
    assert "environment_id" in missing_argument.stderr
