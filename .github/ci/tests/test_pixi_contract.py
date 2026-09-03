import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[3]
PIXI_MANIFEST = ROOT / "pixi.toml"
PIXI_LOCK = ROOT / "pixi.lock"
LEGACY_MANIFEST = ROOT / ".github" / "environment" / "pixi.toml"
LOCAL_CHECK_WRAPPER = ROOT / ".github" / "ci" / "run-local-check.R"
CONTRACT_TEST_WRAPPER = ROOT / ".github" / "ci" / "run-contract-tests.py"

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
    "pandoc",
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
R_TASK_ENV = {
    "LC_ALL": "C",
    "R_ENVIRON_USER": "/dev/null",
    "R_PROFILE_USER": "/dev/null",
    "R_LIBS_USER": "$CONDA_PREFIX/lib/R/library",
}
CI_UNIT_R_TASK_ENV = {
    "LC_ALL": "C",
    "R_ENVIRON_USER": "{% if pixi.is_win %}NUL{% else %}/dev/null{% endif %}",
    "R_PROFILE_USER": "{% if pixi.is_win %}NUL{% else %}/dev/null{% endif %}",
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
        "ci-prepare-source",
        "ci-verify-source",
        "ci-unit",
        "ci-check",
        "ci-summary",
    }

    assert task_args(tasks["ci-unit"]) == [
        "mode",
        "output",
        "r_executable",
        "runner_context",
    ]
    assert task_args(tasks["ci-prepare"]) == [
        "tarball",
        "metadata",
        "source_sha",
        "event_sha",
    ]
    assert task_args(tasks["ci-prepare-source"]) == [
        "repository",
        "output_dir",
        "event_sha",
        "github_output",
    ]
    assert task_args(tasks["ci-verify-source"]) == [
        "source_dir",
        "source_sha",
        "event_sha",
        "github_output",
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
    assert task_args(tasks["ci-summary"]) == [
        "manifest",
        "results_root",
        "expected_source_sha",
        "expected_event_sha",
        "expected_tarball_sha256",
        "summary",
    ]

    commands = {name: task_command(task) for name, task in tasks.items()}
    assert commands["ci-validate"] == (
        "python .github/ci/validate_manifest.py .github/ci/check-matrix.yml"
    )
    assert commands["ci-contract-tests"] == (
        "python -B .github/ci/run-contract-tests.py"
    )
    assert tasks["ci-contract-tests"].get("env") == {
        "PYTHONDONTWRITEBYTECODE": "1"
    }
    assert ".github/ci/artifact_contract.py create" in commands["ci-prepare"]
    assert commands["ci-prepare-source"].startswith(
        "python -B .github/ci/prepare_source.py"
    )
    assert tasks["ci-prepare-source"].get("env") == R_TASK_ENV
    assert commands["ci-verify-source"].startswith(
        "python -B .github/ci/verify_source.py"
    )
    assert ".github/ci/run_unit_driver.py" in commands["ci-unit"]
    assert commands["ci-unit"].startswith("python ")
    assert tasks["ci-unit"].get("env") == CI_UNIT_R_TASK_ENV
    assert ".github/ci/run_driver.py" in commands["ci-check"]
    assert commands["ci-check"].rstrip().endswith("--")
    assert ".github/ci/aggregate_results.py" in commands["ci-summary"]
    assert all(task.get("cwd") == "." for task in tasks.values())


def test_typed_task_values_are_shell_quoted(manifest):
    tasks = manifest["tasks"]
    required_quoted_values = {
        "ci-prepare": ("tarball", "metadata", "source_sha", "event_sha"),
        "ci-prepare-source": (
            "repository",
            "output_dir",
            "event_sha",
            "github_output",
        ),
        "ci-verify-source": (
            "source_dir",
            "source_sha",
            "event_sha",
            "github_output",
        ),
        "ci-unit": ("mode", "output", "r_executable", "runner_context"),
        "ci-check": (
            "environment_id",
            "driver",
            "executable",
            "tarball",
            "metadata",
            "source_sha",
            "event_sha",
            "evidence",
            "r_entrypoint",
        ),
        "ci-summary": (
            "manifest",
            "results_root",
            "expected_source_sha",
            "expected_event_sha",
            "expected_tarball_sha256",
            "summary",
        ),
    }
    quote_filter = " | replace(\"'\", \"'\\\"'\\\"'\")"
    for task_name, values in required_quoted_values.items():
        command = task_command(tasks[task_name])
        for value in values:
            assert f"'{{{{ {value}{quote_filter} }}}}'" in command, (
                task_name,
                value,
            )


def test_top_level_ci_unit_preserves_caller_r_library(tmp_path):
    pixi = shutil.which("pixi") or "/Users/xueweic/.pixi/bin/pixi"
    if not Path(pixi).is_file():
        pytest.skip("Pixi is not installed")

    caller_library = tmp_path / "caller R library"
    caller_library.mkdir()
    stub_bin = tmp_path / "external R bin"
    stub_bin.mkdir()
    if os.name == "nt":
        pytest.skip("the external-R runtime probe uses POSIX executable fixtures")
    r_binary = stub_bin / "R"
    r_binary.write_text("#!/bin/sh\nexit 99\n", encoding="utf-8")
    r_binary.chmod(0o755)
    rscript = stub_bin / "Rscript"
    rscript.write_text(
        "#!/bin/sh\nprintf 'R_LIBS_USER=%s\\n' \"$R_LIBS_USER\"\nexit 19\n",
        encoding="utf-8",
    )
    rscript.chmod(0o755)

    environment = os.environ.copy()
    environment["PATH"] = os.pathsep.join(
        [os.fspath(stub_bin), environment.get("PATH", "")]
    )
    environment["R_LIBS_USER"] = os.fspath(caller_library)
    completed = subprocess.run(
        [
            pixi,
            "run",
            "--locked",
            "--manifest-path",
            os.fspath(PIXI_MANIFEST),
            "ci-unit",
            "source",
            os.fspath(tmp_path / "unused.json"),
            os.fspath(r_binary),
            "external-r-probe",
        ],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )

    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert f"R_LIBS_USER={caller_library}" in output


@pytest.mark.parametrize("feature", ["local-r44", "local-r45"])
def test_local_r_tasks_keep_strict_policy_and_full_check(manifest, feature):
    tasks = manifest["feature"][feature]["tasks"]
    assert set(tasks) == {"test", "check", "ci-r-contract-tests"}

    test_command = task_command(tasks["test"])
    assert test_command.startswith("Rscript --vanilla ")
    assert ".github/ci/run-unit-tests.R" in test_command
    assert f"--context={feature}" in test_command
    assert "--policy=.github/ci/check-policy.yml" in test_command
    assert f".pixi/{feature}/unit-results.json" in test_command

    assert task_command(tasks["ci-r-contract-tests"]) == (
        "Rscript --vanilla .github/ci/tests/test-unit-policy.R"
    )
    assert all(task.get("env") == R_TASK_ENV for task in tasks.values())
    check_command = task_command(tasks["check"])
    assert check_command == "Rscript --vanilla .github/ci/run-local-check.R"
    assert "devtools::check" not in check_command
    assert "--no-manual" not in check_command
    assert "--no-build-vignettes" not in check_command
    assert "--ignore-vignettes" not in check_command


def test_local_check_wrapper_keeps_full_uncompromised_policy():
    assert LOCAL_CHECK_WRAPPER.is_file(), "the local check wrapper is missing"
    source = LOCAL_CHECK_WRAPPER.read_text(encoding="utf-8")
    assert "devtools::check" in source
    assert "--as-cran" in source
    assert "manual = TRUE" in source
    assert "build_args = character()" in source
    assert "--no-manual" not in source
    assert "--no-build-vignettes" not in source
    assert "--ignore-vignettes" not in source


def test_contract_runner_suppresses_pytest_and_bytecode_caches(tmp_path):
    assert CONTRACT_TEST_WRAPPER.is_file(), "the contract-test wrapper is missing"
    project = tmp_path / "contract project;literal$"
    tests = project / ".github" / "ci" / "tests"
    tests.mkdir(parents=True)
    wrapper = project / ".github" / "ci" / CONTRACT_TEST_WRAPPER.name
    shutil.copy2(CONTRACT_TEST_WRAPPER, wrapper)
    (tests / "test_sample.py").write_text(
        "def test_sample():\n    assert True\n",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    completed = subprocess.run(
        [sys.executable, "-B", os.fspath(wrapper)],
        cwd=project,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "1 passed" in completed.stdout
    assert not list(project.rglob(".pytest_cache"))
    assert not list(project.rglob("__pycache__"))
    assert not list(project.rglob("*.pyc"))


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


def test_ci_prepare_preserves_practical_shell_metacharacters_in_paths(tmp_path):
    directory = tmp_path / "paths with spaces;dollar$HOME&ampersand'apostrophe"
    directory.mkdir()
    tarball = directory / "package $HOME;literal&'quote.tar.gz"
    metadata = directory / "metadata $HOME;literal&'quote.json"
    tarball.write_bytes(b"task-6 path quoting fixture")

    rendered = run_pixi_dry(
        "ci-prepare",
        os.fspath(tarball),
        os.fspath(metadata),
        "a" * 40,
        "b" * 40,
    )
    assert rendered.returncode == 0, rendered.stdout + rendered.stderr
    task_line = next(
        line
        for line in (rendered.stdout + rendered.stderr).splitlines()
        if ".github/ci/artifact_contract.py create" in line
    )
    command = task_line.partition(": ")[2]
    argv = shlex.split(command, posix=True)
    assert argv[0] == "python"
    argv[0] = sys.executable

    completed = subprocess.run(
        argv,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    document = json.loads(metadata.read_text(encoding="utf-8"))
    assert document["filename"] == tarball.name


def test_pixi_renders_typed_task_arguments_and_passthrough():
    unit = run_pixi_dry(
        "ci-unit",
        "source",
        "unit.json",
        "/opt/R path's/bin/R",
        "r-cmd-check-installed",
    )
    assert unit.returncode == 0, unit.stderr
    unit_output = unit.stdout + unit.stderr
    assert "--load-package='source'" in unit_output
    assert "--output='unit.json'" in unit_output
    assert "--r-executable='/opt/R path'\"'\"'s/bin/R'" in unit_output
    assert "--context='r-cmd-check-installed'" in unit_output

    prepare = run_pixi_dry(
        "ci-prepare",
        "package.tar.gz",
        "metadata.json",
        "a" * 40,
        "b" * 40,
    )
    assert prepare.returncode == 0, prepare.stderr
    prepare_output = prepare.stdout + prepare.stderr
    assert "--tarball='package.tar.gz'" in prepare_output
    assert "--metadata='metadata.json'" in prepare_output

    prepare_source = run_pixi_dry(
        "ci-prepare-source",
        "repository's path",
        "output directory's path",
        "a" * 40,
        "github output's path",
    )
    assert prepare_source.returncode == 0, prepare_source.stderr
    prepare_source_output = prepare_source.stdout + prepare_source.stderr
    assert "--repository='repository'\"'\"'s path'" in prepare_source_output
    assert "--output-dir='output directory'\"'\"'s path'" in prepare_source_output
    assert "--github-output='github output'\"'\"'s path'" in prepare_source_output

    verify_source = run_pixi_dry(
        "ci-verify-source",
        "downloaded source's path",
        "a" * 40,
        "b" * 40,
        "github output's path",
    )
    assert verify_source.returncode == 0, verify_source.stderr
    verify_source_output = verify_source.stdout + verify_source.stderr
    assert "--source-dir='downloaded source'\"'\"'s path'" in verify_source_output
    assert "--github-output='github output'\"'\"'s path'" in verify_source_output

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
    assert "--r-entrypoint='r-cmd'" in check_output
    assert "-- check --as-cran {tarball}" in check_output

    summary = run_pixi_dry(
        "ci-summary",
        "manifest's.yml",
        "results root's",
        "a" * 40,
        "b" * 40,
        "c" * 64,
        "summary's.md",
    )
    assert summary.returncode == 0, summary.stderr
    summary_output = summary.stdout + summary.stderr
    assert "--manifest='manifest'\"'\"'s.yml'" in summary_output
    assert "--results-root='results root'\"'\"'s'" in summary_output
    assert "--summary='summary'\"'\"'s.md'" in summary_output

    missing_argument = run_pixi_dry(
        "ci-unit", "source", "unit.json", "/opt/R/bin/R"
    )
    assert missing_argument.returncode != 0
    assert "runner_context" in missing_argument.stderr
