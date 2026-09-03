import hashlib
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / ".github" / "ci" / "prepare_source.py"
SOURCE_SHA = "a" * 40


def run_git(repository, *arguments, check=True):
    return subprocess.run(
        ["git", "-C", os.fspath(repository), *arguments],
        check=check,
        capture_output=True,
        text=True,
    )


def create_repository(tmp_path):
    repository = tmp_path / "repository"
    repository.mkdir()
    run_git(repository, "init")
    run_git(repository, "config", "user.name", "CI Contract")
    run_git(repository, "config", "user.email", "ci@example.invalid")
    (repository / "DESCRIPTION").write_text(
        "Package: colocboost\nVersion: 1.2.3\n", encoding="utf-8"
    )
    (repository / "committed-marker.txt").write_text("committed\n", encoding="utf-8")
    run_git(repository, "add", "DESCRIPTION", "committed-marker.txt")
    run_git(repository, "commit", "-m", "fixture")
    head = run_git(repository, "rev-parse", "HEAD").stdout.strip()
    return repository, head


def create_fake_r(tmp_path):
    prefix = tmp_path / "pixi-prefix"
    bin_directory = prefix / "bin"
    r_home = prefix / "lib" / "R"
    bin_directory.mkdir(parents=True)
    r_home.mkdir(parents=True)
    executable = bin_directory / "R"
    pandoc = bin_directory / "pandoc"
    pandoc.write_text("#!/bin/sh\nprintf 'pandoc 3.0\\n'\n", encoding="utf-8")
    pandoc.chmod(0o755)
    executable.write_text(
        "\n".join(
            [
                "#!/bin/sh",
                'if [ "$1" = "RHOME" ]; then',
                f"  printf '%s\\n' {shlex.quote(os.fspath(r_home))}",
                "  exit 0",
                "fi",
                'if [ "$1" != "--vanilla" ] || [ "$2" != "CMD" ] || [ "$3" != "build" ]; then',
                "  exit 81",
                "fi",
                f'test "$RSTUDIO_PANDOC" = {shlex.quote(os.fspath(bin_directory))} || exit 85',
                'source_directory="$4"',
                'test ! -e "$source_directory/.git" || exit 82',
                'test "$(cat "$source_directory/committed-marker.txt")" = "committed" || exit 83',
                'case "${FAKE_R_BUILD_MODE:-one}" in',
                "  one)",
                "    printf 'source tarball bytes' > colocboost_1.2.3.tar.gz",
                "    ;;",
                "  zero)",
                "    ;;",
                "  two)",
                "    printf one > first.tar.gz",
                "    printf two > second.tar.gz",
                "    ;;",
                "  symlink)",
                "    printf payload > payload.bin",
                "    ln -s payload.bin linked.tar.gz",
                "    ;;",
                "  *) exit 84 ;;",
                "esac",
                "exit 0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return prefix


def run_prepare(
    repository,
    output_directory,
    event_sha,
    github_output,
    prefix,
    *,
    mode="one",
):
    environment = os.environ.copy()
    environment["CONDA_PREFIX"] = os.fspath(prefix)
    environment["PATH"] = os.pathsep.join(
        [os.fspath(prefix / "bin"), environment.get("PATH", "")]
    )
    environment["FAKE_R_BUILD_MODE"] = mode
    return subprocess.run(
        [
            sys.executable,
            os.fspath(SCRIPT),
            "--repository",
            os.fspath(repository),
            "--output-dir",
            os.fspath(output_directory),
            "--event-sha",
            event_sha,
            "--github-output",
            os.fspath(github_output),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=15,
    )


def parse_github_outputs(path):
    pairs = [line.split("=", 1) for line in path.read_text(encoding="utf-8").splitlines()]
    return dict(pairs)


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
def test_prepare_builds_once_from_clean_exact_head_and_emits_verified_outputs(tmp_path):
    repository, head = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    output_directory = tmp_path / "source output"
    output_directory.mkdir()
    github_output = tmp_path / "github output.txt"
    github_output.write_text("existing=value\n", encoding="utf-8")

    completed = run_prepare(
        repository, output_directory, head, github_output, prefix
    )

    assert completed.returncode == 0, completed.stderr
    tarball = output_directory / "colocboost_1.2.3.tar.gz"
    metadata_path = output_directory / "source-metadata.json"
    assert sorted(path.name for path in output_directory.iterdir()) == [
        "colocboost_1.2.3.tar.gz",
        "source-metadata.json",
    ]
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata == {
        "source_sha": head,
        "event_sha": head,
        "filename": tarball.name,
        "size": len(b"source tarball bytes"),
        "sha256": hashlib.sha256(b"source tarball bytes").hexdigest(),
    }
    outputs = parse_github_outputs(github_output)
    assert outputs == {
        "existing": "value",
        "source_sha": head,
        "event_sha": head,
        "tarball_filename": tarball.name,
        "tarball_path": os.fspath(tarball.resolve()),
        "metadata_path": os.fspath(metadata_path.resolve()),
        "tarball_size": str(metadata["size"]),
        "tarball_sha256": metadata["sha256"],
    }


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
@pytest.mark.parametrize("dirty_kind", ["tracked", "untracked"])
def test_prepare_rejects_any_dirty_repository(tmp_path, dirty_kind):
    repository, head = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    if dirty_kind == "tracked":
        (repository / "committed-marker.txt").write_text("dirty\n", encoding="utf-8")
    else:
        (repository / "untracked.txt").write_text("dirty\n", encoding="utf-8")
    output_directory = tmp_path / "source-output"
    output_directory.mkdir()

    completed = run_prepare(
        repository,
        output_directory,
        head,
        tmp_path / "github-output",
        prefix,
    )

    assert completed.returncode != 0
    assert "clean" in completed.stderr
    assert not list(output_directory.iterdir())


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
def test_prepare_rejects_event_sha_that_is_not_exact_head(tmp_path):
    repository, _ = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    output_directory = tmp_path / "source-output"
    output_directory.mkdir()

    completed = run_prepare(
        repository,
        output_directory,
        "f" * 40,
        tmp_path / "github-output",
        prefix,
    )

    assert completed.returncode != 0
    assert "HEAD" in completed.stderr
    assert not list(output_directory.iterdir())


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
def test_prepare_requires_an_existing_empty_output_directory_outside_repo(tmp_path):
    repository, head = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    github_output = tmp_path / "github-output"

    missing = tmp_path / "missing"
    missing_result = run_prepare(repository, missing, head, github_output, prefix)
    assert missing_result.returncode != 0
    assert "output" in missing_result.stderr

    internal = repository / "output"
    internal.mkdir()
    internal_result = run_prepare(repository, internal, head, github_output, prefix)
    assert internal_result.returncode != 0
    assert "outside" in internal_result.stderr

    nonempty = tmp_path / "nonempty"
    nonempty.mkdir()
    (nonempty / "keep.txt").write_text("keep", encoding="utf-8")
    nonempty_result = run_prepare(repository, nonempty, head, github_output, prefix)
    assert nonempty_result.returncode != 0
    assert "empty" in nonempty_result.stderr
    assert (nonempty / "keep.txt").read_text(encoding="utf-8") == "keep"


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
@pytest.mark.parametrize("mode", ["zero", "two", "symlink"])
def test_prepare_requires_exactly_one_regular_non_symlink_tarball(tmp_path, mode):
    repository, head = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    output_directory = tmp_path / "source-output"
    output_directory.mkdir()

    completed = run_prepare(
        repository,
        output_directory,
        head,
        tmp_path / "github-output",
        prefix,
        mode=mode,
    )

    assert completed.returncode != 0
    assert "one regular non-symlink" in completed.stderr
    assert not list(output_directory.iterdir())


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
def test_prepare_rejects_r_outside_the_declared_pixi_prefix(tmp_path):
    repository, head = create_repository(tmp_path)
    real_prefix = create_fake_r(tmp_path)
    declared_prefix = tmp_path / "different-prefix"
    declared_prefix.mkdir()
    output_directory = tmp_path / "source-output"
    output_directory.mkdir()
    environment = os.environ.copy()
    environment["CONDA_PREFIX"] = os.fspath(declared_prefix)
    environment["PATH"] = os.pathsep.join(
        [os.fspath(real_prefix / "bin"), environment.get("PATH", "")]
    )

    completed = subprocess.run(
        [
            sys.executable,
            os.fspath(SCRIPT),
            "--repository",
            os.fspath(repository),
            "--output-dir",
            os.fspath(output_directory),
            "--event-sha",
            head,
            "--github-output",
            os.fspath(tmp_path / "github-output"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=15,
    )

    assert completed.returncode != 0
    assert "Pixi" in completed.stderr
    assert not list(output_directory.iterdir())


@pytest.mark.skipif(os.name == "nt", reason="prepare source runs on Ubuntu")
def test_prepare_rejects_pandoc_outside_the_declared_pixi_prefix(tmp_path):
    repository, head = create_repository(tmp_path)
    prefix = create_fake_r(tmp_path)
    (prefix / "bin" / "pandoc").unlink()
    outside = tmp_path / "outside-bin"
    outside.mkdir()
    outside_pandoc = outside / "pandoc"
    outside_pandoc.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    outside_pandoc.chmod(0o755)
    output_directory = tmp_path / "source-output"
    output_directory.mkdir()
    environment = os.environ.copy()
    environment["CONDA_PREFIX"] = os.fspath(prefix)
    environment["PATH"] = os.pathsep.join(
        [os.fspath(prefix / "bin"), os.fspath(outside), environment.get("PATH", "")]
    )

    completed = subprocess.run(
        [
            sys.executable,
            os.fspath(SCRIPT),
            "--repository",
            os.fspath(repository),
            "--output-dir",
            os.fspath(output_directory),
            "--event-sha",
            head,
            "--github-output",
            os.fspath(tmp_path / "github-output"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=15,
    )

    assert completed.returncode != 0
    assert "Pandoc" in completed.stderr
    assert not list(output_directory.iterdir())
