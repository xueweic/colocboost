import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(CI_DIR))

import verify_platform_r as module  # noqa: E402
from verify_platform_r import verify_platform_r  # noqa: E402


def make_r(tmp_path, name=None):
    executable = tmp_path / "setup R with spaces" / "bin" / (name or "R")
    executable.parent.mkdir(parents=True)
    executable.write_bytes(b"fake R executable\n")
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def row(r_binary, **overrides):
    document = {
        "id": "platform-row",
        "driver": "r-binary",
        "runner": "ubuntu-24.04",
        "setup_r_selector": "devel",
        "system_r": os.fspath(r_binary).replace("\\", "/"),
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
    }
    document.update(overrides)
    return document


def identity_output(**overrides):
    fields = {
        "version": "4.7.0",
        "version_string": "R Under development (unstable) (2026-09-01 r99999)",
        "status": "Under development (unstable)",
        "svn_revision": "99999",
        "platform": "x86_64-pc-linux-gnu",
        "os": "Linux",
        "os_release": "6.11.0",
        "arch": "x86_64",
        "r_home": "/opt/R/devel/lib/R",
        "compiler": "gcc 16.1.0",
        "cc": "gcc-16",
        "cc_path": "/usr/bin/gcc-16",
        "cc_version": "gcc (GCC) 16.1.0",
        "cxx": "g++-16",
        "fc": "gfortran-16",
        "locale": "LC_CTYPE=C;LC_COLLATE=C",
    }
    fields.update(overrides)
    return "".join(f"{key}={value}\n" for key, value in fields.items())


def completed(stdout, returncode=0, stderr=""):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


def test_verifies_manifest_setup_r_selection_and_actual_identity(tmp_path, monkeypatch):
    r_binary = make_r(tmp_path)
    evidence = tmp_path / "proof" / "platform-r.json"
    github_output = tmp_path / "github-output.txt"
    github_output.write_text("prior=value\n", encoding="utf-8")
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        home = os.fspath(r_binary.parent.parent)
        return completed(identity_output(r_home=home))

    monkeypatch.setattr(subprocess, "run", fake_run)
    proof = verify_platform_r(
        row(r_binary), environment_id="platform-row", r_executable=r_binary,
        selected_selector="devel", setup_r_version="4.7.0",
        runner_label="ubuntu-24.04", evidence_path=evidence,
        github_output=github_output,
    )

    assert calls[0][0][0] == os.fspath(r_binary)
    assert 'system2(cc_path, "--version"' in calls[0][0][-1]
    assert calls[0][1]["shell"] is False
    assert proof["status"] == "pass"
    assert proof["expected_r_kind"] == "devel"
    assert proof["identity"]["svn_revision"] == "99999"
    assert proof["identity"]["cc"] == "gcc-16"
    assert json.loads(evidence.read_text(encoding="utf-8")) == proof
    assert "r_executable=" in github_output.read_text(encoding="utf-8")


def test_accepts_official_devel_revision_embedded_in_version_string(
    tmp_path, monkeypatch
):
    r_binary = make_r(tmp_path)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: completed(
            identity_output(
                r_home=os.fspath(r_binary.parent.parent),
                svn_revision="<none>",
                version_string="R Under development (unstable) (2026-09-01 r99999)",
            )
        ),
    )
    proof = verify_platform_r(
        row(r_binary), environment_id="platform-row", r_executable=r_binary,
        selected_selector="devel", setup_r_version="4.7.0",
        runner_label="ubuntu-24.04",
    )
    assert proof["identity"]["version_string"].endswith("r99999)")


def test_macos_framework_parent_symlink_binds_resolved_r_home(tmp_path, monkeypatch):
    versioned_home = tmp_path / "Framework" / "Versions" / "4.6" / "Resources"
    real_r = versioned_home / "bin" / "R"
    real_r.parent.mkdir(parents=True)
    real_r.write_bytes(b"fake R executable\n")
    real_r.chmod(real_r.stat().st_mode | stat.S_IXUSR)
    lexical_home = tmp_path / "Framework" / "Resources"
    try:
        lexical_home.symlink_to(versioned_home, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks unavailable")
    lexical_r = lexical_home / "bin" / "R"
    manifest_row = row(
        lexical_r,
        runner="macos-15",
        setup_r_selector="release",
        expected_r_kind="release",
        expected_os="macos",
        expected_architecture="aarch64",
    )
    monkeypatch.setattr(module, "_MACOS_R", os.fspath(lexical_r))
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: completed(
            identity_output(
                version="4.6.2",
                version_string="R version 4.6.2 (2026-06-12)",
                status="<none>",
                svn_revision="<none>",
                os="Darwin",
                arch="aarch64",
                platform="aarch64-apple-darwin20",
                r_home=os.fspath(versioned_home),
                compiler="<none>",
                cc="clang",
                cc_path="/usr/bin/clang",
                cc_version="Apple clang version 17.0.0",
                cxx="clang++",
                fc="gfortran",
            )
        ),
    )
    proof = verify_platform_r(
        manifest_row, environment_id="platform-row", r_executable=lexical_r,
        selected_selector="release", setup_r_version="4.6.2",
        runner_label="macos-15",
    )
    assert proof["r_executable"]["path"] == os.fspath(lexical_r)
    assert proof["r_executable"]["resolved_path"] == os.fspath(real_r)

    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: completed(
            identity_output(
                version="4.6.2",
                version_string="R version 4.6.2 (2026-06-12)",
                status="<none>",
                svn_revision="<none>",
                os="Darwin",
                arch="aarch64",
                platform="aarch64-apple-darwin20",
                r_home=os.fspath(versioned_home),
                compiler="gcc (GCC) 14.2.0",
                cc="gcc",
                cc_path="/usr/bin/gcc",
                cc_version="gcc (GCC) 14.2.0",
                cxx="g++",
                fc="gfortran",
            )
        ),
    )
    with pytest.raises(ValueError, match="Apple/Clang"):
        verify_platform_r(
            manifest_row, environment_id="platform-row", r_executable=lexical_r,
            selected_selector="release", setup_r_version="4.6.2",
            runner_label="macos-15",
        )


def test_windows_identity_requires_actual_rtools_compiler_path(tmp_path, monkeypatch):
    r_binary = make_r(tmp_path, "R.exe")
    manifest_row = row(
        r_binary,
        runner="windows-2022",
        setup_r_selector="release",
        expected_r_kind="release",
        expected_os="windows",
        expected_architecture="x86_64",
    )
    monkeypatch.setattr(
        module, "_WINDOWS_R", os.fspath(r_binary).replace("\\", "/")
    )

    def use_compiler(path):
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *args, **kwargs: completed(
                identity_output(
                    version="4.6.2",
                    version_string="R version 4.6.2 (2026-06-12)",
                    status="<none>",
                    svn_revision="<none>",
                    os="Windows",
                    arch="x86_64",
                    platform="x86_64-w64-mingw32",
                    r_home=os.fspath(r_binary.parent.parent),
                    compiler="gcc.exe (GCC) 14.2.0",
                    cc="gcc",
                    cc_path=path,
                    cc_version="gcc.exe (GCC) 14.2.0",
                    cxx="g++",
                    fc="gfortran",
                )
            ),
        )

    use_compiler("C:/rtools46/x86_64-w64-mingw32.static.posix/bin/gcc.exe")
    proof = verify_platform_r(
        manifest_row, environment_id="platform-row", r_executable=r_binary,
        selected_selector="release", setup_r_version="4.6.2",
        runner_label="windows-2022",
    )
    assert "rtools46" in proof["identity"]["cc_path"]

    use_compiler("C:/unrelated/bin/gcc.exe")
    with pytest.raises(ValueError, match="Rtools"):
        verify_platform_r(
            manifest_row, environment_id="platform-row", r_executable=r_binary,
            selected_selector="release", setup_r_version="4.6.2",
            runner_label="windows-2022",
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"selected_selector": "release"}, "selector"),
        ({"setup_r_version": "4.6.0"}, "setup-r.*version"),
        ({"runner_label": "windows-2022"}, "runner"),
    ],
)
def test_rejects_caller_values_that_disagree_with_manifest(
    tmp_path, monkeypatch, mutation, message
):
    r_binary = make_r(tmp_path)
    monkeypatch.setattr(
        subprocess, "run",
        lambda *args, **kwargs: completed(
            identity_output(r_home=os.fspath(r_binary.parent.parent))
        ),
    )
    arguments = {
        "manifest_row": row(r_binary), "environment_id": "platform-row",
        "r_executable": r_binary, "selected_selector": "devel",
        "setup_r_version": "4.7.0", "runner_label": "ubuntu-24.04",
    }
    arguments.update(mutation)
    with pytest.raises(ValueError, match=message):
        verify_platform_r(**arguments)


@pytest.mark.parametrize(
    ("identity_mutation", "message"),
    [
        ({"status": "<none>", "version_string": "R version 4.7.0"}, "release kind"),
        ({"os": "Windows"}, "OS or architecture"),
        ({"arch": "aarch64"}, "OS or architecture"),
        ({"os_release": "<none>"}, "OS release"),
        ({"locale": "<none>"}, "locale"),
        ({"cc": "<none>"}, "compiler"),
        ({"cc_path": "<none>"}, "compiler executable"),
        ({"cc_version": "<none>"}, "compiler version"),
        (
            {
                "svn_revision": "<none>",
                "version_string": "R Under development (unstable) (2026-09-01)",
            },
            "development revision",
        ),
        (
            {
                "svn_revision": "not-a-revision",
                "version_string": "R Under development (unstable)",
            },
            "development revision",
        ),
        ({"r_home": "/wrong/R"}, "R.home"),
    ],
)
def test_rejects_wrong_or_incomplete_actual_r_identity(
    tmp_path, monkeypatch, identity_mutation, message
):
    r_binary = make_r(tmp_path)
    values = {"r_home": os.fspath(r_binary.parent.parent), **identity_mutation}
    monkeypatch.setattr(
        subprocess, "run", lambda *args, **kwargs: completed(identity_output(**values))
    )
    with pytest.raises(ValueError, match=message):
        verify_platform_r(
            row(r_binary), environment_id="platform-row", r_executable=r_binary,
            selected_selector="devel", setup_r_version="4.7.0",
            runner_label="ubuntu-24.04",
        )


@pytest.mark.parametrize(
    "stdout",
    [
        identity_output(version="4.7.0\ninjected"),
        identity_output().replace("arch=x86_64\n", ""),
        identity_output() + "arch=x86_64\n",
        "not=the protocol\n",
    ],
)
def test_rejects_ambiguous_or_open_identity_protocol(tmp_path, monkeypatch, stdout):
    r_binary = make_r(tmp_path)
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: completed(stdout))
    with pytest.raises(ValueError, match="ambiguous|identity"):
        verify_platform_r(
            row(r_binary), environment_id="platform-row", r_executable=r_binary,
            selected_selector="devel", setup_r_version="4.7.0",
            runner_label="ubuntu-24.04",
        )


def test_rejects_symlink_and_manifest_path_substitution(tmp_path, monkeypatch):
    r_binary = make_r(tmp_path)
    linked = tmp_path / "linked-R"
    try:
        linked.symlink_to(r_binary)
    except OSError:
        pytest.skip("symlinks unavailable")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: completed(identity_output()))
    with pytest.raises(ValueError, match="R executable"):
        verify_platform_r(
            row(linked), environment_id="platform-row", r_executable=linked,
            selected_selector="devel", setup_r_version="4.7.0",
            runner_label="ubuntu-24.04",
        )
    with pytest.raises(ValueError, match="manifest system_r"):
        verify_platform_r(
            row(r_binary, system_r="/different/R"), environment_id="platform-row",
            r_executable=r_binary, selected_selector="devel",
            setup_r_version="4.7.0", runner_label="ubuntu-24.04",
        )


def test_committed_windows_and_macos_paths_are_canonical():
    manifest = module.validate_manifest(module.load_manifest(CI_DIR / "check-matrix.yml"))
    rows = {entry["id"]: entry for entry in manifest["coverage"]}
    windows = [row for row in rows.values() if row.get("expected_os") == "windows"]
    macos = [row for row in rows.values() if row.get("expected_os") == "macos"]
    assert windows and all(row["system_r"] == "C:/R/bin/R.exe" for row in windows)
    assert macos and all(
        row["system_r"] == "/Library/Frameworks/R.framework/Resources/bin/R"
        for row in macos
    )
