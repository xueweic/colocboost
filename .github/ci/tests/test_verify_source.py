import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
CI_DIR = ROOT / ".github" / "ci"
SCRIPT = CI_DIR / "verify_source.py"
SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def write_contract(directory, *, filename="colocboost_1.0.0.tar.gz", content=b"raw"):
    tarball = directory / filename
    tarball.write_bytes(content)
    metadata = directory / "source-metadata.json"
    metadata.write_text(
        json.dumps(
            {
                "source_sha": SOURCE_SHA,
                "event_sha": EVENT_SHA,
                "filename": filename,
                "size": len(content),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
        ),
        encoding="utf-8",
    )
    return tarball, metadata


def run_verify(source_dir, github_output, *, timeout=5):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--source-dir",
            str(source_dir),
            "--source-sha",
            SOURCE_SHA,
            "--event-sha",
            EVENT_SHA,
            "--github-output",
            str(github_output),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def parse_outputs(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines())


def test_verify_source_discovers_exact_contract_and_emits_raw_digest(tmp_path):
    source_dir = tmp_path / "downloaded source"
    source_dir.mkdir()
    tarball, metadata = write_contract(source_dir, content=b"raw downloaded bytes")
    github_output = tmp_path / "github-output"

    completed = run_verify(source_dir, github_output)

    assert completed.returncode == 0, completed.stderr
    outputs = parse_outputs(github_output)
    assert outputs == {
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_filename": tarball.name,
        "tarball_path": str(tarball.resolve()),
        "metadata_path": str(metadata.resolve()),
        "tarball_size": str(len(b"raw downloaded bytes")),
        "tarball_sha256": hashlib.sha256(b"raw downloaded bytes").hexdigest(),
    }


@pytest.mark.parametrize("tarball_count", [0, 2])
def test_verify_source_rejects_zero_or_multiple_tarballs(tmp_path, tarball_count):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    _, metadata = write_contract(source_dir)
    (source_dir / "colocboost_1.0.0.tar.gz").unlink()
    if tarball_count == 2:
        (source_dir / "first.tar.gz").write_bytes(b"first")
        (source_dir / "second.tar.gz").write_bytes(b"second")

    completed = run_verify(source_dir, tmp_path / "github-output")

    assert completed.returncode == 2
    assert "exactly one" in completed.stderr.lower()
    assert metadata.exists()


def test_verify_source_rejects_symlink_tarball(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    outside = tmp_path / "outside.tar.gz"
    outside.write_bytes(b"outside")
    symlink = source_dir / "colocboost_1.0.0.tar.gz"
    symlink.symlink_to(outside)
    (source_dir / "source-metadata.json").write_text("{}", encoding="utf-8")

    completed = run_verify(source_dir, tmp_path / "github-output")

    assert completed.returncode == 2
    assert "regular non-symlink" in completed.stderr.lower()


def test_verify_source_rejects_symlink_metadata(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    tarball = source_dir / "colocboost_1.0.0.tar.gz"
    tarball.write_bytes(b"raw")
    outside = tmp_path / "outside.json"
    outside.write_text("{}", encoding="utf-8")
    (source_dir / "source-metadata.json").symlink_to(outside)

    completed = run_verify(source_dir, tmp_path / "github-output")

    assert completed.returncode == 2
    assert "metadata" in completed.stderr.lower()
    assert "regular non-symlink" in completed.stderr.lower()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFO is unavailable")
def test_verify_source_rejects_fifo_tarball_without_blocking(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    os.mkfifo(source_dir / "package.tar.gz")
    (source_dir / "source-metadata.json").write_text("{}", encoding="utf-8")

    completed = run_verify(source_dir, tmp_path / "github-output", timeout=5)

    assert completed.returncode == 2
    assert "regular non-symlink" in completed.stderr.lower()


def test_verify_source_rejects_extra_downloaded_entries(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    write_contract(source_dir)
    (source_dir / "unexpected.txt").write_text("unexpected", encoding="utf-8")

    completed = run_verify(source_dir, tmp_path / "github-output")

    assert completed.returncode == 2
    assert "exactly" in completed.stderr.lower()
