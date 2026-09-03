import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from artifact_contract import create_metadata, verify_tarball  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def write_tarball(path, content=b"raw tar.gz bytes"):
    path.write_bytes(content)
    return path


def create_valid_metadata(tmp_path, content=b"raw tar.gz bytes"):
    tarball = write_tarball(tmp_path / "colocboost_1.0.0.tar.gz", content)
    metadata_path = tmp_path / "source-metadata.json"
    metadata = create_metadata(
        tarball, metadata_path, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
    )
    return tarball, metadata_path, metadata


def test_create_metadata_hashes_raw_bytes_and_records_event_sha(tmp_path):
    content = b"\x1f\x8b\x08\x00not reconstructed tar members\x00\xff"
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path, content)

    assert metadata == {
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "filename": tarball.name,
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }
    assert json.loads(metadata_path.read_text()) == metadata
    assert verify_tarball(
        tarball,
        metadata_path,
        expected_source_sha=SOURCE_SHA,
        expected_event_sha=EVENT_SHA,
    ) == metadata


def test_verify_rejects_corrupted_tarball_bytes(tmp_path):
    tarball, metadata_path, _ = create_valid_metadata(tmp_path)
    tarball.write_bytes(tarball.read_bytes() + b"corruption")

    with pytest.raises(ValueError, match="size|sha256"):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


def test_verify_rejects_filename_mismatch(tmp_path):
    tarball, metadata_path, _ = create_valid_metadata(tmp_path)
    renamed = tarball.with_name("different.tar.gz")
    tarball.rename(renamed)

    with pytest.raises(ValueError, match="filename"):
        verify_tarball(
            renamed,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


@pytest.mark.parametrize(
    ("expected_source_sha", "expected_event_sha", "match"),
    [
        ("c" * 40, EVENT_SHA, "source_sha"),
        (SOURCE_SHA, "d" * 40, "event_sha"),
    ],
)
def test_verify_rejects_wrong_expected_identity(
    tmp_path, expected_source_sha, expected_event_sha, match
):
    tarball, metadata_path, _ = create_valid_metadata(tmp_path)

    with pytest.raises(ValueError, match=match):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=expected_source_sha,
            expected_event_sha=expected_event_sha,
        )


@pytest.mark.parametrize("field", ["source_sha", "event_sha", "filename", "size", "sha256"])
def test_verify_rejects_missing_metadata_keys(tmp_path, field):
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path)
    metadata.pop(field)
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match=field):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


def test_verify_rejects_unknown_metadata_keys(tmp_path):
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path)
    metadata["actions_artifact_digest"] = "c" * 64
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match="actions_artifact_digest"):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("source_sha", "A" * 40),
        ("event_sha", "b" * 39),
        ("size", True),
        ("size", -1),
        ("sha256", "C" * 64),
        ("sha256", "c" * 63),
    ],
)
def test_verify_rejects_malformed_metadata_values(tmp_path, field, value):
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path)
    metadata[field] = value
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match=field):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


@pytest.mark.parametrize(
    "unsafe_name", ["../package.tar.gz", "sub/package.tar.gz", r"sub\package.tar.gz", ".", ".."]
)
def test_verify_rejects_filename_traversal(tmp_path, unsafe_name):
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path)
    metadata["filename"] = unsafe_name
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match="filename"):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


@pytest.mark.parametrize(("field", "value"), [("size", 1), ("sha256", "c" * 64)])
def test_verify_recomputes_size_and_sha256(tmp_path, field, value):
    tarball, metadata_path, metadata = create_valid_metadata(tmp_path)
    metadata[field] = value
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match=field):
        verify_tarball(
            tarball,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


def test_create_and_verify_reject_symlink_tarballs(tmp_path):
    target = write_tarball(tmp_path / "target.tar.gz")
    symlink = tmp_path / "linked.tar.gz"
    symlink.symlink_to(target)
    metadata_path = tmp_path / "source-metadata.json"

    with pytest.raises(ValueError, match="regular"):
        create_metadata(
            symlink, metadata_path, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
        )

    metadata = {
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "filename": symlink.name,
        "size": target.stat().st_size,
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
    }
    metadata_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="regular"):
        verify_tarball(
            symlink,
            metadata_path,
            expected_source_sha=SOURCE_SHA,
            expected_event_sha=EVENT_SHA,
        )


def test_create_and_verify_reject_non_regular_tarball(tmp_path):
    directory = tmp_path / "not-a-tarball.tar.gz"
    directory.mkdir()
    metadata_path = tmp_path / "source-metadata.json"

    with pytest.raises(ValueError, match="regular"):
        create_metadata(
            directory, metadata_path, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
        )


def test_create_metadata_validates_identity_before_writing(tmp_path):
    tarball = write_tarball(tmp_path / "package.tar.gz")
    metadata_path = tmp_path / "source-metadata.json"

    with pytest.raises(ValueError, match="source_sha"):
        create_metadata(
            tarball, metadata_path, source_sha="A" * 40, event_sha=EVENT_SHA
        )
    assert not metadata_path.exists()


def test_atomic_replace_failure_preserves_existing_metadata(tmp_path, monkeypatch):
    tarball, metadata_path, _ = create_valid_metadata(tmp_path, b"first")
    before = metadata_path.read_bytes()
    tarball.write_bytes(b"second")

    def fail_replace(source, destination):
        raise OSError("simulated replacement failure")

    monkeypatch.setattr(os, "replace", fail_replace)
    with pytest.raises(OSError, match="replacement failure"):
        create_metadata(
            tarball, metadata_path, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
        )

    assert metadata_path.read_bytes() == before
    assert not list(tmp_path.glob(f".{metadata_path.name}.*.tmp"))


def test_artifact_contract_cli_creates_then_verifies_metadata(tmp_path):
    tarball = write_tarball(tmp_path / "package with spaces.tar.gz", b"cli bytes")
    metadata_path = tmp_path / "source metadata.json"
    script = CI_DIR / "artifact_contract.py"

    created = subprocess.run(
        [
            sys.executable,
            str(script),
            "create",
            "--tarball",
            str(tarball),
            "--metadata",
            str(metadata_path),
            "--source-sha",
            SOURCE_SHA,
            "--event-sha",
            EVENT_SHA,
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert created.returncode == 0, created.stderr
    assert json.loads(created.stdout)["sha256"] == hashlib.sha256(b"cli bytes").hexdigest()

    verified = subprocess.run(
        [
            sys.executable,
            str(script),
            "verify",
            "--tarball",
            str(tarball),
            "--metadata",
            str(metadata_path),
            "--source-sha",
            SOURCE_SHA,
            "--event-sha",
            EVENT_SHA,
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert verified.returncode == 0, verified.stderr
    assert json.loads(verified.stdout) == json.loads(created.stdout)
