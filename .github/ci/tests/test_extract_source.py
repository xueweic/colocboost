import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "extract_source.py"
sys.path.insert(0, os.fspath(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def add_bytes(archive, name, data=b"x", kind="file"):
    member = tarfile.TarInfo(name)
    if kind == "symlink":
        member.type = tarfile.SYMTYPE
        member.linkname = "DESCRIPTION"
        archive.addfile(member)
    else:
        member.size = len(data)
        archive.addfile(member, io.BytesIO(data))


def make_contract(tmp_path, members):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    tarball = source_dir / "colocboost_1.0.9.tar.gz"
    with tarfile.open(tarball, "w:gz") as archive:
        for member in members:
            add_bytes(archive, *member)
    metadata = source_dir / "source-metadata.json"
    create_metadata(tarball, metadata, source_sha=SOURCE_SHA, event_sha=EVENT_SHA)
    return tarball, metadata


def invoke(tmp_path, tarball, metadata, destination=None):
    if destination is None:
        destination = tmp_path / "extracted"
    github_output = tmp_path / "github-output"
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            os.fspath(SCRIPT),
            f"--tarball={tarball}",
            f"--metadata={metadata}",
            f"--source-sha={SOURCE_SHA}",
            f"--event-sha={EVENT_SHA}",
            f"--destination={destination}",
            f"--github-output={github_output}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed, destination, github_output


def test_extracts_verified_archive_to_one_safe_top_level_package(tmp_path):
    tarball, metadata = make_contract(
        tmp_path,
        [
            ("colocboost/DESCRIPTION", b"Package: colocboost\n"),
            ("colocboost/R/example.R", b"x <- 1\n"),
        ],
    )

    completed, destination, github_output = invoke(tmp_path, tarball, metadata)

    assert completed.returncode == 0, completed.stderr
    package = destination / "colocboost"
    assert (package / "DESCRIPTION").is_file()
    assert (package / "R" / "example.R").is_file()
    assert github_output.read_text().strip() == f"package_path={package}"


@pytest.mark.parametrize(
    "members",
    [
        [("../escape", b"bad")],
        [("colocboost/DESCRIPTION", b"Package: colocboost\n"), ("other/file", b"x")],
        [("colocboost/DESCRIPTION", b"", "symlink")],
        [("/absolute", b"bad")],
    ],
    ids=("traversal", "two-roots", "symlink", "absolute"),
)
def test_rejects_unsafe_or_ambiguous_archive_members(tmp_path, members):
    tarball, metadata = make_contract(tmp_path, members)

    completed, destination, github_output = invoke(tmp_path, tarball, metadata)

    assert completed.returncode != 0
    assert not destination.exists()
    assert not github_output.exists()


def test_rejects_nonempty_destination_and_corrupt_source(tmp_path):
    tarball, metadata = make_contract(
        tmp_path, [("colocboost/DESCRIPTION", b"Package: colocboost\n")]
    )
    destination = tmp_path / "extracted"
    destination.mkdir()
    (destination / "owned").write_text("keep", encoding="utf-8")

    completed, _, _ = invoke(tmp_path, tarball, metadata, destination)

    assert completed.returncode != 0
    assert (destination / "owned").read_text(encoding="utf-8") == "keep"

    tarball.write_bytes(b"corrupt")
    completed, _, _ = invoke(tmp_path, tarball, metadata, tmp_path / "other")
    assert completed.returncode != 0
