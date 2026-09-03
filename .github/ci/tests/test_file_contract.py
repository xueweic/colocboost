import hashlib
import json
import os
import stat
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(CI_DIR))

from verify_file_contract import verify_bound_file  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
DIGEST = "c" * 64


def executable(path, content=b"#!/bin/sh\nexit 0\n"):
    path.write_bytes(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def invoke(tmp_path, *, purpose="vnu-dispatcher", mutation=None):
    output = tmp_path / "proof.json"
    if purpose == "vnu-dispatcher":
        candidate = executable(tmp_path / "vnu.sh")
        row = {
            "id": "vnu",
            "vnu_path": os.fspath(candidate),
            "vnu_sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
        }
        reference = None
    else:
        reference = tmp_path / "rhub-valgrind.supp"
        reference.write_bytes(b"official suppression\n")
        candidate = tmp_path / "default.supp"
        candidate.write_bytes(b"system prefix\n" + reference.read_bytes())
        row = {
            "id": "valgrind",
            "suppression_path": os.fspath(candidate),
            "suppression_suffix_bytes": reference.stat().st_size,
            "suppression_suffix_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        }
    if mutation:
        mutation(row, candidate, reference)
    proof = verify_bound_file(
        row, environment_id=row["id"], purpose=purpose, path=candidate,
        source_sha=SOURCE_SHA, event_sha=EVENT_SHA,
        tarball_sha256=DIGEST, output=output, reference=reference,
    )
    return proof, output, candidate, reference


def test_exact_dispatcher_contract_is_identity_bound_and_atomic(tmp_path):
    proof, output, candidate, _ = invoke(tmp_path)
    assert proof == json.loads(output.read_text())
    assert proof["purpose"] == "vnu-dispatcher"
    assert proof["path"] == os.fspath(candidate)
    assert proof["executable"] is True
    assert proof["suffix_bytes"] is None
    assert proof["source_sha"] == SOURCE_SHA
    assert proof["event_sha"] == EVENT_SHA
    assert proof["tarball_sha256"] == DIGEST


def test_valgrind_contract_requires_exact_committed_suffix_bytes(tmp_path):
    proof, _, candidate, reference = invoke(tmp_path, purpose="valgrind-suppression")
    assert proof["path"] == os.fspath(candidate)
    assert proof["suffix_bytes"] == reference.stat().st_size
    assert proof["sha256"] == hashlib.sha256(candidate.read_bytes()).hexdigest()
    assert proof["suffix_sha256"] == hashlib.sha256(reference.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda row, candidate, ref: row.__setitem__("vnu_path", "/tmp/other-vnu.sh"),
        lambda row, candidate, ref: row.__setitem__("vnu_sha256", "0" * 64),
        lambda row, candidate, ref: candidate.chmod(stat.S_IRUSR | stat.S_IWUSR),
    ],
)
def test_dispatcher_rejects_wrong_path_hash_or_mode(tmp_path, mutation):
    with pytest.raises(ValueError):
        invoke(tmp_path, mutation=mutation)


@pytest.mark.parametrize("kind", ["mutated", "truncated", "wrong-size", "wrong-hash"])
def test_suppression_rejects_any_suffix_drift(tmp_path, kind):
    def mutate(row, candidate, reference):
        if kind == "mutated":
            candidate.write_bytes(candidate.read_bytes()[:-1] + b"X")
        elif kind == "truncated":
            candidate.write_bytes(candidate.read_bytes()[:-1])
        elif kind == "wrong-size":
            row["suppression_suffix_bytes"] += 1
        else:
            row["suppression_suffix_sha256"] = "0" * 64

    with pytest.raises(ValueError):
        invoke(tmp_path, purpose="valgrind-suppression", mutation=mutate)


@pytest.mark.parametrize("kind", ["symlink", "fifo"])
def test_bound_file_rejects_links_and_special_files(tmp_path, kind):
    def mutate(row, candidate, _):
        candidate.unlink()
        target = tmp_path / "target"
        target.write_bytes(b"target")
        if kind == "symlink":
            candidate.symlink_to(target)
        else:
            os.mkfifo(candidate)

    with pytest.raises(ValueError, match="regular non-symlink"):
        invoke(tmp_path, mutation=mutate)


def test_file_proof_rejects_open_output_document_contract(tmp_path):
    proof, output, candidate, _ = invoke(tmp_path)
    proof["unexpected"] = True
    output.write_text(json.dumps(proof))
    # A new verification replaces stale/open output rather than trusting it.
    new = verify_bound_file(
        {
            "id": "vnu", "vnu_path": os.fspath(candidate),
            "vnu_sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
        },
        environment_id="vnu", purpose="vnu-dispatcher", path=candidate,
        source_sha=SOURCE_SHA, event_sha=EVENT_SHA, tarball_sha256=DIGEST,
        output=output,
    )
    assert set(json.loads(output.read_text())) == set(new)
