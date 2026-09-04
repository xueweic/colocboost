import json
import sys
from pathlib import Path

import pytest

CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))
from compare_arm_pair import ARM_INDEX, CHECK_ARGS, CHILD_IMAGES, compare_pair  # noqa: E402


def child(arch, source="a" * 40, event="b" * 40, tarball="c" * 64):
    machine = "x86_64" if arch == "amd64" else "aarch64"
    return {"environment_id": "linux-arm64", "architecture": arch, "image": CHILD_IMAGES[arch], "index": ARM_INDEX,
            "system_r": "/usr/bin/R", "source_sha": source, "event_sha": event, "tarball_sha256": tarball,
            "status": "pass", "check_status": "pass", "check_exit_code": 0,
            "check_args": CHECK_ARGS,
            "configuration": {"args": CHECK_ARGS, "r": "/usr/bin/R"}, "tarball": "shared-source",
            "native": {"clean": True, "uname_m": machine, "r_platform": machine + "-linux-gnu",
                       "r_home": "/usr/lib/R", "image": CHILD_IMAGES[arch], "index": ARM_INDEX}}


def write(path, doc):
    path.write_text(json.dumps(doc))


def test_pair_requires_two_explicit_native_identities(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    write(a, child("amd64")); write(b, child("arm64"))
    result = tmp_path / "result.json"; evidence = tmp_path / "comparison.json"
    comparison = compare_pair(a, b, result_path=result, evidence_path=evidence, source_sha="a" * 40, event_sha="b" * 40, tarball_sha256="c" * 64)
    assert comparison["errors"] == []


def test_none_or_wrong_child_identity_cannot_compare_green(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    left, right = child("amd64"), child("arm64")
    right["native"] = None
    write(a, left); write(b, right)
    comparison = compare_pair(a, b, result_path=tmp_path / "result.json", evidence_path=tmp_path / "comparison.json", source_sha="a" * 40, event_sha="b" * 40, tarball_sha256="c" * 64)
    assert comparison["errors"]
    assert json.loads((tmp_path / "result.json").read_text())["status"] == "fail"


def test_wrong_architecture_image_is_rejected(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    left, right = child("amd64"), child("arm64")
    right["image"] = left["image"]
    write(a, left); write(b, right)
    comparison = compare_pair(a, b, result_path=tmp_path / "result.json", evidence_path=tmp_path / "comparison.json", source_sha="a" * 40, event_sha="b" * 40, tarball_sha256="c" * 64)
    assert comparison["errors"]


def test_extra_fields_and_nonexact_configuration_are_rejected(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    left, right = child("amd64"), child("arm64")
    left["unexpected"] = True
    right["configuration"]["args"] = ["--as-cran"]
    write(a, left); write(b, right)
    comparison = compare_pair(a, b, result_path=tmp_path / "result.json", evidence_path=tmp_path / "comparison.json", source_sha="a" * 40, event_sha="b" * 40, tarball_sha256="c" * 64)
    assert any("fields differ" in error for error in comparison["errors"])
    assert any("configuration" in error for error in comparison["errors"])


def test_duplicate_json_keys_fail_closed(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text('{"environment_id":"linux-arm64","environment_id":"linux-arm64"}')
    write(b, child("arm64"))
    with pytest.raises(ValueError, match="duplicate JSON key"):
        compare_pair(a, b, result_path=tmp_path / "result.json", evidence_path=tmp_path / "comparison.json", source_sha="a" * 40, event_sha="b" * 40, tarball_sha256="c" * 64)
