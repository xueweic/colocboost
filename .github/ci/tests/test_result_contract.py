import json
import os
import sys
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator


CI_DIR = Path(__file__).resolve().parents[1]
FIXTURE_DIR = Path(__file__).parent / "fixtures" / "results"
SCHEMA_PATH = CI_DIR / "result.schema.json"
SCHEMA = json.loads(SCHEMA_PATH.read_text())
SCHEMA_VALIDATOR = Draft202012Validator(SCHEMA)
DETAIL_VALUES = {
    "tests": {"pass": 0, "fail": 0, "error": 0, "warning": 0, "skip": 0},
    "check": {"errors": 0, "warnings": 0, "notes": 0, "log_path": "00check.log"},
    "predicate": "no-package-owned-native-code",
    "evidence": "No package-owned native source was found.",
    "classification": "runner",
    "message": "The runner stopped.",
}
KIND_FIELDS = {
    "unit": {"tests"},
    "package-check": {"check"},
    "applicability": {"predicate", "evidence"},
    "infrastructure": {"classification", "message"},
}
sys.path.insert(0, str(CI_DIR))

from result_contract import finalize, validate_result, write_provisional  # noqa: E402


def load_fixture(name):
    return json.loads((FIXTURE_DIR / f"{name}.json").read_text())


def assert_schema_invalid(document):
    assert list(SCHEMA_VALIDATOR.iter_errors(document))


@pytest.mark.parametrize(
    "name", ["unit", "package-check", "applicability", "infrastructure"]
)
def test_valid_kind_specific_documents(name):
    document = load_fixture(name)

    SCHEMA_VALIDATOR.validate(document)
    assert validate_result(document) == document


def test_schema_declares_discriminators_and_closes_objects():
    Draft202012Validator.check_schema(SCHEMA)

    assert SCHEMA["properties"]["result_kind"]["enum"] == [
        "unit",
        "package-check",
        "applicability",
        "infrastructure",
    ]
    assert SCHEMA["properties"]["status"]["enum"] == [
        "provisional",
        "pass",
        "fail",
        "not-applicable",
        "missing",
        "cancelled",
        "uncovered",
    ]
    assert SCHEMA["additionalProperties"] is False
    assert SCHEMA["properties"]["tests"]["additionalProperties"] is False
    assert SCHEMA["properties"]["check"]["additionalProperties"] is False


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        (lambda value: value.update(status="unknown"), "status"),
        (lambda value: value.pop("source_sha"), "source_sha"),
        (lambda value: value.update(source_sha="A" * 40), "source_sha"),
        (lambda value: value.update(event_sha="b" * 39), "event_sha"),
        (lambda value: value.update(tarball_sha256="c" * 63), "tarball_sha256"),
        (lambda value: value.update(extra_key=True), "extra_key"),
        (lambda value: value["tests"].pop("skip"), "skip"),
    ],
)
def test_invalid_documents_are_rejected(mutate, match):
    document = load_fixture("unit")
    mutate(document)

    assert_schema_invalid(document)
    with pytest.raises(ValueError, match=match):
        validate_result(document)


@pytest.mark.parametrize(
    ("fixture", "missing"),
    [
        ("unit", "tests"),
        ("package-check", "check"),
        ("applicability", "predicate"),
        ("applicability", "evidence"),
        ("infrastructure", "classification"),
        ("infrastructure", "message"),
    ],
)
def test_final_documents_require_kind_specific_fields(fixture, missing):
    document = load_fixture(fixture)
    document.pop(missing)

    assert_schema_invalid(document)
    with pytest.raises(ValueError, match=missing):
        validate_result(document)


def test_only_provisional_results_may_have_a_null_tarball_digest():
    document = load_fixture("unit")
    document["tarball_sha256"] = None

    assert_schema_invalid(document)
    with pytest.raises(ValueError, match="tarball_sha256"):
        validate_result(document)

    for field in ("tests",):
        document.pop(field)
    document["status"] = "provisional"
    SCHEMA_VALIDATOR.validate(document)
    assert validate_result(document) == document


def test_nested_counts_are_nonnegative_integers_and_objects_are_closed():
    boolean_count = load_fixture("unit")
    boolean_count["tests"]["pass"] = True
    assert_schema_invalid(boolean_count)
    with pytest.raises(ValueError, match="tests.pass"):
        validate_result(boolean_count)

    negative_count = load_fixture("package-check")
    negative_count["check"]["notes"] = -1
    assert_schema_invalid(negative_count)
    with pytest.raises(ValueError, match="check.notes"):
        validate_result(negative_count)

    extra_nested_key = load_fixture("package-check")
    extra_nested_key["check"]["duration"] = 10
    assert_schema_invalid(extra_nested_key)
    with pytest.raises(ValueError, match="duration"):
        validate_result(extra_nested_key)


@pytest.mark.parametrize(
    ("fixture", "foreign_field"),
    [
        (kind, field)
        for kind, allowed_fields in KIND_FIELDS.items()
        for field in sorted(DETAIL_VALUES.keys() - allowed_fields)
    ],
)
def test_result_kind_details_are_exclusive(fixture, foreign_field):
    document = load_fixture(fixture)
    document[foreign_field] = DETAIL_VALUES[foreign_field]

    assert_schema_invalid(document)
    with pytest.raises(ValueError, match=foreign_field):
        validate_result(document)


def test_write_provisional_then_finalize_atomically(tmp_path):
    path = tmp_path / "nested" / "result.json"
    identity = {
        "result_kind": "unit",
        "environment_id": "linux-r-release-openblas",
        "source_sha": "a" * 40,
        "event_sha": "b" * 40,
        "tarball_sha256": None,
    }

    provisional = write_provisional(path, identity)
    assert json.loads(path.read_text()) == provisional
    assert provisional["schema_version"] == 1
    assert provisional["status"] == "provisional"

    final = finalize(
        path,
        "pass",
        {
            "tarball_sha256": "c" * 64,
            "tests": {"pass": 3, "fail": 0, "error": 0, "warning": 0, "skip": 0},
        },
    )
    assert json.loads(path.read_text()) == final
    assert final["status"] == "pass"
    assert not list(path.parent.glob(f".{path.name}.*.tmp"))


def test_failed_finalization_preserves_valid_provisional(tmp_path):
    path = tmp_path / "result.json"
    identity = {
        "result_kind": "unit",
        "environment_id": "linux-r-release-openblas",
        "source_sha": "a" * 40,
        "event_sha": "b" * 40,
        "tarball_sha256": None,
    }
    write_provisional(path, identity)
    before = path.read_bytes()

    with pytest.raises(ValueError, match="tarball_sha256"):
        finalize(
            path,
            "pass",
            {"tests": {"pass": 3, "fail": 0, "error": 0, "warning": 0, "skip": 0}},
        )

    assert path.read_bytes() == before
    assert validate_result(json.loads(path.read_text()))["status"] == "provisional"


@pytest.mark.parametrize(
    ("identity_field", "attempted_value"),
    [
        ("schema_version", 1),
        ("result_kind", "package-check"),
        ("environment_id", "different-environment"),
        ("source_sha", "d" * 40),
        ("event_sha", "e" * 40),
    ],
)
def test_finalize_rejects_identity_fields_in_details(
    tmp_path, identity_field, attempted_value
):
    path = tmp_path / "result.json"
    identity = {
        "result_kind": "unit",
        "environment_id": "linux-r-release-openblas",
        "source_sha": "a" * 40,
        "event_sha": "b" * 40,
        "tarball_sha256": "c" * 64,
    }
    write_provisional(path, identity)
    before = path.read_bytes()
    details = {
        identity_field: attempted_value,
        "tests": {"pass": 3, "fail": 0, "error": 0, "warning": 0, "skip": 0},
    }
    if identity_field == "result_kind":
        details["check"] = DETAIL_VALUES["check"]

    with pytest.raises(ValueError, match=identity_field):
        finalize(path, "pass", details)

    assert path.read_bytes() == before


def test_finalize_rejects_details_for_a_different_result_kind(tmp_path):
    path = tmp_path / "result.json"
    identity = {
        "result_kind": "unit",
        "environment_id": "linux-r-release-openblas",
        "source_sha": "a" * 40,
        "event_sha": "b" * 40,
        "tarball_sha256": "c" * 64,
    }
    write_provisional(path, identity)
    before = path.read_bytes()

    with pytest.raises(ValueError, match="check"):
        finalize(
            path,
            "pass",
            {
                "tests": {
                    "pass": 3,
                    "fail": 0,
                    "error": 0,
                    "warning": 0,
                    "skip": 0,
                },
                "check": {
                    "errors": 0,
                    "warnings": 0,
                    "notes": 0,
                    "log_path": "00check.log",
                },
            },
        )

    assert path.read_bytes() == before


def test_replace_failure_preserves_valid_provisional(tmp_path, monkeypatch):
    path = tmp_path / "result.json"
    identity = {
        "result_kind": "infrastructure",
        "environment_id": "runner-startup",
        "source_sha": "a" * 40,
        "event_sha": "b" * 40,
        "tarball_sha256": "c" * 64,
    }
    write_provisional(path, identity)
    before = path.read_bytes()

    def fail_replace(source, destination):
        raise OSError("simulated replacement failure")

    monkeypatch.setattr(os, "replace", fail_replace)
    with pytest.raises(OSError, match="replacement failure"):
        finalize(
            path,
            "fail",
            {"classification": "runner", "message": "runner stopped"},
        )

    assert path.read_bytes() == before
    assert not list(path.parent.glob(f".{path.name}.*.tmp"))
