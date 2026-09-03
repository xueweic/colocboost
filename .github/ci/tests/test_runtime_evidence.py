import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
SCRIPT = CI_DIR / "verify_runtime_evidence.py"
MANIFEST = CI_DIR / "check-matrix.yml"
SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40
TARBALL_SHA256 = "c" * 64


def valid_document(environment_id="r-devel-linux-x86-64-debian-clang"):
    contracts = {
        "r-devel-linux-x86-64-debian-clang": (
            "clang22", "/opt/R/devel/bin/R", "Under development (unstable)",
            "clang-22", "clang++-22", "flang-new-22", "ubuntu", "22.04",
        ),
        "r-devel-linux-x86-64-debian-gcc": (
            "ubuntu-gcc16", "/opt/R/devel/bin/R", "Under development (unstable)",
            "gcc-16", "g++-16", "gfortran-16", "ubuntu", "24.04",
        ),
        "r-devel-linux-x86-64-fedora-gcc": (
            "gcc16", "/opt/R/devel-gcc16/bin/R", "Under development (unstable)",
            "gcc", "g++", "gfortran", "fedora", "44",
        ),
        "r-release-linux-x86-64": (
            "ubuntu-release", "/opt/R/release/bin/R", "", "gcc", "g++", "gfortran",
            "ubuntu", "24.04",
        ),
        "atlas": (
            "atlas", "/opt/R/devel/bin/R", "Under development (unstable)",
            "gcc", "g++", "gfortran", "fedora", "42",
        ),
    }
    (
        profile, selected_r, r_status, cc, cxx, fc, distribution,
        distribution_version,
    ) = contracts[environment_id]
    compiler_family = "clang" if "clang" in cc else "gcc"
    compiler_major = "22" if compiler_family == "clang" else (
        "16" if profile in {"ubuntu-gcc16", "gcc16"} else "13"
    )
    libraries = ["/usr/lib64/libR.so", "/usr/lib64/atlas/libsatlas.so.3"] if environment_id == "atlas" else ["/opt/R/lib/R/lib/libR.so"]
    return {
        "schema_version": 1,
        "kind": "r-runtime-proof",
        "status": "pass",
        "environment_id": environment_id,
        "runtime_profile": profile,
        "source_sha": SOURCE_SHA,
        "event_sha": EVENT_SHA,
        "tarball_sha256": TARBALL_SHA256,
        "r_executable": selected_r,
        "r_resolved": selected_r,
        "r_version": "R Under development (unstable) (2026-09-01 r99999)" if r_status else "R version 4.6.1 (2026-06-12)",
        "r_status": r_status,
        "r_platform": "x86_64-pc-linux-gnu",
        "os": "linux",
        "os_release": "6.11.0",
        "distribution_id": distribution,
        "distribution_version": distribution_version,
        "architecture": "x86_64",
        "compiler": f"{compiler_family} version {compiler_major}.0.0",
        "cc": cc,
        "cc_path": f"/usr/bin/{cc}",
        "cc_version": f"{compiler_family} version {compiler_major}.0.0",
        "cxx": cxx,
        "cxx_path": f"/usr/bin/{cxx}",
        "cxx_version": f"{compiler_family}++ version {compiler_major}.0.0",
        "fc": fc,
        "fc_path": f"/usr/bin/{fc}",
        "fc_version": (
            f"flang-new version {compiler_major}.0.0"
            if "flang" in fc else f"GNU Fortran version {compiler_major}.0.0"
        ),
        "locale": "C",
        "session_info": "R version evidence",
        "ext_soft_version": {"BLAS": "atlas" if environment_id == "atlas" else "generic"},
        "la_library": "/usr/lib64/atlas/libsatlas.so.3" if environment_id == "atlas" else "/usr/lib/libblas.so",
        "la_version": "3.10.3",
        "blas_libs": "-lsatlas" if environment_id == "atlas" else "-lblas",
        "matrix_dimension": 3,
        "matrix_checksum": 729.0,
        "maps_method": "proc-self-maps",
        "loaded_libraries": libraries,
        "long_double": True,
        "environment": {"PATH": f"{Path(selected_r).parent}:/usr/bin"},
    }


def invoke(tmp_path, document, environment_id=None):
    evidence = tmp_path / "runtime.json"
    evidence.write_text(json.dumps(document) + "\n", encoding="utf-8")
    expected_id = environment_id or document["environment_id"]
    return subprocess.run(
        [
            sys.executable, "-B", os.fspath(SCRIPT), f"--manifest={MANIFEST}",
            f"--environment-id={expected_id}", f"--evidence={evidence}",
            f"--source-sha={SOURCE_SHA}", f"--event-sha={EVENT_SHA}",
            f"--tarball-sha256={TARBALL_SHA256}",
        ], check=False, capture_output=True, text=True,
    )


@pytest.mark.parametrize(
    "environment_id",
    [
        "r-devel-linux-x86-64-debian-clang",
        "r-devel-linux-x86-64-debian-gcc",
        "r-devel-linux-x86-64-fedora-gcc",
        "r-release-linux-x86-64",
        "atlas",
    ],
)
def test_accepts_selected_r_same_process_runtime_proof(tmp_path, environment_id):
    completed = invoke(tmp_path, valid_document(environment_id))
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("status", "fail"),
        ("runtime_profile", "ubuntu-release"),
        ("source_sha", "d" * 40),
        ("event_sha", "d" * 40),
        ("tarball_sha256", "d" * 64),
        ("r_executable", "/tmp/R"),
        ("r_resolved", "/tmp/R"),
        ("r_status", ""),
        ("os", "darwin"),
        ("architecture", "arm64"),
        ("cc", "gcc"),
        ("cc_path", "/usr/bin/clang-21"),
        ("cc_version", "clang version 21.0.0"),
        ("distribution_id", "fedora"),
        ("distribution_version", "44"),
        ("matrix_checksum", 0),
        ("maps_method", "python-process"),
    ],
)
def test_rejects_mismatched_or_nonexecuted_runtime_proof(tmp_path, field, value):
    document = valid_document()
    document[field] = value
    completed = invoke(tmp_path, document, "r-devel-linux-x86-64-debian-clang")
    assert completed.returncode != 0


@pytest.mark.parametrize("library", ["libmkl_core.so", "libopenblas.so", "libblis.so"])
def test_atlas_rejects_forbidden_loaded_blas(tmp_path, library):
    document = valid_document("atlas")
    document["loaded_libraries"].append(f"/usr/lib/{library}")
    completed = invoke(tmp_path, document)
    assert completed.returncode != 0


def test_atlas_requires_lassatlas_from_same_r_process(tmp_path):
    document = valid_document("atlas")
    document["loaded_libraries"] = ["/usr/lib/libR.so"]
    completed = invoke(tmp_path, document)
    assert completed.returncode != 0


def test_rejects_open_document(tmp_path):
    document = valid_document()
    document["untrusted"] = True
    assert invoke(tmp_path, document).returncode != 0
