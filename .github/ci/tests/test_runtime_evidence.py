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
        "clang-asan": (
            "clang-asan", "/opt/R/devel-asan/bin/R", "Under development (unstable)",
            "clang -fsanitize=address,undefined -fno-sanitize=float-divide-by-zero -fno-sanitize=alignment -fno-omit-frame-pointer -fsanitize=pointer-overflow -fsanitize=signed-integer-overflow",
            "clang++ -fsanitize=address,undefined -fno-sanitize=float-divide-by-zero -fno-sanitize=alignment -fno-omit-frame-pointer -fsanitize=pointer-overflow -fsanitize=signed-integer-overflow -frtti",
            "flang-new", "ubuntu", "22.04",
        ),
        "clang-ubsan": (
            "clang-ubsan", "/opt/R/devel-asan/bin/R", "Under development (unstable)",
            "clang -fsanitize=undefined -fno-sanitize=function -fno-omit-frame-pointer",
            "clang++ -fsanitize=undefined -fno-sanitize=function -fno-omit-frame-pointer -frtti",
            "flang-new", "ubuntu", "22.04",
        ),
        "donttest": (
            "donttest", "/opt/R/devel/bin/R", "Under development (unstable)",
            "gcc", "g++", "gfortran", "ubuntu", "22.04",
        ),
        "gcc-asan": (
            "gcc-asan", "/opt/R/devel/bin/R", "Under development (unstable)",
            "gcc -fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer",
            "g++ -fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer",
            "gfortran", "fedora", "42",
        ),
        "gcc-ubsan": (
            "gcc-ubsan", "/opt/R/devel/bin/R", "Under development (unstable)",
            "gcc -fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer",
            "g++ -fsanitize=address,undefined,bounds-strict -fno-omit-frame-pointer",
            "gfortran", "fedora", "42",
        ),
        "nold": (
            "nold", "/opt/R/devel-nold/bin/R", "Under development (unstable)",
            "gcc", "g++", "gfortran", "ubuntu", "22.04",
        ),
        "valgrind": (
            "valgrind", "/opt/R/devel-valgrind/bin/R", "Under development (unstable)",
            "gcc", "g++", "gfortran", "fedora", "42",
        ),
        "vnu": (
            "vnu", "/opt/R/release/bin/R", "",
            "gcc", "g++", "gfortran", "ubuntu", "24.04",
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
    libraries = ["/opt/R/lib/R/lib/libR.so"]
    if environment_id == "atlas":
        libraries.append("/usr/lib64/atlas/libsatlas.so.3")
    elif environment_id == "clang-asan":
        libraries.append("/usr/lib/llvm-22/lib/clang/22/lib/linux/libclang_rt.asan-x86_64.so")
    elif environment_id == "clang-ubsan":
        libraries.append("/usr/lib/llvm-22/lib/clang/22/lib/linux/libclang_rt.ubsan_standalone-x86_64.so")
    elif environment_id in {"gcc-asan", "gcc-ubsan"}:
        libraries.extend(["/usr/lib64/libasan.so.8", "/usr/lib64/libubsan.so.1"])
    environment_names = {
        "PATH", "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE",
        "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "CC", "CXX", "FC", "F77",
        "LANG", "LC_ALL", "LC_CTYPE", "ASAN_OPTIONS", "UBSAN_OPTIONS",
        "LD_PRELOAD", "VALGRIND_OPTS", "CHECK_ARGS",
        "_R_CHECK_DONTTEST_EXAMPLES_",
        "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "R_COMPILE_PKGS",
        "R_JIT_STRATEGY", "R_CHECK_CONSTANTS",
    }
    environment = {name: None for name in environment_names}
    environment["PATH"] = f"{Path(selected_r).parent}:/usr/bin"
    active_args = {
        "atlas": "--no-manual --no-build-vignettes",
        "clang-asan": "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes",
        "clang-ubsan": "--extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes",
        "donttest": "--no-manual --no-build-vignettes",
        "gcc-asan": "--no-manual --no-build-vignettes",
        "gcc-ubsan": "--no-manual --no-build-vignettes",
        "nold": "--no-manual --no-build-vignettes",
        "valgrind": "--use-valgrind --extra-arch --no-stop-on-test-error --no-manual --no-build-vignettes",
        "vnu": "--no-manual --no-build-vignettes",
    }
    if environment_id in active_args:
        environment["CHECK_ARGS"] = active_args[environment_id]
    if environment_id == "clang-asan":
        environment.update({
            "ASAN_OPTIONS": "detect_leaks=0:alloc_dealloc_mismatch=0",
            "UBSAN_OPTIONS": "print_stacktrace=1",
        })
    elif environment_id == "clang-ubsan":
        environment.update({
            "ASAN_OPTIONS": "detect_leaks=0:alloc_dealloc_mismatch=0",
            "UBSAN_OPTIONS": "print_stacktrace=1",
        })
    elif environment_id in {"gcc-asan", "gcc-ubsan"}:
        environment.update({
            "ASAN_OPTIONS": "detect_leaks=0",
            "UBSAN_OPTIONS": "print_stacktrace=1",
            "LD_PRELOAD": "/usr/lib64/libasan.so.8:/usr/lib64/libubsan.so.1",
        })
    elif environment_id == "donttest":
        environment["_R_CHECK_DONTTEST_EXAMPLES_"] = "true"
    elif environment_id == "valgrind":
        environment["VALGRIND_OPTS"] = "--track-origins=yes --leak-check=full"
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
        "cc_path": f"/usr/bin/{cc.split()[0]}",
        "cc_version": f"{compiler_family} version {compiler_major}.0.0",
        "cxx": cxx,
        "cxx_path": f"/usr/bin/{cxx.split()[0]}",
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
        "makeconf": {
            "CFLAGS": "-g -O2 -fsanitize=address" if environment_id in {"gcc-asan", "gcc-ubsan"} else "-g -O2",
            "CXXFLAGS": "-g -O2",
            "FFLAGS": "-g -O2",
            "MAIN_LDFLAGS": "-fsanitize=address,undefined -pthread" if environment_id in {"gcc-asan", "gcc-ubsan"} else "",
            "SAN_LIBS": "-lclang_rt.ubsan_standalone-x86_64" if environment_id == "clang-ubsan" else "",
        },
        "matrix_dimension": 3,
        "matrix_checksum": 729.0,
        "maps_method": "proc-self-maps",
        "loaded_libraries": libraries,
        "long_double": environment_id != "nold",
        "environment": environment,
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
        "clang-asan",
        "clang-ubsan",
        "donttest",
        "gcc-asan",
        "gcc-ubsan",
        "nold",
        "valgrind",
        "vnu",
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


@pytest.mark.parametrize(
    ("environment_id", "mutation"),
    [
        ("clang-asan", ("cc", "clang")),
        ("clang-asan", ("cc_path", "/usr/bin/gcc")),
        ("clang-ubsan", ("cc_version", "gcc (GCC) 22.0.0")),
        ("clang-asan", ("loaded_libraries", ["/opt/R/lib/libR.so"])),
        ("clang-ubsan", ("makeconf", {"CFLAGS": "", "CXXFLAGS": "", "FFLAGS": "", "MAIN_LDFLAGS": "", "SAN_LIBS": ""})),
        ("clang-ubsan", ("loaded_libraries", ["/opt/R/lib/libR.so"])),
        ("donttest", ("environment", {"PATH": "/opt/R/devel/bin:/usr/bin"})),
        ("gcc-asan", ("environment", {"PATH": "/opt/R/devel/bin:/usr/bin", "LD_PRELOAD": "/usr/lib64/libubsan.so.1", "ASAN_OPTIONS": "detect_leaks=0", "UBSAN_OPTIONS": "print_stacktrace=1"})),
        ("gcc-asan", ("cc_path", "/usr/bin/clang")),
        ("gcc-ubsan", ("cc_version", "Apple clang version 17.0.0")),
        ("gcc-ubsan", ("loaded_libraries", ["/opt/R/lib/libR.so", "/usr/lib64/libasan.so.8"])),
        ("nold", ("long_double", True)),
        ("valgrind", ("environment", {"PATH": "/opt/R/devel-valgrind/bin:/usr/bin"})),
    ],
)
def test_active_runtime_profiles_reject_missing_special_identity(tmp_path, environment_id, mutation):
    document = valid_document(environment_id)
    field, value = mutation
    document[field] = value
    assert invoke(tmp_path, document).returncode != 0


def test_vnu_preserves_manifest_lexical_r_while_accepting_versioned_resolution(tmp_path):
    document = valid_document("vnu")
    document["r_resolved"] = "/opt/R/4.6.1/bin/R"
    completed = invoke(tmp_path, document)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "environment_id",
    ["r-devel-linux-x86-64-debian-clang", "clang-asan", "clang-ubsan"],
)
def test_clang_runtime_accepts_resolved_compiler_symlink_targets(
    tmp_path, environment_id
):
    document = valid_document(environment_id)
    document["cc_path"] = "/usr/lib/llvm-22/bin/clang"
    document["cxx_path"] = "/usr/lib/llvm-22/bin/clang"
    document["fc_path"] = "/usr/lib/llvm-22/bin/flang-new"

    completed = invoke(tmp_path, document)

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "resolved_r",
    [
        "/opt/R/not-a-version/bin/R",
        "/opt/R/4.6.1/bin/Rscript",
        "/opt/R/4.6.1/../other/bin/R",
    ],
)
def test_vnu_rejects_unbound_resolved_r_alias(tmp_path, resolved_r):
    document = valid_document("vnu")
    document["r_resolved"] = resolved_r

    completed = invoke(tmp_path, document)

    assert completed.returncode != 0


def test_rejects_open_document(tmp_path):
    document = valid_document()
    document["untrusted"] = True
    assert invoke(tmp_path, document).returncode != 0


def test_runtime_probe_has_no_unprepared_jsonlite_dependency_and_keeps_lexical_r():
    source = (CI_DIR / "probe-runtime.R").read_text(encoding="utf-8")
    assert "requireNamespace(\"jsonlite\"" not in source
    assert "jsonlite::" not in source
    assert "r_executable = lexical_r" in source
    assert "r_resolved = selected_r" in source
    assert '"CHECK_ARGS"' in source
    assert '"VALGRIND_OPTS"' in source
