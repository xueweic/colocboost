#!/usr/bin/env python3
"""Validate runtime evidence captured inside the selected R process."""

from __future__ import annotations

import argparse
import json
import math
import re
import stat
import sys
from collections.abc import Mapping
from pathlib import Path

from validate_manifest import load_manifest, validate_manifest


FIELDS = {
    "schema_version", "kind", "status", "environment_id", "runtime_profile",
    "source_sha", "event_sha", "tarball_sha256", "r_executable", "r_resolved",
    "r_version", "r_status", "r_platform", "os", "os_release",
    "distribution_id", "distribution_version", "architecture", "compiler",
    "cc", "cc_path", "cc_version", "cxx", "cxx_path", "cxx_version",
    "fc", "fc_path", "fc_version", "locale", "session_info", "ext_soft_version",
    "la_library", "la_version", "blas_libs", "makeconf", "matrix_dimension",
    "matrix_checksum", "maps_method", "loaded_libraries", "long_double",
    "environment",
}
ENVIRONMENT_FIELDS = {
    "PATH", "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE",
    "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "CC", "CXX", "FC", "F77",
    "LANG", "LC_ALL", "LC_CTYPE", "ASAN_OPTIONS", "UBSAN_OPTIONS",
    "LD_PRELOAD", "VALGRIND_OPTS", "CHECK_ARGS",
    "_R_CHECK_DONTTEST_EXAMPLES_",
}
MAKECONF_FIELDS = {"CFLAGS", "CXXFLAGS", "FFLAGS", "MAIN_LDFLAGS", "SAN_LIBS"}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _single_line(value, label, *, allow_empty=False):
    if (
        not isinstance(value, str)
        or (not allow_empty and not value)
        or any(character in value for character in "\r\n\x00")
    ):
        raise ValueError(f"{label} must be a single-line string")


def _classify_r(document):
    status = document["r_status"].lower()
    version = document["r_version"].lower()
    if "under development" in status or "under development" in version:
        return "devel"
    if "patched" in status or "patched" in version:
        return "patched"
    return "release"


def _require_active_compiler(document, family, *, major=None):
    if document["compiler"] != document["cc_version"]:
        raise ValueError("runtime compiler summary does not match the executed C compiler")
    cc_name = Path(document["cc_path"]).name.lower()
    cxx_name = Path(document["cxx_path"]).name.lower()
    configured_cc = Path(document["cc"].split()[0]).name.lower()
    configured_cxx = Path(document["cxx"].split()[0]).name.lower()
    cc_version = document["cc_version"].lower()
    cxx_version = document["cxx_version"].lower()
    if family == "clang":
        if (
            re.fullmatch(r"clang(?:-[0-9]+)?", configured_cc) is None
            or re.fullmatch(r"clang\+\+(?:-[0-9]+)?", configured_cxx) is None
            or re.fullmatch(r"clang(?:-[0-9]+)?", cc_name) is None
            or re.fullmatch(r"clang(?:\+\+)?(?:-[0-9]+)?", cxx_name) is None
        ):
            raise ValueError("active Clang compiler paths are invalid")
        if "clang" not in cc_version or "clang" not in cxx_version:
            raise ValueError("active Clang compiler versions are invalid")
    elif family == "gcc":
        if re.fullmatch(r"gcc(?:-[0-9]+)?", cc_name) is None or re.fullmatch(
            r"g\+\+(?:-[0-9]+)?", cxx_name
        ) is None:
            raise ValueError("active GCC compiler paths are invalid")
        if not re.search(r"(?:gcc|gnu)", cc_version) or not re.search(
            r"(?:g\+\+|gcc|gnu)", cxx_version
        ):
            raise ValueError("active GCC compiler versions are invalid")
    else:
        raise ValueError("unsupported active compiler family")
    if major is not None and (
        re.search(rf"\b{major}(?:\.|\b)", cc_version) is None
        or re.search(rf"\b{major}(?:\.|\b)", cxx_version) is None
    ):
        raise ValueError("active compiler major version is invalid")


def validate_runtime_evidence(
    document,
    *,
    manifest,
    environment_id,
    source_sha,
    event_sha,
    tarball_sha256,
):
    if not isinstance(document, Mapping) or set(document) != FIELDS:
        raise ValueError("runtime evidence fields are not closed")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("schema_version must be integer 1")
    if document["kind"] != "r-runtime-proof" or document["status"] != "pass":
        raise ValueError("runtime evidence is not a passing R runtime proof")
    rows = [row for row in manifest["coverage"] if row["id"] == environment_id]
    if len(rows) != 1:
        raise ValueError("runtime environment must select exactly one manifest row")
    row = rows[0]
    expected = {
        "environment_id": environment_id,
        "runtime_profile": row.get("runtime_profile"),
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "r_executable": row.get("system_r"),
        "os": row.get("expected_os"),
        "distribution_id": row.get("expected_distribution"),
        "distribution_version": row.get("expected_distribution_version"),
        "architecture": row.get("expected_architecture"),
    }
    if any(value is None for value in expected.values()):
        raise ValueError("manifest row lacks a closed runtime binding")
    for field, value in expected.items():
        if document[field] != value:
            raise ValueError(f"runtime {field} does not match manifest/caller")
    _single_line(document["r_resolved"], "r_resolved")
    if not document["r_resolved"].startswith("/"):
        raise ValueError("runtime r_resolved must be absolute")
    if document["r_resolved"] != row.get("system_r"):
        if re.fullmatch(
            r"/opt/R/[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9._-]+)?/bin/R",
            document["r_resolved"],
        ) is None:
            raise ValueError("runtime resolved R is not an approved manifest alias")
    for field in (
        "r_executable", "r_resolved", "r_version", "r_platform", "os",
        "os_release", "distribution_id", "distribution_version", "architecture",
        "compiler", "cc", "cc_path", "cc_version", "cxx", "cxx_path",
        "cxx_version", "fc", "fc_path", "fc_version", "locale",
        "session_info", "la_library", "la_version", "blas_libs",
    ):
        _single_line(document[field], field)
    _single_line(document["r_status"], "r_status", allow_empty=True)
    expected_kind = row["expected_r_kind"]
    status_lower = document["r_status"].lower()
    if (
        _classify_r(document) != expected_kind
        or expected_kind == "devel" and "under development" not in status_lower
        or expected_kind == "patched" and "patched" not in status_lower
        or expected_kind == "release" and status_lower
    ):
        raise ValueError("runtime R release kind does not match manifest")

    profile = document["runtime_profile"]
    for field in ("cc_path", "cxx_path", "fc_path"):
        if not document[field].startswith("/"):
            raise ValueError(f"runtime {field} is not an absolute compiler path")
    if profile == "clang22":
        if not re.search(r"clang(?:\+\+)?-?22(?:\b|$)", document["cc"], re.I):
            raise ValueError("clang22 runtime does not prove clang 22")
        if not re.search(r"clang(?:\+\+)?-?22(?:\b|$)", document["cxx"], re.I):
            raise ValueError("clang22 runtime does not prove clang++ 22")
        if not re.search(r"flang(?:-new)?-?22(?:\b|$)", document["fc"], re.I):
            raise ValueError("clang22 runtime does not prove flang 22")
        for field, pattern in (
            ("cc_path", r"clang(?:-?22)?$"),
            ("cxx_path", r"clang(?:\+\+)?(?:-?22)?$"),
            ("fc_path", r"flang(?:-new)?(?:-?22)?$"),
        ):
            if re.search(pattern, document[field], re.I) is None:
                raise ValueError("clang22 runtime compiler executable/version is wrong")
        for field, family in (
            ("cc_version", "clang"), ("cxx_version", "clang"),
            ("fc_version", "flang"),
        ):
            value = document[field]
            if family not in value.lower() or re.search(r"\b22(?:\.|\b)", value) is None:
                raise ValueError("clang22 runtime compiler executable/version is wrong")
    elif profile == "ubuntu-gcc16":
        if not re.search(r"gcc-?16(?:\b|$)", document["cc"], re.I):
            raise ValueError("ubuntu-gcc16 runtime does not prove GCC 16")
        if re.search(r"/gcc-?16$", document["cc_path"], re.I) is None:
            raise ValueError("ubuntu-gcc16 runtime compiler executable is wrong")
        if "gcc" not in document["cc_version"].lower() or re.search(
            r"\b16(?:\.|\b)", document["cc_version"]
        ) is None:
            raise ValueError("ubuntu-gcc16 runtime compiler version is wrong")
    elif profile == "gcc16":
        if "gcc" not in document["cc"].lower() or "gcc" not in document["cc_version"].lower():
            raise ValueError("gcc16 runtime does not prove GCC")
        if re.search(r"\b16(?:\.|\b)", document["cc_version"]) is None:
            raise ValueError("gcc16 runtime compiler version is wrong")

    dimension = document["matrix_dimension"]
    checksum = document["matrix_checksum"]
    if type(dimension) is not int or dimension < 2:
        raise ValueError("matrix operation dimension is invalid")
    if type(checksum) not in {int, float} or isinstance(checksum, bool) or not math.isfinite(checksum) or checksum == 0:
        raise ValueError("matrix operation checksum is invalid")
    if document["maps_method"] != "proc-self-maps":
        raise ValueError("runtime libraries were not captured from this R process")
    libraries = document["loaded_libraries"]
    if (
        not isinstance(libraries, list) or not libraries
        or any(not isinstance(item, str) or not item.startswith("/") or "\x00" in item for item in libraries)
        or len(libraries) != len(set(libraries))
    ):
        raise ValueError("loaded library evidence is invalid")
    if type(document["long_double"]) is not bool:
        raise ValueError("long_double must be boolean")
    if not isinstance(document["ext_soft_version"], Mapping):
        raise ValueError("ext_soft_version must be an object")
    makeconf = document["makeconf"]
    if not isinstance(makeconf, Mapping) or set(makeconf) != MAKECONF_FIELDS or any(
        not isinstance(value, str) or any(character in value for character in "\r\n\x00")
        for value in makeconf.values()
    ):
        raise ValueError("runtime Makeconf evidence is invalid")
    environment = document["environment"]
    if not isinstance(environment, Mapping) or set(environment) != ENVIRONMENT_FIELDS or any(
        not isinstance(key, str)
        or value is not None and (
            not isinstance(value, str) or any(character in value for character in "\r\n\x00")
        )
        for key, value in environment.items()
    ):
        raise ValueError("runtime environment evidence is invalid")

    lowered = "\n".join(libraries).lower()
    if profile == "atlas":
        if re.search(r"lib(?:mkl|openblas|blis)", lowered):
            raise ValueError("ATLAS runtime loaded a forbidden BLAS")
        if re.search(r"lib(?:s?atlas)(?:\.so(?:\.\d+)*)?(?:$|\s|/)", lowered) is None:
            raise ValueError("ATLAS runtime did not load libsatlas")
        if "atlas" not in document["la_library"].lower():
            raise ValueError("ATLAS La_library identity is missing")
        if "atlas" not in str(document["ext_soft_version"].get("BLAS", "")).lower():
            raise ValueError("ATLAS extSoftVersion BLAS identity is missing")
        if "atlas" not in document["blas_libs"].lower():
            raise ValueError("ATLAS BLAS_LIBS identity is missing")
    elif profile in {"clang-asan", "clang-ubsan"}:
        _require_active_compiler(document, "clang", major=22)
        cc = document["cc"]
        cxx = document["cxx"]
        required = (
            {"-fsanitize=address,undefined", "-fno-sanitize=float-divide-by-zero", "-fno-sanitize=alignment", "-fno-omit-frame-pointer", "-fsanitize=pointer-overflow", "-fsanitize=signed-integer-overflow"}
            if profile == "clang-asan"
            else {"-fsanitize=undefined", "-fno-sanitize=function", "-fno-omit-frame-pointer"}
        )
        if not required <= set(cc.split()) or not required <= set(cxx.split()):
            raise ValueError("Clang sanitizer compiler flags do not match the active image")
        runtime_name = "libclang_rt.asan" if profile == "clang-asan" else "libclang_rt.ubsan"
        if runtime_name not in lowered:
            raise ValueError("Clang sanitizer runtime is not loaded in the selected R process")
        if profile == "clang-ubsan" and "ubsan_standalone" not in makeconf["SAN_LIBS"]:
            raise ValueError("Clang UBSAN linker configuration is missing")
        if environment["UBSAN_OPTIONS"] != "print_stacktrace=1":
            raise ValueError("Clang UBSAN runtime options do not match")
        if environment["ASAN_OPTIONS"] != "detect_leaks=0:alloc_dealloc_mismatch=0":
            raise ValueError("Clang ASAN runtime options do not match")
    elif profile in {"gcc-asan", "gcc-ubsan"}:
        _require_active_compiler(document, "gcc")
        required_flags = {"-fsanitize=address,undefined,bounds-strict", "-fno-omit-frame-pointer"}
        if not required_flags <= set(document["cc"].split()) or not required_flags <= set(document["cxx"].split()):
            raise ValueError("GCC sanitizer compiler flags do not match the active image")
        if "-fsanitize=address" not in makeconf["CFLAGS"] or "-fsanitize=address,undefined" not in makeconf["MAIN_LDFLAGS"]:
            raise ValueError("GCC sanitizer compiler/linker configuration is missing")
        preload = (environment["LD_PRELOAD"] or "").split(":")
        if set(preload) != {"/usr/lib64/libasan.so.8", "/usr/lib64/libubsan.so.1"}:
            raise ValueError("GCC sanitizer preload does not match the active image")
        if "libasan.so" not in lowered or "libubsan.so" not in lowered:
            raise ValueError("GCC sanitizer runtimes are not loaded in the selected R process")
        if environment["ASAN_OPTIONS"] != "detect_leaks=0" or environment["UBSAN_OPTIONS"] != "print_stacktrace=1":
            raise ValueError("GCC sanitizer runtime options do not match")
    elif profile == "donttest":
        if environment["_R_CHECK_DONTTEST_EXAMPLES_"] != "true":
            raise ValueError("donttest policy is not enabled")
    elif profile == "nold":
        if document["long_double"] is not False:
            raise ValueError("noLD runtime still reports long-double capability")
    elif profile == "valgrind":
        if environment["VALGRIND_OPTS"] != "--track-origins=yes --leak-check=full":
            raise ValueError("Valgrind runtime options do not match the active image")

    if profile in {"atlas", "clang-asan", "clang-ubsan", "donttest", "gcc-asan", "gcc-ubsan", "nold", "valgrind", "vnu"}:
        if environment["CHECK_ARGS"] != " ".join(row["check_args"]):
            raise ValueError("runtime CHECK_ARGS do not match the manifest")
    return document


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    args = parser.parse_args(argv)
    try:
        path = Path(args.evidence)
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
            raise ValueError("evidence must be a regular non-symlink file")
        manifest = load_manifest(args.manifest)
        validate_manifest(manifest)
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
        validate_runtime_evidence(
            document,
            manifest=manifest,
            environment_id=args.environment_id,
            source_sha=args.source_sha,
            event_sha=args.event_sha,
            tarball_sha256=args.tarball_sha256,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        print(f"runtime evidence error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(document, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
