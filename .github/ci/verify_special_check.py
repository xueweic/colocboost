#!/usr/bin/env python3
"""Fail-closed semantic validation for CRAN check trees and special lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

from artifact_contract import append_github_outputs
from validate_manifest import load_manifest, validate_manifest
from verify_runtime_evidence import validate_runtime_evidence


FULL_DOCUMENTATION = "full-documentation"
NO_DOCUMENTATION = "standard-no-documentation"
ACTIVE_PROFILES = {
    "openblas", "rcnst",
    "atlas", "clang-asan", "clang-ubsan", "donttest", "gcc-asan",
    "gcc-ubsan", "nold", "valgrind", "vnu",
}
PROFILES = {FULL_DOCUMENTATION, NO_DOCUMENTATION, *ACTIVE_PROFILES}
RECONCILABLE = {"clang-asan", "clang-ubsan", "valgrind"}
DOC_HEADINGS = {
    "manual": re.compile(r"^\* checking PDF version of manual \.\.\.(?: \[[^ ]+\])? OK$"),
    "package": re.compile(r"^\* checking package vignettes \.\.\.(?: \[[^ ]+\])? OK$"),
    "rebuilding": re.compile(r"^\* checking re-building of vignette outputs \.\.\.(?: \[[^ ]+\])? OK$"),
}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _read_json(path_value, label):
    path = Path(path_value)
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is missing") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file")
    try:
        document = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is malformed") from error
    if not isinstance(document, Mapping):
        raise ValueError(f"{label} must be an object")
    return document


def _scan_tree(root_value):
    root = Path(root_value)
    if not root.is_absolute():
        raise ValueError("check root must be absolute")
    try:
        root_status = root.lstat()
    except OSError as error:
        raise ValueError("check root is missing") from error
    if stat.S_ISLNK(root_status.st_mode) or not stat.S_ISDIR(root_status.st_mode):
        raise ValueError("check root must be a regular non-symlink directory")
    for current, directories, files in os.walk(root, followlinks=False):
        for name in [*directories, *files]:
            path = Path(current) / name
            status = path.lstat()
            if stat.S_ISLNK(status.st_mode):
                raise ValueError("check tree must not contain symlinks")
            if name in directories and not stat.S_ISDIR(status.st_mode):
                raise ValueError("check tree contains a non-directory entry")
            if name in files and not stat.S_ISREG(status.st_mode):
                raise ValueError("check tree contains a non-regular file")
    candidates = [entry for entry in root.iterdir() if entry.name.endswith(".Rcheck")]
    if len(candidates) != 1:
        raise ValueError("check root must contain exactly one *.Rcheck directory")
    check = candidates[0]
    if not check.is_dir():
        raise ValueError("*.Rcheck must be a directory")
    log = check / "00check.log"
    try:
        status = log.lstat()
    except OSError as error:
        raise ValueError("check tree lacks 00check.log") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("00check.log must be a regular non-symlink file")
    return check, log


def _tree_outputs(check, log):
    outputs = []
    for path in sorted(check.rglob("*")):
        if not path.is_file():
            continue
        name = path.name
        data = path.read_bytes()
        strict_text_output = (
            path == log
            or name == "00install.out"
            or name.endswith("-Ex.Rout")
            or re.search(r"\.Rout(?:\..+)?$", name) is not None
        )
        if strict_text_output:
            try:
                text = data.decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError(f"check output {path.relative_to(check)} is not UTF-8") from error
        else:
            text = data.decode("utf-8", errors="replace")
        outputs.append((path, data, text))
    return outputs


def _heading_count(lines, stem):
    pattern = re.compile(
        rf"^\* {re.escape(stem)} \.\.\.(?: \[[^ ]+\])? OK$"
    )
    return sum(pattern.fullmatch(line) is not None for line in lines)


def _installation_heading_count(lines):
    pattern = re.compile(
        r"^\* checking whether package(?: (?:'[^']+'|‘[^’]+’))? can be installed \.\.\.(?: \[[^ ]+\])? OK$"
    )
    return sum(pattern.fullmatch(line) is not None for line in lines)


def _validate_native(native, row, log):
    raw_exit = native.get("wrapper_exit_code")
    if type(raw_exit) is not int or raw_exit < 0:
        raise ValueError("native wrapper exit is invalid")
    executable = native.get("executable")
    required_r = native.get("required_r_executable")
    if not isinstance(executable, Mapping) or not isinstance(required_r, Mapping):
        raise ValueError("native executable identities are missing")
    expected_hash = row.get("wrapper_sha256")
    if (
        executable.get("path") != row.get("wrapper_path")
        or executable.get("sha256") != expected_hash
        or native.get("wrapper_sha256") != expected_hash
    ):
        raise ValueError("native wrapper identity does not match manifest")
    if (
        required_r.get("path") != row.get("system_r")
        or native.get("path_r_resolution") != row.get("system_r")
    ):
        raise ValueError("native system R identity does not match manifest")
    if native.get("check_args") != " ".join(row.get("check_args", [])):
        raise ValueError("native CHECK_ARGS do not match manifest")
    if native.get("check_log") != os.fspath(log):
        raise ValueError("native check log path does not match the unique check tree")
    return raw_exit


def _scan_sanitizers(outputs):
    pattern = re.compile(
        r"(?:ERROR:\s*(?:Address|Leak)Sanitizer|SUMMARY:\s*(?:Address|Leak)Sanitizer|UndefinedBehaviorSanitizer|runtime error:)",
        re.IGNORECASE,
    )
    for path, _, text in outputs:
        if pattern.search(text):
            raise ValueError(f"sanitizer diagnostic found in {path.name}")


def _scan_valgrind(outputs):
    harmful = re.compile(
        r"(?:Invalid (?:read|write|free)|Conditional jump or move depends on uninitialised|Use of uninitialised|(?:definitely|indirectly|possibly) lost:\s*[1-9][0-9,]*|ERROR SUMMARY:\s*[1-9][0-9,]* errors)",
        re.IGNORECASE,
    )
    zero_summaries = 0
    for path, _, text in outputs:
        if harmful.search(text):
            raise ValueError(f"Valgrind diagnostic found in {path.name}")
        zero_summaries += len(re.findall(r"ERROR SUMMARY:\s*0 errors", text, re.I))
    if zero_summaries == 0:
        raise ValueError("Valgrind check has no clean instrumented transcript")
    return zero_summaries


def _count_donttest_blocks(check):
    source_root = check / "00_pkg_src"
    try:
        source_status = source_root.lstat()
    except OSError as error:
        raise ValueError("donttest check lacks 00_pkg_src source evidence") from error
    if stat.S_ISLNK(source_status.st_mode) or not stat.S_ISDIR(source_status.st_mode):
        raise ValueError("donttest 00_pkg_src must be a regular directory")
    packages = []
    for candidate in source_root.iterdir():
        candidate_status = candidate.lstat()
        if stat.S_ISLNK(candidate_status.st_mode):
            raise ValueError("donttest source package must not be a symlink")
        if stat.S_ISDIR(candidate_status.st_mode):
            packages.append(candidate)
    if len(packages) != 1:
        raise ValueError("donttest check must retain exactly one source package")
    man = packages[0] / "man"
    try:
        man_status = man.lstat()
    except OSError as error:
        raise ValueError("donttest source package lacks man directory") from error
    if stat.S_ISLNK(man_status.st_mode) or not stat.S_ISDIR(man_status.st_mode):
        raise ValueError("donttest man evidence must be a regular directory")
    rd_files = sorted(man.glob("*.Rd"))
    if not rd_files:
        raise ValueError("donttest source package has no Rd files to audit")
    count = 0
    for path in rd_files:
        path_status = path.lstat()
        if stat.S_ISLNK(path_status.st_mode) or not stat.S_ISREG(path_status.st_mode):
            raise ValueError("donttest Rd evidence must be regular files")
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeError as error:
            raise ValueError("donttest Rd evidence must be UTF-8") from error
        count += text.count("\\donttest{")
    return count


def _atomic_write(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(document, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def verify_special_check(
    manifest_row,
    *,
    environment_id,
    profile,
    check_root,
    native_evidence,
    runtime_evidence,
    source_sha,
    event_sha,
    tarball_sha256,
    output,
    github_output=None,
):
    if profile not in PROFILES:
        raise ValueError("unsupported special check profile")
    if manifest_row.get("id") != environment_id:
        raise ValueError("environment ID does not match manifest")
    expected_profile = manifest_row.get("runtime_profile") if environment_id in ACTIVE_PROFILES else (
        FULL_DOCUMENTATION if manifest_row.get("check_args") == [] else NO_DOCUMENTATION
    )
    if profile != expected_profile:
        raise ValueError("check profile does not match manifest check arguments")
    native = _read_json(native_evidence, "native evidence")
    runtime = _read_json(runtime_evidence, "runtime evidence")
    for field, expected in {
        "environment_id": environment_id,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
    }.items():
        if runtime.get(field) != expected:
            raise ValueError(f"runtime evidence {field} does not match")
    if environment_id in ACTIVE_PROFILES:
        validate_runtime_evidence(
            runtime, manifest={"coverage": [manifest_row]},
            environment_id=environment_id, source_sha=source_sha,
            event_sha=event_sha, tarball_sha256=tarball_sha256,
        )

    check, log = _scan_tree(check_root)
    raw_exit = _validate_native(native, manifest_row, log)
    log_bytes = log.read_bytes()
    try:
        text = log_bytes.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("00check.log must be UTF-8") from error
    if not text.endswith("* DONE\nStatus: OK\n"):
        raise ValueError("00check.log does not end in a strict successful footer")
    lines = text.splitlines()
    critical = {
        "installation": _installation_heading_count(lines),
        "tests": _heading_count(lines, "checking tests"),
        "examples": _heading_count(lines, "checking examples"),
    }
    if environment_id in ACTIVE_PROFILES and any(count != 1 for count in critical.values()):
        raise ValueError("critical package check stages were not each completed exactly once")
    matched = {
        name: sum(bool(pattern.fullmatch(line)) for line in lines)
        for name, pattern in DOC_HEADINGS.items()
    }
    documentation_policy = manifest_row.get("documentation_policy") if environment_id in ACTIVE_PROFILES else (
        "required" if profile == FULL_DOCUMENTATION else "forbidden"
    )
    if documentation_policy == "required":
        if any(count != 1 for count in matched.values()):
            raise ValueError("full documentation stages were not each executed exactly once")
    elif documentation_policy == "forbidden" and (
        matched["manual"] or matched["rebuilding"]
    ):
        raise ValueError("no-documentation profile unexpectedly executed a doc stage")

    outputs = _tree_outputs(check, log)
    scan_records = [
        {
            "path": os.fspath(path.relative_to(check)),
            "size": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        for path, data, _ in outputs
    ]
    if profile in {"clang-asan", "clang-ubsan", "gcc-asan", "gcc-ubsan"}:
        _scan_sanitizers(outputs)
    valgrind_zero_summaries = _scan_valgrind(outputs) if profile == "valgrind" else 0

    donttest_policy = False
    donttest_block_count = None
    examples_executed = critical["examples"] == 1
    if profile == "donttest":
        donttest_policy = runtime["environment"].get("_R_CHECK_DONTTEST_EXAMPLES_") == "true"
        if not donttest_policy:
            raise ValueError("donttest policy was not enabled in the selected R process")
        for name in ("colocboost-Ex.R", "colocboost-Ex.Rout"):
            candidate = check / name
            try:
                status = candidate.lstat()
            except OSError as error:
                raise ValueError(f"donttest check lacks nonempty {name}") from error
            if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode) or status.st_size <= 0:
                raise ValueError(f"donttest {name} must be a nonempty regular file")
        if not examples_executed:
            raise ValueError("donttest examples stage was not completed")
        donttest_block_count = _count_donttest_blocks(check)
        if donttest_block_count != 0:
            raise ValueError("donttest source unexpectedly contains donttest blocks")

    reconciliation = None
    effective_exit = raw_exit
    if raw_exit != 0:
        if raw_exit != 1 or profile not in RECONCILABLE:
            raise ValueError("native wrapper exit cannot be reconciled")
        tests = check / "tests"
        try:
            tests_status = tests.lstat()
        except OSError as error:
            raise ValueError("known wrapper defect requires a tests directory") from error
        if stat.S_ISLNK(tests_status.st_mode) or not stat.S_ISDIR(tests_status.st_mode):
            raise ValueError("known wrapper defect requires a regular tests directory")
        rout = [path for path in tests.rglob("*") if path.is_file() and re.search(r"\.Rout(?:\..+)?$", path.name)]
        failed = [path for path in rout if path.name.endswith(".Rout.fail")]
        if not rout or failed:
            raise ValueError("native exit 1 does not match the empty Rout.fail grep defect")
        reconciliation = "rhub-empty-rout-fail-grep"
        effective_exit = 0

    proof = {
        "schema_version": 1,
        "kind": "special-check-proof",
        "status": "pass",
        "environment_id": environment_id,
        "profile": profile,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": tarball_sha256,
        "check_directory": os.fspath(check),
        "check_log_sha256": hashlib.sha256(log_bytes).hexdigest(),
        "manual_executed": matched["manual"] == 1,
        "vignettes_executed": matched["package"] == matched["rebuilding"] == 1,
        "critical_stages": critical,
        "scanned_outputs": scan_records,
        "donttest_policy_enabled": donttest_policy,
        "donttest_block_count": donttest_block_count,
        "examples_executed": examples_executed,
        "valgrind_zero_error_summaries": valgrind_zero_summaries,
        "raw_exit_code": raw_exit,
        "effective_exit_code": effective_exit,
        "reconciliation": reconciliation,
    }
    _atomic_write(output, proof)
    if github_output is not None:
        append_github_outputs(github_output, {"effective_exit_code": str(effective_exit)})
    return proof


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--profile", required=True, choices=sorted(PROFILES))
    parser.add_argument("--check-root", required=True)
    parser.add_argument("--native-evidence", required=True)
    parser.add_argument("--runtime-evidence", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--tarball-sha256", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--github-output")
    args = parser.parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(args.manifest))
        rows = [row for row in manifest["coverage"] if row["id"] == args.environment_id]
        if len(rows) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        verify_special_check(
            rows[0], environment_id=args.environment_id, profile=args.profile,
            check_root=args.check_root, native_evidence=args.native_evidence,
            runtime_evidence=args.runtime_evidence, source_sha=args.source_sha,
            event_sha=args.event_sha, tarball_sha256=args.tarball_sha256,
            output=args.output, github_output=args.github_output,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"special check evidence error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
