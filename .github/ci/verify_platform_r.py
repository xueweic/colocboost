#!/usr/bin/env python3
"""Bind setup-r selection to a committed platform row and actual R identity."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import stat
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from artifact_contract import _atomic_write_json, append_github_outputs
from validate_manifest import load_manifest, validate_manifest


_SELECTOR = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_IDENTITY_FIELDS = (
    "version", "version_string", "status", "svn_revision", "platform", "os",
    "os_release", "arch", "r_home", "compiler", "cc", "cc_path",
    "cc_version", "cxx", "fc", "locale",
)
_MACOS_R = "/Library/Frameworks/R.framework/Resources/bin/R"
_WINDOWS_R = "C:/R/bin/R.exe"
_R_IDENTITY_EXPRESSION = r'''
value <- function(x) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) "<none>" else as.character(x[[1L]])
r_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
config <- function(name) {
  output <- tryCatch(
    suppressWarnings(system2(r_bin, c("CMD", "config", name), stdout = TRUE, stderr = TRUE)),
    error = function(error) character()
  )
  status <- attr(output, "status")
  if (!length(output) || !is.null(status) && status != 0L) "<none>" else paste(output, collapse = " ")
}
rtools_executable <- function(executable) {
  if (.Platform$OS.type != "windows") return("")
  environment_names <- grep("^RTOOLS[0-9]+_HOME$", names(Sys.getenv()), value = TRUE)
  roots <- unique(unname(Sys.getenv(environment_names, unset = "")))
  roots <- roots[nzchar(roots)]
  candidates <- unique(unlist(lapply(roots, function(root) {
    Sys.glob(file.path(root, "*", "bin", executable))
  }), use.names = FALSE))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) "" else normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}
cc <- config("CC")
cc_command <- if (identical(cc, "<none>")) "" else strsplit(trimws(cc), "[[:space:]]+")[[1L]][[1L]]
cc_path <- if (nzchar(cc_command)) unname(Sys.which(cc_command)) else ""
if (!nzchar(cc_path)) cc_path <- rtools_executable("gcc.exe")
if (identical(cc, "<none>") && nzchar(cc_path)) cc <- basename(cc_path)
cc_version <- if (nzchar(cc_path)) {
  output <- system2(cc_path, "--version", stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || !length(output)) "<none>" else output[[1L]]
} else "<none>"
cxx <- config("CXX")
if (identical(cxx, "<none>")) {
  cxx_path <- rtools_executable("g++.exe")
  if (nzchar(cxx_path)) cxx <- basename(cxx_path)
}
fc <- config("FC")
if (identical(fc, "<none>")) {
  fc_path <- rtools_executable("gfortran.exe")
  if (nzchar(fc_path)) fc <- basename(fc_path)
}
fields <- c(
  version = paste(R.version$major, R.version$minor, sep = "."),
  version_string = R.version.string,
  status = value(R.version$status),
  svn_revision = value(R.version$svn.rev),
  platform = R.version$platform,
  os = unname(Sys.info()[["sysname"]]),
  os_release = unname(Sys.info()[["release"]]),
  arch = value(R.version$arch),
  r_home = R.home(),
  compiler = cc_version,
  cc = cc,
  cc_path = value(cc_path),
  cc_version = value(cc_version),
  cxx = cxx,
  fc = fc,
  locale = Sys.getlocale()
)
cat(paste(names(fields), fields, sep = "="), sep = "\n")
'''.strip()


def _single_line(value: str, label: str) -> str:
    if not isinstance(value, str) or not value or any(c in value for c in "\r\n\x00"):
        raise ValueError(f"{label} must be an unambiguous nonempty single-line string")
    return value


def _portable_path(value: str) -> str:
    return value.replace("\\", "/")


def _portable_absolute(value: str) -> bool:
    portable = _portable_path(value)
    return portable.startswith("/") or re.fullmatch(r"[A-Za-z]:/.*", portable) is not None


def _require_r_executable(value: str | os.PathLike[str]) -> tuple[Path, Path]:
    raw = _single_line(os.fspath(value), "R executable path")
    portable = _portable_path(raw)
    if not (Path(raw).is_absolute() or re.fullmatch(r"[A-Za-z]:/.*", portable)):
        raise ValueError("R executable path must be absolute")
    lexical = Path(raw)
    try:
        status = lexical.lstat()
    except OSError as error:
        raise ValueError("R executable must be an existing regular file") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError("R executable must be a regular non-symlink file")
    if not os.access(lexical, os.X_OK):
        raise ValueError("R executable must be executable")
    if lexical.name.lower() not in {"r", "r.exe"}:
        raise ValueError("R executable basename must be R or R.exe")
    resolved = lexical.resolve(strict=True)
    if not stat.S_ISREG(resolved.stat().st_mode) or not os.access(resolved, os.X_OK):
        raise ValueError("resolved R executable must be a regular executable file")
    return lexical, resolved


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_identity(stdout: str) -> dict[str, str]:
    if not isinstance(stdout, str) or any(c in stdout for c in "\r\x00"):
        raise ValueError("R identity output is ambiguous")
    identity: dict[str, str] = {}
    for line in stdout.splitlines():
        if "=" not in line:
            raise ValueError("R identity output does not match the closed protocol")
        key, value = line.split("=", 1)
        if key not in _IDENTITY_FIELDS or key in identity or not value:
            raise ValueError("R identity output does not match the closed protocol")
        identity[key] = value
    if tuple(identity) != _IDENTITY_FIELDS:
        raise ValueError("R identity output fields are missing, reordered, or duplicated")
    return identity


def _normal_os(value: str) -> str:
    aliases = {"linux": "linux", "windows": "windows", "mingw32": "windows", "darwin": "macos", "macos": "macos"}
    try:
        return aliases[value.casefold()]
    except KeyError as error:
        raise ValueError(f"unsupported operating system identity: {value!r}") from error


def _normal_architecture(value: str) -> str:
    aliases = {"amd64": "x86_64", "x86_64": "x86_64", "x64": "x86_64", "arm64": "aarch64", "aarch64": "aarch64"}
    try:
        return aliases[value.casefold()]
    except KeyError as error:
        raise ValueError(f"unsupported R architecture identity: {value!r}") from error


def _r_kind(identity: Mapping[str, str]) -> str:
    combined = f"{identity['status']} {identity['version_string']}".lower()
    if "under development" in combined:
        return "devel"
    if "patched" in combined:
        return "patched"
    return "release"


def _has_development_revision(identity: Mapping[str, str]) -> bool:
    return bool(re.fullmatch(r"[0-9]+", identity["svn_revision"])) or bool(
        re.search(r"\br[0-9]{4,}\b", identity["version_string"])
    )


def verify_platform_r(
    manifest_row: Mapping[str, Any],
    *,
    environment_id: str,
    r_executable: str | os.PathLike[str],
    selected_selector: str,
    setup_r_version: str,
    runner_label: str,
    evidence_path: str | os.PathLike[str] | None = None,
    github_output: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    if manifest_row.get("id") != environment_id or manifest_row.get("driver") != "r-binary":
        raise ValueError("environment ID must select the committed r-binary row")
    selected = _single_line(selected_selector, "selected setup-r selector")
    if _SELECTOR.fullmatch(selected) is None or selected != manifest_row.get("setup_r_selector"):
        raise ValueError("selected setup-r selector does not match manifest")
    runner = _single_line(runner_label, "runner label")
    if runner != manifest_row.get("runner"):
        raise ValueError("runner label does not match manifest")
    setup_version = _single_line(setup_r_version, "setup-r output version")
    lexical, resolved = _require_r_executable(r_executable)
    lexical_portable = _portable_path(os.fspath(lexical))
    declared_r = manifest_row.get("system_r")
    if declared_r is None or lexical_portable.casefold() != declared_r.casefold():
        raise ValueError("R executable path does not match manifest system_r")
    expected_os = manifest_row.get("expected_os")
    expected_arch = manifest_row.get("expected_architecture")
    if expected_os == "macos" and lexical_portable != _MACOS_R:
        raise ValueError(f"macOS R must be {_MACOS_R}")
    if expected_os == "windows" and lexical_portable.casefold() != _WINDOWS_R.casefold():
        raise ValueError(f"Windows R must be {_WINDOWS_R}")

    completed = subprocess.run(
        [os.fspath(lexical), "--vanilla", "--slave", "-e", _R_IDENTITY_EXPRESSION],
        check=False, capture_output=True, text=True, shell=False,
    )
    if completed.returncode != 0:
        detail = " ".join(completed.stderr.split())[:500]
        suffix = f": {detail}" if detail else ""
        raise ValueError(
            f"selected R identity probe failed with exit code "
            f"{completed.returncode}{suffix}"
        )
    identity = _parse_identity(completed.stdout)
    if identity["version"] != setup_version:
        raise ValueError("setup-r output version does not match selected R version")
    actual_os = _normal_os(identity["os"])
    actual_arch = _normal_architecture(identity["arch"])
    if actual_os != expected_os or actual_arch != expected_arch:
        raise ValueError("selected R OS or architecture does not match manifest")
    actual_kind = _r_kind(identity)
    if actual_kind != manifest_row.get("expected_r_kind"):
        raise ValueError("selected R release kind does not match manifest")
    if actual_kind == "devel" and not _has_development_revision(identity):
        raise ValueError("selected R development revision proof is missing")
    if identity["os_release"] == "<none>" or identity["locale"] == "<none>":
        raise ValueError("selected R OS release or locale proof is missing")
    if identity["cc"] == "<none>" or identity["cxx"] == "<none>" or identity["fc"] == "<none>":
        raise ValueError("selected R compiler proof is missing")
    if identity["cc_path"] == "<none>" or not _portable_absolute(identity["cc_path"]):
        raise ValueError("selected R compiler executable proof is missing")
    if identity["cc_version"] == "<none>":
        raise ValueError("selected R compiler version proof is missing")
    cc_path = _portable_path(identity["cc_path"]).casefold()
    cc_version = identity["cc_version"].casefold()
    if expected_os == "windows" and not (
        "rtools" in cc_path
        and re.search(r"/gcc(?:\.exe)?$", cc_path) is not None
        and "gcc" in cc_version
    ):
        raise ValueError("Windows selected R does not prove Rtools GCC")
    if expected_os == "macos" and not (
        re.search(r"/clang$", cc_path) is not None and "clang" in cc_version
    ):
        raise ValueError("macOS selected R does not prove Apple/Clang compiler")
    home = _portable_path(identity["r_home"]).rstrip("/")
    expected_from_home = f"{home}/bin/{'R.exe' if expected_os == 'windows' else 'R'}"
    resolved_portable = _portable_path(os.fspath(resolved))
    if expected_from_home.casefold() not in {
        lexical_portable.casefold(), resolved_portable.casefold(),
    }:
        raise ValueError("R.home does not bind the selected R executable")

    status = lexical.stat()
    proof: dict[str, Any] = {
        "schema_version": 1,
        "kind": "platform-r-identity-proof",
        "status": "pass",
        "environment_id": environment_id,
        "runner": runner,
        "selector": selected,
        "setup_r_version": setup_version,
        "expected_r_kind": manifest_row["expected_r_kind"],
        "expected_os": expected_os,
        "expected_architecture": expected_arch,
        "r_executable": {
            "path": lexical_portable,
            "resolved_path": resolved_portable,
            "size": status.st_size,
            "sha256": _sha256(lexical),
        },
        "identity": identity,
    }
    if evidence_path is not None:
        _atomic_write_json(Path(evidence_path), proof)
    if github_output is not None:
        append_github_outputs(github_output, {
            "r_executable": lexical_portable,
            "r_resolved_executable": resolved_portable,
            "r_version": identity["version"],
            "r_os": actual_os,
            "r_architecture": actual_arch,
        })
    return proof


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--r-executable", required=True)
    parser.add_argument("--selected-selector", required=True)
    parser.add_argument("--setup-r-version", required=True)
    parser.add_argument("--runner-label", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--github-output")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        manifest = validate_manifest(load_manifest(args.manifest))
        rows = [row for row in manifest["coverage"] if row["id"] == args.environment_id]
        if len(rows) != 1:
            raise ValueError("environment ID must select exactly one manifest row")
        verify_platform_r(
            rows[0], environment_id=args.environment_id,
            r_executable=args.r_executable, selected_selector=args.selected_selector,
            setup_r_version=args.setup_r_version, runner_label=args.runner_label,
            evidence_path=args.evidence, github_output=args.github_output,
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"platform R verification error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
