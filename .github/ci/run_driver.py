"""Run an explicitly selected CRAN preflight driver without PATH fallback."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from artifact_contract import verify_tarball
from validate_manifest import load_manifest, validate_manifest


_DRIVERS = frozenset({"r-binary", "native-wrapper"})
_TARBALL_TOKEN = "{tarball}"
_TARBALL_PARENT_TOKEN = "{tarball-parent}"
_EVIDENCE_ENVIRONMENT_FIELDS = (
    "PATH",
    "R_HOME",
    "R_LIBS",
    "R_LIBS_USER",
    "R_LIBS_SITE",
    "LD_LIBRARY_PATH",
    "DYLD_LIBRARY_PATH",
    "LIBPATH",
    "SHLIB_PATH",
    "CC",
    "CXX",
    "FC",
    "F77",
    "LANG",
    "LC_ALL",
    "LC_CTYPE",
)


def _require_executable(
    value: str | os.PathLike[str], *, label: str = "executable"
) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"{label} path must be absolute: {path}")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} must be an existing regular executable file: {path}") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink executable file: {path}")
    if not os.access(path, os.X_OK):
        raise ValueError(f"{label} must be executable: {path}")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _dynamic_library_evidence() -> dict[str, Any]:
    maps_path = Path("/proc/self/maps")
    if not maps_path.is_file():
        return {"method": "unsupported", "libraries": []}
    libraries = set()
    for line in maps_path.read_text(encoding="utf-8", errors="replace").splitlines():
        candidate = line.rsplit(maxsplit=1)[-1]
        if candidate.startswith("/"):
            libraries.add(candidate)
    return {"method": "proc-self-maps", "libraries": sorted(libraries)}


def capture_environment_evidence(
    executable: str | os.PathLike[str],
    *,
    environment: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Capture a fixed allowlist of execution identity fields."""

    executable_path = _require_executable(executable)
    selected_environment = os.environ if environment is None else environment
    status = executable_path.stat()
    return {
        "environment": {
            field: selected_environment.get(field)
            for field in _EVIDENCE_ENVIRONMENT_FIELDS
        },
        "executable": {
            "path": str(executable_path),
            "resolved_path": str(executable_path.resolve(strict=True)),
            "size": status.st_size,
            "sha256": _sha256(executable_path),
        },
        "dynamic_libraries": _dynamic_library_evidence(),
    }


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(document, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _normalize_argv(argv: Sequence[str | os.PathLike[str]]) -> list[str]:
    if isinstance(argv, (str, bytes)) or not isinstance(argv, Sequence):
        raise ValueError("argv must be a sequence of arguments")
    normalized = []
    for argument in argv:
        if not isinstance(argument, (str, os.PathLike)):
            raise ValueError("every argv item must be a string or path")
        normalized.append(os.fspath(argument))
    return normalized


def _require_tarball_parent(tarball: Path) -> Path:
    parent = tarball.parent
    try:
        parent_status = parent.lstat()
    except OSError as error:
        raise ValueError(f"tarball parent does not exist: {parent}") from error
    if stat.S_ISLNK(parent_status.st_mode) or not stat.S_ISDIR(parent_status.st_mode):
        raise ValueError(f"tarball parent must be a regular non-symlink directory: {parent}")

    candidates = []
    try:
        entries = list(parent.iterdir())
    except OSError as error:
        raise ValueError(f"could not inspect tarball parent {parent}: {error}") from error
    for entry in entries:
        if entry.name.endswith(".tar.gz"):
            candidates.append(entry)
    if len(candidates) != 1:
        raise ValueError(
            "tarball parent must contain exactly one regular non-symlink *.tar.gz"
        )
    candidate = candidates[0]
    try:
        candidate_status = candidate.lstat()
    except OSError as error:
        raise ValueError(
            "tarball parent must contain exactly one regular non-symlink *.tar.gz"
        ) from error
    if stat.S_ISLNK(candidate_status.st_mode) or not stat.S_ISREG(candidate_status.st_mode):
        raise ValueError(
            "tarball parent must contain exactly one regular non-symlink *.tar.gz"
        )
    if candidate != tarball:
        raise ValueError("tarball parent does not contain the verified tarball")
    return parent


def _verify_native_r_binding(
    manifest_row: Mapping[str, Any],
    executable: Path,
    required_r_executable: str | os.PathLike[str] | None,
    environment: Mapping[str, str],
) -> tuple[Path, Path]:
    if required_r_executable is None:
        raise ValueError("required_r_executable is mandatory for native-wrapper")
    system_r = _require_executable(required_r_executable, label="required system R")
    if system_r.name not in {"R", "R.exe"}:
        raise ValueError("required system R executable basename must be R or R.exe")

    declared_wrapper = manifest_row.get("wrapper_path")
    if declared_wrapper is not None and os.fspath(executable) != declared_wrapper:
        raise ValueError("native wrapper executable does not match manifest wrapper_path")
    declared_r = manifest_row.get("system_r")
    if declared_r is not None and os.fspath(system_r) != declared_r:
        raise ValueError("required system R does not match manifest system_r")

    declared_wrapper_sha256 = manifest_row.get("wrapper_sha256")
    if declared_wrapper_sha256 is not None:
        if (
            not isinstance(declared_wrapper_sha256, str)
            or len(declared_wrapper_sha256) != 64
            or any(character not in "0123456789abcdef" for character in declared_wrapper_sha256)
        ):
            raise ValueError("manifest wrapper sha256 must be 64 lowercase hex characters")
        if _sha256(executable) != declared_wrapper_sha256:
            raise ValueError("native wrapper sha256 does not match manifest")

    declared_check_args = manifest_row.get("check_args")
    if declared_check_args is not None:
        if (
            not isinstance(declared_check_args, list)
            or any(not isinstance(item, str) or not item or any(c.isspace() for c in item)
                   for item in declared_check_args)
        ):
            raise ValueError("manifest check_args must be a list of single-token strings")
        actual_check_args = environment.get("CHECK_ARGS", "")
        if actual_check_args != " ".join(declared_check_args):
            raise ValueError("CHECK_ARGS does not match manifest check_args")

    path_value = environment.get("PATH")
    if not isinstance(path_value, str) or not path_value:
        raise ValueError("PATH must resolve the required system R for native-wrapper")
    selected = shutil.which("R", path=path_value)
    if selected is None:
        raise ValueError("PATH does not resolve R for native-wrapper")
    path_r = _require_executable(Path(os.path.abspath(selected)), label="PATH R")
    if path_r != system_r:
        raise ValueError(
            f"PATH R {path_r} is not the required system R {system_r}"
        )
    return system_r, path_r


def _build_command(
    driver: str,
    executable: Path,
    argv: Sequence[str | os.PathLike[str]],
    r_entrypoint: str | None,
) -> tuple[list[str], Path]:
    arguments = _normalize_argv(argv)
    if driver == "native-wrapper":
        if r_entrypoint is not None:
            raise ValueError("r_entrypoint is forbidden for native-wrapper")
        return [str(executable), *arguments], executable

    if executable.name not in {"R", "R.exe"}:
        raise ValueError("r-binary executable basename must be R or R.exe")
    if r_entrypoint == "r-cmd":
        return [str(executable), "CMD", *arguments], executable
    if r_entrypoint == "rscript":
        sibling_name = "Rscript.exe" if executable.name == "R.exe" else "Rscript"
        rscript = _require_executable(
            executable.with_name(sibling_name), label="Rscript sibling"
        )
        return [str(rscript), *arguments], rscript
    raise ValueError("r_entrypoint must be rscript or r-cmd for r-binary")


def run_driver(
    manifest_row: Mapping[str, Any],
    *,
    requested_driver: str,
    executable: str | os.PathLike[str],
    argv: Sequence[str | os.PathLike[str]],
    tarball: str | os.PathLike[str],
    metadata: str | os.PathLike[str],
    expected_source_sha: str,
    expected_event_sha: str,
    r_entrypoint: str | None = None,
    required_r_executable: str | os.PathLike[str] | None = None,
    evidence_path: str | os.PathLike[str] | None = None,
    environment: Mapping[str, str] | None = None,
    cwd: str | os.PathLike[str] | None = None,
) -> int:
    """Verify identity and execute exactly one declared external driver."""

    if not isinstance(manifest_row, Mapping):
        raise ValueError("manifest_row must be an object")
    manifest_driver = manifest_row.get("driver")
    if manifest_driver not in _DRIVERS:
        raise ValueError(f"unsupported driver in manifest row: {manifest_driver!r}")
    if requested_driver != manifest_driver:
        raise ValueError(
            f"driver mismatch: manifest declares {manifest_driver!r}, "
            f"caller requested {requested_driver!r}"
        )

    executable_path = _require_executable(executable)
    selected_environment = os.environ if environment is None else environment
    required_system_r = None
    path_r = None
    if manifest_driver == "native-wrapper":
        required_system_r, path_r = _verify_native_r_binding(
            manifest_row,
            executable_path,
            required_r_executable,
            selected_environment,
        )
    elif required_r_executable is not None:
        raise ValueError("required_r_executable is forbidden for r-binary")

    unbound_arguments = _normalize_argv(argv)
    token_count = sum(
        unbound_arguments.count(token)
        for token in (_TARBALL_TOKEN, _TARBALL_PARENT_TOKEN)
    )
    if token_count != 1:
        raise ValueError(
            "argv must contain exactly one literal {tarball} or {tarball-parent} placeholder"
        )
    if manifest_driver == "r-binary" and _TARBALL_PARENT_TOKEN in unbound_arguments:
        raise ValueError("r-binary requires the literal {tarball} placeholder")
    declared_input = manifest_row.get("wrapper_input")
    if declared_input is not None:
        expected_token = {
            "tarball": _TARBALL_TOKEN,
            "tarball-parent": _TARBALL_PARENT_TOKEN,
        }.get(declared_input)
        if expected_token is None or expected_token not in unbound_arguments:
            raise ValueError("argv does not match manifest wrapper_input")
    verified_tarball = Path(os.path.abspath(os.fspath(tarball)))
    verify_tarball(
        verified_tarball,
        metadata,
        expected_source_sha=expected_source_sha,
        expected_event_sha=expected_event_sha,
    )
    tarball_parent = None
    if _TARBALL_PARENT_TOKEN in unbound_arguments:
        tarball_parent = _require_tarball_parent(verified_tarball)
    bound_arguments = []
    for argument in unbound_arguments:
        if argument == _TARBALL_TOKEN:
            bound_arguments.append(str(verified_tarball))
        elif argument == _TARBALL_PARENT_TOKEN:
            bound_arguments.append(str(tarball_parent))
        else:
            bound_arguments.append(argument)
    command, invoked_executable = _build_command(
        manifest_driver, executable_path, bound_arguments, r_entrypoint
    )

    process_environment = None if environment is None else dict(environment)
    evidence = capture_environment_evidence(
        executable_path,
        environment=os.environ if process_environment is None else process_environment,
    )
    if invoked_executable != executable_path:
        invoked_evidence = capture_environment_evidence(
            invoked_executable,
            environment=(
                os.environ if process_environment is None else process_environment
            ),
        )
        evidence["invoked_executable"] = invoked_evidence["executable"]
    if required_system_r is not None:
        r_evidence = capture_environment_evidence(
            required_system_r,
            environment=selected_environment,
        )
        evidence["required_r_executable"] = r_evidence["executable"]
        evidence["path_r_resolution"] = os.fspath(path_r)
        evidence["wrapper_sha256"] = _sha256(executable_path)
        evidence["check_args"] = selected_environment.get("CHECK_ARGS", "")
    if evidence_path is not None:
        _atomic_write_json(Path(evidence_path), evidence)

    completed = subprocess.run(
        command,
        check=False,
        shell=False,
        env=process_environment,
        cwd=cwd,
    )
    return completed.returncode


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify and run one explicit CRAN preflight driver."
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--driver", required=True, choices=sorted(_DRIVERS))
    parser.add_argument("--executable", required=True)
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--r-entrypoint", choices=("rscript", "r-cmd"))
    parser.add_argument("--required-r-executable")
    parser.add_argument("driver_argv", nargs=argparse.REMAINDER)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _argument_parser().parse_args(argv)
    driver_argv = args.driver_argv
    if driver_argv[:1] == ["--"]:
        driver_argv = driver_argv[1:]
    try:
        manifest = load_manifest(args.manifest)
        validate_manifest(manifest)
        matches = [
            row
            for row in manifest["coverage"]
            if row["id"] == args.environment_id
        ]
        if len(matches) != 1:
            raise ValueError(
                f"environment_id must select exactly one manifest row: "
                f"{args.environment_id!r}"
            )
        return run_driver(
            matches[0],
            requested_driver=args.driver,
            executable=args.executable,
            argv=driver_argv,
            tarball=args.tarball,
            metadata=args.metadata,
            expected_source_sha=args.source_sha,
            expected_event_sha=args.event_sha,
            r_entrypoint=args.r_entrypoint,
            required_r_executable=args.required_r_executable,
            evidence_path=args.evidence,
        )
    except ValueError as error:
        print(f"driver configuration error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
