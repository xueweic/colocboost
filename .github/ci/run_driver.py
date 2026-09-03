"""Run an explicitly selected CRAN preflight driver without PATH fallback."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
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
    unbound_arguments = _normalize_argv(argv)
    if unbound_arguments.count(_TARBALL_TOKEN) != 1:
        raise ValueError(
            "argv must contain exactly one literal {tarball} placeholder"
        )
    verified_tarball = Path(os.path.abspath(os.fspath(tarball)))
    verify_tarball(
        verified_tarball,
        metadata,
        expected_source_sha=expected_source_sha,
        expected_event_sha=expected_event_sha,
    )
    bound_arguments = [
        str(verified_tarball) if argument == _TARBALL_TOKEN else argument
        for argument in unbound_arguments
    ]
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
            evidence_path=args.evidence,
        )
    except ValueError as error:
        print(f"driver configuration error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
