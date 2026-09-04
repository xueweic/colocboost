#!/usr/bin/env python3
"""Run a pinned CRAN-like check through Docker from a hosted runner.

This module deliberately owns the Docker boundary.  It never invokes a shell,
never resolves an image tag, and treats a failed identity/proof check as a
failed producer rather than substituting the host's R installation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from artifact_contract import verify_tarball
from result_contract import finalize, write_provisional


BASE_DIGEST = "sha256:be9d65e2344d805cc11114319c685ecaa96b6d9b4350a0a6460cdb931babbd19"
R_SVN_URL = "https://svn.r-project.org/R/trunk"
R_SVN_REVISION = "90483"
BLIS_COMMIT = "e8566eb3e773fb54d11b33e371d13f22d2941e50"
MUSL_IMAGE = "ghcr.io/bastistician/rcheck-musl@sha256:8b33a511897c8025a84072efdb97dc0a58292428700bf90e6628a438134c453d"
ARM_INDEX = "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:74f3781764e90b6c4bf8e300c788b0f853e867289585f020dd3b842180c08bf6"
ARM_AMD64 = "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:ef57b4fbbcc76a173e53622143d59f8d4a97978481237e81153af23b03c66e9c"
ARM_ARM64 = "ghcr.io/r-devel/rcheckserver/ubuntu@sha256:bcfe37df999b716bd41ffc139d242c251c06b79d82f3632b4e8c65b6eb13f050"


def _digest_image(image: str) -> str:
    if isinstance(image, str) and re.fullmatch(r"local/colocboost-(?:blis|noomp)-proxy", image):
        return image
    if not isinstance(image, str) or not re.fullmatch(
        r"[A-Za-z0-9./_-]+@sha256:[0-9a-f]{64}", image
    ):
        raise ValueError(f"Docker image must use an immutable digest: {image!r}")
    return image


def _regular_file(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is not accessible: {path}") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file: {path}")
    return path


def _regular_directory(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    try:
        status = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} is not accessible: {path}") from error
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink directory: {path}")
    return path


def validate_proxy_dockerfile(path: str | os.PathLike[str], kind: str) -> Path:
    """Validate the immutable source/build contract before ``docker build``."""

    dockerfile = _regular_file(path, "Dockerfile")
    text = dockerfile.read_text(encoding="utf-8")
    required = [
        f"FROM docker.io/library/fedora@{BASE_DIGEST}",
        f"svn export -r {R_SVN_REVISION} {R_SVN_URL}",
    ]
    if kind == "blis":
        required += [
            BLIS_COMMIT,
            "./configure --enable-threading=no",
            "BLIS_NUM_THREADS=1",
        ]
    elif kind == "noomp":
        required += [
            "./configure --disable-openmp",
            '"Rfast"',
            "SHLIB_OPENMP_CFLAGS=",
            "SHLIB_OPENMP_LDFLAGS=",
        ]
    else:
        raise ValueError(f"unsupported proxy Dockerfile kind: {kind}")
    missing = [item for item in required if item not in text]
    if missing:
        raise ValueError(f"{kind} Dockerfile is missing contract: {missing[0]}")
    return dockerfile


def docker_argv(
    image: str,
    *,
    mounts: Sequence[tuple[str | os.PathLike[str], str, bool]],
    command: Sequence[str],
    user: str = "0:0",
    workdir: str = "/work",
    environment: Mapping[str, str] | None = None,
    platform: str | None = None,
    network: str = "none",
) -> list[str]:
    """Build an explicit, shell-free Docker invocation."""

    image = _digest_image(image)
    if not command or any(not isinstance(item, str) or not item for item in command):
        raise ValueError("container command must be nonempty strings")
    if not re.fullmatch(r"[0-9]+(?::[0-9]+)?", user):
        raise ValueError("container user must be numeric uid[:gid]")
    if network not in {"none", "default"}:
        raise ValueError("network must be none or default")
    argv = ["docker", "run", "--rm", "--network", network, "--user", user]
    if platform is not None:
        if platform not in {"linux/amd64", "linux/arm64"}:
            raise ValueError("platform must be linux/amd64 or linux/arm64")
        argv += ["--platform", platform]
    argv += ["--workdir", workdir]
    for source, destination, readonly in mounts:
        source_value = Path(source)
        source_path = (
            _regular_file(source_value, "mount source")
            if source_value.is_file()
            else _regular_directory(source_value, "mount source")
        )
        spec = f"type=bind,src={source_path},dst={destination}"
        if readonly:
            spec += ",readonly"
        argv += ["--mount", spec]
    for key, value in sorted((environment or {}).items()):
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise ValueError(f"invalid container environment key: {key}")
        argv += ["--env", f"{key}={value}"]
    return [*argv, image, *command]


def _atomic_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
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


def _check_counts(log: Path) -> dict[str, Any]:
    text = log.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    errors = sum(bool(re.search(r"(?:^|\s)ERROR(?:$|\s)", line)) for line in lines)
    warnings = sum(bool(re.search(r"(?:^|\s)WARNING(?:$|\s)", line)) for line in lines)
    notes = sum(bool(re.search(r"(?:^|\s)NOTE(?:$|\s)", line)) for line in lines)
    ok = text.endswith("* DONE\nStatus: OK\n") and not any((errors, warnings, notes))
    return {"errors": int(errors), "warnings": int(warnings), "notes": int(notes), "log_path": os.fspath(log)}, ok


DEPENDENCY_PACKAGES = ("Rfast", "matrixStats", "testthat", "knitr", "rmarkdown", "ashr", "MASS", "susieR")


def run_prepare(
    *, image: str, r_binary: str, tarball: str | os.PathLike[str], metadata: str | os.PathLike[str],
    source_sha: str, event_sha: str, library: str | os.PathLike[str], work_dir: str | os.PathLike[str],
    evidence: str | os.PathLike[str], environment_id: str, user: str = "0:0",
    platform: str | None = None,
) -> int:
    """Populate a fresh library while networking is enabled, before an isolated check."""
    verified = verify_tarball(tarball, metadata, expected_source_sha=source_sha, expected_event_sha=event_sha)
    if user != "0:0":
        raise ValueError("dependency preparation must run as root")
    library_path = Path(library).absolute()
    if library_path.exists() and any(library_path.iterdir()):
        raise ValueError("dependency library must be fresh and empty")
    library_path.mkdir(parents=True, exist_ok=True)
    work = Path(work_dir).absolute(); work.mkdir(parents=True, exist_ok=True)
    input_dir = work / "prep-input"; input_dir.mkdir(exist_ok=True)
    source = Path(tarball).absolute(); copied = input_dir / verified["filename"]
    copied.write_bytes(source.read_bytes())
    package_vector = ",".join(repr(x) for x in DEPENDENCY_PACKAGES)
    r_code = (
        "options(repos=c(CRAN='https://cloud.r-project.org')); "
        f"wanted <- c({package_vector}); "
        "install.packages(wanted, lib='/library', "
        "dependencies=c('Depends','Imports','LinkingTo')); "
        "missing <- setdiff(wanted, rownames(installed.packages(lib.loc='/library'))); "
        "if (length(missing)) stop('dependency installation incomplete: ', paste(missing, collapse=',')); "
        "status <- system2(command=Sys.getenv('CB_R_BINARY'), args=c('CMD','INSTALL','--library=/library', "
        "paste0('/input/', Sys.getenv('CB_TARBALL')))); "
        "if (!identical(status, 0L)) stop('package installation failed with status ', status); "
        "installed <- find.package('colocboost', lib.loc='/library'); "
        "stopifnot(startsWith(normalizePath(installed), '/library/'));"
    )
    command = docker_argv(
        image, mounts=[(input_dir, "/input", True), (library_path, "/library", False), (work, "/work", False)],
        command=[r_binary, "--vanilla", "--slave", "-e", r_code],
        user=user, platform=platform, network="default",
        environment={
            "CB_R_BINARY": r_binary,
            "CB_TARBALL": verified["filename"],
            "R_LIBS_USER": "/library",
            "R_PROFILE_USER": "/dev/null",
            "R_ENVIRON_USER": "/dev/null",
        },
    )
    try:
        completed = subprocess.run(command, check=False, shell=False, timeout=2700)
        exit_code = completed.returncode
    except subprocess.TimeoutExpired:
        exit_code = 124
    proof = {"schema_version": 1, "kind": "dependency-prep", "image": image,
             "environment_id": environment_id, "r_binary": r_binary,
             "source_sha": source_sha, "event_sha": event_sha, "tarball_sha256": verified["sha256"],
             "library": "/library", "packages": list(DEPENDENCY_PACKAGES), "network": "default",
             "exit_code": exit_code}
    _atomic_json(Path(evidence), proof)
    return exit_code


def run_check(
    *,
    image: str,
    wrapper: str | None = None,
    r_binary: str | None = None,
    direct: bool = False,
    tarball: str | os.PathLike[str],
    metadata: str | os.PathLike[str],
    source_sha: str,
    event_sha: str,
    environment_id: str,
    work_dir: str | os.PathLike[str],
    result: str | os.PathLike[str],
    evidence: str | os.PathLike[str],
    user: str = "0:0",
    platform: str | None = None,
    readonly_library: str | os.PathLike[str] | None = None,
    library: str | os.PathLike[str] | None = None,
    environment: Mapping[str, str] | None = None,
    runtime_proof: str | os.PathLike[str] | None = None,
) -> int:
    """Run a wrapper or direct R CMD check and emit an aggregate package-check result."""
    if direct and not r_binary:
        raise ValueError("direct checks require r_binary")
    if not direct and not wrapper:
        raise ValueError("wrapper checks require wrapper")

    verified = verify_tarball(tarball, metadata, expected_source_sha=source_sha, expected_event_sha=event_sha)
    work = Path(work_dir).absolute()
    work.mkdir(parents=True, exist_ok=True)
    input_dir = work / "input"
    input_dir.mkdir(exist_ok=True)
    source = Path(tarball).absolute()
    copied = input_dir / verified["filename"]
    copied.write_bytes(source.read_bytes())
    mounts: list[tuple[str | os.PathLike[str], str, bool]] = [
        (input_dir, "/input", direct),
        (work, "/work", False),
    ]
    if library is not None and readonly_library is not None:
        raise ValueError("library and readonly_library are mutually exclusive")
    if library is not None:
        mounts.append((library, "/library", False))
    if readonly_library is not None:
        mounts.append((readonly_library, "/library", True))
    check_command = ([r_binary, "CMD", "check", "--as-cran", "--no-manual", "--no-build-vignettes",
                      f"/input/{verified['filename']}"] if direct else [wrapper, "/input"])
    check_environment = {
        "R_PROFILE_USER": "/dev/null",
        "R_ENVIRON_USER": "/dev/null",
        "_R_CHECK_FORCE_SUGGESTS_": "true",
        **(environment or {}),
        **({"R_LIBS_USER": "/library"} if (readonly_library is not None or library is not None) else {}),
    }
    if environment_id == "musl":
        check_environment.update({"LC_ALL": "C.UTF-8", "LANG": "C.UTF-8", "TZ": "Europe/Berlin"})
    command = docker_argv(
        image,
        mounts=mounts,
        command=check_command,
        user=user,
        environment=check_environment,
        platform=platform,
    )
    readonly_proof = None
    if readonly_library is not None:
        if user.startswith("0"):
            raise ValueError("read-only library checks must run as a non-root uid")
        probe = docker_argv(
            image,
            mounts=mounts,
            command=["/bin/sh", "-c", "grep -F /library /proc/self/mountinfo | grep -E '(^|,)ro(,|$)' && touch /library/.colocboost-write-probe"],
            user=user,
            environment={"R_LIBS_USER": "/library"},
            platform=platform,
        )
        probe_run = subprocess.run(probe, check=False, shell=False, capture_output=True, text=True)
        if probe_run.returncode == 0:
            raise ValueError("read-only bind mount accepted a write")
        mountinfo = probe_run.stdout or probe_run.stderr
        if "/library" not in mountinfo or not re.search(r"(?:^|,)ro(?:,|$)", mountinfo):
            raise ValueError("mountinfo did not prove /library is read-only")
        readonly_proof = {"write_probe_exit": probe_run.returncode, "mount": "/library:ro", "mountinfo": mountinfo, "nonroot_user": user}
        if direct:
            identity_probe = docker_argv(
                image, mounts=mounts,
                command=[r_binary, "--vanilla", "--slave", "-e",
                         "stopifnot(any(grepl('/library', .libPaths(), fixed=TRUE))); "
                         "p <- find.package('colocboost'); stopifnot(startsWith(normalizePath(p), '/library')); "
                         "writeLines(c(paste(.libPaths(), collapse=';'), p), '/work/library-identity.txt')"],
                user=user, environment={"R_LIBS_USER": "/library"}, platform=platform,
            )
            identity_run = subprocess.run(identity_probe, check=False, shell=False, capture_output=True, text=True)
            if identity_run.returncode != 0:
                raise ValueError("R library identity did not resolve to read-only /library")
            identity_path = work / "library-identity.txt"
            if not identity_path.is_file() or "/library" not in identity_path.read_text(encoding="utf-8"):
                raise ValueError("R library identity proof is missing")
            readonly_proof["r_library_identity"] = identity_path.read_text(encoding="utf-8").splitlines()
    completed = subprocess.run(command, check=False, shell=False)
    # Direct R checks write below /work.  Native wrappers differ: some change
    # into /input while the purpose-built proxies retain /work.
    logs = list(work.glob("*.Rcheck/00check.log")) + list(input_dir.glob("*.Rcheck/00check.log"))
    if len(logs) != 1 or logs[0].is_symlink() or not logs[0].is_file():
        raise ValueError("Docker check must produce exactly one regular *.Rcheck/00check.log")
    check, clean = _check_counts(logs[0])
    proof = {
        "schema_version": 1,
        "environment_id": environment_id,
        "source_sha": source_sha,
        "event_sha": event_sha,
        "tarball_sha256": verified["sha256"],
        "image": image,
        "command": command,
        "user": user,
        "platform": platform,
        "readonly_bind_mount": readonly_proof,
        "check_exit_code": completed.returncode,
        "mode": "direct-r-cmd-check" if direct else "wrapper",
        "r_binary": r_binary,
        "wrapper": wrapper,
        "check_log": os.fspath(logs[0]),
    }
    if runtime_proof is not None:
        runtime = _regular_file(runtime_proof, "runtime proof")
        proof["runtime_proof"] = os.fspath(runtime)
    _atomic_json(Path(evidence), proof)
    identity = {"result_kind": "package-check", "environment_id": environment_id, "source_sha": source_sha, "event_sha": event_sha, "tarball_sha256": verified["sha256"]}
    write_provisional(result, identity)
    finalize(result, "pass" if completed.returncode == 0 and clean else "fail", {"check": check})
    if completed.returncode != 0:
        return completed.returncode
    return 0 if clean else 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--wrapper")
    parser.add_argument("--r-binary")
    parser.add_argument("--direct", action="store_true")
    parser.add_argument("--tarball", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True)
    parser.add_argument("--event-sha", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--result")
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--user", default="0:0")
    parser.add_argument("--platform")
    parser.add_argument("--readonly-library")
    parser.add_argument("--runtime-proof")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--library")
    parser.add_argument("--check-library")
    args = parser.parse_args(argv)
    try:
        values = vars(args)
        if values.pop("prepare"):
            library = values.pop("library")
            values.pop("result")
            values.pop("wrapper")
            values.pop("direct")
            values.pop("readonly_library")
            values.pop("runtime_proof")
            values.pop("check_library")
            if not values.get("r_binary") or not library:
                raise ValueError("prepare requires --r-binary and --library")
            return run_prepare(library=library, **values)
        values.pop("library")
        values["library"] = values.pop("check_library")
        if not values.get("result"):
            raise ValueError("check requires --result")
        return run_check(**values)
    except (OSError, TypeError, ValueError, subprocess.SubprocessError) as error:
        print(f"host Docker check failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
