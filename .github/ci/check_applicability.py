#!/usr/bin/env python3
"""Fail-closed, tarball-only applicability predicates for CRAN-like lanes."""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
from collections.abc import Iterable
from pathlib import PurePosixPath

from artifact_contract import verify_tarball
from extract_source import _validated_members
from result_contract import finalize, write_provisional


_NATIVE_SUFFIXES = {".c", ".cc", ".cp", ".cpp", ".cxx", ".f", ".f90", ".f95", ".for", ".m", ".mm"}
_OBJECT_SUFFIXES = {".o", ".obj", ".a", ".so", ".dylib", ".dll", ".lib"}
_CPP_SUFFIXES = {".cc", ".cp", ".cpp", ".cxx", ".c++", ".mm"}
_DIRECT_CALL = re.compile(r"\.(?:C|Call|Fortran|External)\s*\(")
_OPENMP = re.compile(r"(?:^|[^A-Za-z0-9_])(?:openmp|fopenmp|libgomp|libomp)(?:$|[^A-Za-z0-9_])", re.I)
_RULES = frozenset({
    "no-native-source-linkingto-compilation-or-direct-call", "no-c-source",
    "no-package-owned-native-source", "no-native-object",
    "no-native-compilation-target", "no-compiled-source", "no-cpp-source",
    "no-direct-native-call", "no-c-cpp-object",
})


def _description(members: Iterable[tarfile.TarInfo], archive: tarfile.TarFile, root: str) -> str:
    target = f"{root}/DESCRIPTION"
    for member in members:
        if member.name == target:
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError("could not read package DESCRIPTION")
            with stream:
                return stream.read().decode("utf-8", errors="strict")
    raise ValueError("package DESCRIPTION is missing")


def _description_field(text: str, field: str) -> str | None:
    match = re.search(rf"(?im)^{re.escape(field)}\s*:\s*(.*(?:\n[ \t]+.*)*)", text)
    return None if match is None else " ".join(match.group(1).split())


def inspect_tarball(tarball: str, metadata: str, source_sha: str, event_sha: str) -> dict[str, bool]:
    """Verify the artifact identity and inspect only safe archive members."""
    verify_tarball(tarball, metadata, expected_source_sha=source_sha, expected_event_sha=event_sha)
    try:
        with tarfile.open(tarball, mode="r:gz") as archive:
            members, root = _validated_members(archive)
            description = _description(members, archive, root)
            names = [PurePosixPath(member.name) for member in members if member.isfile()]
            text_members = []
            for member in members:
                path = PurePosixPath(member.name)
                if not member.isfile():
                    continue
                suffix = path.suffix.lower()
                if suffix in {".r", ".rmd", ".rd", ".in", ".ac", ".am", ".mk", ".make", ".sh", ".txt"} or path.name in {"configure", "configure.ac", "makevars", "makevars.in"}:
                    stream = archive.extractfile(member)
                    if stream is None:
                        raise ValueError(f"could not read {member.name}")
                    with stream:
                        text_members.append((path, stream.read().decode("utf-8", errors="replace")))
    except (OSError, tarfile.TarError) as error:
        raise ValueError(f"could not inspect source tarball: {error}") from error

    relative = [path.relative_to(root) for path in names]
    native_source = any(path.suffix.lower() in _NATIVE_SUFFIXES for path in relative)
    c_source = any(path.suffix.lower() == ".c" for path in relative)
    cpp_source = any(path.suffix.lower() in _CPP_SUFFIXES for path in relative)
    native_object = any(path.suffix.lower() in _OBJECT_SUFFIXES for path in relative)
    compilation_target = any(path.parts and path.parts[0].lower() == "src" for path in relative) or any(path.name.lower().startswith("makevars") or path.name in {"configure", "configure.ac", "configure.in"} for path in relative)
    text = "\n".join(value for _, value in text_members)
    linking_to = _description_field(description, "LinkingTo")
    needs = _description_field(description, "NeedsCompilation")
    invalid_needs = needs is not None and needs.lower() != "no"
    direct_call = _DIRECT_CALL.search(text) is not None
    openmp = _OPENMP.search(text) is not None
    return {
        "native_source": native_source, "c_source": c_source, "cpp_source": cpp_source,
        "native_object": native_object, "compilation_target": compilation_target,
        "linking_to": linking_to is not None, "needs_compilation": invalid_needs,
        "direct_call": direct_call, "openmp": openmp,
    }


def rule_holds(rule: str, features: dict[str, bool]) -> bool:
    if rule not in _RULES:
        raise ValueError(f"unknown applicability rule: {rule}")
    if rule == "no-native-source-linkingto-compilation-or-direct-call":
        return not any(features[key] for key in ("native_source", "linking_to", "needs_compilation", "compilation_target", "direct_call", "openmp"))
    if rule == "no-c-source": return not features["c_source"]
    if rule == "no-package-owned-native-source": return not features["native_source"]
    if rule == "no-native-object": return not any(features[key] for key in ("native_object", "compilation_target", "needs_compilation"))
    if rule == "no-native-compilation-target": return not any(features[key] for key in ("compilation_target", "native_source", "needs_compilation", "openmp"))
    if rule == "no-compiled-source": return not features["native_source"]
    if rule == "no-cpp-source": return not features["cpp_source"]
    if rule == "no-direct-native-call": return not features["direct_call"]
    return not any(features[key] for key in ("c_source", "cpp_source", "native_object", "needs_compilation"))


def run(arguments: argparse.Namespace) -> int:
    identity = {"result_kind": "applicability", "environment_id": arguments.environment_id, "source_sha": arguments.source_sha, "event_sha": arguments.event_sha, "tarball_sha256": None}
    write_provisional(arguments.output, identity)
    try:
        features = inspect_tarball(arguments.tarball, arguments.metadata, arguments.source_sha, arguments.event_sha)
        verified = verify_tarball(arguments.tarball, arguments.metadata, expected_source_sha=arguments.source_sha, expected_event_sha=arguments.event_sha)
        holds = rule_holds(arguments.rule, features)
        finalize(arguments.output, "not-applicable" if holds else "uncovered", {"tarball_sha256": verified["sha256"], "predicate": arguments.rule, "evidence": ",".join(key for key, value in sorted(features.items()) if value) or "no-native-features"})
        return 0 if holds else 1
    except (OSError, ValueError, tarfile.TarError) as error:
        # The final result remains non-green and carries a stable diagnostic.
        finalize(arguments.output, "uncovered", {"tarball_sha256": "0" * 64, "predicate": arguments.rule, "evidence": f"verification-failed:{error}"})
        return 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tarball", required=True); parser.add_argument("--metadata", required=True)
    parser.add_argument("--source-sha", required=True); parser.add_argument("--event-sha", required=True)
    parser.add_argument("--environment-id", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--rule", required=True, choices=sorted(_RULES))
    arguments = parser.parse_args(argv)
    try:
        return run(arguments)
    except (OSError, ValueError) as error:
        print(f"applicability check failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
