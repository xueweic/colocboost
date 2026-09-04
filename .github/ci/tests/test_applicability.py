import io
import sys
import tarfile
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402
from check_applicability import inspect_tarball, rule_holds  # noqa: E402


SHA = "a" * 40


def make_artifact(tmp_path, files=None, description="Package: example\nVersion: 0.0.1\nNeedsCompilation: no\n"):
    files = files or {}
    tarball = tmp_path / "example_0.0.1.tar.gz"
    with tarfile.open(tarball, "w:gz") as archive:
        payload = description.encode()
        info = tarfile.TarInfo("example/DESCRIPTION"); info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
        for name, content in files.items():
            payload = content.encode() if isinstance(content, str) else content
            info = tarfile.TarInfo(f"example/{name}"); info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    metadata = tmp_path / "source-metadata.json"
    create_metadata(tarball, metadata, source_sha=SHA, event_sha=SHA)
    return tarball, metadata


def features(tmp_path, files=None, description="Package: example\nVersion: 0.0.1\nNeedsCompilation: no\n"):
    tarball, metadata = make_artifact(tmp_path, files, description)
    return inspect_tarball(str(tarball), str(metadata), SHA, SHA)


def test_needs_compilation_no_is_a_valid_non_native_package(tmp_path):
    observed = features(tmp_path)
    assert not observed["needs_compilation"]
    assert rule_holds("no-native-source-linkingto-compilation-or-direct-call", observed)


@pytest.mark.parametrize(
    ("files", "description", "rule"),
    [
        ({"src/x.c": "int x;"}, None, "no-native-source-linkingto-compilation-or-direct-call"),
        ({"src/x.c": "int x;"}, None, "no-c-source"),
        ({"src/x.f90": "end"}, None, "no-package-owned-native-source"),
        ({"src/x.o": b"object"}, None, "no-native-object"),
        ({"configure": "#!/bin/sh"}, None, "no-native-compilation-target"),
        ({"src/x.c": "int x;"}, None, "no-compiled-source"),
        ({"src/x.cpp": "int x;"}, None, "no-cpp-source"),
        ({"R/x.R": ".Call('x')"}, None, "no-direct-native-call"),
        ({"src/x.cpp": "int x;"}, None, "no-c-cpp-object"),
        ({"src/Makevars": "PKG_CFLAGS=-fopenmp"}, None, "no-native-compilation-target"),
        ({"R/x.R": "x <- 1"}, "Package: example\nVersion: 0.0.1\nNeedsCompilation: yes\n", "no-native-source-linkingto-compilation-or-direct-call"),
        ({"R/x.R": "x <- 1"}, "Package: example\nVersion: 0.0.1\nLinkingTo: Rcpp\nNeedsCompilation: no\n", "no-native-source-linkingto-compilation-or-direct-call"),
        ({"src/Makevars": "PKG_LIBS=-lgomp"}, None, "no-native-source-linkingto-compilation-or-direct-call"),
    ],
)
def test_native_mutations_invalidate_the_relevant_predicate(tmp_path, files, description, rule):
    observed = features(tmp_path, files, description or "Package: example\nVersion: 0.0.1\nNeedsCompilation: no\n")
    assert not rule_holds(rule, observed)


def test_rejects_tarball_with_duplicate_members(tmp_path):
    tarball = tmp_path / "duplicate.tar.gz"
    with tarfile.open(tarball, "w:gz") as archive:
        for _ in range(2):
            payload = b"Package: example\nVersion: 0.0.1\n"
            info = tarfile.TarInfo("example/DESCRIPTION"); info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    metadata = tmp_path / "source-metadata.json"
    create_metadata(tarball, metadata, source_sha=SHA, event_sha=SHA)
    with pytest.raises(ValueError, match="duplicate"):
        inspect_tarball(str(tarball), str(metadata), SHA, SHA)
