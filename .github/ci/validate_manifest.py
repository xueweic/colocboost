#!/usr/bin/env python3
"""Validate the committed CRAN preflight coverage manifest."""

from __future__ import annotations

import re
import sys
from collections import Counter
from collections.abc import Mapping, Sequence
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import yaml


PRIMARY_STATES = {
    "r-devel-linux-x86_64-debian-clang": "proxy",
    "r-devel-linux-x86_64-debian-gcc": "proxy",
    "r-devel-linux-x86_64-fedora-clang": "proxy",
    "r-devel-linux-x86_64-fedora-gcc": "direct",
    "r-devel-windows-x86_64": "proxy",
    "r-patched-linux-x86_64": "proxy",
    "r-release-linux-x86_64": "proxy",
    "r-release-macos-arm64": "proxy",
    "r-release-macos-x86_64": "proxy",
    "r-release-windows-x86_64": "proxy",
    "r-oldrel-macos-arm64": "proxy",
    "r-oldrel-macos-x86_64": "proxy",
    "r-oldrel-windows-x86_64": "proxy",
}
ADDITIONAL_STATES = {
    "ATLAS": "direct",
    "BLIS": "proxy",
    "BLAS": "not-applicable",
    "C23": "not-applicable",
    "Intel": "not-applicable",
    "LTO": "not-applicable",
    "M1mac": "proxy",
    "MKL": "direct",
    "OpenBLAS": "proxy",
    "Strict": "not-applicable",
    "clang-ASAN": "direct",
    "clang-UBSAN": "direct",
    "donttest": "direct",
    "gcc-ASAN": "direct",
    "gcc-UBSAN": "direct",
    "gcc": "not-applicable",
    "gcc15": "not-applicable",
    "noLD": "direct",
    "noOMP": "proxy",
    "noRemap": "not-applicable",
    "noSuggests": "direct",
    "valgrind": "direct",
    "0len": "not-applicable",
    "rchk": "not-applicable",
    "rcnst": "proxy",
    "rlibro": "proxy",
    "musl": "direct",
    "linux-arm64": "direct",
    "vnu": "direct",
}
EXPECTED_STATES = {"primary": PRIMARY_STATES, "additional": ADDITIONAL_STATES}
EXPECTED_ROW_CORE = {
    "r-devel-linux-x86_64-debian-clang": ("primary", "r-devel-linux-x86-64-debian-clang", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e"),
    "r-devel-linux-x86_64-debian-gcc": ("primary", "r-devel-linux-x86-64-debian-gcc", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-gcc16@sha256:2e9576e51ad17a706887b7e06fc4057388765226ec795f3c7af0e7348fb8fbf1"),
    "r-devel-linux-x86_64-fedora-clang": ("primary", "r-devel-linux-x86-64-fedora-clang", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang22@sha256:f4193769412c461365849dd664b3d42bd2a0aebe520eab3ae6cd426a19d8b71e"),
    "r-devel-linux-x86_64-fedora-gcc": ("primary", "r-devel-linux-x86-64-fedora-gcc", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2"),
    "r-devel-windows-x86_64": ("primary", "r-devel-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "r-patched-linux-x86_64": ("primary", "r-patched-linux-x86-64", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-next@sha256:1c29eed93b0aa05fe147e476e6b3510a09fb0476460eeefbbe99a4a1eb2c3bd9"),
    "r-release-linux-x86_64": ("primary", "r-release-linux-x86-64", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-release@sha256:714722b7ecb4307fbf88a707f83fb045181482c0c8ce6086a0e26860416ca9c2"),
    "r-release-macos-arm64": ("primary", "r-release-macos-arm64", "proxy", "r-binary", "runner", "macos-15"),
    "r-release-macos-x86_64": ("primary", "r-release-macos-x86-64", "proxy", "r-binary", "runner", "macos-15-intel"),
    "r-release-windows-x86_64": ("primary", "r-release-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "r-oldrel-macos-arm64": ("primary", "r-oldrel-macos-arm64", "proxy", "r-binary", "runner", "macos-15"),
    "r-oldrel-macos-x86_64": ("primary", "r-oldrel-macos-x86-64", "proxy", "r-binary", "runner", "macos-15-intel"),
    "r-oldrel-windows-x86_64": ("primary", "r-oldrel-windows-x86-64", "proxy", "r-binary", "runner", "windows-2022"),
    "ATLAS": ("additional", "atlas", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/atlas@sha256:7597f7d0b6b2f009ae7bb425391523d8f4388223238db50d7dfb1572add63a88"),
    "BLIS": ("additional", "blis", "proxy", "native-wrapper", "image", "local/colocboost-blis-proxy"),
    "BLAS": ("additional", "blas", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "C23": ("additional", "c23", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "Intel": ("additional", "intel", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "LTO": ("additional", "lto", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "M1mac": ("additional", "m1mac", "proxy", "r-binary", "runner", "macos-15"),
    "MKL": ("additional", "mkl", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/mkl@sha256:d84847130b3ae0b9b0402208ee17360ab527df1c448c44ba045b44524b17620a"),
    "OpenBLAS": ("additional", "openblas", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc16@sha256:1127418efe3938f0e72fc29c55596a8cdde03293928b60aeb0ac44b6586960d2"),
    "Strict": ("additional", "strict", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "clang-ASAN": ("additional", "clang-asan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang-asan@sha256:dfab3d2274151577eb705d2be9acd3391798ba4b56d384ca9df84544e5b6be96"),
    "clang-UBSAN": ("additional", "clang-ubsan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/clang-ubsan@sha256:a58b00b52e9c4b210c3474eac4c28dc18bf70c912d75bb45e466410452a694e6"),
    "donttest": ("additional", "donttest", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/donttest@sha256:fd1942c8b8627d7e1d80582a023b2792acfa158a6edfeeeba3edd27a52c57967"),
    "gcc-ASAN": ("additional", "gcc-asan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f"),
    "gcc-UBSAN": ("additional", "gcc-ubsan", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/gcc-asan@sha256:32c9ad423cbca983bae893c2f0f44021ddbd5ac236c29fba8078b97acd80d78f"),
    "gcc": ("additional", "gcc", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "gcc15": ("additional", "gcc15", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "noLD": ("additional", "nold", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/nold@sha256:9dac23acd0e610b5fb5f82957c8136ffc6acf216a328a067795241542b3dc091"),
    "noOMP": ("additional", "noomp", "proxy", "native-wrapper", "image", "local/colocboost-noomp-proxy"),
    "noRemap": ("additional", "noremap", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "noSuggests": ("additional", "nosuggests", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/nosuggests@sha256:588bd73470d657d0d2560a90e5fe036d191a89fde394b2c67c6423b71d7a12df"),
    "valgrind": ("additional", "valgrind", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/valgrind@sha256:98dfda5016513c269e33c8155ea57745d977dad4b1a703b28db82fcb1237ef9c"),
    "0len": ("additional", "0len", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "rchk": ("additional", "rchk", "not-applicable", "r-binary", "runner", "ubuntu-24.04"),
    "rcnst": ("additional", "rcnst", "proxy", "native-wrapper", "image", "ghcr.io/r-hub/containers/ubuntu-clang@sha256:b66e5f86ce6f8fa3e6afd497fabfdb3aeac79c74ee3525e8a8814de79aa2bc82"),
    "rlibro": ("additional", "rlibro", "proxy", "native-wrapper", "image", "local/colocboost-rlibro-proxy"),
    "musl": ("additional", "musl", "direct", "native-wrapper", "image", "cran-linked-public-musl-reproduction"),
    "linux-arm64": ("additional", "linux-arm64", "direct", "r-binary", "runner", ("ubuntu-24.04", "ubuntu-24.04-arm")),
    "vnu": ("additional", "vnu", "direct", "native-wrapper", "image", "ghcr.io/r-hub/containers/vnu@sha256:03f45d5944fc092627cae2e944f16f037fed9795f8768c90d9511799098a7884"),
}
ARTIFACT_ID = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
RHUB_IMAGE = re.compile(
    r"^ghcr\.io/r-hub/containers/(?P<name>[a-z0-9-]+)@sha256:[0-9a-f]{64}$"
)
DEPRECATED_RHUB_IMAGES = {"c23", "gcc15", "intel", "noremap", "ubuntu-gcc12"}
PRIMARY_PROOF = {
    "source-tarball-sha256",
    "r-executable",
    "r-version",
    "operating-system",
    "compiler",
    "architecture",
    "locale",
}


def _proofs(tokens: str) -> frozenset[str]:
    return frozenset(tokens.split())


REQUIRED_PROOFS = {
    "r-devel-linux-x86_64-debian-clang": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler clang-22 flang-22 libcxx architecture locale"),
    "r-devel-linux-x86_64-debian-gcc": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler gcc-16 architecture locale"),
    "r-devel-linux-x86_64-fedora-clang": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system compiler clang-22 flang-22 libcxx architecture locale"),
    "r-devel-linux-x86_64-fedora-gcc": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-devel-revision operating-system fedora-44 compiler gcc-version architecture locale"),
    "r-devel-windows-x86_64": _proofs("source-tarball-sha256 r-executable r-version r-devel-revision operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "r-patched-linux-x86_64": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-revision r-patched operating-system compiler architecture locale"),
    "r-release-linux-x86_64": _proofs("source-tarball-sha256 container-system-r r-executable r-version r-release operating-system compiler architecture locale"),
    "r-release-macos-arm64": _proofs("source-tarball-sha256 r-executable r-version r-release operating-system macos-version compiler apple-and-gnu-compilers architecture arm64 locale"),
    "r-release-macos-x86_64": _proofs("source-tarball-sha256 r-executable r-version r-release operating-system macos-version compiler apple-and-gnu-compilers architecture x86-64 locale"),
    "r-release-windows-x86_64": _proofs("source-tarball-sha256 r-executable r-version r-release operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "r-oldrel-macos-arm64": _proofs("source-tarball-sha256 r-executable r-version r-oldrel operating-system macos-version compiler apple-and-gnu-compilers architecture arm64 locale"),
    "r-oldrel-macos-x86_64": _proofs("source-tarball-sha256 r-executable r-version r-oldrel operating-system macos-version compiler apple-and-gnu-compilers architecture x86-64 locale"),
    "r-oldrel-windows-x86_64": _proofs("source-tarball-sha256 r-executable r-version r-oldrel operating-system windows-server-2022 compiler rtools-gcc architecture locale"),
    "ATLAS": _proofs("source-tarball-sha256 container-system-r r-executable extsoftversion-blas la-library process-mappings loaded-atlas"),
    "BLIS": _proofs("source-tarball-sha256 r-executable r-devel-revision operating-system compiler architecture loaded-blis blis-version single-thread"),
    "BLAS": _proofs("source-tarball-sha256 no-package-owned-native-source no-linkingto no-compilation-marker no-direct-native-call"),
    "C23": _proofs("source-tarball-sha256 no-package-c-source obsolete-special-area"),
    "Intel": _proofs("source-tarball-sha256 no-package-owned-native-source obsolete-special-area"),
    "LTO": _proofs("source-tarball-sha256 no-package-owned-native-object no-compilation-target"),
    "M1mac": _proofs("source-tarball-sha256 r-executable r-version r-devel-revision operating-system macos-version compiler architecture arm64"),
    "MKL": _proofs("source-tarball-sha256 container-system-r r-executable r-version operating-system compiler architecture session-info blas-identity lapack-identity la-library la-version process-mappings loaded-mkl blas-operation mkl-verbose single-thread"),
    "OpenBLAS": _proofs("source-tarball-sha256 container-system-r r-executable r-devel-revision operating-system compiler loaded-openblas single-thread"),
    "Strict": _proofs("source-tarball-sha256 no-package-owned-compilation-target strict-r-headers-not-applicable"),
    "clang-ASAN": _proofs("source-tarball-sha256 container-system-r r-executable sanitizer-compiler-flags sanitizer-linker-flags asan-runtime sanitizer-output-scan"),
    "clang-UBSAN": _proofs("source-tarball-sha256 container-system-r r-executable sanitizer-compiler-flags sanitizer-linker-flags ubsan-runtime undefined-behavior-output-scan"),
    "donttest": _proofs("source-tarball-sha256 container-system-r r-executable check-donttest-examples-true expanded-examples-executed"),
    "gcc-ASAN": _proofs("source-tarball-sha256 container-system-r r-executable gcc-asan-flags asan-preload asan-runtime sanitizer-output-scan"),
    "gcc-UBSAN": _proofs("source-tarball-sha256 container-system-r r-executable gcc-ubsan-flags ubsan-preload ubsan-runtime sanitizer-output-scan"),
    "gcc": _proofs("source-tarball-sha256 no-package-owned-compiled-source"),
    "gcc15": _proofs("source-tarball-sha256 no-package-owned-native-source obsolete-special-area"),
    "noLD": _proofs("source-tarball-sha256 container-system-r r-executable opt-r-devel-nold long-double-disabled"),
    "noOMP": _proofs("source-tarball-sha256 r-executable r-devel-revision no-openmp-compile-flags no-openmp-link-flags no-loaded-openmp-runtime dependency-openmp-audit"),
    "noRemap": _proofs("source-tarball-sha256 no-package-cpp-source obsolete-special-area"),
    "noSuggests": _proofs("source-tarball-sha256 container-system-r r-executable depends-only-policy allowed-test-frameworks allowed-vignette-builders ashr-absent susier-absent nonzero-installed-test-count"),
    "valgrind": _proofs("source-tarball-sha256 container-system-r r-executable opt-r-devel-valgrind use-valgrind valgrind-runtime suppression-and-error-scan"),
    "0len": _proofs("source-tarball-sha256 no-direct-native-interface removed-official-experimental-support"),
    "rchk": _proofs("source-tarball-sha256 no-package-c-cpp-object compiled-dependencies-visible-limitation"),
    "rcnst": _proofs("source-tarball-sha256 container-system-r r-executable r-devel-revision r-compile-pkgs-1 r-jit-strategy-4 r-check-constants-5 constant-corruption-diagnostics"),
    "rlibro": _proofs("source-tarball-sha256 r-executable nonroot-uid readonly-bind-mount write-failure installed-library-path package-check-executed"),
    "musl": _proofs("source-tarball-sha256 r-executable r-version musl-libc alpine-version locale architecture"),
    "linux-arm64": _proofs("source-tarball-sha256 r-executable native-amd64 native-arm64 identical-rcheckserver-image identical-tarball-and-configuration architecture-only-comparison"),
    "vnu": _proofs("source-tarball-sha256 container-system-r r-executable vnu-special-dispatch nu-validator-executed zero-bad-entries validator-output"),
}

PREDICATE_PREFIX = (
    "python",
    ".github/ci/evaluate_native_features.py",
    "--tarball",
    "{tarball}",
    "--rule",
)
EXPECTED_PREDICATES = {
    "BLAS": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-source-linkingto-compilation-or-direct-call",)),
    "C23": ("built-source-tarball", PREDICATE_PREFIX + ("no-c-source",)),
    "Intel": ("built-source-tarball", PREDICATE_PREFIX + ("no-package-owned-native-source",)),
    "LTO": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-object",)),
    "Strict": ("built-source-tarball", PREDICATE_PREFIX + ("no-native-compilation-target",)),
    "gcc": ("built-source-tarball", PREDICATE_PREFIX + ("no-compiled-source",)),
    "gcc15": ("built-source-tarball", PREDICATE_PREFIX + ("no-package-owned-native-source",)),
    "noRemap": ("built-source-tarball", PREDICATE_PREFIX + ("no-cpp-source",)),
    "0len": ("built-source-tarball", PREDICATE_PREFIX + ("no-direct-native-call",)),
    "rchk": ("built-source-tarball", PREDICATE_PREFIX + ("no-c-cpp-object",)),
}
ACCEPTANCE_SPIKE_PROOFS = {
    "BLIS": {"loaded-blis", "blis-version", "single-thread"},
    "noOMP": {
        "no-openmp-compile-flags",
        "no-openmp-link-flags",
        "no-loaded-openmp-runtime",
        "dependency-openmp-audit",
    },
}
SPECIAL_IMAGE_PROOFS = {
    "MKL": {
        "container-system-r",
        "session-info",
        "blas-identity",
        "lapack-identity",
        "la-library",
        "la-version",
        "process-mappings",
        "loaded-mkl",
        "blas-operation",
        "mkl-verbose",
        "single-thread",
    },
    "noSuggests": {
        "container-system-r",
        "depends-only-policy",
        "allowed-test-frameworks",
        "allowed-vignette-builders",
        "ashr-absent",
        "susier-absent",
        "nonzero-installed-test-count",
    },
    "valgrind": {
        "container-system-r",
        "opt-r-devel-valgrind",
        "use-valgrind",
        "valgrind-runtime",
        "suppression-and-error-scan",
    },
    "vnu": {
        "container-system-r",
        "vnu-special-dispatch",
        "nu-validator-executed",
        "zero-bad-entries",
        "validator-output",
    },
}
APPROVED_RHUB_IMAGE = {
    "r-devel-linux-x86_64-debian-clang": "clang22",
    "r-devel-linux-x86_64-debian-gcc": "ubuntu-gcc16",
    "r-devel-linux-x86_64-fedora-clang": "clang22",
    "r-devel-linux-x86_64-fedora-gcc": "gcc16",
    "r-patched-linux-x86_64": "ubuntu-next",
    "r-release-linux-x86_64": "ubuntu-release",
    "ATLAS": "atlas",
    "MKL": "mkl",
    "OpenBLAS": "gcc16",
    "clang-ASAN": "clang-asan",
    "clang-UBSAN": "clang-ubsan",
    "donttest": "donttest",
    "gcc-ASAN": "gcc-asan",
    "gcc-UBSAN": "gcc-asan",
    "noLD": "nold",
    "noSuggests": "nosuggests",
    "valgrind": "valgrind",
    "rcnst": "ubuntu-clang",
    "vnu": "vnu",
}
DRIVER_ENDPOINT = {"native-wrapper": "image", "r-binary": "runner"}
SPECIAL_NATIVE_BINDINGS = {
    "atlas": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/devel/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "atlas",
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "fedora",
        "expected_distribution_version": "42",
    },
    "r-devel-linux-x86-64-debian-clang": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/devel/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "clang22",
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "ubuntu",
        "expected_distribution_version": "22.04",
    },
    "r-devel-linux-x86-64-debian-gcc": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/devel/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "ubuntu-gcc16",
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "ubuntu",
        "expected_distribution_version": "24.04",
    },
    "r-devel-linux-x86-64-fedora-clang": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/devel/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "clang22",
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "ubuntu",
        "expected_distribution_version": "22.04",
    },
    "r-devel-linux-x86-64-fedora-gcc": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/devel-gcc16/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "gcc16",
        "expected_r_kind": "devel",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "fedora",
        "expected_distribution_version": "44",
    },
    "r-patched-linux-x86-64": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/next/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": ["--no-manual", "--no-build-vignettes"],
        "runtime_profile": "ubuntu-next",
        "expected_r_kind": "patched",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "ubuntu",
        "expected_distribution_version": "24.04",
    },
    "r-release-linux-x86-64": {
        "wrapper_path": "/usr/local/bin/r-check",
        "wrapper_sha256": "a42092f0de63c4a9c1bed3c1c9b341b32c51f72335169d02732318c102646090",
        "system_r": "/opt/R/release/bin/R",
        "wrapper_input": "tarball-parent",
        "check_args": [],
        "runtime_profile": "ubuntu-release",
        "expected_r_kind": "release",
        "expected_os": "linux",
        "expected_architecture": "x86_64",
        "expected_distribution": "ubuntu",
        "expected_distribution_version": "24.04",
    },
    "mkl": {
        "wrapper_path": "/usr/local/bin/r-check",
        "system_r": "/opt/R/devel-mkl/bin/R",
        "wrapper_input": "tarball-parent",
    },
    "nosuggests": {
        "wrapper_path": "/usr/local/bin/r-check",
        "system_r": "/opt/R/devel/bin/R",
        "wrapper_input": "tarball-parent",
    },
}
PLATFORM_R_BINDINGS = {
    "r-devel-windows-x86-64": {
        "setup_r_selector": "devel",
        "system_r": "C:/R/bin/R.exe",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "devel",
        "expected_os": "windows",
        "expected_architecture": "x86_64",
    },
    "r-release-macos-arm64": {
        "setup_r_selector": "release",
        "system_r": "/Library/Frameworks/R.framework/Resources/bin/R",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "macos",
        "expected_architecture": "aarch64",
    },
    "r-release-macos-x86-64": {
        "setup_r_selector": "release",
        "system_r": "/Library/Frameworks/R.framework/Resources/bin/R",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "macos",
        "expected_architecture": "x86_64",
    },
    "r-release-windows-x86-64": {
        "setup_r_selector": "release",
        "system_r": "C:/R/bin/R.exe",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "windows",
        "expected_architecture": "x86_64",
    },
    "r-oldrel-macos-arm64": {
        "setup_r_selector": "oldrel-1",
        "system_r": "/Library/Frameworks/R.framework/Resources/bin/R",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "macos",
        "expected_architecture": "aarch64",
    },
    "r-oldrel-macos-x86-64": {
        "setup_r_selector": "oldrel-1",
        "system_r": "/Library/Frameworks/R.framework/Resources/bin/R",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "macos",
        "expected_architecture": "x86_64",
    },
    "r-oldrel-windows-x86-64": {
        "setup_r_selector": "oldrel-1",
        "system_r": "C:/R/bin/R.exe",
        "check_args": ["--as-cran", "--no-manual", "--no-build-vignettes"],
        "expected_r_kind": "release",
        "expected_os": "windows",
        "expected_architecture": "x86_64",
    },
}
EXECUTION_BINDING_FIELDS = {
    field
    for binding in (*SPECIAL_NATIVE_BINDINGS.values(), *PLATFORM_R_BINDINGS.values())
    for field in binding
}
UNIT_LANE_FIELDS = {
    "environment_id",
    "coverage_id",
    "mode",
    "suite",
    "runner_context",
    "required_sidecar_tests",
}
UNIT_COVERAGE_IDS = [
    "r-devel-linux-x86-64-debian-clang",
    "r-devel-linux-x86-64-debian-gcc",
    "r-devel-linux-x86-64-fedora-clang",
    "r-devel-linux-x86-64-fedora-gcc",
    "r-devel-windows-x86-64",
    "r-patched-linux-x86-64",
    "r-release-linux-x86-64",
    "r-release-macos-arm64",
    "r-release-macos-x86-64",
    "r-release-windows-x86-64",
    "r-oldrel-macos-arm64",
    "r-oldrel-macos-x86-64",
    "r-oldrel-windows-x86-64",
    "atlas",
    "blis",
    "mkl",
    "openblas",
    "nosuggests",
]
EXPECTED_UNIT_LANES = [
    (
        coverage_id,
        coverage_id,
        "installed" if coverage_id == "nosuggests" else "source",
        "full",
        "r-cmd-check-installed" if coverage_id == "nosuggests" else coverage_id,
        ("test_utils.R", "test_Xref.R") if coverage_id == "mkl" else (),
    )
    for coverage_id in UNIT_COVERAGE_IDS
]
WAIVER_FIELDS = {
    "id",
    "context",
    "file",
    "test_title",
    "reason",
    "expected_count",
    "rationale",
    "expires",
}
EXPECTED_WAIVERS = {
    "installed-model-init-data": ("r-cmd-check-installed", "test_model.R", "colocboost_init_data correctly initializes data", "colocboost_init_data not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-dictionary-mapping": ("r-cmd-check-installed", "test_model.R", "colocboost correctly maps focal outcome to keep_variables with dict_keep_variables", "colocboost_init_data not directly accessible for integration test", 1, "The legacy installed-package integration branch cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-assemble": ("r-cmd-check-installed", "test_model.R", "colocboost_assemble processes model results", "colocboost_assemble not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-model-workhorse": ("r-cmd-check-installed", "test_model.R", "colocboost_workhorse performs boosting iterations", "colocboost_workhorse not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
    "installed-utils-dictionary-mapping": ("r-cmd-check-installed", "test_utils.R", "colocboost_init_data handles complex dictionary mappings", "colocboost_init_data not directly accessible", 1, "The legacy installed-package test cannot access this unexported internal; source-mode tests exercise it.", "2027-03-03"),
}
EXPECTED_MKL_POLICY = {
    "required_library_patterns": [
        r"(?i)(?:^|/)libmkl_intel_lp64(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libmkl_core(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libmkl_sequential(?:\.so(?:\.\d+)*)?(?:$|\s)",
    ],
    "forbidden_library_patterns": [
        r"(?i)(?:^|/)libmkl_(?:intel|gnu|tbb)_thread(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)libopenblas",
        r"(?i)(?:^|/)lib(?:s?atlas)(?:\.so(?:\.\d+)*)?(?:$|\s)",
        r"(?i)(?:^|/)libblis(?:\.so(?:\.\d+)*)?(?:$|\s)",
    ],
}


def _require_nonempty_string(value: Any, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")


def _require_string_list(value: Any, label: str) -> None:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be a non-empty list of unique strings")


def _valid_endpoint(value: Any) -> bool:
    if isinstance(value, str):
        return bool(value.strip())
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(item, str) and bool(item.strip()) for item in value)
        and len(value) == len(set(value))
    )


def _expected_id(group: str, cran_name: str) -> str:
    del group
    return cran_name.lower().replace("_", "-")


def _validate_predicate(row: Mapping[str, Any]) -> None:
    predicate = row.get("predicate")
    if not isinstance(predicate, Mapping):
        raise ValueError(f"{row['cran_name']} predicate must be an object")
    if set(predicate) != {"input", "command"}:
        raise ValueError(
            f"{row['cran_name']} predicate must contain only input and command"
        )
    if predicate["input"] != "built-source-tarball":
        raise ValueError(f"{row['cran_name']} predicate must inspect built-source-tarball")
    command = predicate["command"]
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(part, str) or not part for part in command)
        or "{tarball}" not in command
    ):
        raise ValueError(
            f"{row['cran_name']} predicate command must be an executable argv list using {{tarball}}"
        )
    expected = EXPECTED_PREDICATES.get(row["cran_name"])
    actual = (predicate["input"], tuple(command))
    if actual != expected:
        raise ValueError(f"{row['cran_name']} predicate does not match approved contract")


def _validate_row(row: Any, index: int) -> None:
    if not isinstance(row, Mapping):
        raise ValueError(f"coverage row {index} must be an object")
    required = {"id", "group", "cran_name", "state", "driver", "proof"}
    missing = required - set(row)
    if missing:
        raise ValueError(f"coverage row {index} is missing {sorted(missing)[0]}")

    for field in ("id", "group", "cran_name", "state", "driver"):
        _require_nonempty_string(row[field], f"coverage row {index}.{field}")
    if ARTIFACT_ID.fullmatch(row["id"]) is None:
        raise ValueError(f"{row['id']!r} is not a lowercase artifact-safe id")
    if row["group"] not in EXPECTED_STATES:
        raise ValueError(f"{row['id']} has invalid group {row['group']!r}")
    if row["state"] not in {"direct", "proxy", "not-applicable", "uncovered"}:
        raise ValueError(f"{row['id']} has invalid state {row['state']!r}")
    if row["state"] == "uncovered":
        raise ValueError(f"{row['id']} is uncovered")

    expected_id = _expected_id(row["group"], row["cran_name"])
    if row["id"] != expected_id:
        raise ValueError(
            f"{row['cran_name']} id must remain the stable artifact-safe id {expected_id}"
        )

    _require_string_list(row["proof"], f"{row['id']}.proof")
    proof = set(row["proof"])
    if "source-tarball-sha256" not in proof:
        raise ValueError(f"{row['id']}.proof must include source-tarball-sha256")

    driver = row["driver"]
    if driver not in DRIVER_ENDPOINT:
        raise ValueError(f"{row['id']} has invalid driver {driver!r}")
    endpoint = DRIVER_ENDPOINT[driver]
    other_endpoint = "runner" if endpoint == "image" else "image"
    if endpoint not in row or not _valid_endpoint(row[endpoint]):
        raise ValueError(f"{driver} row {row['id']} requires a valid {endpoint}")
    if other_endpoint in row:
        raise ValueError(f"{driver} row {row['id']} must not declare {other_endpoint}")

    if endpoint == "image" and isinstance(row["image"], str):
        image = row["image"]
        if image.startswith("ghcr.io/r-hub/containers/"):
            image_name = image.split("/", maxsplit=4)[-1].split("@", maxsplit=1)[0]
            if image_name in DEPRECATED_RHUB_IMAGES:
                raise ValueError(f"{row['id']} uses deprecated R-hub image {image_name}")
            match = RHUB_IMAGE.fullmatch(image)
            if match is None:
                raise ValueError(f"{row['id']} R-hub image must use an immutable digest")
            approved = APPROVED_RHUB_IMAGE.get(row["cran_name"])
            if approved != match.group("name"):
                raise ValueError(
                    f"{row['id']} must use approved R-hub image {approved}, "
                    f"not {match.group('name')}"
                )
            if "container-system-r" not in proof:
                raise ValueError(
                    f"{row['id']}.proof must include container-system-r"
                )

    endpoint_value = row[endpoint]
    if isinstance(endpoint_value, list):
        endpoint_value = tuple(endpoint_value)
    actual_core = (
        row["group"],
        row["id"],
        row["state"],
        row["driver"],
        endpoint,
        endpoint_value,
    )
    if actual_core != EXPECTED_ROW_CORE.get(row["cran_name"]):
        raise ValueError(f"{row['id']} does not match approved row contract")

    expected_binding = SPECIAL_NATIVE_BINDINGS.get(row["id"])
    expected_platform_binding = PLATFORM_R_BINDINGS.get(row["id"])
    if expected_binding is not None:
        actual_binding = {field: row.get(field) for field in expected_binding}
        if actual_binding != expected_binding:
            raise ValueError(
                f"{row['id']} native wrapper binding does not match approved contract"
            )
    elif expected_platform_binding is not None:
        actual_binding = {
            field: row.get(field) for field in expected_platform_binding
        }
        if actual_binding != expected_platform_binding:
            raise ValueError(
                f"{row['id']} platform R binding does not match approved contract"
            )
    elif EXECUTION_BINDING_FIELDS & set(row):
        raise ValueError(
            f"{row['id']} must not declare an execution binding"
        )

    if row["group"] == "primary":
        missing_proof = PRIMARY_PROOF - proof
        if missing_proof:
            raise ValueError(
                f"{row['id']} primary proof is missing {sorted(missing_proof)[0]}"
            )
    if row["state"] == "direct" and "r-executable" not in proof:
        raise ValueError(f"{row['id']} direct row lacks an identity assertion")
    if row["state"] == "proxy":
        if not proof:
            raise ValueError(f"{row['id']}.proof is required for proxy coverage")
        _require_nonempty_string(row.get("limitation"), f"{row['id']}.limitation")
    if row["state"] == "not-applicable":
        _validate_predicate(row)

    required_spike_proof = ACCEPTANCE_SPIKE_PROOFS.get(row["cran_name"])
    if required_spike_proof and not required_spike_proof <= proof:
        missing_spike = sorted(required_spike_proof - proof)[0]
        raise ValueError(
            f"{row['cran_name']} acceptance spike proof is missing {missing_spike}"
        )
    required_special_proof = SPECIAL_IMAGE_PROOFS.get(row["cran_name"])
    if required_special_proof and not required_special_proof <= proof:
        missing_special = sorted(required_special_proof - proof)[0]
        raise ValueError(
            f"{row['cran_name']} special image proof is missing {missing_special}"
        )
    required_proof = REQUIRED_PROOFS.get(row["cran_name"])
    if required_proof is None or not required_proof <= proof:
        missing_proof = sorted((required_proof or set()) - proof)
        raise ValueError(
            f"{row['id']} proof contract is missing "
            f"{missing_proof[0] if missing_proof else 'an approved inventory entry'}"
        )


def _validate_unit_lanes(lanes: Any, rows: Sequence[Mapping[str, Any]]) -> None:
    if not isinstance(lanes, list):
        raise ValueError("unit_lanes must be a list")

    environment_ids = [
        lane.get("environment_id")
        for lane in lanes
        if isinstance(lane, Mapping)
        and isinstance(lane.get("environment_id"), str)
    ]
    duplicates = sorted(
        environment_id
        for environment_id, count in Counter(environment_ids).items()
        if isinstance(environment_id, str) and count > 1
    )
    if duplicates:
        raise ValueError(f"duplicate unit environment_id {duplicates[0]}")

    coverage_by_id = {row["id"]: row for row in rows}
    actual = []
    for index, lane in enumerate(lanes):
        if not isinstance(lane, Mapping) or set(lane) != UNIT_LANE_FIELDS:
            raise ValueError(
                f"unit lane {index} fields must be exactly {sorted(UNIT_LANE_FIELDS)}"
            )

        for field in (
            "environment_id",
            "coverage_id",
            "mode",
            "suite",
            "runner_context",
        ):
            _require_nonempty_string(lane[field], f"unit lane {index}.{field}")

        environment_id = lane["environment_id"]
        if ARTIFACT_ID.fullmatch(environment_id) is None:
            raise ValueError(
                f"unit environment_id {environment_id!r} is not artifact-safe"
            )

        coverage = coverage_by_id.get(lane["coverage_id"])
        if coverage is None or coverage["state"] not in {"direct", "proxy"}:
            raise ValueError(
                f"unit lane {environment_id} must reference direct or proxy coverage"
            )

        sidecars = lane["required_sidecar_tests"]
        if (
            not isinstance(sidecars, list)
            or any(not isinstance(item, str) or not item for item in sidecars)
            or len(sidecars) != len(set(sidecars))
        ):
            raise ValueError(
                f"unit lane {environment_id}.required_sidecar_tests must be a list "
                "of unique strings"
            )

        actual.append(
            (
                environment_id,
                lane["coverage_id"],
                lane["mode"],
                lane["suite"],
                lane["runner_context"],
                tuple(sidecars),
            )
        )

    if actual != EXPECTED_UNIT_LANES:
        raise ValueError("unit lanes do not match approved contract")


def _coverage_result_keys(
    rows: Sequence[Mapping[str, Any]],
) -> list[tuple[str, str]]:
    return [
        (
            "applicability" if row["state"] == "not-applicable" else "package-check",
            row["id"],
        )
        for row in rows
    ]


def _unit_result_keys(lanes: Sequence[Mapping[str, Any]]) -> list[tuple[str, str]]:
    return [("unit", lane["environment_id"]) for lane in lanes]


def load_manifest(path: str | Path) -> Mapping[str, Any]:
    """Load a YAML manifest without accepting an empty document."""

    manifest_path = Path(path)
    try:
        data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise ValueError(f"could not load manifest {manifest_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ValueError("manifest must be an object")
    return data


def validate_manifest(data: Mapping[str, Any]) -> Mapping[str, Any]:
    """Fail closed unless *data* declares the exact approved inventory."""

    if not isinstance(data, Mapping):
        raise ValueError("manifest must be an object")
    if set(data) != {"version", "coverage", "unit_lanes"}:
        raise ValueError(
            "manifest top-level fields must be only version, coverage, and unit_lanes"
        )
    if isinstance(data["version"], bool) or data["version"] != 1:
        raise ValueError("manifest version must be 1")
    rows = data["coverage"]
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        raise ValueError("coverage must be a list")

    ids = [row.get("id") for row in rows if isinstance(row, Mapping)]
    duplicates = sorted(
        name
        for name, count in Counter(ids).items()
        if isinstance(name, str) and count > 1
    )
    if duplicates:
        raise ValueError(f"duplicate id {duplicates[0]}")

    for index, row in enumerate(rows):
        _validate_row(row, index)

    names = [(row["group"], row["cran_name"]) for row in rows]
    duplicate_names = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicate_names:
        raise ValueError(f"duplicate inventory entry {duplicate_names[0]}")

    actual_by_group = {
        group: {row["cran_name"]: row["state"] for row in rows if row["group"] == group}
        for group in EXPECTED_STATES
    }
    for group, expected in EXPECTED_STATES.items():
        actual = actual_by_group[group]
        if actual != expected:
            missing = sorted(set(expected) - set(actual))
            extra = sorted(set(actual) - set(expected))
            wrong = sorted(
                name
                for name in set(actual) & set(expected)
                if actual[name] != expected[name]
            )
            detail = (
                f"missing={missing}, extra={extra}, wrong_state={wrong}"
            )
            raise ValueError(f"{group} inventory does not match approved mapping: {detail}")

    expected_order = list(EXPECTED_ROW_CORE)
    actual_order = [row["cran_name"] for row in rows]
    if actual_order != expected_order:
        raise ValueError("coverage order does not match approved aggregate contract")

    _validate_unit_lanes(data["unit_lanes"], rows)
    result_keys = _coverage_result_keys(rows) + _unit_result_keys(data["unit_lanes"])
    if len(result_keys) != len(set(result_keys)):
        raise ValueError("duplicate composite result key")

    return data


def expected_result_ids(data: Mapping[str, Any]) -> list[str]:
    """Return all stable result IDs after validating the complete inventory."""

    validate_manifest(data)
    return [row["id"] for row in data["coverage"]]


def expected_coverage_result_keys(
    data: Mapping[str, Any],
) -> list[tuple[str, str]]:
    """Return ordered package-check and applicability aggregate keys."""

    validate_manifest(data)
    return _coverage_result_keys(data["coverage"])


def expected_unit_result_keys(data: Mapping[str, Any]) -> list[tuple[str, str]]:
    """Return ordered unit aggregate keys without sidecar-only targeted runs."""

    validate_manifest(data)
    return _unit_result_keys(data["unit_lanes"])


def expected_result_keys(data: Mapping[str, Any]) -> list[tuple[str, str]]:
    """Return the complete ordered composite aggregate-result inventory."""

    validate_manifest(data)
    return _coverage_result_keys(data["coverage"]) + _unit_result_keys(
        data["unit_lanes"]
    )


def validate_policy(
    data: Mapping[str, Any], today: date | None = None
) -> Mapping[str, Any]:
    """Validate exact, unexpired check waivers and numerical-backend rules."""

    if not isinstance(data, Mapping):
        raise ValueError("policy must be an object")
    if set(data) != {"version", "unit_tests", "r_cmd_check", "numerical_backends"}:
        raise ValueError("policy has unexpected or missing top-level fields")
    if isinstance(data["version"], bool) or data["version"] != 1:
        raise ValueError("policy version must be 1")
    if today is None:
        today = datetime.now(timezone.utc).date()
    if not isinstance(today, date):
        raise ValueError("today must be a date")

    unit_tests = data["unit_tests"]
    if not isinstance(unit_tests, Mapping) or set(unit_tests) != {"allowed_skips"}:
        raise ValueError("unit_tests must contain only allowed_skips")
    waivers = unit_tests["allowed_skips"]
    if not isinstance(waivers, list):
        raise ValueError("allowed_skips must be a list")

    actual_waivers = {}
    for index, waiver in enumerate(waivers):
        if not isinstance(waiver, Mapping) or set(waiver) != WAIVER_FIELDS:
            raise ValueError(f"waiver {index} has unexpected or missing fields")
        expires = waiver["expires"]
        if not isinstance(expires, str):
            raise ValueError(f"waiver {index} expires must be an ISO expiry date")
        try:
            expiry = date.fromisoformat(expires)
        except ValueError as error:
            raise ValueError(
                f"waiver {index} expires must be an ISO expiry date"
            ) from error
        if expiry.isoformat() != expires:
            raise ValueError(f"waiver {index} expires must be an ISO expiry date")
        if expiry < today:
            raise ValueError(f"waiver {waiver['id']} expired on {expires}")

        waiver_id = waiver["id"]
        _require_nonempty_string(waiver_id, f"waiver {index}.id")
        if waiver_id in actual_waivers:
            raise ValueError(f"duplicate waiver id {waiver_id}")
        actual_waivers[waiver_id] = (
            waiver["context"],
            waiver["file"],
            waiver["test_title"],
            waiver["reason"],
            waiver["expected_count"],
            waiver["rationale"],
            expires,
        )

    if actual_waivers != EXPECTED_WAIVERS:
        raise ValueError("installed skip waivers do not match approved policy contract")
    if data["r_cmd_check"] != {"allowed_notes": []}:
        raise ValueError("r_cmd_check must have an empty allowed_notes list")
    if data["numerical_backends"] != {"mkl": EXPECTED_MKL_POLICY}:
        raise ValueError("MKL policy does not match approved library patterns")
    for pattern in (
        EXPECTED_MKL_POLICY["required_library_patterns"]
        + EXPECTED_MKL_POLICY["forbidden_library_patterns"]
    ):
        re.compile(pattern)
    return data


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {Path(argv[0]).name} MANIFEST", file=sys.stderr)
        return 2
    try:
        data = load_manifest(argv[1])
        validate_manifest(data)
        policy_path = Path(argv[1]).with_name("check-policy.yml")
        policy = load_manifest(policy_path)
        validate_policy(policy)
    except ValueError as error:
        print(f"manifest invalid: {error}", file=sys.stderr)
        return 1

    rows = data["coverage"]
    primary = sum(row["group"] == "primary" for row in rows)
    additional = sum(row["group"] == "additional" for row in rows)
    uncovered = sum(row["state"] == "uncovered" for row in rows)
    unit = len(data["unit_lanes"])
    results = len(_coverage_result_keys(rows) + _unit_result_keys(data["unit_lanes"]))
    print(
        f"manifest valid: primary={primary}, additional={additional}, "
        f"uncovered={uncovered}, unit={unit}, results={results}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
