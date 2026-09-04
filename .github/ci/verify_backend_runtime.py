#!/usr/bin/env python3
"""Fail-closed validation for the structured BLIS/noOMP runtime proof."""
from __future__ import annotations

import json
import sys
from pathlib import Path


def validate(path: str | Path, backend: str) -> dict:
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("kind") != "backend-runtime":
        raise ValueError("runtime proof is not structured backend evidence")
    if document.get("backend") != backend:
        raise ValueError("runtime proof identity is incomplete")
    if backend == "blis":
        if document.get("same_process_maps") is not True:
            raise ValueError("BLIS runtime proof lacks same-process mappings")
        required = {"matrix_operation", "loaded_blis", "loaded_openblas", "loaded_mkl", "loaded_atlas", "serial"}
        if not required <= document.keys() or not document["matrix_operation"] or not document["loaded_blis"] or document["serial"] is not True:
            raise ValueError("BLIS runtime proof lacks real matrix/serial evidence")
        if any(document[name] for name in ("loaded_openblas", "loaded_mkl", "loaded_atlas")):
            raise ValueError("BLIS proof permits a fallback BLAS")
    elif backend == "noomp":
        if document.get("same_process_maps") is not True:
            raise ValueError("noOMP runtime proof lacks same-process mappings")
        flags = document.get("shlib_openmp_flags")
        if not isinstance(flags, dict) or any(flags.get(name) != "" for name in ("C", "CXX", "FFLAGS")):
            raise ValueError("noOMP Makeconf OpenMP flags are not empty")
        if document.get("rfast_crossprod") is not True or document.get("forbidden_elf_libraries") or document.get("forbidden_maps"):
            raise ValueError("noOMP runtime proof permits OpenMP fallback")
        if not document.get("hard_dependency_elf"):
            raise ValueError("noOMP hard-dependency ELF evidence is missing")
    elif backend == "musl":
        expected = {
            "r_binary": "/R/bin/R",
            "alpine_version": "3.24.1",
            "musl": True,
            "architecture": "x86_64",
            "locale": "C.UTF-8",
        }
        if any(document.get(key) != value for key, value in expected.items()):
            raise ValueError("musl runtime proof does not match the pinned image contract")
        if not str(document.get("r_home", "")).startswith("/R"):
            raise ValueError("musl runtime proof has the wrong R home")
        if not str(document.get("r_version", "")).startswith("R version "):
            raise ValueError("musl runtime proof lacks an R version")
    else:
        raise ValueError(f"unsupported backend: {backend}")
    return document


if __name__ == "__main__":
    try:
        validate(sys.argv[1], sys.argv[2])
    except (IndexError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"backend runtime proof rejected: {error}", file=sys.stderr)
        raise SystemExit(1)
