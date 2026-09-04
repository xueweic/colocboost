#!/usr/bin/env python3
"""Execute the complete unit suite inside a pinned host-Docker image."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from run_host_docker import docker_argv


def run_unit(*, image: str, package: str, output: str, environment_id: str, policy: str, user: str = "0:0", platform: str | None = None, runtime_proof: str | None = None) -> int:
    if not runtime_proof or not Path(runtime_proof).is_file():
        raise ValueError("BLIS unit requires the verified runtime proof")
    proof = json.loads(Path(runtime_proof).read_text(encoding="utf-8"))
    if proof.get("backend") != "blis" or proof.get("matrix_operation") is not True or proof.get("same_process_maps") is not True or proof.get("loaded_blis") is not True:
        raise ValueError("BLIS unit runtime proof is incomplete or permits fallback")
    package_path = Path(package).absolute()
    output_path = Path(output).absolute()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    command = docker_argv(
        image,
        mounts=[(package_path, "/package", True), (output_path.parent, "/output", False)],
        command=[
            "/opt/R/devel-blis/bin/Rscript" if "blis" in image else "/opt/R/devel-noomp/bin/Rscript",
            "--vanilla",
            "/package/.github/ci/run-unit-tests.R",
            "--package=/package",
            "--load-package=source",
            f"--context={environment_id}",
            "--policy=/package/.github/ci/check-policy.yml",
            f"--output=/output/{output_path.name}",
        ],
        user=user,
        platform=platform,
        environment={"R_PROFILE_USER": "/dev/null", "R_ENVIRON_USER": "/dev/null"},
    )
    completed = subprocess.run(command, check=False, shell=False)
    if not output_path.is_file():
        raise ValueError("Docker unit run did not produce its diagnostic")
    return completed.returncode


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--package", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--environment-id", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--user", default="0:0")
    parser.add_argument("--platform")
    parser.add_argument("--runtime-proof")
    args = parser.parse_args(argv)
    try:
        return run_unit(**vars(args))
    except (OSError, TypeError, ValueError, subprocess.SubprocessError) as error:
        print(f"host Docker unit run failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
