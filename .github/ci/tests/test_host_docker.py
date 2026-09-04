import json
import sys
from pathlib import Path

import pytest

CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from run_host_docker import _check_counts, docker_argv, validate_proxy_dockerfile  # noqa: E402
from verify_backend_runtime import validate as validate_backend_runtime  # noqa: E402


def test_direct_r_cmd_check_is_network_isolated_and_explicit_mounts(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    argv = docker_argv(
        "ghcr.io/bastistician/rcheck-musl@sha256:" + "a" * 64,
        mounts=[(source, "/input", True), (tmp_path, "/work", False)],
        command=["/R/bin/R", "CMD", "check", "--as-cran", "/input/pkg.tar.gz"],
        user="1000:1000",
    )
    assert argv[:8] == ["docker", "run", "--rm", "--network", "none", "--user", "1000:1000", "--workdir"]
    assert "/R/bin/R" in argv and "/usr/local/bin/r-check" not in argv
    assert any("dst=/input,readonly" in value for value in argv)


def test_network_prep_is_explicitly_distinct_from_check(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    none = docker_argv("local/colocboost-blis-proxy", mounts=[(source, "/input", True)], command=["R"])
    enabled = docker_argv("local/colocboost-blis-proxy", mounts=[(source, "/input", True)], command=["R"], network="default")
    assert none[4] == "none"
    assert enabled[4] == "default"


def test_proxy_identity_drift_fails_closed(tmp_path):
    dockerfile = tmp_path / "Dockerfile"
    dockerfile.write_text("FROM docker.io/library/fedora@sha256:" + "0" * 64)
    with pytest.raises(ValueError, match="missing contract"):
        validate_proxy_dockerfile(dockerfile, "blis")


def test_external_tag_and_root_readonly_identity_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="immutable digest"):
        docker_argv("ghcr.io/example/rcheck:latest", mounts=[(tmp_path, "/work", False)], command=["R"])


def test_check_log_requires_exact_clean_ending_and_counts_diagnostics(tmp_path):
    log = tmp_path / "00check.log"
    log.write_text("* checking ... WARNING\n* DONE\nStatus: 1 WARNING\n")
    counts, clean = _check_counts(log)
    assert counts["warnings"] == 2
    assert clean is False
    log.write_text("* checking ... OK\n* DONE\nStatus: OK\n")
    assert _check_counts(log)[1] is True


def test_musl_runtime_proof_is_pinned_and_structured(tmp_path):
    proof = tmp_path / "runtime.json"
    document = {
        "kind": "backend-runtime", "backend": "musl", "r_binary": "/R/bin/R",
        "r_home": "/R/lib/R", "r_version": "R version 4.6.1",
        "alpine_version": "3.24.1", "musl": True,
        "architecture": "x86_64", "locale": "C.UTF-8",
    }
    proof.write_text(json.dumps(document))
    assert validate_backend_runtime(proof, "musl")["musl"] is True
    document["alpine_version"] = "latest"
    proof.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="pinned image contract"):
        validate_backend_runtime(proof, "musl")
