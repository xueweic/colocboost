import hashlib
import io
import json
import os
import stat
import sys
import tarfile
from pathlib import Path

import pytest


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(CI_DIR))

from artifact_contract import create_metadata  # noqa: E402
from run_vnu import run_vnu  # noqa: E402
from verify_file_contract import verify_bound_file  # noqa: E402


SOURCE_SHA = "a" * 40
EVENT_SHA = "b" * 40


def executable(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def contract(tmp_path):
    source = tmp_path / "contract"
    source.mkdir()
    tarball = source / "colocboost_1.0.9.tar.gz"
    description = b"Package: colocboost\nVersion: 1.0.9\n"
    with tarfile.open(tarball, "w:gz") as archive:
        root = tarfile.TarInfo("colocboost")
        root.type = tarfile.DIRTYPE
        archive.addfile(root)
        info = tarfile.TarInfo("colocboost/DESCRIPTION")
        info.size = len(description)
        archive.addfile(info, io.BytesIO(description))
    metadata = source / "source-metadata.json"
    document = create_metadata(
        tarball, metadata, source_sha=SOURCE_SHA, event_sha=EVENT_SHA
    )
    return tarball, metadata, document["sha256"]


def invoke(tmp_path, *, dispatcher_body=None, mutation=None, environment=None):
    tarball, metadata, digest = contract(tmp_path)
    r = executable(tmp_path / "bin" / "R", "raise SystemExit(0)\n")
    if dispatcher_body is None:
        dispatcher_body = (
            "import pathlib\n"
            "pathlib.Path('pkg.html').write_text('<html>valid</html>')\n"
            "print('No validation errors')\n"
        )
    dispatcher = executable(tmp_path / "bin" / "vnu.sh", dispatcher_body)
    row = {
        "id": "vnu", "driver": "native-wrapper",
        "system_r": os.fspath(r), "vnu_path": os.fspath(dispatcher),
        "vnu_sha256": hashlib.sha256(dispatcher.read_bytes()).hexdigest(),
    }
    library = tmp_path / "check-library"
    library.mkdir()
    (library / "colocboost").mkdir()
    file_evidence = tmp_path / "dispatcher-proof.json"
    verify_bound_file(
        row, environment_id="vnu", purpose="vnu-dispatcher", path=dispatcher,
        source_sha=SOURCE_SHA, event_sha=EVENT_SHA, tarball_sha256=digest,
        output=file_evidence,
    )
    if mutation:
        mutation(row, dispatcher, file_evidence)
    selected_environment = dict(os.environ if environment is None else environment)
    if environment is None:
        selected_environment["PATH"] = os.fspath(r.parent)
    selected_environment.setdefault("R_LIBS_USER", os.fspath(library))
    work = tmp_path / "vnu-source"
    output = tmp_path / "vnu-proof.json"
    result = run_vnu(
        row, environment_id="vnu", dispatcher=dispatcher,
        r_executable=r, tarball=tarball, metadata=metadata,
        source_sha=SOURCE_SHA, event_sha=EVENT_SHA,
        tarball_sha256=digest, work_dir=work, library=library,
        file_evidence=file_evidence, output=output,
        environment=selected_environment,
    )
    return result, output, work


def test_vnu_executes_checksum_bound_dispatcher_without_arguments_in_verified_source(tmp_path):
    proof, output, work = invoke(tmp_path)
    assert proof == json.loads(output.read_text())
    assert proof["status"] == "pass"
    assert proof["dispatcher_exit_code"] == 0
    assert proof["zero_bad_entries"] is True
    assert proof["stdout"].strip() == "No validation errors"
    assert Path(proof["package_root"]).parent == work
    assert proof["pkg_html_size"] > 0


@pytest.mark.parametrize("kind", ["wrong-path", "wrong-hash", "nonzero", "no-output", "no-html", "empty-html", "symlink-html", "stale-file-proof"])
def test_vnu_fails_closed_on_dispatch_or_identity_ambiguity(tmp_path, kind):
    body = None
    mutation = None
    if kind == "wrong-path":
        mutation = lambda row, dispatcher, evidence: row.__setitem__("vnu_path", "/tmp/vnu.sh")
    elif kind == "wrong-hash":
        mutation = lambda row, dispatcher, evidence: row.__setitem__("vnu_sha256", "0" * 64)
    elif kind == "nonzero":
        body = "print('validator failed')\nraise SystemExit(7)\n"
    elif kind == "no-output":
        body = "import pathlib\npathlib.Path('pkg.html').write_text('<html/>')\n"
    elif kind == "no-html":
        body = "print('validator ran')\n"
    elif kind == "empty-html":
        body = "import pathlib\npathlib.Path('pkg.html').write_text('')\nprint('validator ran')\n"
    elif kind == "symlink-html":
        body = "import pathlib\np=pathlib.Path('target'); p.write_text('x'); pathlib.Path('pkg.html').symlink_to(p)\nprint('validator ran')\n"
    else:
        def mutation(row, dispatcher, evidence):
            document = json.loads(evidence.read_text())
            document["source_sha"] = "d" * 40
            evidence.write_text(json.dumps(document))

    with pytest.raises(ValueError):
        invoke(tmp_path, dispatcher_body=body, mutation=mutation)


def test_vnu_rejects_path_r_mismatch_and_reused_workdir(tmp_path):
    first = tmp_path / "first"
    first.mkdir()
    wrong_r = executable(first / "wrong" / "R", "raise SystemExit(0)\n")
    with pytest.raises(ValueError, match="PATH"):
        invoke(first, environment={"PATH": os.fspath(wrong_r.parent)})

    second = tmp_path / "second"
    second.mkdir()
    stale = second / "vnu-source"
    stale.mkdir()
    (stale / "stale").write_text("old")
    with pytest.raises(ValueError, match="initially empty"):
        invoke(second)
