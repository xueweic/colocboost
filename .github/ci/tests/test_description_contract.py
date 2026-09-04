from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def test_description_omits_stale_optional_date_field():
    lines = (ROOT / "DESCRIPTION").read_text(encoding="utf-8").splitlines()
    assert not any(line.startswith("Date:") for line in lines)
