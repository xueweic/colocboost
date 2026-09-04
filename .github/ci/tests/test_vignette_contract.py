from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def test_animation_lookup_is_not_executed_when_r_cmd_check_sources_vignette_code():
    source = (ROOT / "vignettes" / "ColocBoost_Update.Rmd").read_text(
        encoding="utf-8"
    )
    assert "if (knitr::is_html_output()) {" in source
    assert 'knitr::include_graphics("../man/figures/ColocBoost_update.gif")' in source
