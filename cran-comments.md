## colocboost 1.0.9 release comments

This is a CRAN-requested patch update to colocboost 1.0.8.

This patch includes:

* Fixed the CRAN-reported macOS arm64 test issue in the uCoS robustness tests.
  The tests now use a stronger simulation setting and aligned inputs for
  `get_robust_ucos()` and `get_ucos_evidence()`, so they no longer depend on
  weak or platform-sensitive simulated signals.

No R source code or package dependency changes were made for this CRAN-requested
patch.

## CRAN-requested macOS arm64 fix

CRAN reported test failures for colocboost 1.0.8 on macOS arm64 and requested a
correction before 2026-06-21. The failing test checked robust trait-specific
uncolocalized event filtering and evidence calculation. The issue has been
addressed by strengthening the simulation used in the test suite and by using
matched `cb_obj` and `cb_res` inputs for the evidence check.

CRAN also pointed to the M1mac check service for arm64 issues:
https://www.stats.ox.ac.uk/pub/bdr/M1mac/README.txt. CRAN noted that this
service runs a much older OS/toolchain and that toolchain differences often
matter. This submission therefore fixes the test design itself, rather than
relying on a platform-specific workaround.

## R CMD check results

There is one NOTE about installed package size:

* checking installed package size ... NOTE
  installed size is 5.0 MB
  sub-directories of 1Mb or more:
    data 2.0 MB
    doc  1.9 MB

This NOTE is expected. The installed size is mainly due to reduced example datasets and rendered vignettes with figures. These files are kept to make the tutorials reproducible and self-contained for multi-trait colocalization workflows. No external data are downloaded during examples or vignette rendering.

## Previous comments

* This package implements methods described in our paper "ColocBoost"
  (Cao et al., 2025), now cited in DESCRIPTION.
* Previous CRAN-requested fixes addressed tarball size, user option handling,
  LICENSE metadata, and accepted package/domain terms in inst/WORDLIST.
* The examples and vignettes use small datasets to avoid long check times.
