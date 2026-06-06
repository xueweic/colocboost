## colocboost 1.0.8 release comments

This is an update to colocboost 1.0.7.

This release includes:

* Added X_ref support as a memory-efficient alternative to precomputed LD matrices for summary-statistics workflows.
* Added and refined robust post-filtering for colocalization and trait-specific uncolocalized events.
* Improved computational efficiency for repeated matrix products in reference-panel workflows.
* Improved plotting robustness for extreme association signals and coefficient or z-score displays.
* Updated documentation and vignettes, including the bioinformatics pipeline vignette.

## R CMD check results

There is one NOTE about installed package size:

* checking installed package size ... NOTE
  installed size is 5.0 MB
  sub-directories of 1Mb or more:
    data 2.0 MB
    doc  1.9 MB

This NOTE is expected. The installed size is mainly due to reduced example datasets and rendered vignettes with figures. These files are kept to make the tutorials reproducible and self-contained for multi-trait colocalization workflows. No external data are downloaded during examples or vignette rendering.

## Previous comments

* This is an important update for the latest release colocboost_1.0.5. We resolved an computational issue caused by the latest update.
* This package implements methods described in our paper "ColocBoost" (Cao et al., 2025), added in DESCRIPTION
* Fixed issues requested by CRAN in previous submission:
  - Reduced tarball less than 5 MB
  - Fixed reset users' options issues
  - Added proper COPYRIGHT HOLDER and ORGANIZATION to LICENSE
  - Added explanation of acronyms used in this package to inst/WORDLIST
* The examples and vignettes use small datasets to avoid long check times
