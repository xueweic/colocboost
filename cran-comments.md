## colocboost 1.0.10 release comments

This is a CRAN-requested patch update to colocboost 1.0.10.

This patch includes:

* Fixed the CRAN-reported MKL test issue in the CoS robustness tests.
  The tests now use a stronger simulation setting and aligned inputs for
  `get_robust_ucos()` and `get_ucos_evidence()`, so they no longer depend on
  weak or platform-sensitive simulated signals.

* R source code changes were limited to targeted computational optimizations and robustness fixes. No package dependencies were changed.


## Advanced documentation and tutorials

* Added conceptual and advanced-scenario vignettes to help users understand multi-trait colocalization events and interpret CoS, VCP, and NPC.
* Added practical illustrations of multiple causal variants and weaker disease GWAS signals, together with supporting figures and guidance.

## R CMD check results

There is one NOTE about installed package size:

* checking installed package size ... NOTE
    installed size is  7.1Mb
    sub-directories of 1Mb or more:
      data   2.0Mb
      doc    4.1Mb

This NOTE is expected. The increase from the previous release is primarily due to two substantive documentation vignettes that introduce the conceptual framework for multi-trait colocalization and provide tutorials on advanced analysis scenarios in ColocBoost. These vignettes include six high-level summary figures: three explain the key conceptual definitions and their relationships to existing methods, and three provide empirical illustrations based on representative simulation results. The figures are embedded in the rendered, self-contained vignettes in the `doc` directory to help users understand the methodology and interpret ColocBoost results. The `data` directory contains reduced example datasets used in reproducible tutorials. No external data are downloaded during examples or vignette rendering.
