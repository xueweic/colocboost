# ColocBoost for multi-context colocalization in molecular QTL and GWAS studies

This R package implements ColocBoost --- motivated and designed for colocalization analysis ([first formulated here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)) of multiple genetic association studies --- as a multi-task learning approach to variable selection regression with highly correlated predictors and sparse effects, based on frequentist statistical inference. It provides statistical evidence to identify which subsets of predictors have non-zero effects on which subsets of response variables.

Temporary usage before packaging (v0.1.0-alpha)
- clone the repo to your local folder, then
```bash
cd colocboost
R --slave -e "devtools::install()"
```
- To run FineBoost, you need `colocboost(X=X, Y=y)`, where X and y are the same as `susie(X,y)`
- To run ColocBoost we suggest using [this pipeline wrapper](https://github.com/StatFunGen/pecotmr/blob/main/R/colocboost_pipeline.R) to manage multiple data-sets mixing individual level and summary statistics data.

