# ColocBoost for multi-context colocalization in molecular QTL and GWAS studies

This R package implements ColocBoost --- motivated and designed for colocalization analysis ([first formulated here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)) of multiple genetic association studies --- as a multi-task learning approach to variable selection regression with highly correlated predictors and sparse effects, based on frequentist statistical inference. It provides statistical evidence to identify which subsets of predictors have non-zero effects on which subsets of response variables.

Temporary usage before packaging (v0.1.0-alpha)
- Download the source code in `R/` folder. 
- In R, you should source code using `for(file in list.files("R/", full.names = T)) {source(file)}`.
- To run FineBoost, you need `colocboost(X=X, Y=y)`, where X and y are the same as `susie(X,y)`
