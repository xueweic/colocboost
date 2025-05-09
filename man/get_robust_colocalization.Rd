% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocboost_output.R
\name{get_robust_colocalization}
\alias{get_robust_colocalization}
\title{Recalibrate and summarize robust colocalization events.}
\source{
See detailed instructions in our tutorial portal:
\url{https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html}
}
\usage{
get_robust_colocalization(
  cb_output,
  cos_npc_cutoff = 0.5,
  npc_outcome_cutoff = 0.2,
  pvalue_cutoff = NULL,
  weight_fudge_factor = 1.5,
  coverage = 0.95
)
}
\arguments{
\item{cb_output}{Output object from \code{colocboost} analysis}

\item{cos_npc_cutoff}{Minimum threshold of normalized probability of colocalization (NPC) for CoS.}

\item{npc_outcome_cutoff}{Minimum threshold of normalized probability of colocalized traits in each CoS.}

\item{pvalue_cutoff}{Maximum threshold of marginal p-values of colocalized variants on colocalized traits in each CoS.}

\item{weight_fudge_factor}{The strength to integrate weight from different outcomes, default is 1.5}

\item{coverage}{A number between 0 and 1 specifying the \dQuote{coverage} of the estimated colocalization confidence sets (CoS) (default is 0.95).}
}
\value{
A \code{"colocboost"} object with some or all of the following elements:

\item{cos_summary}{A summary table for colocalization events.}
\item{vcp}{The variable colocalized probability for each variable.}
\item{cos_details}{A object with all information for colocalization results.}
\item{data_info}{A object with detailed information from input data}
\item{model_info}{A object with detailed information for colocboost model}
\item{ucos_from_cos}{A object with information for trait-specific effects if exists after removing weaker signals.}
}
\description{
\code{get_robust_colocalization} get the colocalization by discarding the weaker colocalization events or colocalized outcomes
}
\examples{
# colocboost example
set.seed(1)
N <- 1000
P <- 100
# Generate X with LD structure
sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
X <- MASS::mvrnorm(N, rep(0, P), sigma)
colnames(X) <- paste0("SNP", 1:P)
L <- 3
true_beta <- matrix(0, P, L)
true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
Y <- matrix(0, N, L)
for (l in 1:L) {
  Y[, l] <- X \%*\% true_beta[, l] + rnorm(N, 0, 1)
}
res <- colocboost(X = X, Y = Y)
res$cos_details$cos$cos_index
filter_res <- get_robust_colocalization(res, cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2)
filter_res$cos_details$cos$cos_index

}
\seealso{
Other colocboost_inference: 
\code{\link{get_ambiguous_colocalization}()},
\code{\link{get_colocboost_summary}()}
}
\concept{colocboost_inference}
