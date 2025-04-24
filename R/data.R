#' Individual level data for 5 traits
#'
#' An example dataset with simulated genotypes and traits for 5 traits
#'
#' @format ## `Ind_5traits`
#' A list with 3 elements
#' \describe{
#'   \item{X}{List of genotype matrices}
#'   \item{Y}{List of traits}
#'   \item{true_effect_variants}{List of causal variants}
#' }
#' @source The Ind_5traits dataset contains 5 simulated phenotypes alongside corresponding genotype matrices.
#' The dataset is specifically designed for evaluating and demonstrating the capabilities of ColocBoost in multi-trait colocalization analysis 
#' with individual-level data. See Cao etc. 2025 for details.
#' Due to the file size limitation of CRAN release, this is a subset of simulated data. 
#' See full dataset in colocboost paper repo \url{https://github.com/StatFunGen/colocboost-paper}.
#' 
#' @family colocboost_data
"Ind_5traits"

#' Summary level data for 5 traits
#'
#' An example dataset with simulated statistics for 5 traits
#'
#' @format ## `Sumstat_5traits`
#' A list with 2 elements
#' \describe{
#'   \item{sumstat}{Summary statistics for 5 traits}
#'   \item{true_effect_variants}{List of causal variants}
#' }
#' @source The Sumstat_5traits dataset contains 5 simulated summary statistics, 
#' where it is directly derived from the Ind_5traits dataset using marginal association. 
#' The dataset is specifically designed for evaluating and demonstrating the capabilities of ColocBoost 
#' in multi-trait colocalization analysis with summary association data. See Cao etc. 2025 for details. 
#' Due to the file size limitation of CRAN release, this is a subset of simulated data. 
#' See full dataset in colocboost paper repo \url{https://github.com/StatFunGen/colocboost-paper}.
#' 
#' @family colocboost_data
"Sumstat_5traits"


#' Individual level data for 2 traits and 2 causal variants with heterogeneous effects
#'
#' An example dataset with simulated genotypes and traits for 2 traits and 2 common causal variants with heterogeneous effects
#'
#' @format ## `Heterogeneous_Effect`
#' A list with 3 elements
#' \describe{
#'   \item{X}{List of genotype matrices}
#'   \item{Y}{List of traits}
#'   \item{variant}{indices of two causal variants}
#' }
#' @source The Heterogeneous_Effect dataset contains 2 simulated phenotypes alongside corresponding genotype matrices.
#' There are two causal variants, both of which have heterogeneous effects on two traits.
#' Due to the file size limitation of CRAN release, this is a subset of simulated data to generate Figure 2b in Cao etc. 2025. 
#' See full dataset in colocboost paper repo \url{https://github.com/StatFunGen/colocboost-paper}.
#' 
#' @family colocboost_data
"Heterogeneous_Effect"


#' Individual level data for 2 traits and 2 causal variants with weaker effects for focal trait
#'
#' An example dataset with simulated genotypes and traits for 2 traits and 2 common causal variants with heterogeneous effects
#'
#' @format ## `Weaker_GWAS_Effect`
#' A list with 3 elements
#' \describe{
#'   \item{X}{List of genotype matrices}
#'   \item{Y}{List of traits}
#'   \item{variant}{indices of two causal variants}
#' }
#' @source The Weaker_GWAS_Effect dataset contains 2 simulated phenotypes alongside corresponding genotype matrices.
#' There are two causal variants, one of which has a weaker effect on the focal trait compared to the other trait.
#' Due to the file size limitation of CRAN release, this is a subset of simulated data to generate Figure 2b in Cao etc. 2025. 
#' See full dataset in colocboost paper repo \url{https://github.com/StatFunGen/colocboost-paper}.
#' 
#' @family colocboost_data
"Weaker_GWAS_Effect"


#' Individual level data for 2 traits and 2 causal variants, but the strongest marginal association is not causal
#'
#' An example dataset with simulated genotypes and traits for 2 traits and 2 common causal variants, but the strongest marginal association is not causal variant.
#'
#' @format ## `Non_Causal_Strongest_Marginal`
#' A list with 3 elements
#' \describe{
#'   \item{X}{List of genotype matrices}
#'   \item{Y}{List of traits}
#'   \item{variant}{indices of two causal variants}
#' }
#' @source The Non_Causal_Strongest_Marginal dataset contains 2 simulated phenotypes alongside corresponding genotype matrices.
#' There are two causal variants, but the strongest marginal association is not a causal variant.
#' Due to the file size limitation of CRAN release, this is a subset of simulated data to generate Figure 2b in Cao etc. 2025. 
#' See full dataset in colocboost paper repo \url{https://github.com/StatFunGen/colocboost-paper}.
#' 
#' @family colocboost_data
"Non_Causal_Strongest_Marginal"



#' A real data example includes an ambiguous colocalization between eQTL and GWAS
#'
#' An example result from one of our real data applications, which shows an ambiguous colocalization between eQTL and GWAS.
#'
#' @format ## `Ambiguous_Colocalization`
#' A list with 2 elements
#' \describe{
#'   \item{ColocBoost_Results}{A `colocboost` output object}
#'   \item{SuSiE_Results}{Two `susie` output object for eQTL and GWAS}
#'   \item{COLOC_V5_Results}{A `coloc` output object}
#' }
#' @source The Ambiguous_Colocalization dataset contains a real data example from one of our real data applications,
#' which shows an ambiguous colocalization between eQTL and GWAS.
#' The dataset is specifically designed for evaluating and demonstrating the capabilities of ColocBoost in real data applications.
#' See details in tutorial vignette \url{https://statfungen.github.io/colocboost/articles/index.html}.
#' 
#' @family colocboost_data
"Ambiguous_Colocalization"
