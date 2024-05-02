#' VEMIRT: A package for high-dimensional IRT models
#'
#' VEMIRT is created to assist researchers to conduct exploratory and confirmatory
#' multidimensional item response theory (MIRT) analysis and cooresponding item
#' differential functioning (DIF) analysis. The core computation
#' engine of VEMIRT is a family of Gaussian Variational EM algorithms that are
#' considerably more efficient than currently available algorithms in other
#' software packages, especially when the number of latent factors exceeds four.
#'
#' @section Identifying the number of factors:
#' \code{\link{pa_poly}} identifies the number of factors via parallel analysis.
#' @section Exploratory factor analysis:
#' \itemize{
#'   \item \code{\link{gvem_2PLEFA_rot}} conducts M2PL Analysis with post-hoc rotation (Promax & CF-Quartimax)
#'   \item \code{\link{gvem_2PLEFA_lasso}} conducts M2PL Analysis with Lasso penalty
#'   \item \code{\link{gvem_2PLEFA_adaptlasso}} conducts M2PL Analysis with adaptive Lasso penalty
#'   \item \code{\link{sgvem_3PLEFA_rot}} conducts stochastic GVEM to futher imporve the computational effficiency for exploratory M3PL analysis
#'   \item \code{\link{sgvem_3PLEFA_lasso}} conducts M3PL Analysis with Lasso penalty
#'   \item \code{\link{sgvem_3PLEFA_adaptlasso}} conducts M3PL Analysis with adaptive Lasso penalty
#' }
#' @section Confirmatory factor analysis:
#' \itemize{
#'   \item \code{\link{gvem_2PLCFA}} conducts GVEM for confirmatory M2PL analysis
#'   \item \code{\link{sgvem_3PLCFA}} conducts stochastic GVEM for confirmatory M3PL analysis
#'   \item \code{\link{bs_2PLCFA}} conducts bootstrap sampling to correct bias and produce standard errors for confirmatory M2PL analysis
#'   \item \code{\link{importanceSampling}} conducts importance sampling to correct bias for M2PL analysis
#' }
#' @section Differential item functioning analysis:
#' \itemize{
#'   \item \code{\link{em_DIF}} conducts DIF analysis for M2PL models using EM algorithms
#'   \item \code{\link{gvemm_DIF}} conducts DIF analysis for M2PL models using GVEMM algorithms
#'   \item \code{\link{lrt_DIF}} conducts DIF analysis for M2PL models using the likelihood ratio test
#' }
# #' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib VEMIRT
#' @import Rcpp RcppArmadillo
#' @importFrom psych Promax
#' @importFrom GPArotation cfQ
#' @importFrom polycor polychor
#' @import tibble
#' @import MASS
#' @import Matrix
#' @import testit
#' @importFrom abind abind
#' @import rlang
#' @import mirt
#' @import torch
## usethis namespace: end
NULL
