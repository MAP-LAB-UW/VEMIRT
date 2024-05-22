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
#'   \item \code{\link{E2PL_gvem_rot}} conducts M2PL Analysis with post-hoc rotation (Promax & CF-Quartimax)
#'   \item \code{\link{E2PL_gvem_lasso}} conducts M2PL Analysis with Lasso penalty
#'   \item \code{\link{E2PL_gvem_adaptlasso}} conducts M2PL Analysis with adaptive Lasso penalty
#'   \item \code{\link{E3PL_sgvem_rot}} conducts stochastic GVEM to futher imporve the computational effficiency for exploratory M3PL analysis
#'   \item \code{\link{E3PL_sgvem_lasso}} conducts M3PL Analysis with Lasso penalty
#'   \item \code{\link{E3PL_sgvem_adaptlasso}} conducts M3PL Analysis with adaptive Lasso penalty
#' }
#' @section Confirmatory factor analysis:
#' \itemize{
#'   \item \code{\link{C2PL_gvem}} conducts GVEM for confirmatory M2PL analysis
#'   \item \code{\link{C3PL_sgvem}} conducts stochastic GVEM for confirmatory M3PL analysis
#'   \item \code{\link{C2PL_bs}} conducts bootstrap sampling to correct bias and produce standard errors for confirmatory M2PL analysis
#'   \item \code{\link{importanceSampling}} conducts importance sampling to correct bias for M2PL analysis
#' }
#' @section Differential item functioning analysis:
#' \itemize{
#'   \item \code{\link{DIF_em}} conducts DIF analysis for M2PL models using EM algorithms
#'   \item \code{\link{DIF_gvem}} conducts DIF analysis for M2PL models using GVEMM algorithms
#'   \item \code{\link{DIF_lrt}} conducts DIF analysis for M2PL models using the likelihood ratio test
#' }
# #' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib VEMIRT
#' @import Rcpp RcppArmadillo
#' @importFrom psych Promax
#' @importFrom GPArotation cfQ
#' @importFrom polycor polychor
#' @importFrom tibble tibble
#' @importFrom tibble lst
#' @importFrom mvnfast dmvn
#' @import MASS
#' @import Matrix
#' @import testit
#' @importFrom abind abind
#' @import mirt
#' @import torch
## usethis namespace: end
NULL
