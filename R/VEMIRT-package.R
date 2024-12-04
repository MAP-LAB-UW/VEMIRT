#' VEMIRT: A Package for High-Dimensional IRT Models
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
#'   \item \code{\link{E2PL_iw}} conducts importance sampling to correct bias for M2PL analysis
#'   \item \code{\link{E3PL_sgvem_rot}} conducts stochastic GVEM to futher imporve the computational effficiency for exploratory M3PL analysis
#'   \item \code{\link{E3PL_sgvem_lasso}} conducts M3PL Analysis with Lasso penalty
#'   \item \code{\link{E3PL_sgvem_adaptlasso}} conducts M3PL Analysis with adaptive Lasso penalty
#' }
#' @section Confirmatory factor analysis:
#' \itemize{
#'   \item \code{\link{C2PL_gvem}} conducts GVEM for confirmatory M2PL analysis
#'   \item \code{\link{C2PL_bs}} conducts bootstrap sampling to correct bias and produce standard errors for confirmatory M2PL analysis
#'   \item \code{\link{C2PL_iw}} conducts importance sampling to correct bias for M2PL analysis
#'   \item \code{\link{C2PL_iw2}} conducts IW-GVEM for confirmatory M2PL analysis (alternative implementation to \code{\link{C2PL_iw}})
#'   \item \code{\link{C3PL_sgvem}} conducts stochastic GVEM for confirmatory M3PL analysis
#' }
#' @section Differential item functioning analysis:
#' \itemize{
#'   \item \code{\link{D2PL_em}} conducts DIF analysis for M2PL models using EM algorithms
#'   \item \code{\link{D2PL_pair_em}} conducts DIF analysis for 2PL models using EM algorithms with group pairwise truncated \eqn{L_1} penalty
#'   \item \code{\link{D2PL_gvem}} conducts DIF analysis for M2PL models using GVEM algorithms
#'   \item \code{\link{D2PL_lrt}} conducts DIF analysis for M2PL models using the likelihood ratio test
#' }
#' @section Shiny app for VEMIRT:
#' \itemize{
#'   \item \code{\link{shinyVEMIRT}} Run the shiny app for VEMIRT
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
#' @import mvQuad
## usethis namespace: end
NULL
