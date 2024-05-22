#' Response data set for M2PL
#'
#' The response data set is simulated based on a
#' between-item M2PL model with 5 factors. The true factor correlations are set
#' as 0.1.
#' @format A data frame with 2000 respondents and 75 items
"exampleData_2pl"

#' True item parameters for M2PL
'exampleItem_2pl'

#' Factor-loading indicator matrix for M2PL-CFA
#'
#' The factor-loading indicator matrix can be used as an input for
#' confirmatory factor analysis.
#' @format A data frame with 75 items and 5 factors
"exampleIndic_cfa2pl"

#' Factor-loading indicator matrix for M2PL-EFA with lasso/ adaptive penalty under constraint 1
#'
#' The factor-loading indicator matrix can be used as an input for
#' exploratory factor analysis with lasso/ adaptive lasso penalty under constraint 1.
#' @format A data frame with 75 items and 5 factors. Items 1, 16, 31, 46 and 61 can
#' be combined as an identity matrix to satisfy constraint 1
"exampleIndic_efa2pl_c1"

#' Factor-loading indicator matrix for M2PL-EFA with lasso/ adaptive penalty under constraint 2
#'
#' The factor-loading indicator matrix can be used as an input for
#' exploratory factor analysis with lasso/ adaptive lasso penalty for constraint 1.
#' @format A data frame with 75 items and 5 factors. Items 1, 16, 31, 46 and 61 can
#' be combined as a triangular matrix to satisfy constraint 2
"exampleIndic_efa2pl_c2"

#' Response data set for M3PL
#'
#' The response data set is simulated based on a
#' within-item M3PL model with 3 factors. The true factor correlations are set
#' as 0.1.
#' @format A data frame with 2000 respondents and 45 items
"exampleData_3pl"

#' True item parameters for M3PL
'exampleItem_3pl'

#' Factor-loading indicator matrix for M3PL-CFA
#'
#' The factor-loading indicator matrix can be used as an input for
#' confirmatory factor analysis.
#' @format A data frame with 45 items and 3 factors
"exampleIndic_cfa3pl"

#' Factor-loading indicator matrix for M3PL-EFA with lasso/ adaptive penalty under constraint 1
#'
#' The factor-loading indicator matrix can be used as an input for
#' exploratory factor analysis with lasso/ adaptive lasso penalty under constraint 1.
#' @format A data frame with 45 items and 3 factors. Items 1, 16, and 19 can
#' be combined as an identity matrix to satisfy constraint 1
"exampleIndic_efa3pl_c1"

#' Factor-loading indicator matrix for M3PL-EFA with lasso/ adaptive penalty under constraint 2
#'
#' The factor-loading indicator matrix can be used as an input for
#' exploratory factor analysis with lasso/ adaptive lasso penalty for constraint 1.
#' @format A data frame with 45 items and 3 factors. Items 1, 16, and 19 can
#' be combined as a triangular matrix to satisfy constraint 2
"exampleIndic_efa3pl_c2"

#' Simulated Data Set for DIF Analysis
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators\cr
#' \tab\cr
#' \code{group}\tab Group indicators\cr
#' \tab\cr
#' \code{j}\tab Number of DIF items (the first \code{j} items have DIF)\cr
#' \tab\cr
#' \code{params}\tab A list of true parameters used for generating the item responses:\cr
#' \tab\cr
#' \code{ ...$a}\tab Slopes\cr
#' \tab\cr
#' \code{ ...$b}\tab Negated intercepts\cr
#' \tab\cr
#' \code{ ...$theta}\tab Latent traits\cr
#' }
'exampleDIF'
