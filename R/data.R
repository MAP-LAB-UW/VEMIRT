#' Simulated Data Set for Confirmatory M2PL Analysis
#'
#' Responses are simulated based on a between-item M2PL model with 5 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'C2PL_data'

#' Simulated Data Set for Confirmatory M3PL Analysis
#'
#' Responses are simulated based on a within-item M3PL model with 3 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'C3PL_data'

#' Simulated Data Set for Exploratory M2PL Analysis Under C1 Constraints
#'
#' Responses are simulated based on a between-item M2PL model with 5 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C1'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors (\code{NULL} under C1)
#' }
'E2PL_data_C1'

#' Simulated Data Set for Exploratory M2PL Analysis Under C2 Constraints
#'
#' Responses are simulated based on a between-item M2PL model with 5 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C2'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors
#' }
'E2PL_data_C2'

#' Simulated Data Set for Exploratory M3PL Analysis Under C1 Constraints
#'
#' Responses are simulated based on a within-item M3PL model with 3 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C1'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors (\code{NULL} under C1)
#' }
'E3PL_data_C1'

#' Simulated Data Set for Exploratory M3PL Analysis Under C2 Constraints
#'
#' Responses are simulated based on a within-item M3PL model with 3 factors. The true factor correlations are set as 0.1.
#'
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C2'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors
#' }
'E3PL_data_C2'

#' Simulated Data Set for DIF M2PL Analysis
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
'D2PL_data'
