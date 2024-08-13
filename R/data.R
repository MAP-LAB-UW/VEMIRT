#' Simulated Data Set for Confirmatory M2PL Analysis
#'
#' Responses are simulated based on an M2PL model with 2 factors. The true factor correlations are set as 0.8.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
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
#' Responses are simulated based on an M3PL model with 2 factors. The true factor correlations are set as 0.8.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'C3PL_data'

#' Simulated Data Set for Exploratory M2PL Analysis Under C1 Constraint
#'
#' Responses are simulated based on an M2PL model with 3 factors. The true factor correlations are set as 0.5.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C1'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors (\code{NULL} under \code{C1})\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'E2PL_data_C1'

#' Simulated Data Set for Exploratory M2PL Analysis Under C2 Constraint
#'
#' Responses are simulated based on an M2PL model with 3 factors. The true factor correlations are set as 0.5.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C2'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'E2PL_data_C2'

#' Simulated Data Set for Exploratory M3PL Analysis Under C1 Constraint
#'
#' Responses are simulated based on an M3PL model with 3 factors. The true factor correlations are set as 0.5.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C1'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors (\code{NULL} under \code{C1})\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'E3PL_data_C1'

#' Simulated Data Set for Exploratory M3PL Analysis Under C2 Constraint
#'
#' Responses are simulated based on an M3PL model with 3 factors. The true factor correlations are set as 0.5.
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @format A list of components of the data set:
#' \tabular{ll}{
#' \code{ data}\tab Item responses\cr
#' \tab\cr
#' \code{model}\tab Loading indicators for (adaptive) lasso penalty\cr
#' \tab\cr
#' \code{constrain}\tab Constraint for model identification (\code{'C2'})\cr
#' \tab\cr
#' \code{non_pen}\tab Index of an item that is associated with all the factors\cr
#' \tab\cr
#' \code{params}\tab True parameters used for generating the item responses
#' }
'E3PL_data_C2'

#' Simulated Data Set for DIF M2PL Analysis
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
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
