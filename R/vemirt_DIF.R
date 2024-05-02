# Written by Weicong Lyu

best_vemirt_DIF <- function(all, N, criterion) {
  ic <- if (is.character(criterion))
    sapply(all, `[[`, criterion)
  else
    sapply(all, function(fit) {
      with(fit, -2 * ll + criterion * l0 * log(N) * log(log(N)))
    })
  best <- which.min(ic)
  if (length(ic) > 1)
    if (best == 1)
      warning('Optimal lambda0 may be less than ', all[[best]]$lambda0, '.')
    else if (best == length(ic))
      warning('Optimal lambda0 may be greater than ', all[[best]]$lambda0, '.')
  best
}

new_vemirt_DIF <- function(all, N, criterion) {
  best <- best_vemirt_DIF(all, N, criterion)
  structure(list(N = N, fit = all[[best]], best = best, all = all), class = 'vemirt_DIF')
}

#' Extract Parameter Estimates from DIF Analysis
#'
#' @param object An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage coef(object, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
coef.vemirt_DIF <- function(object, criterion = NULL, ...) {
  if (length(criterion) != 1)
    object$fit
  else
    object$all[[best_vemirt_DIF(object$all, object$N, criterion)]]
}

#' Print DIF Items
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}
#' @export
print.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  print(do.call(cbind, lapply(2:nrow(fit$beta), function(k) {
    dif <- as.data.frame(ifelse(cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0, 'X', ''))
    colnames(dif) <- paste0(c(paste0('a', 1:(ncol(dif) - 1), ','), 'b'), k)
    rownames(dif) <- paste0(1:nrow(dif))
    dif
  })))
  cat(sep = '', '* lambda0 = ', fit$lambda0, '\n')
  invisible(fit)
}
