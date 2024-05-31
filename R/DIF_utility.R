# Written by Weicong Lyu

#' @export
`%*%.torch_tensor` <- function(x, y) {
  torch_matmul(x, y)
}

#' @export
t.torch_tensor <- function(x) {
  torch_transpose(x, -1, -2)
}

keep <- function(...) {
  env <- parent.frame()
  env$ls.keep <- names(lst(...))
  with(parent.frame(), {
    rm(list = setdiff(ls(), ls.keep))
  })
}

NULL.tensor <- function() {
  torch_tensor(NULL, dtype = torch_float())
}

diagonal <- function(e) {
  torch_diagonal(e, dim1 = -1, dim2 = -2)
}

sym <- function(e) {
  (e + t(e)) / 2
}

prox <- function(x, lambda) {
  sign(x) * (abs(x) - lambda)$maximum(0)
}

groupsum <- function(n, x) {
  torch_stack(lapply(n, function(n) {
    x[n]$sum(1)
  }))
}

groupmean <- function(n, x) {
  torch_stack(lapply(n, function(n) {
    x[n]$mean(1)
  }))
}

distance <- function(x, y) {
  mapply(function(x, y) tryCatch(as.array(max(abs(x - y))), error = function(e) 0), x, y)
}

IC <- function(ll, l0, N, c) {
  -2 * ll + l0 * c(AIC = 2, BIC = log(N), GIC = c * log(N) * log(log(N)))
}

init.new <- function() {
  with(parent.frame(), {
    vars <- ls()
    tensors <- list()
    arrays <- lst(...)
    for (var in vars) {
      value <- get(var)
      if (inherits(value, 'torch_tensor'))
        tensors[[var]] <- value
      else
        arrays[[var]] <- value
    }
    function(var = NULL) {
      if (is.null(var))
        c(arrays, with_no_grad(lapply(tensors, function(tensor) {
          tensor$clone()$requires_grad_(tensor$requires_grad)
        })))
      else
        get(var)
    }
  })
}

check.vemirt_DIF <- function(all, fit, var) {
  var.range <- range(sapply(all, `[[`, var))
  if (var.range[1] < var.range[2]) {
    x <- fit[[var]]
    if (x > 0 && x == var.range[1])
      warning('Optimal ', var, ' may be less than ', x, '.')
    if (fit$l0 > 0 && x == var.range[2])
      warning('Optimal ', var, ' may be greater than ', x, '.')
  }
}

best.vemirt_DIF <- function(all, N, criterion) {
  ic <- if (is.character(criterion))
    sapply(all, `[[`, criterion)
  else
    sapply(all, function(fit) {
      with(fit, -2 * ll + criterion * l0 * log(N) * log(log(N)))
    })
  best <- which.min(ic)
  fit <- all[[best]]
  if (!is.null(fit$lambda0))
    check.vemirt_DIF(all, fit, 'lambda0')
  else {
    check.vemirt_DIF(all, fit, 'lambda')
    check.vemirt_DIF(all, fit, 'tau')
  }
  best
}

new.vemirt_DIF <- function(niter0, all, N, criterion) {
  best <- best.vemirt_DIF(all, N, criterion)
  structure(lst(niter0, all, N, best, fit = all[[best]]), class = 'vemirt_DIF')
}

#' Extract Parameter Estimates from DIF Analysis
#'
#' @param object An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage coef(object, criterion = NULL)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
coef.vemirt_DIF <- function(object, criterion = NULL, ...) {
  if (length(criterion) != 1)
    object$fit
  else
    object$all[[best.vemirt_DIF(object$all, object$N, criterion)]]
}

#' Print DIF Items by Group
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{DIF_em}}, \code{\link{DIF_gvem}}, \code{\link{DIF_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
print.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  if (is.null(fit$tau)) {
    K <- nrow(fit$beta)
    if (K == 1) {
      dat <- as.data.frame(t(cbind(fit$a, as.matrix(fit$b))))
      rownames(dat) <- c(paste0('a', 1:(nrow(dat) - 1)), 'b')
      colnames(dat) <- paste0(1:ncol(dat))
      print(dat, digits = 3)
    } else {
      dat <- do.call(rbind, lapply(2:K, function(k) {
        dif <- as.data.frame(t(ifelse(cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0, 'X', '')))
        rownames(dif) <- paste0(c(paste0('a', 1:(nrow(dif) - 1), ','), 'b'), k)
        colnames(dif) <- paste0(1:ncol(dif))
        dif
      }))
      print(dat, max = 100000L)
    }
    if (K > 1 || length(x$all) != 1 || fit$lambda0 != 0)
      cat(sep = '', '* lambda0 = ', fit$lambda0, '\n')
  }
  invisible(fit)
}

#' Summarize DIF Items
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @usage print(x, criterion = NULL)
#'
#' @seealso \code{\link{DIF_em}}, \code{\link{DIF_gvem}}, \code{\link{DIF_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
summary.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  if (is.null(fit$tau)) {
    K <- nrow(fit$beta)
    if (K == 1) {
      dat <- as.data.frame(t(cbind(fit$a, as.matrix(fit$b))))
      rownames(dat) <- c(paste0('a', 1:(nrow(dat) - 1)), 'b')
      colnames(dat) <- paste0(1:ncol(dat))
      print(dat, digits = 3)
    } else {
      dat <- (Reduce(`+`, lapply(2:nrow(fit$beta), function(k) {
        dif <- t((cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0) + 0)
        rownames(dif) <- c(paste0('a', 1:(nrow(dif) - 1)), 'b')
        colnames(dif) <- paste0(1:ncol(dif))
        dif
      })) != 0) + 0
      print(as.data.frame(ifelse(dat, 'X', '')), max = 100000L)
    }
    if (K > 1 || length(x$all) != 1 || fit$lambda0 != 0)
      cat(sep = '', '* lambda0 = ', fit$lambda0, '\n')
  }
  invisible(dat)
}
