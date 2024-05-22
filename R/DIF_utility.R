# Written by Weicong Lyu

#' @export
`%*%.torch_tensor` <- function(x, y) {
  torch_matmul(x, y)
}

#' @export
t.torch_tensor <- function(x) {
  torch_transpose(x, -1, -2)
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

proximal.adam <- function(x, state, params, lambda) {
  betas.t <- params$betas ^ state$step
  lr <- params$lr / (1 - betas.t[1]) / (sqrt(state$exp_avg_sq / (1 - betas.t[2])) + params$eps)
  x$set_data(prox(x, lr * lambda))
}

distance <- function(x, y) {
  mapply(function(x, y) as.array(max(abs(x - y))), x, y)
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

best_vemirt_DIF <- function(all, N, criterion) {
  ic <- if (is.character(criterion))
    sapply(all, `[[`, criterion)
  else
    sapply(all, function(fit) {
      with(fit, -2 * ll + criterion * l0 * log(N) * log(log(N)))
    })
  best <- which.min(ic)
  fit <- all[[best]]
  lambda0 <- range(sapply(all, `[[`, 'lambda0'))
  if (length(ic) > 1)
    if (fit$lambda0 == lambda0[1] && fit$lambda0 > 0)
      warning('Optimal lambda0 may be less than ', all[[best]]$lambda0, '.')
    else if (fit$lambda0 == lambda0[2] && fit$l0 > 0)
      warning('Optimal lambda0 may be greater than ', all[[best]]$lambda0, '.')
  best
}

new_vemirt_DIF <- function(niter0, all, N, criterion) {
  best <- best_vemirt_DIF(all, N, criterion)
  structure(lst(niter0, all, N, best, fit = all[[best]]), class = 'vemirt_DIF')
}

#' Extract Parameter Estimates from DIF Analysis
#'
#' @param object An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage coef(object, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
coef.vemirt_DIF <- function(object, criterion = NULL, ...) {
  if (length(criterion) != 1)
    object$fit
  else
    object$all[[best_vemirt_DIF(object$all, object$N, criterion)]]
}

#' Print DIF Items by Group
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
print.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
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
  invisible(fit)
}

#' Summarize DIF Items
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
summary.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
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
  invisible(dat)
}
