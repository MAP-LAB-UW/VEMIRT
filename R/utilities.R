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

proximal.adam <- function(x, state, params, lambda) {
  betas.t <- params$betas ^ state$step
  lr <- params$lr / (1 - betas.t[1]) / (sqrt(state$exp_avg_sq / (1 - betas.t[2])) + params$eps)
  x$set_data(prox(x, lr * lambda))
}

groupsum <- function(n, x) {
  torch_stack(lapply(n, function(n) {
    x[n]$nansum(1)
  }))
}

groupmean <- function(n, x) {
  torch_stack(lapply(n, function(n) {
    x[n]$nanmean(1)
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

#' Extract Parameter Estimates from DIF 2PL Analysis
#'
#' @param object An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage coef(object, criterion = NULL)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{em_D2PL}}, \code{\link{gvemm_D2PL}}, \code{\link{lrt_D2PL}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
coef.vemirt_DIF <- function(object, criterion = NULL, ...) {
  if (length(criterion) != 1)
    object$fit
  else
    object$all[[best.vemirt_DIF(object$all, object$N, criterion)]]
}

#' Print DIF 2PL Items by Group
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL, max = 99999L, digits = 3, ...)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{D2PL_em}}, \code{\link{D2PL_gvem}}, \code{\link{D2PL_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
print.vemirt_DIF <- function(x, criterion = NULL, max = 99999L, digits = 3, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  if (is.null(fit$tau)) {
    K <- nrow(fit$beta)
    if (K == 1) {
      dat <- as.data.frame(t(cbind(fit$a, as.matrix(fit$b))))
      rownames(dat) <- c(paste0('a', 1:(nrow(dat) - 1)), 'b')
      colnames(dat) <- paste0(1:ncol(dat))
      print(dat, max = max, digits = digits, ...)
    } else {
      dat <- do.call(rbind, lapply(2:K, function(k) {
        dif <- as.data.frame(t(ifelse(cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0, 'X', '')))
        rownames(dif) <- paste0(k, ':', c(paste0('a', 1:(nrow(dif) - 1)), 'b'))
        colnames(dif) <- paste0(1:ncol(dif))
        dif
      }))
      print(dat, max = max, ...)
    }
    if (K > 1 || length(x$all) != 1 || fit$lambda0 != 0)
      cat(sep = '', '* lambda0 = ', fit$lambda0, '\n')
  }
  invisible(fit)
}

#' Summarize DIF 2PL Items
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'}, \code{'GIC'}, or the constant for computing GIC, otherwise use the criterion specified when fitting the model(s)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @usage summary(x, criterion = NULL, max = 99999L, digits = 3, ...)
#'
#' @seealso \code{\link{D2PL_em}}, \code{\link{D2PL_gvem}}, \code{\link{D2PL_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
summary.vemirt_DIF <- function(x, criterion = NULL, max = 99999L, digits = 3, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  if (is.null(fit$tau)) {
    K <- nrow(fit$beta)
    if (K == 1) {
      dat <- as.data.frame(t(cbind(fit$a, as.matrix(fit$b))))
      rownames(dat) <- c(paste0('a', 1:(nrow(dat) - 1)), 'b')
      colnames(dat) <- paste0(1:ncol(dat))
      print(dat, max = max, digits = digits, ...)
    } else {
      dat <- (Reduce(`+`, lapply(2:nrow(fit$beta), function(k) {
        dif <- t((cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0) + 0)
        rownames(dif) <- c(paste0('a', 1:(nrow(dif) - 1)), 'b')
        colnames(dif) <- paste0(1:ncol(dif))
        dif
      })) != 0) + 0
      print(as.data.frame(ifelse(dat, 'X', '')), max = max, ...)
    }
    if (K > 1 || length(x$all) != 1 || fit$lambda0 != 0)
      cat(sep = '', '* lambda0 = ', fit$lambda0, '\n')
  }
  invisible(dat)
}

new.vemirt_FA <- function(raw) {
  structure(raw, class = 'vemirt_FA')
}

#' Extract Parameter Estimates from Explanatory or Confirmatory Analysis
#'
#' @param object An object of class \code{vemirt_FA}
#'
#' @usage coef(object)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{C2PL_gvem}}, \code{\link{C2PL_bs}}, \code{\link{C2PL_iw}}, \code{\link{C3PL_sgvem}}, \code{\link{E2PL_gvem_adaptlasso}}, \code{\link{E2PL_gvem_lasso}}, \code{\link{E2PL_gvem_rot}}, \code{\link{E2PL_IS}}, \code{\link{E3PL_sgvem_adaptlasso}}, \code{\link{E3PL_sgvem_lasso}}, \code{\link{E3PL_sgvem_rot}}, \code{\link{print.vemirt_FA}}
#' @export
coef.vemirt_FA <- function(object, ...) {
  coefs <- as.data.frame(if (!is.null(object$boots_a))
    cbind(object$boots_a, object$boots_b)
  else if (!is.null(object$new_a))
    cbind(object$new_a, object$new_b)
  else
    cbind(object$ra, object$rb)
  )
  rownames(coefs) <- paste0('Item', 1:nrow(coefs))
  colnames(coefs) <- c(paste0('a', 1:(ncol(coefs) - 1)), 'b')
  if (!is.null(object$rc))
    coefs <- cbind(coefs, c = object$rc)
  coefs
}

#' Print Parameter Estimates from Explanatory or Confirmatory Analysis
#'
#' @param x An object of class \code{vemirt_FA}
#'
#' @usage print(x)
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{C2PL_gvem}}, \code{\link{C2PL_bs}}, \code{\link{C2PL_iw}}, \code{\link{C3PL_sgvem}}, \code{\link{E2PL_gvem_adaptlasso}}, \code{\link{E2PL_gvem_lasso}}, \code{\link{E2PL_gvem_rot}}, \code{\link{E2PL_IS}}, \code{\link{E3PL_sgvem_adaptlasso}}, \code{\link{E3PL_sgvem_lasso}}, \code{\link{E3PL_sgvem_rot}}, \code{\link{coef.vemirt_FA}}
#' @export
print.vemirt_FA <- function(x, max = 99999L, digits = 3, ...) {
  coefs <- coef.vemirt_FA(x)
  print(coefs, max = max, digits = digits, ...)
  invisible(coefs)
}
