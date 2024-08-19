# Written by Weicong Lyu

final.C2PL_iw2 <- function(e, SE.level) {
  list2env(e, environment())
  with_no_grad({
    niter <- niter.init[2]
    ll <- as.array(Q)
    l0 <- sum(D)
    Mu <- Mu[1]
    Sigma <- (Sigma.L %*% t(Sigma.L))[1]
  })
  if (is.null(SE.level)) {
    SE.a <- NULL
    SE.b <- NULL
  } else {
    if (!is.numeric(SE.level)) {
      warning('No accuracy level provided; use 10 by default.')
      SE.level <- 10
    }
    grid <- suppressMessages(createNIGrid(K, 'nHN', SE.level, 'sparse'))
    rescale(grid, m = as.array(Mu), C = as.array(Sigma), dec.type = 1)
    z <- getNodes(grid)
    w <- as.vector(getWeights(grid)) / as.array(Sigma$det()$sqrt())
    a <- torch_zeros(a.mask$shape)$masked_scatter(a.mask, a.v)
    xi <- (z %*% t(a) - b)$masked_fill(!Y.mask, NaN)
    logP <- torch_where(Y, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3)
    Pw <- with_no_grad(exp(logP) * w)
    l <- (Pw * logP)$sum(2) / Pw$sum(2)
    I <- Reduce(`+`, lapply(l$unbind(), function(l) {
      d <- torch_cat(autograd_grad(l, list(a.v, b), retain_graph = T))
      d$unsqueeze(2) %*% d$unsqueeze(1)
    }))
    SE <- I$pinverse()$diag()$sqrt()
    SE.a <- as.array(torch_zeros(a.mask$shape)$masked_scatter(a.mask, SE[1:-(J + 1)]))
    SE.b <- as.array(SE[-J:-1])
  }
  with_no_grad(c(lst(niter, ll, l0), lapply(lst(SIGMA, MU, Sigma, Mu, a, b), as.array), lst(SE.a, SE.b), IC(ll, l0, N, c)))
}

#' IW-GVEM Algorithm for Confirmatory M2PL Analysis
#'
#' @param data An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})
#' @param model A \eqn{J\times K} binary matrix of loading indicators  (all items load on the only dimension by default)
#' @param criterion Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'BIC'}, or \code{'AIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param S Sample size for approximating the expected lower bound
#' @param M Sample size for approximating a tighter lower bound
#' @param lr Learning rate for the Adam optimizer
#' @param SE.level Accuracy level of Gaussian quadrature for \code{mvQuad} to compute standard errors (SEs are not computed if \code{SE.level} is \code{NULL})
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:
#'   \item{N}{Number of respondents}
#'   \item{niter0}{Number(s) of iterations for initialization}
#'   \item{fit}{The only element of \code{all}}
#'   \item{best}{Equal to \code{1}}
#'   \item{all}{A list of model which has one element:}
#'   \item{ ...$niter}{Number(s) of iterations}
#'   \item{ ...$SIGMA}{Person-level posterior covariance matrices}
#'   \item{ ...$MU}{Person-level posterior mean vectors}
#'   \item{ ...$Sigma}{Population covariance matrix}
#'   \item{ ...$Mu}{Population mean vector}
#'   \item{ ...$a}{Slopes}
#'   \item{ ...$b}{Intercepts}
#'   \item{ ...$SE.a}{Standard errors of \code{a}}
#'   \item{ ...$SE.b}{Standard errors of \code{b}}
#'   \item{ ...$ll}{Estimated lower bound of log-likelihood}
#'   \item{ ...$l0}{Number of nonzero elements in \code{model}}
#'   \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
#'   \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{C2PL_gvem}}, \code{\link{C2PL_iw}}, \code{\link{D2PL_gvem}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(C2PL_data, C2PL_iw2(data, model, SE = TRUE))}
C2PL_iw2 <- function(data, model = matrix(1, ncol(data)), criterion = 'BIC', iter = 200, eps = 1e-3, c = 1, S = 10, M = 10, lr = 0.1, SE.level = NULL) {
  Y <- as.matrix(data)
  D <- as.matrix(model)
  N <- nrow(Y)

  init <- init.D2PL_iwgvemm(Y, D, rep(1, N), iter, eps, S, M, lr, c)
  results <- list(final.C2PL_iw2(init(), SE.level))
  new.vemirt_DIF(init('niter.init')[1], results, N, criterion)
}
