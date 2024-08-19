# Written by Weicong Lyu

init.D2PL_pair_em <- function(Y, X, level, iter, eps, ...) {
  init <- do.call(init.D2PL_gvem, c(as.list(environment()), list(D = matrix(1, ncol(Y)))))
  e <- init()
  e$Y <- array(Y, append(dim(Y), 1, 1))
  list2env(e, environment())
  Y.mask <- torch_tensor(!is.na(Y))
  Y <- torch_tensor(Y)
  U <- Y$bool()
  a <- (a + gamma)$squeeze(-1)
  b <- (b - beta)$squeeze(-1)
  Mu$squeeze_(-1)
  Sigma <- Sigma$view(-1)

  grid <- suppressMessages(createNIGrid(K, 'nHN', level, 'sparse'))
  mask <- torch_stack(replicate(J, torch_tensor(1:G)$view(c(-1, 1)) < torch_tensor(1:G)$view(c(1, -1))), dim = 3) * group.mask$unsqueeze(1) * group.mask$unsqueeze(2)
  a$mul_(group.mask)
  b$mul_(group.mask)
  d.a <- a$unsqueeze(2) - a$unsqueeze(1)
  u.a <- torch_zeros_like(d.a)
  d.b <- b$unsqueeze(2) - b$unsqueeze(1)
  u.b <- torch_zeros_like(d.b)
  keep(Y, X, level, iter, eps, grid, Y.mask, group.mask, mask, N, n, J, G, U, Sigma, Mu, a, b, d.a, u.a, d.b, u.b, niter.init)
  init.new()
}

estep.D2PL_pair_em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      sigma <- as.array(Sigma[g])
      rescale(grid, m = as.array(Mu[g]), C = sigma, dec.type = 1)
      lapply(list(z = getNodes(grid), w = getWeights(grid) / sqrt(sigma)), as.vector)
    })
    z <- torch_stack(lapply(gd, `[[`, 'z'))
    xi <- ((a$unsqueeze(2) * z$unsqueeze(3) - b$unsqueeze(2))[X])$masked_fill(!Y.mask, NaN)
    P <- torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3)$exp() * torch_stack(lapply(gd, `[[`, 'w'))[X]
    W <- P / P$sum(2)$unsqueeze(2)
    w <- groupsum(n, W)
    w$div_(w$sum(2)$unsqueeze(2))

    Mu <- (z * w)$sum(2)
    Sigma <- (z ^ 2 * w)$sum(2) - Mu ^ 2
    x <- if (lambda == 0) 1:G else 1
    Mu[x] <- 0
    Sigma[x] <- 1
  })
}

penalty.D2PL_pair_em <- function(d.u, a, mask) {
  sum((d.u - (a$unsqueeze(2) - a$unsqueeze(1)))$masked_fill(!mask, 0) ^ 2)
}

proximal.D2PL_pair_em <- function(d, u, a, tau, eta) {
  x <- a$unsqueeze(2) - a$unsqueeze(1)
  y <- x - u
  d$set_data(torch_where(abs(d) >= tau, y, prox(y, eta)))
  u$add_(d - x)
}

mstep.D2PL_pair_em <- function() {
  with(parent.frame(), {
    eta <- lambda / rho
    a$requires_grad_(T)
    b$requires_grad_(T)
    opt <- optim_lbfgs(lst(a, b))
    loss <- function() {
      xi <- ((a$unsqueeze(2) * z$unsqueeze(3) - b$unsqueeze(2))[X])$masked_fill(!Y.mask, NaN)
      Q <- (torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3) * W)$sum() - rho / 2 * (penalty.D2PL_pair_em(d.a + u.a, a, mask) + penalty.D2PL_pair_em(d.b + u.b, b, mask))
      opt$zero_grad()
      (-Q)$backward()
      -Q
    }
    pars.old <- NULL
    for (ii in 1:iter) {
      opt$step(loss)
      with_no_grad({
        proximal.D2PL_pair_em(d.a, u.a, a, tau, eta)
        proximal.D2PL_pair_em(d.b, u.b, b, tau, eta)
        pars <- lapply(lst(a, b, d.a, d.b), torch_clone)
      })
      if (!is.null(pars.old) && all(distance(pars, pars.old) < eps))
        break
      pars.old <- pars
    }
    a$requires_grad_(F)
    b$requires_grad_(F)
    ii
  })
}

count.D2PL_pair_em <- function(d, mask, G) {
  find <- function(x) {
    if (x != p[x])
      p[x] <<- find(p[x])
    else
      x
  }
  g <- which(mask)
  p <- (1:G) * mask
  for (n in g)
    for (m in g)
      if (m < n && d[m, n] == 0) {
        i <- find(m)
        j <- find(n)
        if (i != j)
          p[i] <- j
      }
  for (n in g)
    for (m in g)
      if (find(m) == find(n))
        d[m, n] <- 0
  lst(d, l0 = sum(p == 1:G) - 1)
}

em.D2PL_pair_em <- function() {
  with(parent.frame(), {
    j <- c()
    params.old <- NULL
    for (i in 1:iter) {
      estep.D2PL_pair_em()
      j <- c(j, mstep.D2PL_pair_em())
      params <- lapply(lst(Sigma, Mu, a, b, d.a, d.b), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
    }
    j
  })
}

init.lambda.D2PL_pair_em <- function(e, lambda, ...) {
  list2env(e, environment())
  tau <- Inf
  niter.init <- em.D2PL_pair_em()
  keep(Y, X, rho, level, iter, eps, c, grid, Y.mask, group.mask, mask, N, n, J, G, U, Sigma, Mu, a, b, d.a, u.a, d.b, u.b, P, niter.init)
  init.new()
}

est.D2PL_pair_em <- function(e, pb, lambda, Tau, ...) {
  init <- init.lambda.D2PL_pair_em(e, lambda)
  if (!is.null(pb)) setTxtProgressBar(pb, pb$getVal() + 1)
  lapply(Tau, function(tau) {
    list2env(init(), environment())
    niter <- if (tau == Inf)
      0
    else
      em.D2PL_pair_em()
    d.a$mul_(mask)
    d.b$mul_(mask)
    ll <- as.array(P$sum(2)$log()$sum())
    l0 <- sum(sapply(1:J, function(j) {
      result.a <- count.D2PL_pair_em(as.array(d.a[.., j]), as.array(group.mask[.., j]), G)
      d.a[.., j] <- result.a$d
      result.b <- count.D2PL_pair_em(as.array(d.b[.., j]), as.array(group.mask[.., j]), G)
      d.b[.., j] <- result.b$d
      result.a$l0 + result.b$l0
    }))
    if (!is.null(pb)) setTxtProgressBar(pb, pb$getVal() + 1)
    c(lst(..., rho, lambda, tau, niter0 = init('niter.init'), niter, ll, l0), lapply(lst(Sigma, Mu, a, b, d.a, d.b, u.a, u.b), as.array), IC(ll, l0, N, c))
  })
}

#' EM Algorithm with ADMM for DIF Detection Using Group Pairwise Truncated \eqn{L_1} Penalty in 2PL Models
#'
#' @param data An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})
#' @param group An \eqn{N} dimensional vector of group indicators from \code{1} to \code{G} (all respondents are in the same group by default)
#' @param Lambda0 A vector of \code{lambda0} values for truncated \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) / G * lambda0})
#' @param Tau A vector of \code{tau} values for truncated \eqn{L_1} penalty (becomes \eqn{L_1} penalty when \code{tau} equals \code{Inf})
#' @param rho0 A value of \code{rho} for augmented Lagrangian in ADMM (\code{tau} equals \code{sqrt(N) / G * tau0})
#' @param level Accuracy level of Gaussian quadrature for \code{mvQuad}
#' @param criterion Information criterion for model selection, one of \code{'BIC'} (recommended), \code{'AIC'}, or \code{'GIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param verbose Whether to show the progress
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:
#'   \item{N}{Number of respondents}
#'   \item{niter0}{Number(s) of iterations for initialization}
#'   \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
#'   \item{best}{The index of \code{fit} in \code{all}}
#'   \item{all}{A list of models which has the same length as \code{Lambda0}:}
#'   \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
#'   \item{ ...$lambda}{\code{sqrt(N) / G * lambda0}}
#'   \item{ ...$tau}{Corresponding element in \code{Tau}}
#'   \item{ ...$rho0}{Same as \code{rho0} in input}
#'   \item{ ...$rho}{\code{sqrt(N) / G * rho0}}
#'   \item{ ...$niter}{Number(s) of iterations}
#'   \item{ ...$Sigma}{Group-level covariance matrices}
#'   \item{ ...$Mu}{Group-level mean vectors}
#'   \item{ ...$a}{Slopes}
#'   \item{ ...$b}{Intercepts}
#'   \item{ ...$d.a}{Group pairwise differences of slopes}
#'   \item{ ...$d.b}{Group pairwise differences of intercepts}
#'   \item{ ...$u.a}{Lagrangian multipliers of corresponding elements in \code{d.a}}
#'   \item{ ...$u.b}{Lagrangian multipliers of corresponding elements in \code{d.b}}
#'   \item{ ...$ll}{Log-likelihood}
#'   \item{ ...$l0}{Number of nonzero D2PL parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
#'   \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{D2PL_em}}, \code{\link{D2PL_gvem}}, \code{\link{D2PL_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(D2PL_data, D2PL_pair_em(data, group, Tau = c(Inf, seq(0.01, 0.05, by = 0.01))))}
D2PL_pair_em <- function(data, group = rep(1, nrow(data)), Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.5, 1.5, by = 0.1), Tau = if (length(unique(group)) == 1) 0 else c(Inf, seq(0.05, 0.3, by = 0.05)), rho0 = 0.5, level = 10, criterion = 'BIC', iter = 200, eps = 0.001, c = 1, verbose = TRUE) {
  Y <- as.matrix(data)
  N <- nrow(Y)
  X <- group
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  if (min(X) == 0)
    X <- X + 1
  m <- sqrt(N) / max(X)
  output <- if (verbose) stderr() else nullfile()

  cat(file = output, 'Fitting the model without regularization for initial values...\n')
  rho <- m * rho0
  init <- init.D2PL_pair_em(Y, X, level, iter, eps, c, rho)
  cat(file = output, 'Fitting the model with different lambdas and taus...\n')
  pb <- if (verbose)
    txtProgressBar(file = output, 0, length(Lambda0) * (length(Tau) + 1), style = 3)
  else
    NULL
  results <- unlist(lapply(Lambda0, function(lambda0) {
    lambda <- m * lambda0
    lapply(est.D2PL_pair_em(init(), pb, lambda, Tau, rho0, lambda0), c)
  }), F)
  if (verbose) close(pb)
  new.vemirt_DIF(init('niter.init'), results, N, criterion)
}
