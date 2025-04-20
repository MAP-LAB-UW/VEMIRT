# Written by Weicong Lyu

init.D2PL_gvem <- function(Y, D, X, iter, eps, ...) {
  N <- nrow(Y)
  n <- tapply(1:N, X, identity)
  J <- ncol(Y)
  K <- ncol(D)
  G <- max(X)
  Y.na <- is.na(Y)
  Y[Y.na] <- 0.5
  Y <- torch_tensor(Y)
  Y.mask <- torch_tensor(!Y.na)
  p <- groupmean(n, Y$masked_fill(!Y.mask, NaN)$float()$view(c(N, J)))
  nn <- ceiling(N / G)
  z0 <- torch_tensor(array(rnorm(G * nn * K), c(G, nn, K, 1)))

  eta <- torch_where(Y.mask, 0.125, 0)
  Sigma <- torch_stack(replicate(G, diag(K), simplify = F))
  Mu <- torch_zeros(G, K)
  a.mask <- torch_tensor(D != 0)
  a <- a.mask * 1
  b <- torch_zeros(J)
  group.mask <- groupsum(n, Y.mask$float()) > 0
  gamma.mask <- torch_stack(replicate(G, a.mask)) * group.mask$unsqueeze(3)
  gamma.mask[1] <- 0
  gamma <- gamma.mask * 0
  beta.mask <- group.mask$clone()
  beta.mask[1] <- 0
  beta <- beta.mask * 0

  lambda <- 0
  niter.init <- mstep.D2PL_gvem()
  keep(Y, D, X, iter, eps, N, n, J, K, G, MU, SIGMA, Mu, Sigma, a, b, gamma, beta, Y.mask, a.mask, gamma.mask, beta.mask, eta, niter.init, p, z0)
  init.new()
}

mstep.D2PL_gvem <- function() {
  with(parent.frame(), {
    params.old <- NULL
    for (i in 1:iter) {
      aG <- (a + gamma)$unsqueeze(4)
      aG.t <- t(aG)
      AG <- aG[X]
      AG.t <- t(AG)
      BB <- (b - beta)[X]

      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * (aG %*% aG.t)[X])$sum(2)
      SIGMA <- sym(SIGMA.inv$inverse())
      MU <- (SIGMA %*% (((Y - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + (Sigma.inv %*% Mu$unsqueeze(3))[X]))$squeeze(3)
      Mu <- groupmean(n, MU)
      mu <- MU - Mu[X]
      sigma.mu <- SIGMA + mu$unsqueeze(3) %*% mu$unsqueeze(2)
      Sigma <- sym(groupmean(n, sigma.mu))

      mu <- (a.mask * MU$unsqueeze(2))$unsqueeze(4)
      sigma.mu <- SIGMA$unsqueeze(2) * a.mask$unsqueeze(2) * a.mask$unsqueeze(3) + mu %*% t(mu)
      xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% mu)$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
      eta <- torch_where(abs(xi) < 1e-3, 0.125, (torch_sigmoid(xi) - 0.5) / (2 * xi))$masked_fill(!Y.mask, 0)
      a <- ((2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$view(c(N, -1, 1, 1)) * mu + 2 * eta$view(c(N, -1, 1, 1)) * (BB$view(c(N, -1, 1, 1)) * mu - sigma.mu %*% gamma[X]$unsqueeze(4)))$sum(1))$squeeze(3)$masked_fill(!a.mask, 0)
      b <- (0.5 - Y + 2 * eta * (beta[X] + (AG.t %*% mu)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))
      gamma.beta  <- torch_stack(lapply(n, function(n) {
        Y <- Y[n]
        eta <- eta[n]
        mu <- mu[n]
        sigma.mu <- sigma.mu[n]
        if (lambda == 0)
          torch_cat(list(((2 * eta$view(c(-1, J, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$unsqueeze(3) * mu$squeeze(4) + 2 * eta$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu$squeeze(4) - (sigma.mu %*% a$unsqueeze(3))$squeeze(4)))$sum(1)$unsqueeze(3))$squeeze(3),
                         (((Y - 0.5) + 2 * eta * (b - (AG.t[n] %*% mu)$view(c(-1, J))))$sum(1) / (2 * eta$sum(1)))$unsqueeze(2)), 2)
        else
          torch_cat(list(prox(((Y - 0.5)$unsqueeze(3) * mu$squeeze(4) + 2 * eta$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu$squeeze(4) - (sigma.mu %*% a$unsqueeze(3))$squeeze(4)))$sum(1), lambda) / diagonal((2 * eta$view(c(-1, J, 1, 1)) * sigma.mu)$sum(1)),
                         (prox(((Y - 0.5) + 2 * eta * (b - (AG.t[n] %*% mu)$view(c(-1, J))))$sum(1), lambda) / (2 * eta$sum(1)))$unsqueeze(2)), 2)
      }))
      gamma <- gamma.beta[.., 1:-2]$masked_fill(!gamma.mask, 0)
      beta <- gamma.beta[.., -1]$masked_fill(!beta.mask, 0)

      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(a %*% mu)
      beta$add_(gamma %*% mu)
      sigma <- Sigma[1]$diag()$sqrt()
      SIGMA <- sym(SIGMA / sigma / sigma$unsqueeze(2))
      Sigma <- sym(Sigma / sigma / sigma$unsqueeze(2))
      a$mul_(sigma)
      gamma$mul_(sigma)

      params <- lapply(lst(SIGMA, MU, Sigma, Mu, a, b, gamma, beta), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
    }
    i
  })
}

est.D2PL_gvem <- function(e, lambda) {
  list2env(e, environment())
  niter <- if (lambda == 0)
    c(0, 0)
  else {
    j <- mstep.D2PL_gvem()
    lambda <- 0
    gamma.mask <- gamma != 0
    beta.mask <- beta != 0
    c(j, mstep.D2PL_gvem())
  }
  AG <- (a + gamma)$unsqueeze(4)[X]
  AG.t <- t(AG)
  BB <- (b - beta)[X]
  mu <- (a.mask * MU$unsqueeze(2))$unsqueeze(4)
  sigma.mu <- SIGMA$unsqueeze(2) * a.mask$unsqueeze(2) * a.mask$unsqueeze(3) + mu %*% t(mu)
  xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% mu)$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
  mu <- (MU - Mu[X])$unsqueeze(3)
  ll <- as.array((nnf_logsigmoid(xi)$masked_fill(!Y.mask, 0) + (0.5 - Y) * (BB - (AG.t %*% MU$view(c(N, 1, -1, 1)))$view(c(N, -1))) - xi / 2)$masked_fill(!Y.mask, 0)$sum() - (Sigma$logdet()[X]$sum() + diagonal(linalg_solve(Sigma[X], SIGMA + mu %*% t(mu)))$sum()) / 2)
  final.D2PL_gvem()
}

init.D2PL_iwgvemm <- function(Y, D, X, iter, eps, S, M, lr, ...) {
  init <- do.call(init.D2PL_gvem, as.list(environment()))
  e <- init()
  e$Y <- array(Y, append(dim(Y), 1, 1))
  list2env(e, environment())

  Y.mask <- torch_tensor(!is.na(Y))
  Y <- torch_tensor(Y)$bool()
  X.rank <- rank(X, ties.method = 'first')
  Sigma.L <- linalg_cholesky(Sigma)
  Sigma1.L.mask <- torch_tensor(lower.tri(matrix(0, K, K)))$bool()
  Sigma1.L.v <- (Sigma.L[1] / torch_cat(list(torch_ones(K, 1), sqrt(1 - Sigma.L[1]$square()$cumsum(2))[, 1:-2]), 2))$masked_select(Sigma1.L.mask)$arctanh()$requires_grad_(T)
  if (G == 1) {
    Sigma2.L.mask <- torch_empty(c(0, K, K), dtype = torch_bool())
    Sigma2.L.v <- torch_empty(0, requires_grad = T)
  } else {
    Sigma2.L.mask <- torch_stack(replicate(G - 1, lower.tri(matrix(0, K, K), T), F))$bool()
    Sigma2.L.v <- Sigma.L[2:G]$masked_select(Sigma2.L.mask)$requires_grad_(T)
  }
  Mu.mask <- torch_cat(list(torch_zeros(1, K), torch_ones(G - 1, K)))$bool()
  Mu.v <- torch_tensor(Mu)$masked_select(Mu.mask)$requires_grad_(T)
  a.mask <- torch_tensor(D != 0)
  a.v <- torch_tensor(a)$masked_select(a.mask)$requires_grad_(T)
  b <- torch_tensor(b)$requires_grad_(T)
  gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
  gamma.v <- torch_tensor(gamma)$masked_select(gamma.mask)$requires_grad_(T)
  gamma.v.mask <- torch_ones_like(gamma.v)$bool()
  beta.mask <- torch_cat(list(torch_zeros(1, J), torch_ones(G - 1, J)))$bool()
  beta.v <- torch_tensor(beta)$masked_select(beta.mask)$requires_grad_(T)
  beta.v.mask <- torch_ones_like(beta.v)$bool()

  z <- array(rnorm(N * S * M * K), c(N, S * M, K))
  theta <- (linalg_cholesky(SIGMA)$unsqueeze(2) %*% torch_tensor(z)$unsqueeze(4))$squeeze(4) + torch_tensor(MU)$unsqueeze(2)
  theta.logd <- torch_tensor(rowSums(dnorm(z, log = T), dims = 2))
  lambda <- 0
  niter.init <- c(init('niter.init'), mstep.D2PL_iwgvemm())
  keep(Y, D, X, iter, eps, S, M, lr, Y.mask, X.rank, N, n, J, K, G, SIGMA, MU, Sigma.L, Sigma1.L.mask, Sigma1.L.v, Sigma2.L.mask, Sigma2.L.v, Mu, Mu.mask, Mu.v, a, a.mask, a.v, b, gamma, gamma.mask, gamma.v, gamma.v.mask, beta, beta.mask, beta.v, beta.v.mask, theta, theta.logd, Q, niter.init, p, z0)
  init.new()
}

assemble.D2PL_iwgvemm <- function() {
  with(parent.frame(), {
    z <- Sigma1.L.v$tanh()
    Sigma1.L <- torch_eye(K)$masked_scatter(Sigma1.L.mask, z) * torch_cat(c(torch_ones(K, 1), (1 - torch_zeros(K, K)$masked_scatter(Sigma1.L.mask, z)[, 1:-2]$square())$cumprod(2)$sqrt()), 2)
    Sigma2.L <- torch_zeros(Sigma2.L.mask$shape)$masked_scatter(Sigma2.L.mask, Sigma2.L.v)
    Sigma.L <- torch_cat(c(Sigma1.L$unsqueeze(1), Sigma2.L))
    Mu <- torch_zeros(Mu.mask$shape)$masked_scatter(Mu.mask, Mu.v)
    a <- torch_zeros(a.mask$shape)$masked_scatter(a.mask, a.v)
    gamma <- torch_zeros(gamma.mask$shape)$masked_scatter(gamma.mask, gamma.v)
    beta <- torch_zeros(beta.mask$shape)$masked_scatter(beta.mask, beta.v)
  })
}

mstep.D2PL_iwgvemm <- function() {
  with(parent.frame(), {
    params <- lst(Sigma1.L.v, Sigma2.L.v, Mu.v, a.v, b, gamma.v, beta.v)
    opt <- optim_adam(params, lr)
    params.old <- NULL
    for (i in 1:iter) {
      assemble.D2PL_iwgvemm()
      xi <- (((a + gamma)[X]$unsqueeze(2) * theta$unsqueeze(3))$sum(4) - (b - beta)[X]$unsqueeze(2))$masked_fill(!Y.mask, NaN)
      log.w <- torch_where(Y, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3) +
        torch_cat(lapply(1:G, function(g) {
          distr_multivariate_normal(Mu[g], scale_tril = Sigma.L[g])$log_prob(theta[n[[g]]])
        }))[X.rank] - theta.logd
      Q <- (log.w$view(c(N, S, M))$logsumexp(3) - log(M))$mean(2)$sum()
      with_no_grad({
        if (!is.null(params.old) && all(distance(params, params.old) < eps))
          return(i - 1)
        params.old <- lapply(params, torch_clone)
      })
      opt$zero_grad()
      (-Q)$backward()
      opt$step()
      with_no_grad(
        if (lambda == 0) {
          gamma.v$masked_fill_(!gamma.v.mask, 0)
          beta.v$masked_fill_(!beta.v.mask, 0)
        } else {
          proximal.adam(gamma.v, opt$state$get(gamma.v), opt$param_groups[[1]], lambda)
          proximal.adam(beta.v, opt$state$get(beta.v), opt$param_groups[[1]], lambda)
        }
      )
    }
    i
  })
}

est.D2PL_iwgvemm <- function(e, lambda) {
  list2env(e, environment())
  niter <- if (lambda == 0)
    c(0, 0)
  else {
    j <- mstep.D2PL_iwgvemm()
    lambda <- 0
    gamma.v.mask <- gamma.v != 0
    beta.v.mask <- beta.v != 0
    c(j, mstep.D2PL_iwgvemm())
  }
  with_no_grad({
    ll <- as.array(Q)
    Sigma <- Sigma.L %*% t(Sigma.L)
    final.D2PL_gvem()
  })
}

final.D2PL_gvem <- function() {
  with(parent.frame(), {
    theta <- (linalg_cholesky(Sigma)$unsqueeze(2) %*% z0)$squeeze(4) + torch_tensor(Mu)$unsqueeze(2)
    xi <- ((a + gamma)$view(c(G, 1, J, 1, K)) %*% theta$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2)
    q <- (xi$sigmoid())$mean(2)
    RMSE <- sqrt(mean((p - q) ^ 2))
    l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
    c(lst(niter, ll, l0), lapply(lst(SIGMA, MU, Sigma, Mu, a, b, gamma, beta, RMSE), as.array), IC(ll, l0, N, c))
  })
}

#' GVEM Algorithms for DIF Detection in M2PL Models
#'
#' @param data An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})
#' @param model A \eqn{J\times K} binary matrix of loading indicators  (all items load on the only dimension by default)
#' @param group An \eqn{N} dimensional vector of group indicators from \code{1} to \code{G} (all respondents are in the same group by default)
#' @param method Estimation algorithm, one of \code{'GVEM'} or \code{'IWGVEMM'}
#' @param Lambda0 A vector of \code{lambda0} values for \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) * lambda0})
#' @param criterion Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'BIC'}, or \code{'AIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param S Sample size for approximating the expected lower bound (\code{'IWGVEMM'} only)
#' @param M Sample size for approximating a tighter lower bound (\code{'IWGVEMM'} only)
#' @param lr Learning rate for the Adam optimizer (\code{'IWGVEMM'} only)
#' @param verbose Whether to show the progress
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:
#'   \item{N}{Number of respondents}
#'   \item{niter0}{Number(s) of iterations for initialization}
#'   \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
#'   \item{best}{The index of \code{fit} in \code{all}}
#'   \item{all}{A list of models which has the same length as \code{Lambda0}:}
#'   \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
#'   \item{ ...$lambda}{\code{sqrt(N) * lambda0}}
#'   \item{ ...$niter}{Number(s) of iterations}
#'   \item{ ...$SIGMA}{Person-level posterior covariance matrices}
#'   \item{ ...$MU}{Person-level posterior mean vectors}
#'   \item{ ...$Sigma}{Group-level covariance matrices}
#'   \item{ ...$Mu}{Group-level mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$gamma}{D2PL parameters for the slopes}
#'   \item{ ...$beta}{D2PL parameters for the intercepts}
#'   \item{ ...$RMSE}{Root mean square error of fitted probability of each item for each group}
#'   \item{ ...$ll}{Estimated lower bound of log-likelihood}
#'   \item{ ...$l0}{Number of nonzero D2PL parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
#'   \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
#'
#' @author Weicong Lyu <weiconglyu@um.edu.mo>
#' @seealso \code{\link{D2PL_pair_em}}, \code{\link{D2PL_em}}, \code{\link{D2PL_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(D2PL_data, D2PL_gvem(data, model, group))}
D2PL_gvem <- function(data, model = matrix(1, ncol(data)), group = rep(1, nrow(data)), method = 'IWGVEMM', Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.1, 1, by = 0.1), criterion = 'GIC', iter = 200, eps = 1e-3, c = 1, S = 10, M = 10, lr = 0.1, verbose = TRUE) {
  fn <- switch(method,
               GVEM = list(init.D2PL_gvem, est.D2PL_gvem),
               IWGVEMM = list(init.D2PL_iwgvemm, est.D2PL_iwgvemm),
               stop(paste0("Method '", method, "' not supported.")))
  Y <- as.matrix(data)
  D <- as.matrix(model)
  N <- nrow(Y)
  X <- group
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  if (min(X) == 0)
    X <- X + 1
  output <- if (verbose) stderr() else nullfile()

  cat(file = output, 'Fitting the model without regularization for initial values...\n')
  init <- fn[[1]](Y, D, X, iter, eps, S, M, lr, c)
  cat(file = output, 'Fitting the model with different lambdas...\n')
  if (verbose) pb <- txtProgressBar(file = output, 0, length(Lambda0), style = 3)
  results <- lapply(Lambda0, function(lambda0) {
    lambda <- sqrt(N) * lambda0
    result <- fn[[2]](init(), lambda)
    if (verbose) setTxtProgressBar(pb, pb$getVal() + 1)
    c(lst(lambda0, lambda), result)
  })
  if (verbose) close(pb)
  new.vemirt_DIF(init('niter.init'), results, N, criterion)
}
