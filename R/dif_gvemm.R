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

distance <- function(x, y) {
  mapply(function(x, y) {
    max(abs(x - y))
  }, x, y)
}

IC <- function(ll, l0, N, c) {
  -2 * ll + l0 * c(AIC = 2, BIC = log(N), GIC = c * log(N) * log(log(N)))
}

init.gvemm <- function(Y, D, X, ...) {
  env(!!!lst(...)) |> with({
    N <- nrow(Y)
    J <- ncol(Y)
    K <- ncol(D)
    G <- max(X)

    Y.na <- is.na(Y)
    Y[Y.na] <- 0.5
    Y <- torch_tensor(Y)
    Y.mask <- torch_tensor(!Y.na)
    eta <- torch_where(Y.mask, 0.125, 0)
    Sigma <- torch_stack(replicate(G, diag(K), simplify = F))
    Mu <- torch_zeros(G, K)
    a.mask <- torch_tensor(D != 0)
    a <- a.mask * 1
    a.mask.diag <- a$diag_embed()
    b <- torch_zeros(J)
    gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
    gamma <- gamma.mask * 0
    beta.mask <- torch_cat(list(torch_zeros(1, J), torch_ones(G - 1, J)))$bool()
    beta <- beta.mask * 0

    parameters <- function() {
      lapply(lst(SIGMA, MU, Sigma, Mu, a, b, gamma, beta), as.array)
    }

    update <- function() {
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
          Mu <- torch_stack(tapply(1:N, X, function(n) {
            MU[n]$mean(1)
          }))
          mu <- MU - Mu[X]
          sigma.mu <- SIGMA + mu$unsqueeze(3) %*% mu$unsqueeze(2)
          Sigma <- sym(torch_stack(tapply(1:N, X, function(n) {
            sigma.mu[n]$mean(1)
          })))

          mu <- a.mask.diag %*% MU$view(c(N, 1, -1, 1))
          sigma.mu <- a.mask.diag %*% SIGMA$unsqueeze(2) %*% a.mask.diag + mu %*% t(mu)
          xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% mu)$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
          eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))$masked_fill(!Y.mask, 0)
          a <- ((2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$view(c(N, -1, 1, 1)) * mu + 2 * eta$view(c(N, -1, 1, 1)) * (BB$view(c(N, -1, 1, 1)) * mu - sigma.mu %*% gamma[X]$unsqueeze(4)))$sum(1))$squeeze(3)$masked_fill(!a.mask, 0)
          b <- (0.5 - Y + 2 * eta * (beta[X] + (AG.t %*% mu)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))

          gamma.beta  <- torch_stack(tapply(1:N, X, function(n) {
            N <- length(n)
            torch_cat(list(prox(((Y[n] - 0.5)$unsqueeze(3) * mu[n]$squeeze(4) + 2 * eta[n]$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu[n]$squeeze(4) - (sigma.mu[n] %*% a$unsqueeze(3))$squeeze(4)))$sum(1), lambda) / diagonal((2 * eta[n]$view(c(N, -1, 1, 1)) * sigma.mu[n])$sum(1)),
                           (prox(((Y[n] - 0.5) + 2 * eta[n] * (b - (AG.t[n] %*% mu[n])$view(c(N, -1))))$sum(1), lambda) / (2 * eta[n]$sum(1)))$unsqueeze(2)), 2)
          }))
          gamma$set_data(gamma.beta[, , 1:K]$masked_fill(!gamma.mask, 0))
          beta$set_data(gamma.beta[, , (K + 1)]$masked_fill(!beta.mask, 0))

          mu <- Mu[1]$clone()
          MU$sub_(mu)
          Mu$sub_(mu)
          b$sub_(a %*% mu)
          beta$add_(gamma %*% mu)
          sigma <- Sigma[1]$diag()$sqrt()
          a$mul_(sigma)
          gamma$mul_(sigma)
          sigma.inv <- (1 / sigma)$diag()
          SIGMA$set_data(sym(sigma.inv %*% SIGMA %*% sigma.inv))
          Sigma$set_data(sym(sigma.inv %*% Sigma %*% sigma.inv))

          params <- parameters()
          if (!is.null(params.old) && all(distance(params, params.old) < eps))
            break
          params.old <- params
        }
        i
      })
    }

    current_env()
  })
}

gvemm <- function(e, lambda) {
  e$lambda <- lambda
  with(e, {
    if (lambda == 0) {
      niter <- update()
      pars <- parameters()
    } else {
      pars.old <- NULL
      niter <- matrix(nrow = 0, ncol = 2)
      lambda0 <- lambda
      gamma.mask0 <- gamma.mask
      beta.mask0 <- beta.mask
      for (i0 in 1:iter) {
        i1 <- update()
        lambda <- 0
        gamma.mask <- gamma != 0
        beta.mask <- beta != 0
        i2 <- update()
        niter <- rbind(niter, c(i1, i2))
        pars <- parameters()
        if (!is.null(pars.old) && all(distance(pars, pars.old) < eps))
          break
        pars.old <- pars
        lambda <- lambda0
        gamma.mask <- gamma.mask0
        beta.mask <- beta.mask0
      }
    }
    mu <- (MU - Mu[X])$unsqueeze(3)
    ll <- as.array((nnf_logsigmoid(xi) + (0.5 - Y) * (BB - (AG.t %*% MU$view(c(N, 1, -1, 1)))$view(c(N, -1))) - xi / 2)$masked_fill(!Y.mask, 0)$sum() - (Sigma$logdet()[X]$sum() + diagonal(linalg_solve(Sigma[X], SIGMA + mu %*% t(mu)))$sum()) / 2)
    l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
    c(lst(niter, ll, l0), pars, IC(ll, l0, N, c))
  })
}

init.iwgvemm <- function(Y, D, X, iter, eps, c, S, M, ...) {
  e <- do.call(init.gvemm, c(as.list(environment()), lst(...)))
  e$niter.init <- gvemm(e, 0)$niter
  e$Y <- array(Y, append(dim(Y), 1, 1))
  with(e, {
    Y.mask <- torch_tensor(!is.na(Y))
    Y <- torch_tensor(Y)$bool()
    X.rank <- rank(X, ties.method = 'first')
    Sigma.L <- linalg_cholesky(Sigma)
    Sigma1.L.mask <- torch_tensor(lower.tri(matrix(0, K, K)))$bool()
    Sigma1.L.v <- (Sigma.L[1] / torch_cat(list(torch_ones(K, 1), sqrt(1 - Sigma.L[1]$square()$cumsum(2))[, 1:(K - 1)]), 2))$masked_select(Sigma1.L.mask)$arctanh()$requires_grad_(T)
    if (G > 1) {
      Sigma2.L.mask <- torch_stack(replicate(G - 1, lower.tri(matrix(0, K, K), T), F))$bool()
      Sigma2.L.v <- Sigma.L[2:G]$masked_select(Sigma2.L.mask)$requires_grad_(T)
    } else {
      Sigma2.L.mask <- torch_tensor(NULL)
      Sigma2.L.v <- torch_tensor(NULL)
    }
    Mu.mask <- torch_cat(list(torch_zeros(1, K), torch_ones(G - 1, K)))$bool()
    Mu.v <- torch_tensor(Mu)$masked_select(Mu.mask)$requires_grad_(T)
    a.mask <- torch_tensor(D != 0)
    a.v <- torch_tensor(a)$masked_select(a.mask)$requires_grad_(T)
    b <- torch_tensor(b)$requires_grad_(T)
    gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
    gamma.v <- torch_tensor(gamma)$masked_select(gamma.mask)$requires_grad_(T)
    beta.mask <- torch_cat(list(torch_zeros(1, J), torch_ones(G - 1, J)))$bool()
    beta.v <- torch_tensor(beta)$masked_select(beta.mask)$requires_grad_(T)

    z <- array(rnorm(N * S * M * K), c(N, S * M, K))
    theta <- (linalg_cholesky(SIGMA)$unsqueeze(2) %*% torch_tensor(z)$unsqueeze(4))$squeeze(4) + torch_tensor(MU)$unsqueeze(2)
    theta.logd <- torch_tensor(rowSums(dnorm(z, log = T), dims = 2))

    assemble <- function() {
      with(parent.frame(), {
        z <- Sigma1.L.v$tanh()
        Sigma1.L <- torch_eye(K)$masked_scatter(Sigma1.L.mask, z) * torch_cat(c(torch_ones(K, 1), (1 - torch_zeros(K, K)$masked_scatter(Sigma1.L.mask, z)[, 1:(K- 1)]$square())$cumprod(2)$sqrt()), 2)
        Sigma2.L <- if (G > 1)
          torch_zeros(Sigma2.L.mask$shape)$masked_scatter(Sigma2.L.mask, Sigma2.L.v)
        else
          torch_tensor(NULL)
        Sigma.L <- torch_cat(c(Sigma1.L$unsqueeze(1), Sigma2.L))
        Mu <- torch_zeros(Mu.mask$shape)$masked_scatter(Mu.mask, Mu.v)
        a <- torch_zeros(a.mask$shape)$masked_scatter(a.mask, a.v)
        gamma <- torch_zeros(gamma.mask$shape)$masked_scatter(gamma.mask, gamma.v)
        beta <- torch_zeros(beta.mask$shape)$masked_scatter(beta.mask, beta.v)
      })
    }

    parameters <- function() {
      with_no_grad({
        assemble()
        lapply(c(Sigma = Sigma.L %*% t(Sigma.L), lst(Sigma.L, Mu, a, b, gamma, beta)), as.array)
      })
    }

    proximal <- function(x, state, params) {
      beta.t <- params$betas ^ state$step
      lr <- params$lr / (1 - beta.t[1]) / (sqrt(state$exp_avg_sq / (1 - beta.t[2])) + params$eps)
      x$set_data(prox(x, lr * lambda))
    }

    update <- function() {
      with(parent.frame(), {
        opt <- optim_adam(lst(Sigma1.L.v, Sigma2.L.v, Mu.v, a.v, b, gamma.v, beta.v), lr)
        params.old <- parameters()
        for (i in 1:iter) {
          assemble()
          xi <- (((a + gamma)[X]$unsqueeze(2) * theta$unsqueeze(3))$sum(4) - (b - beta)[X]$unsqueeze(2))$masked_fill(!Y.mask, NaN)
          log.w <- torch_where(Y, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3) +
            torch_cat(lapply(1:G, function(g) {
              distr_multivariate_normal(Mu[g], scale_tril = Sigma.L[g])$log_prob(theta[X == g])
            }))[X.rank] - theta.logd
          Q <- (log.w$view(c(N, S, M))$logsumexp(3) - log(M))$mean(2)$sum()
          opt$zero_grad()
          (-Q)$backward()
          opt$step()
          with_no_grad({
            if (lambda == 0) {
              gamma.v$masked_fill_(!gamma.v.mask, 0)
              beta.v$masked_fill_(!beta.v.mask, 0)
            } else {
              proximal(gamma.v, opt$state$get(gamma.v), opt$param_groups[[1]])
              proximal(beta.v, opt$state$get(beta.v), opt$param_groups[[1]])
            }
          })
          params <- parameters()
          if (all(distance(params, params.old)[-1] < eps))
            break
          params.old <- params
        }
        i
      })
    }

    current_env()
  })
}

iwgvemm <- function(e, lambda) {
  e$lambda <- lambda
  with(e, {
    niter <- if (lambda == 0) {
      gamma.v.mask <- torch_ones_like(gamma.v)$bool()
      beta.v.mask <- torch_ones_like(beta.v)$bool()
      update()
    } else {
      i1 <- update()
      lambda <- 0
      gamma.v.mask <- gamma.v != 0
      beta.v.mask <- beta.v != 0
      c(i1, update())
    }
    with_no_grad({
      ll <- as.array(Q)
      l0 <- as.array(sum(gamma.v.mask) + sum(beta.v.mask))
    })
    c(lst(niter, ll, l0), params, IC(ll, l0, N, c))
  })
}

#' GVEMM Algorithms for DIF Detection in 2PL Models
#'
#' @param Y An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})
#' @param D A \eqn{J\times K} binary matrix of loading indicators
#' @param X An \eqn{N} dimensional vector of group indicators (integers from \code{1} to \code{G})
#' @param method Estimation algorithm, one of \code{'GVEMM'} or \code{'IWGVEMM'}
#' @param Lambda0 A vector of \code{lambda0} values (duplicate values removed automatically) for \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) * lambda0})
#' @param criterion Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'BIC'}, or \code{'AIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param S Sample size for approximating the expected lower bound (\code{'IWGVEMM'} only)
#' @param M Sample size for approximating a tighter lower bound (\code{'IWGVEMM'} only)
#' @param lr Learning rate for the Adam optimizer (\code{'IWGVEMM'} only)
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing three elements:
#'   \item{N}{Number of respondents}
#'   \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
#'   \item{best}{The index of \code{fit} in \code{all}}
#'   \item{all}{A list of models whose length is equal to \code{Lambda0}:}
#'   \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
#'   \item{ ...$lambda}{\code{sqrt(N) * lambda0}}
#'   \item{ ...$niter}{Number(s) of iterations}
#'   \item{ ...$SIGMA}{Person-level posterior covariance matrices}
#'   \item{ ...$MU}{Person-level posterior mean vectors}
#'   \item{ ...$Sigma}{Group-level posterior covariance matrices}
#'   \item{ ...$Mu}{Group-level posterior mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$gamma}{DIF parameters for the slopes}
#'   \item{ ...$beta}{DIF parameters for the intercepts}
#'   \item{ ...$ll}{Log-likelihood}
#'   \item{ ...$l0}{Number of nonzero parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion}
#'   \item{ ...$BIC}{Bayesian Information Criterion}
#'   \item{ ...$GIC}{Generalized Information Criterion}
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(exampleDIF, gvemm_DIF(Y, D, X))}
gvemm_DIF <- function(Y, D, X, method = 'IWGVEMM', Lambda0 = seq(0.2, 0.7, by = 0.1), criterion = 'GIC', iter = 1000, eps = 1e-3, c = 0.75, S = 10, M = 10, lr = 0.1) {
  Lambda0 <- unique(sort(Lambda0))
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  if (min(X) == 0)
    X <- X + 1
  N <- nrow(Y)

  est <- if (method == 'GVEMM') {
    cat('Fitting the model using different lambdas...\n')
    e <- init.gvemm(Y, D, X, iter, eps, c)
    gvemm
  } else if (method == 'IWGVEMM') {
    cat('Running GVEMM for initial values...\n')
    e <- init.iwgvemm(Y, D, X, iter, eps, c, S, M, lr)
    cat('Fitting the model using different lambdas...\n')
    iwgvemm
  } else
    stop(paste0("Method '", method, "' not supported."))
  pb <- txtProgressBar(0, length(Lambda0), style = 3)
  result <- lapply(Lambda0, function(lambda0) {
    lambda <- sqrt(N) * lambda0
    result <- est(e, lambda)
    setTxtProgressBar(pb, pb$getVal() + 1)
    c(lst(lambda0, lambda), result)
  })
  close(pb)
  result[[1]]$niter <- c(e$niter.init, result[[1]]$niter)
  new_vemirt_DIF(result, N, criterion)
}
