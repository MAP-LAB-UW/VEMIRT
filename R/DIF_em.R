# Written by Weicong Lyu

groupsum <- function(n, x) {
  torch_stack(lapply(n, function(n) {
    x[n]$sum(1)
  }))
}

estep.em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      if (length(level) == 1) {
        rescale(grid, m = as.array(Mu[g]), C = as.array(Sigma[g]), dec.type = 1)
        list(z = getNodes(grid), w = as.vector(getWeights(grid)))
      } else {
        w <- dmvn(grid, as.array(Mu[g]), as.array(Sigma[g]))
        list(z = torch_tensor(grid), w = torch_tensor(w))
      }
    })
    z <- torch_stack(lapply(gd, `[[`, 'z'))
    z.a <- z$view(c(G, -1, 1, 1, K)) * a.mask$unsqueeze(3)
    ZZ <- sym(z.a %*% t(z.a))[X]
    w <- torch_stack(lapply(gd, `[[`, 'w'))

    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    W <- torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3)$exp() * w[X]
    W <- W / W$sum(2)$unsqueeze(2)
    w <- groupsum(n, W)
    w <- (w / w$sum(2)$unsqueeze(2))$unsqueeze(3)
    W$unsqueeze_(3)
    Mu <- (w * z)$sum(2)
    Sigma <- sym((z$unsqueeze(4) %*% t(z$unsqueeze(4)) * w$unsqueeze(4))$sum(2) - Mu$unsqueeze(3) %*% t(Mu$unsqueeze(3)))

    mu <- Mu[1]$clone()
    z$sub_(mu)
    Mu$sub_(mu)
    b$sub_(a %*% mu)
    beta$add_(gamma %*% mu)
    sigma <- Sigma[1]$diag()$sqrt()
    a$mul_(sigma)
    gamma$mul_(sigma)
    sigma.inv <- (1 / sigma)$diag()
    Sigma <- sym(sigma.inv %*% Sigma %*% sigma.inv)
  })
}

mstep.em <- function() {
  with(parent.frame(), {
    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    p <- torch_sigmoid(xi)
    w.d <- (Y - p) * W * lr
    w.dd <- -(p * (1 - p)) * W

    D.gamma <- (w.d$unsqueeze(4) * z[X]$unsqueeze(3))$sum(2)
    DD.gamma <- (ZZ * w.dd$view(c(N, -1, J, 1, 1)))$sum(2)
    d.gamma <- groupsum(n, D.gamma)$masked_fill(!gamma.mask, 0)
    dd.gamma <- groupsum(n, DD.gamma)
    d.a <- D.gamma$sum(1)
    dd.a <- DD.gamma$sum(1)
    D.beta <- t(w.d$unsqueeze(3))$sum(2)$squeeze(3)
    DD.beta <- w.dd$sum(2)
    d.beta <- groupsum(n, D.beta)$masked_fill(!beta.mask, 0)
    dd.beta <- groupsum(n, DD.beta)
    d.b <- -D.beta$sum(1)
    dd.b <- DD.beta$sum(1)

    a$sub_((dd.a$pinverse() %*% d.a$unsqueeze(3))$squeeze(3))
    b$sub_(d.b / dd.b)
    if (get0('lambda', ifnotfound = 0) == 0) {
      gamma$sub_((dd.gamma$pinverse() %*% d.gamma$unsqueeze(4))$squeeze(4))
      beta$sub_(d.beta / dd.beta)
    } else {
      dd.gamma <- diagonal(dd.gamma)
      gamma <- (-prox(d.gamma - dd.gamma * gamma, lambda) / dd.gamma)$masked_fill(!gamma.mask, 0)
      beta <- -prox(d.beta - dd.beta * beta, lambda) / dd.beta
    }
  })
}

init.em <- function(Y, D, X, level, iter, eps, lr, ...) {
  e <- as.list(environment())
  e$lr <- 0.1
  init <- do.call(init.iwgvemm, c(e, S = 10, M = 10))
  e <- init()
  e$Y <- Y
  e$lr <- lr
  list2env(e, environment())
  with_no_grad(assemble.iwgvemm())
  b$requires_grad_(F)
  rm(Q, MU, SIGMA, Mu.v, Mu.mask, Sigma.L, Sigma1.L, Sigma1.L.v, Sigma1.L.mask, Sigma2.L, Sigma2.L.v, Sigma2.L.mask, a.v, gamma.v, gamma.v.mask, beta.v, beta.v.mask, z, theta, theta.logd, X.rank, S, M, init, e)
  Y <- torch_tensor(Y)$unsqueeze(2)
  U <- Y$bool()
  grid <- if (length(level) == 1)
    createNIGrid(K, 'nHN', level, 'sparse')
  else
    as.matrix(do.call(expand.grid, replicate(K, level, simplify = F)))

  params.old <- NULL
  for (i in 1:iter) {
    estep.em()
    mstep.em()
    params <- lapply(lst(Mu, Sigma), torch_clone)
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  niter.init <- c(niter.init, i)
  rm(sigma, sigma.inv, mu, p, z, z.a, ZZ, w, W, w.d, w.dd, d.a, d.b, d.gamma, d.beta, dd.a, dd.b, dd.gamma, dd.beta, D.gamma, D.beta, DD.gamma, DD.beta, i, params.old)
  init.new()
}

final.em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      if (length(level) == 1) {
        rescale(grid, m = as.array(Mu[g]), C = as.array(Sigma[g]), dec.type = 1)
        list(z = getNodes(grid), w = as.vector(getWeights(grid)))
      } else {
        z <- as.matrix(do.call(expand.grid, replicate(K, (level[-1] + level[-length(level)]) / 2, simplify = F)))
        d <- do.call(expand.grid, replicate(K, diff(level), simplify = F))
        w <- dmvn(z, as.array(Mu[g]), as.array(Sigma[g])) * apply(d, 1, prod)
        lapply(lst(z, w), torch_tensor)
      }
    })
    z <- torch_stack(lapply(gd, `[[`, 'z'))
    w <- torch_stack(lapply(gd, `[[`, 'w'))
    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    ll <- as.array((torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3)$exp() * w[X])$sum(2)$log()$sum())
    l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
    c(lst(niter, ll, l0), lapply(lst(Mu, Sigma, a, b, gamma, beta), as.array), IC(ll, l0, N, c))
  })
}

est.em <- function(e, lambda) {
  list2env(e, environment())
  niter <- if (lambda == 0)
    0
  else {
    params.old <- NULL
    for (i1 in 1:iter) {
      estep.em()
      mstep.em()
      params <- lapply(lst(Mu, Sigma, a, b, gamma, beta), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
    }
    lambda <- 0
    gamma.mask <- gamma != 0
    beta.mask <- beta != 0
    for (i2 in 1:iter) {
      estep.em()
      mstep.em()
      params <- lapply(lst(Mu, Sigma, a, b, gamma, beta), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
    }
    c(i1, i2)
  }
  final.em()
}

est.emm <- function(e, lambda) {
  list2env(e, environment())
  if (lambda == 0)
    niter <- 0
  else {
    lambda0 <- lambda
    gamma.mask0 <- gamma.mask
    beta.mask0 <- beta.mask
    params.old <- NULL
    for (niter in 1:iter) {
      estep.em()
      mstep.em()
      lambda <- 0
      gamma.mask <- gamma != 0
      beta.mask <- beta != 0
      mstep.em()
      params <- lapply(lst(Mu, Sigma, a, b, gamma, beta), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
      lambda <- lambda0
      gamma.mask <- gamma.mask0
      beta.mask <- beta.mask0
    }
  }
  final.em()
}

#' EM Algorithms for DIF Detection in 2PL Models
#'
#' @param data An \eqn{N\times J} binary matrix of item responses
#' @param model A \eqn{J\times K} binary matrix of loading indicators  (all items load on the only dimension by default)
#' @param group An \eqn{N} dimensional vector of group indicators from \code{1} to \code{G} (all respondents are in the same group by default)
#' @param method Estimation algorithm, one of \code{'EM'} or \code{'EMM'}
#' @param Lambda0 A vector of \code{lambda0} values for \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) * lambda0})
#' @param level Accuracy level, either a number for \code{mvQuad} or a vector indicating the grid for each latent dimension
#' @param criterion Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'AIC'}, or \code{'BIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param lr Learning rate for the Newton-Raphson method
#' @param c Constant for computing GIC
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:
#'   \item{N}{Number of respondents}
#'   \item{niter0}{Number of iterations for initialization (GVEM and EM with \code{lambda} equal to \code{0})}
#'   \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
#'   \item{best}{The index of \code{fit} in \code{all}}
#'   \item{all}{A list of models whose length is equal to \code{Lambda0}:}
#'   \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
#'   \item{ ...$lambda}{\code{sqrt(N) * lambda0}}
#'   \item{ ...$niter}{Number(s) of iterations}
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
#' @seealso \code{\link{DIF_gvem}}, \code{\link{DIF_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(exampleDIF, DIF_em(data, model, group))}
DIF_em <- function(data, model = matrix(1, ncol(data)), group = rep(1, nrow(data)), method = 'EMM', Lambda0 = seq(0.1, 0.8, by = 0.1), level = 18, criterion = 'GIC', iter = 200, eps = 1e-3, lr = 0.7, c = 1.7) {
  fn <- switch(method,
               EM = list(init.em, est.em),
               EMM = list(init.em, est.emm),
               stop(paste0("Method '", method, "' not supported.")))
  Y <- as.matrix(data)
  D <- as.matrix(model)
  X <- group
  N <- nrow(Y)
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  if (min(X) == 0)
    X <- X + 1

  cat('Fitting the model with lambda = 0 for initial values...\n')
  init <- fn[[1]](Y, D, X, level, iter, eps, lr, c)
  cat('Fitting the model with different lambdas...\n')
  pb <- txtProgressBar(0, length(Lambda0), style = 3)
  result <- lapply(Lambda0, function(lambda0) {
    lambda <- sqrt(N) * lambda0
    result <- fn[[2]](init(), lambda)
    setTxtProgressBar(pb, pb$getVal() + 1)
    c(lst(lambda0, lambda), result)
  })
  close(pb)
  new_vemirt_DIF(init('niter.init'), result, N, criterion)
}
