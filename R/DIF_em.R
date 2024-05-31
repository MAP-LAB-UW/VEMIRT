# Written by Weicong Lyu

estep.em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      if (inherits(grid, 'NIGrid')) {
        rescale(grid, m = as.array(Mu[g]), C = as.array(Sigma[g]), dec.type = 1)
        z <- getNodes(grid)
        w <- as.vector(getWeights(grid))
      } else {
        z <- grid$z
        w <- dmvn(z, as.array(Mu[g]), as.array(Sigma[g])) * grid$s
      }
      lapply(lst(z, w), torch_tensor)
    })
    z <- torch_stack(lapply(gd, `[[`, 'z'))
    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    W <- torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3)$exp() * torch_stack(lapply(gd, `[[`, 'w'))[X]
    W$div_(W$sum(2)$unsqueeze(2))$unsqueeze_(3)
    w <- groupsum(n, W)
    w$div_(w$sum(2)$unsqueeze(2))
    Mu <- (z * w)$sum(2)
    Sigma <- sym(((z$unsqueeze(4) %*% z$unsqueeze(3)) * w$unsqueeze(4))$sum(2) - Mu$unsqueeze(3) %*% Mu$unsqueeze(2))

    mu <- Mu[1]$clone()
    z$sub_(mu)
    Mu$sub_(mu)
    b$sub_(a %*% mu)
    beta$add_(gamma %*% mu)
    sigma <- Sigma[1]$diag()$sqrt()
    z$div_(sigma)
    Sigma <- sym(Sigma / sigma / sigma$unsqueeze(2))
    a$mul_(sigma)
    gamma$mul_(sigma)
    if (lambda == 0 && exists('init'))
      for (g in 2:G) {
        mu <- Mu[g]$clone()
        z[g]$sub_(mu)
        Mu[g]$sub_(mu)
        beta[g]$add_(gamma[g] %*% mu)
        sigma <- Sigma[g]$diag()$sqrt()
        z[g]$div_(sigma)
        Sigma[g] <- sym(Sigma[g] / sigma / sigma$unsqueeze(2))
        gamma[g] <- (a + gamma[g]) * sigma - a
      }
    Z <- z$unsqueeze(3)[X]
    ZZ <- (z$unsqueeze(4) %*% z$unsqueeze(3))$unsqueeze(3)[X]
  })
}

mstep.em <- function() {
  with(parent.frame(), {
    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    p <- torch_sigmoid(xi)
    q <- torch_sigmoid(-xi)
    W.d <- ((Y - p) * W)$unsqueeze(4) * lr
    W.dd <- (-(p * q) * W)$unsqueeze(4)
    d.a <- (Z * W.d)$sum(1:2)
    dd.a <- (ZZ * W.dd$unsqueeze(5))$sum(1:2)
    d.b <- -W.d$sum(1:2)
    dd.b <- W.dd$sum(1:2)$unsqueeze(3)
    dd.ab <- -(Z * W.dd)$sum(1:2)
    dd <- torch_cat(list(torch_cat(list(dd.a, dd.ab$unsqueeze(3)), 3), torch_cat(list(dd.ab$unsqueeze(2), dd.b), 3)), 2)$masked_fill(!ab.mask, 0)
    ab <- torch_hstack(list(a, b$unsqueeze(2))) - (dd$pinverse()$masked_fill(!ab.mask, 0) %*% torch_cat(list(d.a, d.b), 2)$unsqueeze(3))$squeeze(3)
    a <- ab[, 1:-2]$masked_fill(!a.mask, 0)
    b <- ab[, -1]

    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    p <- torch_sigmoid(xi)
    q <- torch_sigmoid(-xi)
    W.d <- ((Y - p) * W)$unsqueeze(4) * lr
    W.dd <- (-(p * q) * W)$unsqueeze(4)
    d.gamma <- groupsum(n, (Z * W.d)$sum(2))$masked_fill(!gamma.mask, 0)
    dd.gamma <- groupsum(n, (ZZ * W.dd$unsqueeze(5))$sum(2))
    d.beta <- groupsum(n, W.d$sum(2))$masked_fill(!beta.mask$unsqueeze(3), 0)
    dd.beta <- groupsum(n, W.dd$sum(2))$unsqueeze(4)
    if (lambda == 0) {
      mask <- torch_cat(list(gamma.mask, beta.mask$unsqueeze(3)), 3)$unsqueeze(4)$float()
      mask <- (mask %*% t(mask))$bool()
      dd.gammabeta <- groupsum(n, (Z * W.dd)$sum(2))
      dd <- torch_cat(list(torch_cat(list(dd.gamma, dd.gammabeta$unsqueeze(4)), 4), torch_cat(list(dd.gammabeta$unsqueeze(3), dd.beta), 4)), 3)$masked_fill(!mask, 0)
      gammabeta <- torch_cat(list(gamma, beta$unsqueeze(3)), 3) - (dd$pinverse()$masked_fill(!mask, 0) %*% torch_cat(list(d.gamma, d.beta), 3)$unsqueeze(4))$squeeze(4)
      gamma <- gammabeta[.., 1:-2]
      beta <- gammabeta[.., -1]
    } else {
      dd.gamma <- diagonal(dd.gamma)
      dd.beta <- dd.beta$view(c(G, -1))
      gamma <- (-prox(d.gamma - dd.gamma * gamma, lambda) / dd.gamma)
      beta <- (-prox(d.beta$squeeze(3) - dd.beta * beta, lambda) / dd.beta)
    }
    gamma$masked_fill_(!gamma.mask, 0)
    beta$masked_fill_(!beta.mask, 0)
  })
}

init.em <- function(Y, D, X, level, iter, eps, lr, ...) {
  init <- do.call(init.gvem, as.list(environment()))
  e <- init()
  e$Y <- torch_tensor(Y)$unsqueeze(2)
  list2env(e, environment())
  U <- Y$bool()
  ab.mask <- torch_hstack(list(a.mask, torch_ones(J, 1)))$unsqueeze(3)
  ab.mask <- (ab.mask %*% t(ab.mask))$bool()
  grid <- if (length(level) == 1)
    createNIGrid(K, 'nHN', level, 'sparse')
  else {
    z <- as.matrix(do.call(expand.grid, replicate(K, (level[-1] + level[-length(level)]) / 2, simplify = F)))
    s <- apply(do.call(expand.grid, replicate(K, diff(level), simplify = F)), 1, prod)
    lst(z, s)
  }

  lambda <- 0
  params.old <- NULL
  for (i in 1:iter) {
    estep.em()
    mstep.em()
    params <- lapply(lst(Sigma, Mu, a, b, gamma, beta), torch_clone)
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  niter.init <- c(niter.init, i)
  keep(Y, D, X, level, iter, eps, lr, grid, N, n, J, K, G, U, Sigma, Mu, a, a.mask, b, gamma, gamma.mask, beta, beta.mask, ab.mask, niter.init)
  init.new()
}

final.em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      if (inherits(grid, 'NIGrid')) {
        rescale(grid, m = as.array(Mu[g]), C = as.array(Sigma[g]), dec.type = 1)
        z <- getNodes(grid)
        w <- as.vector(getWeights(grid))
      } else {
        z <- grid$z
        w <- dmvn(z, as.array(Mu[g]), as.array(Sigma[g])) * grid$s
      }
      lapply(lst(z, w), torch_tensor)
    })
    xi <- (((a + gamma)$view(c(G, 1, J, 1, K)) %*% torch_stack(lapply(gd, `[[`, 'z'))$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X]
    ll <- as.array((torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3)$exp() * torch_stack(lapply(gd, `[[`, 'w'))[X])$sum(2)$log()$sum())
    l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
    c(lst(niter, ll, l0), lapply(lst(Sigma, Mu, a, b, gamma, beta), as.array), IC(ll, l0, N, c))
  })
}

est.em <- function(e, lambda) {
  list2env(e, environment())
  if (lambda == 0)
    niter <- 0
  else {
    params.old <- NULL
    for (i1 in 1:iter) {
      estep.em()
      mstep.em()
      params <- lapply(lst(Sigma, Mu, a, b, gamma, beta), torch_clone)
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
      params <- lapply(lst(Sigma, Mu, a, b, gamma, beta), torch_clone)
      if (!is.null(params.old) && all(distance(params, params.old) < eps))
        break
      params.old <- params
    }
    niter <- c(i1, i2)
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
      params <- lapply(lst(Sigma, Mu, a, b, gamma, beta), torch_clone)
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
#' @param criterion Information criterion for model selection, one of \code{'BIC'} (recommended), \code{'AIC'}, or \code{'GIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param lr Learning rate for the Newton-Raphson method
#' @param c Constant for computing GIC
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
#'   \item{ ...$Sigma}{Group-level posterior covariance matrices}
#'   \item{ ...$Mu}{Group-level posterior mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$gamma}{DIF parameters for the slopes}
#'   \item{ ...$beta}{DIF parameters for the intercepts}
#'   \item{ ...$ll}{Log-likelihood}
#'   \item{ ...$l0}{Number of nonzero DIF parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
#'   \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
#'
#' @author Weicong Lyu <wlyu4@uw.edu>
#' @seealso \code{\link{DIF_gvem}}, \code{\link{DIF_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(exampleDIF, DIF_em(data, model, group))}
DIF_em <- function(data, model = matrix(1, ncol(data)), group = rep(1, nrow(data)), method = 'EMM', Lambda0 = seq(0.1, 0.8, by = 0.1), level = 12, criterion = 'BIC', iter = 200, eps = 1e-3, lr = 1, c = 1) {
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

  cat('Fitting the model without regularization for initial values...\n')
  init <- fn[[1]](Y, D, X, level, iter, eps, lr, c)
  cat('Fitting the model with different lambdas...\n')
  pb <- txtProgressBar(0, length(Lambda0), style = 3)
  results <- lapply(Lambda0, function(lambda0) {
    lambda <- sqrt(N) * lambda0
    result <- fn[[2]](init(), lambda)
    setTxtProgressBar(pb, pb$getVal() + 1)
    c(lst(lambda0, lambda), result)
  })
  close(pb)
  new.vemirt_DIF(init('niter.init'), results, N, criterion)
}
