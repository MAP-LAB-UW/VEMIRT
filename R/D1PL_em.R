# Written by Weicong Lyu

init.D1PL_em <- function(Y, D, X, a, level, iter, eps, ...) {
  init <- do.call(init.D1PL_gvem, as.list(environment()))
  e <- init()
  e$Y <- array(Y, append(dim(Y), 1, 1))
  list2env(e, environment())
  Y.mask <- torch_tensor(!is.na(Y))
  Y <- torch_tensor(Y)
  U <- Y$bool()
  grid <- if (length(level) == 1)
    suppressMessages(createNIGrid(K, 'nHN', level, 'sparse'))
  else {
    z <- as.matrix(do.call(expand.grid, replicate(K, (level[-1] + level[-length(level)]) / 2, simplify = F)))
    s <- apply(do.call(expand.grid, replicate(K, diff(level), simplify = F)), 1, prod)
    lst(z, s)
  }

  keep(Y, D, X, level, iter, eps, grid, N, n, J, K, G, U, Sigma, Mu, a, b, gamma, Y.mask, beta, beta.mask, niter.init)
  init.new()
}

estep.D1PL_em <- function() {
  with(parent.frame(), {
    gd <- lapply(1:G, function(g) {
      if (inherits(grid, 'NIGrid')) {
        rescale(grid, m = as.array(Mu[g]), C = as.array(Sigma[g]), dec.type = 1)
        z <- getNodes(grid)
        w <- as.vector(getWeights(grid)) / as.array(Sigma[g]$det()$sqrt())
      } else {
        z <- grid$z
        w <- dmvn(z, as.array(Mu[g]), as.array(Sigma[g])) * grid$s
      }
      lapply(lst(z, w), torch_tensor)
    })
    z <- torch_stack(lapply(gd, `[[`, 'z'))
    Z <- z$unsqueeze(3)[X]
    ZZ <- (z$unsqueeze(4) %*% z$unsqueeze(3))$unsqueeze(3)[X]
    xi <- (((a$unsqueeze(2) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X])$masked_fill(!Y.mask, NaN)
    P <- torch_where(U, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$nansum(3)$exp() * torch_stack(lapply(gd, `[[`, 'w'))[X]
    W <- (P / P$sum(2)$unsqueeze(2))$unsqueeze(3)
    w <- groupsum(n, W)
    w$div_(w$sum(2)$unsqueeze(2))

    Mu <- (z * w)$sum(2)
    Sigma <- sym(((z$unsqueeze(4) %*% z$unsqueeze(3)) * w$unsqueeze(4))$sum(2) - Mu$unsqueeze(3) %*% Mu$unsqueeze(2))
    x <- if (lambda.bak == 0) 1:G else 1
    Mu[x] <- 0
    sigma <- Sigma[x]$diagonal(dim1 = -1, dim2 = -2)$sqrt()$view(c(-1, K, 1))
    Sigma[x] <- sym(Sigma[x] / sigma / t(sigma))
  })
}

mstep.D1PL_em <- function() {
  with(parent.frame(), {
    pars.old <- NULL
    for (ii in 1:iter) {
      p <- torch_sigmoid(xi)
      q <- torch_sigmoid(-xi)
      W.d <- ((Y - p) * W)
      W.dd <- (-(p * q) * W)
      d.b <- -W.d$nansum(1:2)
      dd.b <- W.dd$nansum(1:2)
      b <- b - d.b / dd.b

      xi <- (((a$unsqueeze(2) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X])$masked_fill(!Y.mask, NaN)
      p <- torch_sigmoid(xi)
      q <- torch_sigmoid(-xi)
      W.d <- (Y - p) * W
      W.dd <- -(p * q) * W
      d.beta <- groupsum(n, W.d$sum(2))$masked_fill(!beta.mask, 0)
      dd.beta <- groupsum(n, W.dd$sum(2))
      if (lambda == 0)
        beta <- beta - d.beta / dd.beta
      else
        beta <- -prox(d.beta - dd.beta * beta, lambda) / dd.beta
      beta$masked_fill_(!beta.mask, 0)

      xi <- (((a$unsqueeze(2) %*% z$view(c(G, -1, 1, K, 1)))$view(c(G, -1, J)) - (b - beta)$unsqueeze(2))[X])$masked_fill(!Y.mask, NaN)
      pars <- lapply(lst(b, beta), torch_clone)
      if (!is.null(pars.old) && all(distance(pars, pars.old) < eps))
        break
      pars.old <- pars
    }
    ii
  })
}

est.D1PL_em <- function(e, lambda) {
  list2env(e, environment())
  niter <- c()
  lambda.bak <- lambda
  params.old <- NULL
  for (i in 1:iter) {
    estep.D1PL_em()
    niter <- c(niter, mstep.D1PL_em())
    params <- lapply(lst(Sigma, Mu, b, beta), torch_clone)
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  lambda <- 0
  beta.mask <- beta != 0
  for (i in 1:iter) {
    params.old <- params
    estep.D1PL_em()
    niter <- c(niter, mstep.D1PL_em())
    params <- lapply(lst(Sigma, Mu, b, beta), torch_clone)
    if (all(distance(params, params.old) < eps))
      break
  }
  final.D1PL_em()
}

est.D1PL_emm <- function(e, lambda) {
  list2env(e, environment())
  niter <- matrix(nrow = 0, ncol = 2)
  lambda.bak <- lambda
  beta.mask.bak <- beta.mask
  params.old <- NULL
  for (i in 1:iter) {
    estep.D1PL_em()
    j <- mstep.D1PL_em()
    lambda <- 0
    beta.mask <- beta != 0
    niter <- rbind(niter, c(j, mstep.D1PL_em()))
    params <- lapply(lst(Sigma, Mu, b, beta), torch_clone)
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
    lambda <- lambda.bak
    beta.mask <- beta.mask.bak
  }
  final.D1PL_em()
}

final.D1PL_em <- function() {
  with(parent.frame(), {
    ll <- as.array(P$nansum(2)$log()$sum())
    l0 <- as.array(sum(beta != 0))
    c(lst(niter, ll, l0), lapply(lst(Sigma, Mu, a, b, gamma, beta), as.array), IC(ll, l0, N, c))
  })
}

#' EM Algorithms for DIF Detection in M1PL Models
#'
#' @param data An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})
#' @param model A \eqn{J\times K} binary matrix of loading indicators  (all items load on the only dimension by default)
#' @param group An \eqn{N} dimensional vector of group indicators from \code{1} to \code{G} (all respondents are in the same group by default)
#' @param a A scalar indicating the common discrimination parameter for all the dimensions of all the items (takes \code{1} by default)
#' @param method Estimation algorithm, one of \code{'EM'} or \code{'EMM'}
#' @param Lambda0 A vector of \code{lambda0} values for \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) * lambda0})
#' @param level Accuracy level, either a number for \code{mvQuad} or a vector indicating the grid for each latent dimension
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
#'   \item{ ...$lambda}{\code{sqrt(N) * lambda0}}
#'   \item{ ...$niter}{Number(s) of iterations}
#'   \item{ ...$Sigma}{Group-level covariance matrices}
#'   \item{ ...$Mu}{Group-level mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$gamma}{D1PL parameters for the slopes (all elements are zero)}
#'   \item{ ...$beta}{D1PL parameters for the intercepts}
#'   \item{ ...$ll}{Log-likelihood}
#'   \item{ ...$l0}{Number of nonzero D1PL parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
#'   \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
#'
#' @author Weicong Lyu <weiconglyu@um.edu.mo>
#' @seealso \code{\link{D1PL_gvem}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(D1PL_data, D1PL_em(data, model, group))}
D1PL_em <- function(data, model = matrix(1, ncol(data)), group = rep(1, nrow(data)), a = 1, method = 'EMM', Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.1, 1, by = 0.1), level = 10, criterion = 'BIC', iter = 200, eps = 1e-3, c = 1, verbose = TRUE) {
  fn <- switch(method,
               EM = list(init.D1PL_em, est.D1PL_em),
               EMM = list(init.D1PL_em, est.D1PL_emm),
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
  output <- if (verbose) stderr() else nullfile()

  cat(file = output, 'Fitting the model without regularization for initial values...\n')
  init <- fn[[1]](Y, D, X, a, level, iter, eps, c)
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
