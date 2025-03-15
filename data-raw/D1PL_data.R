library(abind)
library(mvtnorm)
set.seed(86420)

Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
J <- 20
D <- cbind(rep(1:0, J / 2), rep(0:1, J / 2))
j <- J * 0.3
n <- 500
X <- rep(1:3, each = n)

a <- matrix(1, nrow = J, ncol = 2) * D
a <- unname(abind(a, a, a, along = 0))
b <- rnorm(J)
b <- unname(rbind(b, b, b))
b[-1, 1:(j / 2)] <- b[-1, 1:(j / 2)] - c(0.5, 1)
b[-1, (j / 2 + 1):j] <- b[-1, (j / 2 + 1):j] + c(0.5, 1)
theta <- rmvnorm(n * 3, rep(0, 2), Sigma)
Y <- t(sapply(1:(n * 3), function(n) {
  rbinom(J, 1, plogis(a[X[n], , ] %*% theta[n, ] - b[X[n], ]))
}))

D1PL_data <- list(data = Y, model = D, group = X, j = j, params = list(a = a, b = b, theta = theta))
usethis::use_data(D1PL_data, overwrite = TRUE)
