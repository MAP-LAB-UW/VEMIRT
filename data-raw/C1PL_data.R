library(abind)
library(mvtnorm)
set.seed(246)

Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
J <- 20
D <- cbind(rep(1:0, J / 2), rep(0:1, J / 2))
N <- 1000

a <- 1.0 * D
b <- rnorm(J)
theta <- rmvnorm(N, rep(0, 2), Sigma)
Y <- t(matrix(rbinom(N * J, 1, plogis(a %*% t(theta) - b)), nrow = J))

C1PL_data <- list(data = Y, model = D, params = list(a = a, b = b, theta = theta))
usethis::use_data(C1PL_data, overwrite = TRUE)
