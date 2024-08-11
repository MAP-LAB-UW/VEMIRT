library(abind)
library(mvtnorm)
set.seed(13579)

Sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), 3)
J <- 30
D <- cbind(rep(c(1, 0, 0), J / 3), rep(c(0, 1, 0), J / 3), rep(c(0, 0, 1), J / 3))
D[1:3, 1:3] <- c(1, 1, 1, 0, 1, 1, 0, 0, 1)
N <- 1000

a <- matrix(runif(J * 3, 1.5, 2.5), ncol = 3) * D
b <- rnorm(J)
theta <- rmvnorm(N, rep(0, 3), Sigma)
Y <- t(matrix(rbinom(N * J, 1, plogis(a %*% t(theta) - b)), nrow = J))
D[-(1:3), ] <- 1

E2PL_data_C2 <- list(data = Y, model = D, constrain = 'C2', non_pen = 3, params = list(a = a, b = b, theta = theta))
usethis::use_data(E2PL_data_C2, overwrite = TRUE)


fit <- VEMIRT::E2PL_gvem_lasso(Y, D, constrain = 'C2', non_pen = 3, max.iter = 10000)
