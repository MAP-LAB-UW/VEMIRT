library(abind)
library(mvtnorm)
set.seed(12345)

Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
J <- 20
D <- cbind(rep(1:0, J / 2), rep(0:1, J / 2))
N <- 2000

a <- matrix(runif(J * 2, 1.5, 2.5), ncol = 2) * D
b <- rnorm(J)
c <- rbeta(J, 5, 20)
theta <- rmvnorm(N, rep(0, 2), Sigma)
Y <- t(matrix(rbinom(N * J, 1, c + (1 - c) * plogis(a %*% t(theta) - b)), nrow = J))

C3PL_data <- list(data = Y, model = D, params = list(a = a, b = b, c = c, theta = theta))
usethis::use_data(C3PL_data, overwrite = TRUE)

#fit <- VEMIRT::C3PL_sgvem(Y, D, mu_b = 0, sigma2_b = 1, Alpha = 5, Beta = 20)
