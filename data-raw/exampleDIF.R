library(abind)
library(mvtnorm)
set.seed(987)

Sigma <- matrix(c(1, 0.85, 0.85, 1), 2)
J <- 20
D <- cbind(rep(1:0, J / 2), rep(0:1, J / 2))
j <- J * 0.6
n <- 1000
X <- rep(1:3, each = n)

a <- matrix(runif(J * 2, 1.5, 2.5), ncol = 2) * D
a <- unname(abind(a, a, a, along = 0))
a[-1, 1:(j / 2), ] <- a[-1, 1:(j / 2), ] + c(0.4, 0.8)
a[-1, (j / 2 + 1):j, ] <- a[-1, (j / 2 + 1):j, ] - c(0.4, 0.8)
a[-1, , ] <- a[-1, , ] * abind(D, D, along = 0)
b <- rnorm(J)
b <- unname(rbind(b, b, b))
b[-1, 1:(j / 2)] <- b[-1, 1:(j / 2)] - c(0.5, 1)
b[-1, (j / 2 + 1):j] <- b[-1, (j / 2 + 1):j] + c(0.5, 1)
theta <- rmvnorm(n * 3, rep(0, 2), Sigma)
Y <- t(sapply(1:(n * 3), function(n) {
  rbinom(J, 1, plogis(a[X[n], , ] %*% theta[n, ] - b[X[n], ]))
}))
exampleDIF <- list(data = Y, model = D, group = X, j = j, params = list(a = a, b = b, theta = theta))

usethis::use_data(exampleDIF, overwrite = TRUE)
