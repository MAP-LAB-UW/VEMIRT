#################### IRT ###################
library(mirt)

set.seed(123)

N <- 2000
J <- 12
K <- 4   # f1, f2, f3, general factor f4

generate_a_bi <- function(loadings) {
  # loadings length = 24
  # 1:12  = specific discriminations
  # 13:24 = general discriminations

  a <- matrix(0, nrow = J, ncol = K)
  colnames(a) <- c("f1", "f2", "f3", "general")
  rownames(a) <- paste0("x", 1:J)

  # specific factors
  a[1:4, 1]  <- loadings[1:4]    # f1
  a[5:8, 2]  <- loadings[5:8]    # f2
  a[9:12, 3] <- loadings[9:12]   # f3

  # general factor
  a[1:12, 4] <- loadings[13:24]

  return(a)
}
loadings1 <- runif(24, 1.5, 2.5)

a1 <- generate_a_bi(loadings1)

# item intercepts
# You can also use rnorm(J, 0, 1)
d1 <- rnorm(J, mean = 0, sd = 1)

# latent abilities: f1, f2, f3, general
Theta1 <- MASS::mvrnorm(
  n = N,
  mu = rep(0, K),
  Sigma = diag(K)
)

dat1 <- simdata(
  a = a1,
  d = d1,
  itemtype = rep("2PL", J),
  Theta = Theta1
)

loadings2 <- loadings1

# change general factor loadings for item 1 and item 3
loadings2[13] <- loadings1[13] + 0.5
loadings2[15] <- loadings1[15] - 0.8
a2 <- generate_a_bi(loadings2)

# keep same d if you only want discrimination DIF
# d2 <- d1
# change some d
d2 <- d1
d2[1] <- d1[1] + 0.5
d2[3] <- d1[3] - 0.8

Theta2 <- MASS::mvrnorm(
  n = N,
  mu = rep(0, K),
  Sigma = diag(K)
)

dat2 <- simdata(
  a = a2,
  d = d2,
  itemtype = rep("2PL", J),
  Theta = Theta2
)

dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dat1$group <- 1
dat2$group <- 2

dat <- rbind(dat1, dat2)

a <- array(dim = c(2, dim(a1)))
a[1, , ] <- a1[, c(4, 1:3)]
a[2, , ] <- a2[, c(4, 1:3)]
D2PL_bifactor_data <- list(data = dat[, 1:12], model = (a1[, c(4, 1:3)] != 0) + 0, ndim = 1, impact = FALSE, group = dat[, 13], j = c(1, 3),
                           params = list(a = a, b = -rbind(d1, d2), theta = rbind(Theta1, Theta2)))
usethis::use_data(D2PL_bifactor_data, overwrite = TRUE)
