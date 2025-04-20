library(MASS)

generate_para_grm <- function(N, J, D, K_list, true_Mod, sL = 0.2, sU = 0.7,aL = 0.5 , aU = 1, bL = 0.3, bU = 1){
  a <- matrix(runif(J*D,aL,aU), J, D)
  true_A <- a * true_Mod
  true_A <- round(true_A,3)

  mat <- matrix(0, J, (max(K_list)))
  for(j in 1:J)
  {
    mat[j,] <-  c(rep(1, (K_list[j])), rep(0, (max(K_list)- K_list[j])))
  }
  if(max(K_list) != 1){
    # a vector including the number of graded score for all items
    true_gam <- matrix(data=runif(J*(max(K_list)),bL,bU), nrow=J, ncol=max(K_list))
    true_B   <- t(scale(apply(true_gam, 1, cumsum), center=TRUE, scale=FALSE))
    attr(true_B, "scaled:center") <- NULL
    true_B <- round(true_B,3)
    true_B <- true_B* mat
  }else
    true_B <- rnorm(J,0,1)

  s <- runif(D*D,sL, sU)
  true_S <- matrix(s,D,D)
  true_S <- .5*true_S + .5*t(true_S)
  diag(true_S) <- 1
  true_S <- round(true_S,3)
  return(list(a_true = true_A, b_true = true_B, Sig_true = true_S))

}
generate_response_grm <- function(N,A,B,Sig,K)
{
  output <- generate_data(N,A,B,Sig,R=K)
  x <- output$x
  y <- output$y
  return(list(y = y, x = x))
}
generate_data <- function(N,A,B,Sig,R=NULL){

  R <- R + 1
  # Laixu3107, 2024.04.16
  # Generate graded responses given true parameters.
  # N:   A positive integer, sample size.
  # A:   J*K mat, slope parameters for each item.
  # B:   J*max(R) mat, threshold parameters for each item.
  # Sig: K*K mat, correlation matrix of latent traits.
  # R:   J*1 vec, number of grade for each item.

  J <- nrow(A)
  K <- ncol(A)

  #if(is.null(R))  R <- rowSums(!is.na(B)) + 1

  x <- mvrnorm(n=N, mu=rep(0,K), Sigma = Sig)
  y <- matrix(NA, nrow=N, ncol=J)
  f <- matrix(0, nrow=N, ncol=max(R)+1)

  min_N_prob <- .02  # the minimum proportion of subjects who get the same graded score.

  for(j in 1:J){

    N_prob <- 0
    while(N_prob < min_N_prob){
      Rj <- R[j]
      Aj <- A[j,,drop=FALSE]
      if(max(R) == 2) Bj <- B[j]
      else Bj <- B[j,1:(Rj-1),drop=FALSE]

      logitf <- matrix(x%*%t(Aj),nrow=N,ncol=Rj-1,byrow=FALSE) - matrix(Bj,nrow=N,ncol=Rj-1,byrow=TRUE)
      f[,]     <- NA
      f[,1]    <- 1
      f[,Rj+1] <- 0
      f[,2:Rj] <- plogis(logitf)
      p <- f[,1:Rj] - f[,(1:Rj)+1]

      yj__ <- apply(X=p, MARGIN=1, FUN=rmultinom, n=1, size=1)
      yj   <- apply(X=yj__, MARGIN=2, FUN=which.max) - 1

      N_prob <- min(table(yj))/N
    }
    y[,j] <- yj
  }

  storage.mode(y) <- "integer"

  return(list(x = x, y = y))

}

N <- 1000; #sample size
J <- 30; # number of item
K_list <- c(rep(3,15), rep(5,15)); # graded responses e.g. 2: 0,1,2; 1: 0,1
D <- 4; # latent dimension
true_Mod <- matrix(0,J,D)
true_Mod[ 1:10,1] <- 1
true_Mod[ 8:20,2] <- 1
true_Mod[ 15:25,3] <- 1
true_Mod[ 20:30,4] <- 1
#demo
set.seed(243)
res3 <- generate_para_grm(N, J, D, K_list, true_Mod,  sL = 0.1, sU = 0.3,aL = 0.5 , aU = 2, bL = 0.3, bU = 1) # generate parameters
# a: Uniformly distributed between aL (lower bound) and aU (upper bound).
# b: Matrix elements are initialized with uniform values in [bL, bU], then transformed via row-wise cumulative summation.
# σ (sigma): Off-diagonal elements drawn from a uniform distribution with bounds sL (lower) and sU (upper).
data3 <- generate_response_grm(N, res3$a_true, res3$b_true, res3$Sig_true, K = K_list)$y # generate responses
MGRM_data <- list(data = data3, model = true_Mod, params = list(a = res3$a_true, b = res3$b_true, Sigma = res3$Sig_true))
usethis::use_data(MGRM_data, overwrite = TRUE)
