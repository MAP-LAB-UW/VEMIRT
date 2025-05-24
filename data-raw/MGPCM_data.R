###########################################################################
# date : Jan 2024
# Auxiliary functions of pGVEM

###########################################################################
generate_para_pcm<-function(N, J, D, model, K_list, sL = 0, sU = 0, aL = 0.5 , aU = 1)
{
  # K_list is a vector indicating K_j-1: K_j is the number of different partial credit scores
  if(length(K_list) != J)
    return(0)
  K=max(K_list);

  #g enerate a_j
  aj <- array(0, dim = c(J, D))
  for(j in 1:J)
    for(d in 1:D)
      aj[j,d] <- runif(1, aL, aU)

  # generate b_j
  # parameter b is stored in a J x K dimension array with K+1 being the largest number of categories of one item
  bet <- array(0, dim = c(J, K))
  for(j in 1:J)
    for(k in 1:K_list[j])
      bet[j,k] <- rnorm(1, 0, 1)

  # generate theta
  thet <- array(0, dim = c(N, D))
  sigt <- diag(D);
  if(D>1)
  {
    for(i in 1:D)
      for(j in 1:D)
        if(i<j){sigt[i,j] <- runif(1,sL,sU); sigt[j,i] <- sigt[i,j];}
  }

  for(i in 1:N)
    thet[i,] <- mvrnorm(1, rep(0, D), sigt);
  aj <- aj * model
  return(list(thet_true = thet, a_true = aj, b_true = bet, Sig_true = sigt))
}

generate_response_pcm <- function(thet, aj, bet, K_list){
  fskip <- 0;
  I <- dim(thet)[1];
  J <- dim(aj)[1];
  D <- dim(thet)[2];
  K <- max(K_list)
  #generate a set of responses such that there does not exist a category for an item which no one scores
  while(!fskip)
  {
    response <- array(0, dim = c(I, J))
    Scorepro <- array(0, dim = c(I, J, (K + 1)))
    for(i in 1:I){
      for(j in 1:J){
        exparr=c();
        for(k in 1:K_list[j])
          exparr[k]<-exp(k*aj[j,]%*%thet[i,]-bet[j,k]);
        s <- sum(exparr)+1;
        Scorepro[i,j,1] <- 1.0 / s;
        for(k in 1:K_list[j])
          Scorepro[i,j,(k+1)] <- exparr[k] / s;
      }
    }
    for(i in 1:I)
      for(j in 1:J)
      {
        x <- runif(1,0,1)
        for(k in 1:(K_list[j]+1))
          if(x>sum(Scorepro[i,j,1:k])) {response[i,j] <- k;}
      }

    fskip <- ifResp(response,K_list);
  }

  return(response);
}

ifResp<-function(Resp, K_list)
{
  fskip <- 1

  for (j in 1:J) {
    for (k in 0:K_list[j]) {
      if (sum(Resp[, j] == k) == 0) {
        fskip <- 0
        break
      }
    }
    if (!fskip) {
      break
    }
  }

  return(fskip);
}



# eta_function
quaeta <- function(x){if(x==0) return(1/8);
  return ((1-exp(-x))/(4*x*(1+exp(-x))));}

#compute evidence lower bound
elbo <- function(I, J, D, K_list, ha, hbet, hSigm, hmu, hxi, response){

  if(length(K_list) != J)
    return(0)

  K <- max(K_list);

  qhxi <- array(1 / 8,dim=c(I,J,(K + 1)))
  for(i in 1:I)
    for(j in 1:J)
      for(k in 1:(K_list[j] + 1))
        if(k != (response[i,j] + 1))
          qhxi[i,j,k] <- quaeta(hxi[i,j,k]);

  llk <- 0;
  hbet <- cbind(0,hbet);
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      llk <- llk + log(2);
      q <- ha[j,] %*% hSigm[i,,] %*% ha[j,];
      p <- ha[j,] %*% hmu[i,];
      for(k in 1:(K_list[j]+1))
      {
        x <- qhxi[i,j,k]
        act <- response[i,j] + 1;
        dw <- (k - act) * p;db <- hbet[j,k] - hbet[j,act];
        vec <- dw - db;
        xs <- vec^2 + (k - act)^2 * q;
        llk <- llk - x * xs - 0.5 * vec + x * hxi[i,j,k]^2 + 0.5 * hxi[i,j,k] - log(1 + exp(hxi[i,j,k]));
      }
    }
    llk <- llk - 0.5 * sum(diag(hSigm[i, , ])) -0.5 * hmu[i,] %*% hmu[i,];
  }
  return(llk);
}
average_mean <- function(estimation, true)
{
  len = length(estimation)
  bias = 0
  y <- estimation
  y.na <- is.na(y)
  y[y.na] <- 0
  estimation <- y
  for(i in c(1:len))
  {
    bias = bias + (estimation[i] - true[i])
  }
  bias/len
}
rmse <- function(estimation, true)
{
  len = length(estimation)
  bias = 0
  y <- estimation
  y.na <- is.na(y)
  y[y.na] <- 0
  estimation <- y
  for(i in c(1:len))
  {
    bias = bias + (estimation[i] - true[i])^2
  }
  sqrt(bias/len)
}
permute <- function(est, real)
{
  D <- ncol(est)
  permn=gtools::permutations(D,D)
  error_last <- 10^6
  permun <- 0
  for(k in 1:(D*(D-1)))
  {
    error <- sum((est[,permn[k,]] - real[1:J, ])^2)
    if(error_last > error)
    {
      permun <- k
      error_last <- error
    }
  }
  return(est[,permn[permun,]])
}

library(MASS)

K_list <- c(rep(3,15), rep(5,15)); # graded responses
D <- 3; N <- 1000; J <- 30;
set.seed(235)
true_Mod <- matrix(1,J,D)
true_Mod[ 2:4,1] <- 0
true_Mod[ 20,2] <- 0
true_Mod[ 10:11,3] <- 0
res1=generate_para_pcm(N, J, D, true_Mod, K_list,sL = 0.2, sU = 0.7,aL = 0.5 , aU = 1)
# a: Uniformly distributed between aL (lower bound) and aU (upper bound).
# b: Normally distributed with mean 0 and standard deviation 1 (standard normal).
# σ (sigma): Off-diagonal elements drawn from a uniform distribution with bounds sL (lower) and sU (upper).
data1=generate_response_pcm(res1$thet_true, res1$a_true, res1$b_true, K_list)
MGPCM_data <- list(data = data1, model = true_Mod, params = list(a = res1$a_true, b = res1$b_true, theta = res1$thet_true))
usethis::use_data(MGPCM_data, overwrite = TRUE)
