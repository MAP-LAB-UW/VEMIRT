generate_para<-function(I, J, D, Model, K_list, sL = 0, sU = 0, aL = 0.5 , aU = 1)
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
  aj <- aj*Model 
  # generate b_j
  # parameter b is stored in a J x K dimension array with K+1 being the largest number of categories of one item
  bet <- array(0, dim = c(J, K))
  for(j in 1:J)
    for(k in 1:K_list[j])
      bet[j,k] <- rnorm(1, 0, 1)
  
  # generate theta
  thet <- array(0, dim = c(I, D))
  sigt <- diag(D);
  if(D>1)
  {
    for(i in 1:D)
      for(j in 1:D)
        if(i<j){sigt[i,j] <- runif(1,sL,sU); sigt[j,i] <- sigt[i,j];}
  }
  
  for(i in 1:I)
    thet[i,] <- mvrnorm(1, rep(0, D), sigt);
  # if(D>1)
  #   aj=(promax(aj,m=4))$loadings
  
  return(list(thet_true = thet, a_true = aj, b_true = bet, Sig_true = sigt))
}

generate_response <- function(thet, aj, bet, K_list){
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
library(MASS)
library(mvtnorm)
set.seed(234)

## Comparison against Vemirt -----
N <- 1000; #sample size
J <- 30; # number of item 
K_list <- c(rep(3,15), rep(5,15)); # graded responses e.g. 2: 0,1,2; 1: 0,1
D <- 3; # latent dimension
true_Mod <- matrix(0,J,D)
true_Mod[ 1:12,1] <- 1
true_Mod[ 8:22,2] <- 1
true_Mod[ 15:30,3] <- 1
res1=generate_para(N, J, D, true_Mod, K_list ,sL = 0.2, sU = 0.7,aL = 0.5 , aU = 1)
res2=generate_response(res1$thet_true, res1$a_true, res1$b_true, K_list)
MGPCM_data <- list(data = res2, model = true_Mod, params = list(a = res1$a_true, b = res1$b_true, theta = res1$thet_true))
usethis::use_data(MGPCM_data, overwrite = TRUE)