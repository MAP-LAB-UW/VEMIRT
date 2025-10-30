options(scipen = 999)
library(mirt)
library(dplyr)
library(tidyr)
library(parallel)
library(rstan)
rstan_options(auto_write = T)

# ---- auxiliary functions ----
## ---- start value ----
start_value <- function(resp,J=20,dim,group,invariance,dif_type,impact_type){
  start_mirt <- mirt::multipleGroup(resp[,1:J],1,itemtype = '2PL',
                                    group = group,
                                    invariance = invariance,TOL=5e-2)
  start <- coef(start_mirt,simplify = T)
  start_mus <- unlist(lapply(start, function(x) as.numeric(x$means)))
  b_j = -as.vector(start[[1]]$items[,2]) # b*a
  beta_jv <- matrix(0,nrow = J, ncol = (dim+1))
  beta_jv[,1] <- b_j
  
  bar_a_j = as.vector(start[[1]]$items[,1])
  gamma_jv <- matrix(0,nrow = J, ncol = (dim+1))
  gamma_jv[,1] <- bar_a_j
  
  if(dif_type == "udif"){
    if(impact_type == "fixed"){
      return(list(a_j = as.vector(start[[1]]$items[,1]),
                  sig2_theta = 1,
                  alpha_1v = rep(0,dim),
                  beta_jv = beta_jv,
                  sig2_b_j = rep(1,J)))
    }else if(impact_type == "random"){
      return(list(a_j = as.vector(start[[1]]$items[,1]),
                  sig2_theta = 1,
                  alpha_1v = rep(0,dim),
                  beta_jv = beta_jv,
                  sig2_b_j = rep(1,J),
                  sig2_alpha_0 = 1))
    }
  }else if(dif_type == "nudif"){
    if(impact_type == "fixed"){
      return(list(gamma_jv = gamma_jv,
                  sig2_theta = 1,
                  alpha_1v = rep(0,dim),
                  beta_jv = beta_jv,
                  sig2_b_j = rep(1,J),
                  sig2_bar_a_j = rep(1,J)))
    }else if(impact_type == "random"){
      return(list(gamma_jv = gamma_jv,
                  sig2_theta = 1,
                  alpha_1v = rep(0,dim),
                  beta_jv = beta_jv,
                  sig2_b_j = rep(1,J),
                  sig2_bar_a_j = rep(1,J),
                  sig2_alpha_0 = 1))
    }
  }
}

## ---- eta function ----
eta <- function(y){
  x <- sqrt(y)
  return((1/(2*x))*((1/(1+exp(-x)))-1/2))
}

## ---- a loss function ----
a_loss_code = "
data {
  int J;
  int S;
  int M; //Dimension of main effect
  real lambda;
  array[J] vector<lower = 0>[S] q_sig2_a_js;
  array[J] vector[S] q_mu_a_js;
  array[S] vector[M] X_tilde_s;
  //real threshold;
}

parameters {
  array[J] real<lower = 0> sig2_bar_a_j;
  array[J] vector[M] gamma_jv;//vector<lower = 0>[J] gamma_j;
}

model {
  for (j in 1:J) {
    vector[S] gamma_jv_T_X;
    vector[S] Phi_arg;
    real sqrt_sig2_bar_a_j = sqrt(sig2_bar_a_j[j]);
    
    // Compute gamma_jv^T * X_tilde_s for each sample
    for (s in 1:S)
      gamma_jv_T_X[s] = dot_product(gamma_jv[j], X_tilde_s[s]);
    
    Phi_arg = gamma_jv_T_X / sqrt_sig2_bar_a_j;
    
    // Negative log-likelihood components
    target += -(S + 2 * lambda) * log(sig2_bar_a_j[j]);
    // Instead of the if/else, you can do:
    // target += -(S + 2 * lambda) * log(fmax(threshold, sig2_bar_a_j[j]));

    target += -sum(
      (q_sig2_a_js[j] + square(q_mu_a_js[j] - gamma_jv_T_X)) / sig2_bar_a_j[j]
      + 2 * log(Phi_approx(Phi_arg))
    );
  }
}
"

# ---- GVEM functions ----
## ---- GVEM bonly fixed impact ----
GVEM_bonly_nomultitheta <- function(resp, group_matrix, all_start, sig_start, lambda, c, iter_criteria = 5e2, tau_criteria = 1e-3, rho_N, rho_N2, debias = F, anchor){
  ## Ysij: a list with S elements, each with J columns, Ns rows
  Ysij <- split(resp%>%select(-starts_with("group")), f = resp$groupfull)
  Ysij <- lapply(Ysij, as.matrix)
  Ns <- as.vector(unlist(lapply(Ysij, nrow)))
  S <- length(unique(resp$groupfull))
  J <- ncol(resp%>%select(-starts_with("group")))
  d <- ncol(group_matrix) #2*5*5： 1+4+4=9 | d <- 9
  group_matrix_tilde <- cbind(1,group_matrix)
  dim <- ncol(group_matrix)
  ksi_sij2 <- list()
  a_j <- list()
  sig2_theta <- list()
  alpha_1v <- list()
  beta_jv <- list()
  sig2_b_j <- list()
  
  q_sig2_theta_is <- list()
  q_mu_theta_is <- list()
  q_sig2_b_js <- list()
  q_mu_b_js <- list()
  
  # start value----
  ksi_sij2[[1]] <- lapply(Ns, function(s) {matrix(0.1, nrow = s, ncol = J)})
  a_j[[1]] <-  all_start$a_j
  sig2_theta[[1]] <- all_start$sig2_theta
  alpha_1v[[1]] <- all_start$alpha_1v
  beta_jv[[1]] <- all_start$beta_jv
  sig2_b_j[[1]] <- sig_start$sig2_b_j
  
  q_sig2_theta_is[[1]] <- list()
  q_mu_theta_is[[1]] <- list()
  for (s in 1:S) {
    q_sig2_theta_is[[1]][[s]] <- rep(1,Ns[s])
    q_mu_theta_is[[1]][[s]] <- rep(1,Ns[s])
  }
  q_sig2_b_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_b_js[[1]] <- matrix(0,nrow = J,ncol = S)
  
  iter <- 1
  tau <- 1
  # iter_criteria = 5e2
  # tau_criteria = 1e-3
  
  rescale <- list()
  out_sig2_b0 <- which(sig2_b_j[[1]]!=0)
  in_sig2_b0 <- which(sig2_b_j[[1]]==0)
  sig2_b_j[[1]][out_sig2_b0] <- 1
  
  q_sig2_b_js[[1]][in_sig2_b0,] <- 0
  q_mu_b_js[[1]][in_sig2_b0,] <- (beta_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_b0,]
  
  while((tau >= tau_criteria)&&(iter<=iter_criteria)){
    ## initialize for E-step
    q_sig2_theta_is[[iter+1]] <- list()
    q_mu_theta_is[[iter+1]] <- list()
    q_sig2_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    
    # E-step----
    for (s in 1:S) {
      eta_temp_ij <- eta(ksi_sij2[[iter]][[s]]) #I*J: eta ksij[[s]]
      alpha_multi_temp <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
      beta_multi_temp <- beta_jv[[iter]]%*%t(group_matrix_tilde)
      
      temp2 <- (1+2*sig2_theta[[iter]]*as.vector(eta_temp_ij%*%(a_j[[iter]]^2)))
      q_sig2_theta_is[[iter+1]][[s]] <- sig2_theta[[iter]]/temp2
      q_sig2_theta_is[[iter]][[s]] <- q_sig2_theta_is[[iter+1]][[s]]
      
      temp3 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*") # I*J matrix
      temp4 <- as.vector((Ysij[[s]]-1/2+2*temp3)%*%(a_j[[iter]]))*sig2_theta[[iter]]+
        alpha_multi_temp
      q_mu_theta_is[[iter+1]][[s]] <- temp4/temp2
      q_mu_theta_is[[iter]][[s]] <- q_mu_theta_is[[iter+1]][[s]]
      
      temp5 <- 1+2*sig2_b_j[[iter]][out_sig2_b0]*colSums(eta_temp_ij[,out_sig2_b0,drop=F])
      q_sig2_b_js[[iter+1]][,s][out_sig2_b0] <- sig2_b_j[[iter]][out_sig2_b0]/temp5
      q_sig2_b_js[[iter]][,s][out_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][out_sig2_b0]
      
      temp6 <- sweep(eta_temp_ij,2,a_j[[iter]],"*") # I*J matrix
      temp6 <- temp6[,out_sig2_b0,drop=F]
      temp7 <- colSums(Ysij[[s]][,out_sig2_b0,drop=F] - 1/2 - 2*temp6*q_mu_theta_is[[iter]][[s]])
      temp8 <- beta_multi_temp[,s][out_sig2_b0] - temp7*sig2_b_j[[iter]][out_sig2_b0]
      q_mu_b_js[[iter+1]][,s][out_sig2_b0] <- temp8/temp5
      q_mu_b_js[[iter]][,s][out_sig2_b0] <- q_mu_b_js[[iter+1]][,s][out_sig2_b0]
      
      q_sig2_b_js[[iter+1]][,s][in_sig2_b0] <- 0
      q_sig2_b_js[[iter]][,s][in_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][in_sig2_b0]
      q_mu_b_js[[iter+1]][,s][in_sig2_b0] <- beta_multi_temp[,s][in_sig2_b0]
      q_mu_b_js[[iter]][,s][in_sig2_b0] <- q_mu_b_js[[iter+1]][,s][in_sig2_b0]
    }
    
    ## initialize for M-step
    tau_M <- 1
    sig2_theta[[iter+1]] <- sig2_theta[[iter]]
    alpha_1v[[iter+1]] <- alpha_1v[[iter]]
    beta_jv[[iter+1]] <- beta_jv[[iter]]
    sig2_b_j[[iter+1]] <- sig2_b_j[[iter]]
    a_j[[iter+1]] <- a_j[[iter]]
    ksi_sij2[[iter+1]] <- ksi_sij2[[iter]]
    
    while (tau_M > tau_criteria) {
      sig2_theta_old <- sig2_theta[[iter+1]]
      alpha_1v_old <- alpha_1v[[iter+1]]
      b_j_old <- beta_jv[[iter+1]]
      sig2_b_j_old<-sig2_b_j[[iter+1]]
      a_j_old <- a_j[[iter+1]]
      ksi_sij2_old <- ksi_sij2[[iter+1]]
      
      a_j_num_temp <- list()
      a_j_deno_temp <- list()
      alpha_num_temp <- list()
      alpha_deno_temp <- list()
      beta_in_deno_temp <- list()
      beta_in_num_temp <- list()
      beta_deno_temp <- list()
      beta_num_temp <- list()
      
      # M-step----
      for (s in 1:S) {
        #alpha_multi_temp_M <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
        eta_temp_ij_M <- eta(ksi_sij2[[iter]][[s]])
        
        outer_temp1 <- outer(q_mu_theta_is[[iter+1]][[s]]^2 + q_sig2_theta_is[[iter+1]][[s]],
                             a_j[[iter]]^2,"*")
        outer_temp2 <- outer(q_mu_theta_is[[iter+1]][[s]],
                             a_j[[iter]]*q_mu_b_js[[iter+1]][,s],'*')
        temp_1_M <- outer_temp1 - 2*outer_temp2
        temp_2_M <- q_mu_b_js[[iter+1]][,s]^2 + q_sig2_b_js[[iter+1]][,s]
        ksi_sij2[[iter+1]][[s]] <- sweep(temp_1_M,2,temp_2_M,"+")
        
        temp_3_M <- sweep(eta_temp_ij_M,2,q_mu_b_js[[iter+1]][,s],"*")
        temp_4_M <- Ysij[[s]]-1/2+2*temp_3_M
        temp_5_M <- temp_4_M*q_mu_theta_is[[iter+1]][[s]] #I*J * I
        a_j_num_temp[[s]] <- colSums(temp_5_M)
        
        temp_6_M <- (q_mu_theta_is[[iter+1]][[s]])^2 + q_sig2_theta_is[[iter+1]][[s]]
        a_j_deno_temp[[s]] <- colSums(2*eta_temp_ij_M*temp_6_M) #I*J * I vector
        
        alpha_deno_temp[[s]] <- Ns[s]*(as.matrix(group_matrix[s,])%*%(group_matrix[s,]))
        
        alpha_num_temp[[s]] <- as.matrix(group_matrix[s,])*(sum(q_mu_theta_is[[iter+1]][[s]]))
        
        beta_deno_temp[[s]] <- (as.matrix(group_matrix_tilde[s,])%*%(group_matrix_tilde[s,]))
        
        deno_temp <- array(beta_deno_temp[[s]], dim = c((1+dim),(1+dim),J))
        
        beta_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M), FUN = "*")
        
        beta_in_num_temp[[s]] <- matrix(colSums(1/2 - Ysij[[s]] +
                                                  sweep((2*eta_temp_ij_M*q_mu_theta_is[[iter+1]][[s]]),
                                                        2,a_j[[iter]],"*")),ncol = 1)%*%group_matrix_tilde[s,]
      }
      
      a_j[[iter+1]] <- as.vector(Reduce('+',a_j_num_temp)/Reduce('+',a_j_deno_temp))
      
      sig2_theta[[iter+1]] <- 1 #Reduce("+",sig2_theta_temp)/S
      
      # rescale[[iter]] <- sqrt(1/sig2_theta[[iter+1]])
      # a_j[[iter+1]] <- a_j[[iter+1]]/rescale[[iter]]
      # sig2_theta[[iter+1]] <- sig2_theta[[iter+1]]*rescale[[iter]]^2
      
      alpha_1v[[iter+1]] <- as.vector(solve(Reduce('+',alpha_deno_temp))%*%Reduce('+',alpha_num_temp))
      
      beta_num_temp <- q_mu_b_js[[iter+1]]%*%group_matrix_tilde
      
      beta_jv[[iter+1]] <- t(solve(Reduce('+',beta_deno_temp))%*%t(beta_num_temp))#matrix(,nrow = 20,ncol = (dim+1))
      
      all_deno <- Reduce('+',beta_in_deno_temp)
      if(length(in_sig2_b0)!=0){
        for (inj in 1:length(in_sig2_b0)) {
          beta_jv[[iter+1]][in_sig2_b0[inj],] <- as.vector(solve(all_deno[,,in_sig2_b0[inj]])%*%matrix((Reduce('+',beta_in_num_temp))[in_sig2_b0[inj],],ncol = 1))
        }
      }
      
      #beta_jv[[iter+1]][in_sig2_b0,1] <- (Reduce('+',beta_in_num_temp))[in_sig2_b0]/(Reduce('+',beta_in_deno_temp))[in_sig2_b0]
      #beta_jv[[iter+1]][out_sig2_b0, ] <- solve(Reduce('+',beta_deno_temp))%*%t(beta_num_temp)
      
      sig2_b_j[[iter+1]][out_sig2_b0] <- rowSums(q_sig2_b_js[[iter]][out_sig2_b0,,drop=F] +
                                                   ((beta_jv[[iter+1]]%*%t(group_matrix_tilde))[out_sig2_b0,] -
                                                      q_mu_b_js[[iter]][out_sig2_b0,,drop=F])^2)/(S+2*lambda)
      
      sig2_b_j[[iter+1]][in_sig2_b0] <- 0
      
      beta_jv[[iter+1]][anchor,-1] <- 0
      
      tau_M <- max(#abs(ksi_sij2[[iter+1]] - ksi_sij2_old),
        abs(sig2_theta[[iter+1]] - sig2_theta_old),
        abs(alpha_1v[[iter+1]] - alpha_1v_old),
        abs(beta_jv[[iter+1]] - b_j_old),
        abs(sig2_b_j[[iter+1]]-sig2_b_j_old),
        abs(a_j[[iter+1]] - a_j_old))
    }
    ## update criteria
    tau <- max(abs(a_j[[iter+1]] - a_j[[iter]]),
               abs(sig2_theta[[iter+1]] - sig2_theta[[iter]]),
               abs(alpha_1v[[iter+1]] - alpha_1v[[iter]]),
               as.vector(abs(beta_jv[[iter+1]] - beta_jv[[iter]])),
               abs(sig2_b_j[[iter+1]]-sig2_b_j[[iter]]))
    
    iter <- iter+1
    #print(tau)
  }
  
  if(!debias){
    sig2_b_j[[iter]][sig2_b_j[[iter]]< rho_N] <- 0
    
    out_sig2_b0 <- which(sig2_b_j[[iter]]!=0)
    in_sig2_b0 <- which(sig2_b_j[[iter]]==0)
  }
  
  sig2_b_j2 <- sig2_b_j[[iter]]
  which_temp <- which((sig2_b_j2 < rho_N2)&(sig2_b_j2 != 0))
  sig2_b_j2[which_temp] <- rho_N2
  
  temp1_like <- list()
  temp2_like <- list()
  temp3_like <- list()
  
  for (s in 1:S) {
    part1 <- -log(1+exp(-sqrt(ksi_sij2[[iter]][[s]]))) # IJ
    part2 <- (Ysij[[s]]-1/2)*
      sweep(outer(q_mu_theta_is[[iter]][[s]], # IJ
                  a_j[[iter]],"*"),2, q_mu_b_js[[iter]][,s],"-") - #IJ-J
      (1/2)*sqrt(ksi_sij2[[iter]][[s]]) #IJ
    
    part31 <- q_mu_b_js[[iter]][,s]^2 + q_sig2_b_js[[iter]][,s] # J
    part32 <- -2*outer(q_mu_theta_is[[iter]][[s]], a_j[[iter]]*q_mu_b_js[[iter]][,s],'*') + #IJ
      outer(q_mu_theta_is[[iter]][[s]]^2 + (q_sig2_theta_is[[iter]][[s]]),#IJ
            a_j[[iter]]^2,"*") - ksi_sij2[[iter]][[s]] #IJ
    part3 <- eta(ksi_sij2[[iter]][[s]])*sweep(part32,2,part31, "+")
    temp1_like[[s]] <- sum(part1 + part2 - part3)
    
    part4 <-  (Ns[s]+length(out_sig2_b0))*log(2*pi)
    part5 <- sum(log(sig2_theta[[iter]]) + 
                   (q_sig2_theta_is[[iter]][[s]] + 
                      (q_mu_theta_is[[iter]][[s]] - as.numeric(alpha_1v[[iter]]%*%group_matrix[s,]))^2)/sig2_theta[[iter]])
    part6 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0])
    part7 <- sum(log(sig2_b_j2[out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j2[out_sig2_b0])
    
    # part7 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
    #                (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
    #                   (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
    #                      b_j[[iter]][out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0]) +
    #   length(in_sig2_b0)*(log(rho_N) + 1)
    
    temp2_like[[s]] <- part4+part5+part6
    temp3_like[[s]] <- part4+part5+part7
  }
  
  likelihood <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)
  likelihood2 <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)
  
  BIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    length(out_sig2_b0)*log(sum(Ns))
  BIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    length(out_sig2_b0)*log(sum(Ns))
  GIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    length(out_sig2_b0)*c*log(sum(Ns))*log(log(sum(Ns)))
  GIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    length(out_sig2_b0)*c*log(sum(Ns))*log(log(sum(Ns)))
  # EBIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
  #   length(out_sig2_b0)*log(sum(Ns)) + 2*g*log(choose(J,length(out_sig2_b0)))
  
  all_return <- list(a_j = a_j[[iter]],
                     sig2_theta = sig2_theta[[iter]],
                     alpha_1v = alpha_1v[[iter]],
                     b_j = beta_jv[[iter]],
                     sig2_b_j = sig2_b_j[[iter]],
                     likelihood = likelihood,
                     likelihood2 = likelihood2,
                     BIC = BIC,
                     BIC2 = BIC2,
                     GIC = GIC,
                     GIC2 = GIC2,
                     iter = iter)
  return(all_return)
}

## ---- GVEM bonly random impact ----
GVEM_bonly_multitheta <- function(resp, group_matrix, all_start, sig_start, lambda, c, iter_criteria = 5e2, tau_criteria = 1e-3, rho_N, rho_N2, debias = F, anchor){
  ## Ysij: a list with S elements, each with J columns, Ns rows
  Ysij <- split(resp%>%select(-starts_with("group")), f = resp$groupfull)
  Ysij <- lapply(Ysij, as.matrix)
  Ns <- as.vector(unlist(lapply(Ysij, nrow)))
  S <- length(unique(resp$groupfull))
  J <- ncol(resp%>%select(-starts_with("group")))
  d <- ncol(group_matrix) #2*5*5： 1+4+4=9 | d <- 9
  group_matrix_tilde <- cbind(1,group_matrix)
  dim <- ncol(group_matrix)
  ksi_sij2 <- list()
  a_j <- list()
  sig2_alpha_0 <- list()
  sig2_theta <- list()
  alpha_1v <- list()
  beta_jv <- list()
  sig2_b_j <- list()
  
  q_sig2_alpha_0s <- list()
  q_mu_alpha_0s <- list()
  q_sig2_theta_is <- list()
  q_mu_theta_is <- list()
  q_sig2_b_js <- list()
  q_mu_b_js <- list()
  
  # start value----
  ksi_sij2[[1]] <- lapply(Ns, function(s) {matrix(0.1, nrow = s, ncol = J)})
  a_j[[1]] <-  all_start$a_j
  sig2_alpha_0[[1]] <- all_start$sig2_alpha_0
  sig2_theta[[1]] <- all_start$sig2_theta
  alpha_1v[[1]] <- all_start$alpha_1v
  beta_jv[[1]] <- all_start$beta_jv
  sig2_b_j[[1]] <- sig_start$sig2_b_j
  
  q_sig2_alpha_0s[[1]] <- rep(1,S)
  q_mu_alpha_0s[[1]] <- rep(0,S)
  q_sig2_theta_is[[1]] <- list()
  q_mu_theta_is[[1]] <- list()
  for (s in 1:S) {
    q_sig2_theta_is[[1]][[s]] <- rep(1,Ns[s])
    q_mu_theta_is[[1]][[s]] <- rep(1,Ns[s])
  }
  q_sig2_b_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_b_js[[1]] <- matrix(0,nrow = J,ncol = S)
  
  iter <- 1
  tau <- 1
  # iter_criteria = 5e2
  # tau_criteria = 1e-3
  
  rescale <- list()
  out_sig2_b0 <- which(sig2_b_j[[1]]!=0)
  in_sig2_b0 <- which(sig2_b_j[[1]]==0)
  sig2_b_j[[1]][out_sig2_b0] <- 1
  
  q_sig2_b_js[[1]][in_sig2_b0,] <- 0
  q_mu_b_js[[1]][in_sig2_b0,] <- (beta_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_b0,]
  
  
  while((tau >= tau_criteria)&&(iter<=iter_criteria)){
    ## initialize for E-step
    q_sig2_alpha_0s[[iter+1]] <- rep(NA,S)
    q_mu_alpha_0s[[iter+1]] <- rep(NA,S)
    q_sig2_theta_is[[iter+1]] <- list()
    q_mu_theta_is[[iter+1]] <- list()
    q_sig2_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    
    # E-step----
    for (s in 1:S) {
      eta_temp_ij <- eta(ksi_sij2[[iter]][[s]]) #I*J: eta ksij[[s]]
      alpha_multi_temp <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
      beta_multi_temp <- beta_jv[[iter]]%*%t(group_matrix_tilde)
      
      temp1 <- Ns[s]*sig2_alpha_0[[iter]]+sig2_theta[[iter]]
      q_sig2_alpha_0s[[iter+1]][s] <- (sig2_alpha_0[[iter]]*sig2_theta[[iter]])/temp1
      q_sig2_alpha_0s[[iter]][s] <- q_sig2_alpha_0s[[iter+1]][s]
      
      q_mu_alpha_0s[[iter+1]][s] <- (sig2_alpha_0[[iter]]/temp1)*
        (sum(q_mu_theta_is[[iter]][[s]]) -  Ns[s]*alpha_multi_temp)
      q_mu_alpha_0s[[iter]][s] <- q_mu_alpha_0s[[iter+1]][s]
      
      temp2 <- (1+2*sig2_theta[[iter]]*as.vector(eta_temp_ij%*%(a_j[[iter]]^2)))
      q_sig2_theta_is[[iter+1]][[s]] <- sig2_theta[[iter]]/temp2
      q_sig2_theta_is[[iter]][[s]] <- q_sig2_theta_is[[iter+1]][[s]]
      
      temp3 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*") # I*J matrix
      temp4 <- as.vector((Ysij[[s]]-1/2+2*temp3)%*%(a_j[[iter]]))*sig2_theta[[iter]]+
        q_mu_alpha_0s[[iter]][s]+alpha_multi_temp
      q_mu_theta_is[[iter+1]][[s]] <- temp4/temp2
      q_mu_theta_is[[iter]][[s]] <- q_mu_theta_is[[iter+1]][[s]]
      
      temp5 <- 1+2*sig2_b_j[[iter]][out_sig2_b0]*colSums(eta_temp_ij[,out_sig2_b0,drop=F])
      q_sig2_b_js[[iter+1]][,s][out_sig2_b0] <- sig2_b_j[[iter]][out_sig2_b0]/temp5
      q_sig2_b_js[[iter]][,s][out_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][out_sig2_b0]
      
      temp6 <- sweep(eta_temp_ij,2,a_j[[iter]],"*") # I*J matrix
      temp6 <- temp6[,out_sig2_b0,drop = F]
      temp7 <- colSums(Ysij[[s]][,out_sig2_b0,drop=F] - 1/2 - 2*temp6*q_mu_theta_is[[iter]][[s]])
      temp8 <- beta_multi_temp[,s][out_sig2_b0] - temp7*sig2_b_j[[iter]][out_sig2_b0]
      q_mu_b_js[[iter+1]][,s][out_sig2_b0] <- temp8/temp5
      q_mu_b_js[[iter]][,s][out_sig2_b0] <- q_mu_b_js[[iter+1]][,s][out_sig2_b0]
      
      q_sig2_b_js[[iter+1]][,s][in_sig2_b0] <- 0
      q_sig2_b_js[[iter]][,s][in_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][in_sig2_b0]
      q_mu_b_js[[iter+1]][,s][in_sig2_b0] <- beta_multi_temp[,s][in_sig2_b0]
      q_mu_b_js[[iter]][,s][in_sig2_b0] <- q_mu_b_js[[iter+1]][,s][in_sig2_b0]
    }
    
    ## initialize for M-step
    tau_M <- 1
    ## initialize for M-step
    sig2_theta[[iter+1]] <- sig2_theta[[iter]]
    alpha_1v[[iter+1]] <- alpha_1v[[iter]]
    beta_jv[[iter+1]] <- beta_jv[[iter]]
    sig2_b_j[[iter+1]] <- sig2_b_j[[iter]]
    a_j[[iter+1]] <- a_j[[iter]]
    ksi_sij2[[iter+1]] <- ksi_sij2[[iter]]
    sig2_alpha_0[[iter+1]] <- sig2_alpha_0[[iter]]
    
    while (tau_M > tau_criteria) {
      
      sig2_theta_old <- sig2_theta[[iter+1]]
      alpha_1v_old <- alpha_1v[[iter+1]]
      b_j_old <- beta_jv[[iter+1]]
      sig2_b_j_old<-sig2_b_j[[iter+1]]
      a_j_old <- a_j[[iter+1]]
      ksi_sij2_old <- ksi_sij2[[iter+1]]
      sig2_alpha_0_old <- sig2_alpha_0[[iter+1]]
      
      a_j_num_temp <- list()
      a_j_deno_temp <- list()
      alpha_num_temp <- list()
      alpha_deno_temp <- list()
      beta_in_deno_temp <- list()
      beta_in_num_temp <- list()
      beta_deno_temp <- list()
      beta_num_temp <- list()
      
      # M-step----
      for (s in 1:S) {
        outer_temp1 <- outer(q_mu_theta_is[[iter+1]][[s]]^2 + q_sig2_theta_is[[iter+1]][[s]],
                             a_j[[iter+1]]^2,"*")
        outer_temp2 <- outer(q_mu_theta_is[[iter+1]][[s]],
                             a_j[[iter+1]]*q_mu_b_js[[iter+1]][,s],'*')
        temp_1_M <- outer_temp1 - 2*outer_temp2
        temp_2_M <- q_mu_b_js[[iter+1]][,s]^2 + q_sig2_b_js[[iter+1]][,s]
        ksi_sij2[[iter+1]][[s]] <- sweep(temp_1_M,2,temp_2_M,"+")
        eta_temp_ij_M <- eta(ksi_sij2[[iter+1]][[s]])
        
        temp_3_M <- sweep(eta_temp_ij_M,2,q_mu_b_js[[iter+1]][,s],"*")
        temp_4_M <- Ysij[[s]]-1/2+2*temp_3_M
        temp_5_M <- temp_4_M*q_mu_theta_is[[iter+1]][[s]] #I*J * I
        a_j_num_temp[[s]] <- colSums(temp_5_M)
        
        temp_6_M <- (q_mu_theta_is[[iter+1]][[s]])^2 + q_sig2_theta_is[[iter+1]][[s]]
        a_j_deno_temp[[s]] <- colSums(2*eta_temp_ij_M*temp_6_M) #I*J * I vector
        
        alpha_deno_temp[[s]] <- Ns[s]*(as.matrix(group_matrix[s,])%*%(group_matrix[s,]))
        
        alpha_num_temp[[s]] <- as.matrix(group_matrix[s,])*(sum(q_mu_theta_is[[iter+1]][[s]]) - 
                                                              Ns[s]*q_mu_alpha_0s[[iter+1]][s])
        
        beta_deno_temp[[s]] <- as.matrix(group_matrix_tilde[s,])%*%group_matrix_tilde[s,]
        
        deno_temp <- array(beta_deno_temp[[s]], dim = c((1+dim),(1+dim),J))
        
        beta_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M), FUN = "*")
        
        beta_in_num_temp[[s]] <- matrix(colSums(1/2 - Ysij[[s]] +
                                                  sweep((2*eta_temp_ij_M*q_mu_theta_is[[iter+1]][[s]]),
                                                        2,a_j[[iter+1]],"*")))%*%group_matrix_tilde[s,]
      }
      
      a_j[[iter+1]] <- as.vector(Reduce('+',a_j_num_temp)/Reduce('+',a_j_deno_temp))
      
      sig2_alpha_0[[iter+1]] <- mean(q_sig2_alpha_0s[[iter+1]] + q_mu_alpha_0s[[iter+1]]^2)
      sig2_theta[[iter+1]] <- 1 # Reduce('+',sig2_theta_temp)/S
      
      alpha_1v[[iter+1]] <- as.vector(solve(Reduce('+',alpha_deno_temp))%*%Reduce('+',alpha_num_temp))
      
      beta_num_temp <- q_mu_b_js[[iter+1]]%*%group_matrix_tilde
      
      beta_jv[[iter+1]] <- t(solve(Reduce('+',beta_deno_temp))%*%t(beta_num_temp))
      all_deno <- Reduce('+',beta_in_deno_temp)
      if(length(in_sig2_b0)!=0){
        for (inj in 1:length(in_sig2_b0)) {
          beta_jv[[iter+1]][in_sig2_b0[inj], ] <- as.vector(solve(all_deno[,,in_sig2_b0[inj]])%*%matrix((Reduce('+',beta_in_num_temp))[in_sig2_b0[inj],],ncol = 1))
        }
      }
      
      # b_j[[iter+1]] <- rowMeans(q_mu_b_js[[iter+1]])
      # 
      # b_j[[iter+1]][in_sig2_b0] <- (Reduce('+',beta_in_num_temp))[in_sig2_b0]/(Reduce('+',beta_in_deno_temp))[in_sig2_b0]
      
      sig2_b_j[[iter+1]][out_sig2_b0] <- rowSums(q_sig2_b_js[[iter+1]][out_sig2_b0,,drop=F] +
                                                   ((beta_jv[[iter]]%*%t(group_matrix_tilde))[out_sig2_b0,] -
                                                      q_mu_b_js[[iter+1]][out_sig2_b0,,drop=F])^2)/(S+2*lambda)
      sig2_b_j[[iter+1]][in_sig2_b0] <- 0
      
      beta_jv[[iter+1]][anchor,-1] <- 0
      
      tau_M <- max(abs(sig2_theta[[iter+1]] - sig2_theta_old),
                   abs(alpha_1v[[iter+1]] - alpha_1v_old),
                   abs(beta_jv[[iter+1]] - b_j_old),
                   abs(sig2_b_j[[iter+1]]-sig2_b_j_old),
                   abs(a_j[[iter+1]] - a_j_old),
                   abs(sig2_alpha_0[[iter+1]] - sig2_alpha_0_old))
    }
    ## update criteria
    tau <- max(abs(a_j[[iter+1]] - a_j[[iter]]),
               abs(sig2_alpha_0[[iter+1]] - sig2_alpha_0[[iter]]),
               abs(sig2_theta[[iter+1]] - sig2_theta[[iter]]),
               abs(alpha_1v[[iter+1]] - alpha_1v[[iter]]),
               as.vector(abs(beta_jv[[iter+1]] - beta_jv[[iter]])),
               abs(sig2_b_j[[iter+1]]-sig2_b_j[[iter]]))
    
    iter <- iter+1
    #print(tau)
  }
  
  if(!debias){
    sig2_b_j[[iter]][sig2_b_j[[iter]]< rho_N] <- 0
    
    out_sig2_b0 <- which(sig2_b_j[[iter]]!=0)
    in_sig2_b0 <- which(sig2_b_j[[iter]]==0)
  }
  
  sig2_b_j2 <- sig2_b_j[[iter]]
  which_temp <- which((sig2_b_j2 < rho_N2)&(sig2_b_j2 != 0))
  sig2_b_j2[which_temp] <- rho_N2
  
  temp1_like <- list()
  temp2_like <- list()
  temp3_like <- list()
  
  for (s in 1:S) {
    part1 <- -log(1+exp(-sqrt(ksi_sij2[[iter]][[s]]))) # IJ
    part2 <- (Ysij[[s]]-1/2)*
      sweep(outer(q_mu_theta_is[[iter]][[s]], # IJ
                  a_j[[iter]],"*"),2, q_mu_b_js[[iter]][,s],"-") - #IJ-J
      (1/2)*sqrt(ksi_sij2[[iter]][[s]]) #IJ
    
    part31 <- q_mu_b_js[[iter]][,s]^2 + q_sig2_b_js[[iter]][,s] # J
    part32 <- -2*outer(q_mu_theta_is[[iter]][[s]], a_j[[iter]]*q_mu_b_js[[iter]][,s],'*') + #IJ
      outer(q_mu_theta_is[[iter]][[s]]^2 + (q_sig2_theta_is[[iter]][[s]]),#IJ
            a_j[[iter]]^2,"*") - ksi_sij2[[iter]][[s]] #IJ
    part3 <- eta(ksi_sij2[[iter]][[s]])*sweep(part32,2,part31, "+")
    temp1_like[[s]] <- sum(part1 + part2 - part3)
    
    part4 <-  (Ns[s]+length(out_sig2_b0)+1)*log(2*pi) + log(sig2_alpha_0[[iter]]) +
      (q_sig2_alpha_0s[[iter]][[s]] + q_mu_alpha_0s[[iter]][[s]]^2)/sig2_alpha_0[[iter]]
    
    part5 <- sum(log(sig2_theta[[iter]]) + 
                   (q_sig2_theta_is[[iter]][[s]] + q_sig2_alpha_0s[[iter]][[s]] +
                      (q_mu_theta_is[[iter]][[s]] - q_mu_alpha_0s[[iter]][[s]] -
                         as.numeric(alpha_1v[[iter]]%*%group_matrix[s,]))^2)/sig2_theta[[iter]])
    part6 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde)[,s])[out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0])
    part7 <- sum(log(sig2_b_j2[out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde)[,s])[out_sig2_b0])^2)/sig2_b_j2[out_sig2_b0])
    
    # part7 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
    #                (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
    #                   (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
    #                      (b_j[[iter]][out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0]) +
    #   length(in_sig2_b0)*(log(rho_N) + 1)
    
    temp2_like[[s]] <- part4+part5+part6
    temp3_like[[s]] <- part4+part5+part7
  }
  
  likelihood <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)
  likelihood2 <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)
  
  BIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    length(out_sig2_b0)*log(sum(Ns))
  BIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    length(out_sig2_b0)*log(sum(Ns))
  GIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    length(out_sig2_b0)*c*log(sum(Ns))*log(log(sum(Ns)))
  GIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    length(out_sig2_b0)*c*log(sum(Ns))*log(log(sum(Ns)))
  # EBIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
  #   length(out_sig2_b0)*log(sum(Ns)) + 2*g*log(choose(J,length(out_sig2_b0)))
  
  all_return <- list(a_j = a_j[[iter]],
                     sig2_theta = sig2_theta[[iter]],
                     alpha_1v = alpha_1v[[iter]],
                     b_j = beta_jv[[iter]],
                     sig2_b_j = sig2_b_j[[iter]],
                     sig2_alpha_0 = sig2_alpha_0[[iter]],
                     likelihood = likelihood,
                     likelihood2 = likelihood2,
                     BIC = BIC,
                     BIC2 = BIC2,
                     GIC = GIC,
                     GIC2 = GIC2,
                     iter = iter,
                     tau = tau)
  return(all_return)
}

## ---- GVEM ab fixed impact ----
GVEM_ab_nomultitheta <- function(resp, group_matrix, all_start, sig_start, lambda, c, iter_criteria = 5e2, tau_criteria = 1e-3,rho_N, rho_N2, rho_Na, rho_Na2, debias = F, anchor){
  
  a_loss <- stan_model(model_code = a_loss_code)
  
  ## Ysij: a list with S elements, each with J columns, Ns rows
  Ysij <- split(resp%>%select(-starts_with("group")), f = resp$groupfull)
  Ysij <- lapply(Ysij, as.matrix)
  Ns <- as.vector(unlist(lapply(Ysij, nrow)))
  S <- length(unique(resp$groupfull))
  J <- ncol(resp%>%select(-starts_with("group")))
  d <- ncol(group_matrix) #2*5*5： 1+4+4=9 | d <- 9
  group_matrix_tilde <- cbind(1,group_matrix)
  dim <- ncol(group_matrix)
  ksi_sij2 <- list()
  sig2_theta <- list()
  alpha_1v <- list()
  beta_jv <- list()
  sig2_b_j <- list()
  gamma_jv <- list()
  sig2_bar_a_j <- list()
  
  q_sig2_theta_is <- list()
  q_mu_theta_is <- list()
  q_sig2_b_js <- list()
  q_mu_b_js <- list()
  q_sig2_bar_a_js <- list()
  q_mu_bar_a_js <- list()
  q_sig2_a_js <- list()
  q_mu_a_js <- list()
  
  # start value----
  ksi_sij2[[1]] <- lapply(Ns, function(s) {matrix(0.1, nrow = s, ncol = J)})
  sig2_theta[[1]] <- all_start$sig2_theta
  alpha_1v[[1]] <- all_start$alpha_1v
  beta_jv[[1]] <- all_start$beta_jv
  sig2_b_j[[1]] <- sig_start$sig2_b_j
  gamma_jv[[1]] <- all_start$gamma_jv
  sig2_bar_a_j[[1]] <- sig_start$sig2_bar_a_j
  
  q_sig2_theta_is[[1]] <- list()
  q_mu_theta_is[[1]] <- list()
  for (s in 1:S) {
    q_sig2_theta_is[[1]][[s]] <- rep(1,Ns[s])
    q_mu_theta_is[[1]][[s]] <- rep(1,Ns[s])
  }
  q_sig2_b_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_b_js[[1]] <- matrix(0,nrow = J,ncol = S)
  q_sig2_bar_a_js[[1]] <- matrix(NA,nrow = J,ncol = S)
  q_mu_bar_a_js[[1]] <- matrix(NA,nrow = J,ncol = S)
  q_sig2_a_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_a_js[[1]] <- matrix(1,nrow = J,ncol = S)
  
  iter <- 1
  tau <- 1
  tau_list <- list()
  # iter_criteria = 5e2
  # tau_criteria = 1e-3
  
  rescale <- list()
  out_sig2_b0 <- which(sig2_b_j[[1]]!=0)
  in_sig2_b0 <- which(sig2_b_j[[1]]==0)
  out_sig2_a0 <- which(sig2_bar_a_j[[1]]!=0)
  in_sig2_a0 <- which(sig2_bar_a_j[[1]]==0)
  sig2_b_j[[1]][out_sig2_b0] <- 1
  sig2_bar_a_j[[1]][out_sig2_a0] <- 1
  
  q_sig2_a_js[[1]][in_sig2_a0,] <- 0
  q_mu_a_js[[1]][in_sig2_a0,] <- (gamma_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_a0,]
  q_sig2_b_js[[1]][in_sig2_b0,] <- 0
  q_mu_b_js[[1]][in_sig2_b0,] <- (beta_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_b0,]
  
  while((tau >= tau_criteria)&&(iter<=iter_criteria)){
    ## initialize for E-step
    q_sig2_theta_is[[iter+1]] <- list()
    q_mu_theta_is[[iter+1]] <- list()
    q_sig2_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_sig2_bar_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_bar_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_sig2_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    
    # E-step----
    for (s in 1:S) {
      eta_temp_ij <- eta(ksi_sij2[[iter]][[s]]) #I*J: eta ksij[[s]]
      alpha_multi_temp <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
      beta_multi_temp <- beta_jv[[iter]]%*%t(group_matrix_tilde)
      gamma_multi_temp <- gamma_jv[[iter]]%*%t(group_matrix_tilde)
      
      temp2 <- (1+2*sig2_theta[[iter]]*
                  as.vector(eta_temp_ij%*%(q_mu_a_js[[iter]][,s]^2+q_sig2_a_js[[iter]][,s])))
      q_sig2_theta_is[[iter+1]][[s]] <- sig2_theta[[iter]]/temp2
      q_sig2_theta_is[[iter]][[s]] <- q_sig2_theta_is[[iter+1]][[s]]
      
      temp3 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*") # I*J matrix
      temp4 <- as.vector((Ysij[[s]]-1/2+2*temp3)%*%(q_mu_a_js[[iter]][,s]))*sig2_theta[[iter]]+
        alpha_multi_temp
      q_mu_theta_is[[iter+1]][[s]] <- temp4/temp2
      q_mu_theta_is[[iter]][[s]] <- q_mu_theta_is[[iter+1]][[s]]
      
      temp9 <- eta_temp_ij*(q_mu_theta_is[[iter]][[s]]^2 + q_sig2_theta_is[[iter]][[s]])
      temp9 <- temp9[,out_sig2_a0,drop=F]
      temp10 <- 1+2*sig2_bar_a_j[[iter]][out_sig2_a0]*colSums(temp9)
      q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0] <- sig2_bar_a_j[[iter]][out_sig2_a0]/temp10
      q_sig2_bar_a_js[[iter]][,s][out_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]
      
      temp11 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*")
      temp11 <- temp11[,out_sig2_a0,drop=F]
      temp12 <- colSums((Ysij[[s]][,out_sig2_a0,drop=F] - 1/2 + 2*temp11)*q_mu_theta_is[[iter]][[s]])
      temp13 <- gamma_multi_temp[,s][out_sig2_a0] + temp12*sig2_bar_a_j[[iter]][out_sig2_a0]
      q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0] <- temp13/temp10
      q_mu_bar_a_js[[iter]][,s][out_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_sig2_bar_a_js[[iter+1]][,s][in_sig2_a0] <- 0
      q_sig2_bar_a_js[[iter]][,s][in_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][in_sig2_a0]
      
      q_mu_bar_a_js[[iter+1]][,s][in_sig2_a0] <- gamma_multi_temp[,s][in_sig2_a0]
      q_mu_bar_a_js[[iter]][,s][in_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][in_sig2_a0]
      
      bar2_ratio <- (q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0]^2)/q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]
      temp14 <- (sqrt(bar2_ratio)*exp(-(bar2_ratio/2)))/
        (sqrt(2*pi)*pnorm(sqrt(bar2_ratio)))
      temp15 <- exp(-bar2_ratio)/(2*pi*(pnorm(sqrt(bar2_ratio)))^2)
      q_sig2_a_js[[iter+1]][,s][out_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]*
        (1-temp14-temp15)
      q_sig2_a_js[[iter]][,s][out_sig2_a0] <- q_sig2_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_mu_a_js[[iter+1]][,s][out_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0] +
        (sqrt(q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0])*exp(-(bar2_ratio/2)))/(sqrt(2*pi)*pnorm(sqrt(bar2_ratio)))
      q_mu_a_js[[iter]][,s][out_sig2_a0] <- q_mu_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_sig2_a_js[[iter+1]][,s][in_sig2_a0] <- 0
      q_sig2_a_js[[iter]][,s][in_sig2_a0] <- q_sig2_a_js[[iter+1]][,s][in_sig2_a0]
      
      q_mu_a_js[[iter+1]][,s][in_sig2_a0] <- gamma_multi_temp[,s][in_sig2_a0]
      q_mu_a_js[[iter]][,s][in_sig2_a0] <- q_mu_a_js[[iter+1]][,s][in_sig2_a0]
      
      temp5 <- 1+2*sig2_b_j[[iter]][out_sig2_b0]*colSums(eta_temp_ij[,out_sig2_b0,drop=F])
      q_sig2_b_js[[iter+1]][,s][out_sig2_b0] <- sig2_b_j[[iter]][out_sig2_b0]/temp5
      q_sig2_b_js[[iter]][,s][out_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][out_sig2_b0]
      
      temp6 <- sweep(eta_temp_ij,2,q_mu_a_js[[iter]][,s],"*") # I*J matrix
      temp6 <- temp6[,out_sig2_b0,drop=F]
      temp7 <- colSums(Ysij[[s]][,out_sig2_b0,drop=F] - 1/2 - 2*temp6*q_mu_theta_is[[iter]][[s]])
      temp8 <- beta_multi_temp[,s][out_sig2_b0] - temp7*sig2_b_j[[iter]][out_sig2_b0]
      q_mu_b_js[[iter+1]][,s][out_sig2_b0] <- temp8/temp5
      q_mu_b_js[[iter]][,s][out_sig2_b0] <- q_mu_b_js[[iter+1]][,s][out_sig2_b0]
      
      q_sig2_b_js[[iter+1]][,s][in_sig2_b0] <- 0
      q_sig2_b_js[[iter]][,s][in_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][in_sig2_b0]
      q_mu_b_js[[iter+1]][,s][in_sig2_b0] <- beta_multi_temp[,s][in_sig2_b0]
      q_mu_b_js[[iter]][,s][in_sig2_b0] <- q_mu_b_js[[iter+1]][,s][in_sig2_b0]
    }
    
    tau_M <- 1
    ## initialize for M-step
    sig2_theta[[iter+1]] <- sig2_theta[[iter]]
    alpha_1v[[iter+1]] <- alpha_1v[[iter]]
    beta_jv[[iter+1]] <- beta_jv[[iter]]
    sig2_b_j[[iter+1]] <- sig2_b_j[[iter]]
    gamma_jv[[iter+1]] <- gamma_jv[[iter]]
    sig2_bar_a_j[[iter+1]] <- sig2_bar_a_j[[iter]]
    ksi_sij2[[iter+1]] <- ksi_sij2[[iter]]
    
    while (tau_M > tau_criteria) {
      sig2_theta_old <- sig2_theta[[iter+1]]
      alpha_1v_old <- alpha_1v[[iter+1]]
      b_j_old <- beta_jv[[iter+1]]
      sig2_b_j_old<-sig2_b_j[[iter+1]]
      bar_a_j_old <- gamma_jv[[iter+1]]
      sig2_bar_a_j_old<-sig2_bar_a_j[[iter+1]]
      ksi_sij2_old <- ksi_sij2[[iter+1]]
      
      alpha_num_temp <- list()
      alpha_deno_temp <- list()
      beta_in_deno_temp <- list()
      beta_in_num_temp <- list()
      beta_deno_temp <- list()
      #beta_num_temp <- list()
      gamma_in_num_temp <- list()
      gamma_in_deno_temp <- list()
      gamma_deno_temp <- list()
      gamma_num_temp <- list()
      sig2_theta_temp <- list()
      
      # M-step----
      for (s in 1:S) {
        # alpha_multi_temp_M <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
        # eta_temp_ij_M <- eta(ksi_sij2[[iter]][[s]])
        outer_temp1 <- outer(q_mu_theta_is[[iter+1]][[s]]^2 + q_sig2_theta_is[[iter+1]][[s]],
                             (q_mu_a_js[[iter+1]][,s]^2 + q_sig2_a_js[[iter+1]][,s]),"*")
        outer_temp2 <- outer(q_mu_theta_is[[iter+1]][[s]],
                             q_mu_a_js[[iter+1]][,s]*q_mu_b_js[[iter+1]][,s],'*')
        temp_1_M <- outer_temp1 - 2*outer_temp2
        temp_2_M <- q_mu_b_js[[iter+1]][,s]^2 + q_sig2_b_js[[iter+1]][,s]
        ksi_sij2[[iter+1]][[s]] <- sweep(temp_1_M,2,temp_2_M,"+")
        eta_temp_ij_M <- eta(ksi_sij2[[iter+1]][[s]])
        
        alpha_deno_temp[[s]] <- Ns[s]*(as.matrix(group_matrix[s,])%*%(group_matrix[s,]))
        
        alpha_num_temp[[s]] <- as.matrix(group_matrix[s,])*(sum(q_mu_theta_is[[iter+1]][[s]]))
        
        beta_deno_temp[[s]] <- (as.matrix(group_matrix_tilde[s,])%*%(group_matrix_tilde[s,]))
        
        deno_temp <- array(beta_deno_temp[[s]], dim = c((1+dim),(1+dim),J))
        
        beta_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M), FUN = "*")
        #colSums(2*eta_temp_ij_M)
        
        beta_in_num_temp[[s]] <- matrix(colSums(1/2 - Ysij[[s]] +
                                                  sweep((2*eta_temp_ij_M*q_mu_theta_is[[iter+1]][[s]]),
                                                        2,q_mu_a_js[[iter+1]][,s],"*")),ncol = 1)%*%group_matrix_tilde[s,]
        
        temp_3_M <- sweep(eta_temp_ij_M,2,q_mu_b_js[[iter+1]][,s],"*")
        temp_4_M <- Ysij[[s]]-1/2+2*temp_3_M
        temp_5_M <- temp_4_M*q_mu_theta_is[[iter+1]][[s]] #I*J * I
        gamma_in_num_temp[[s]] <- matrix(colSums(temp_5_M), ncol = 1)%*%group_matrix_tilde[s,]
        #J*1 * 1*M
        temp_6_M <- (q_mu_theta_is[[iter+1]][[s]])^2 + q_sig2_theta_is[[iter+1]][[s]]
        gamma_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M*temp_6_M), FUN = "*")
        #I*J * I vector
      }
      
      a_optimal <- optimizing(a_loss,
                              dat = list(J = length(out_sig2_a0), S = S, M = dim+1, lambda = lambda,
                                         q_sig2_a_js = q_sig2_a_js[[iter+1]][out_sig2_a0,,drop=F], 
                                         q_mu_a_js = q_mu_a_js[[iter+1]][out_sig2_a0,,drop=F],
                                         X_tilde_s = group_matrix_tilde,
                                         threshold = rho_Na2),
                              init = list(sig2_bar_a_j = array(sig2_bar_a_j[[iter+1]][out_sig2_a0]),
                                          gamma_jv = as.matrix(gamma_jv[[iter+1]][out_sig2_a0,,drop =F])),
                              as_vector = FALSE)
      
      gamma_jv[[iter+1]][out_sig2_a0,] <- a_optimal$par$gamma_jv
      
      all_deno_g <- Reduce('+',gamma_in_deno_temp)
      if(length(in_sig2_a0)!=0){
        for (inj in 1:length(in_sig2_a0)) {
          gamma_jv[[iter+1]][in_sig2_a0[inj],] <- as.vector(solve(all_deno_g[,,in_sig2_a0[inj]])%*%matrix((Reduce('+',gamma_in_num_temp))[in_sig2_a0[inj],],ncol = 1))
        }
      }
      # gamma_jv[[iter+1]][in_sig2_a0, ] <- 
      #   Reduce('+',gamma_in_num_temp)[in_sig2_a0]/Reduce('+',gamma_in_deno_temp)[in_sig2_a0]
      
      sig2_bar_a_j[[iter+1]][out_sig2_a0] <- a_optimal$par$sig2_bar_a_j
      sig2_bar_a_j[[iter+1]][in_sig2_a0] <- 0
      
      gamma_jv[[iter+1]][anchor,-1] <- 0
      
      sig2_theta[[iter+1]] <- 1#Reduce("+",sig2_theta_temp)/S
      
      alpha_1v[[iter+1]] <- as.vector(solve(Reduce('+',alpha_deno_temp))%*%Reduce('+',alpha_num_temp))
      
      beta_num_temp <- q_mu_b_js[[iter+1]]%*%group_matrix_tilde
      beta_jv[[iter+1]] <- t(solve(Reduce('+',beta_deno_temp))%*%t(beta_num_temp))#rowMeans(q_mu_b_js[[iter+1]])
      all_deno <- Reduce('+',beta_in_deno_temp)
      if(length(in_sig2_b0)!=0){
        for (inj in 1:length(in_sig2_b0)) {
          beta_jv[[iter+1]][in_sig2_b0[inj],] <- as.vector(solve(all_deno[,,in_sig2_b0[inj]])%*%matrix((Reduce('+',beta_in_num_temp))[in_sig2_b0[inj],],ncol = 1))
        }
      }
      #b_j[[iter+1]][in_sig2_b0] <- (Reduce('+',beta_in_num_temp))[in_sig2_b0]/(Reduce('+',beta_in_deno_temp))[in_sig2_b0]
      
      sig2_b_j[[iter+1]][out_sig2_b0] <- # sig2_b_j[[iter]]
        rowSums(q_sig2_b_js[[iter]][out_sig2_b0,,drop=F] +
                  ((beta_jv[[iter+1]]%*%t(group_matrix_tilde))[out_sig2_b0,] -
                     q_mu_b_js[[iter]][out_sig2_b0,,drop=F])^2)/(S+2*lambda)
      
      sig2_b_j[[iter+1]][in_sig2_b0] <- 0
      
      beta_jv[[iter+1]][anchor,-1] <- 0
      
      tau_M <- max(#abs(ksi_sij2[[iter+1]] - ksi_sij2_old),
        abs(sig2_theta[[iter+1]] - sig2_theta_old),
        abs(alpha_1v[[iter+1]] - alpha_1v_old),
        abs(beta_jv[[iter+1]] - b_j_old),
        abs(sig2_b_j[[iter+1]]-sig2_b_j_old),
        abs(gamma_jv[[iter+1]] - bar_a_j_old),
        abs(sig2_bar_a_j[[iter+1]]-sig2_bar_a_j_old))
    }
    ## update criteria
    tau <- max(abs(sig2_theta[[iter+1]] - sig2_theta[[iter]]),
               abs(alpha_1v[[iter+1]] - alpha_1v[[iter]]),
               abs(beta_jv[[iter+1]] - beta_jv[[iter]]),
               abs(sig2_b_j[[iter+1]]-sig2_b_j[[iter]]),
               abs(gamma_jv[[iter+1]] - gamma_jv[[iter]]),
               abs(sig2_bar_a_j[[iter+1]]-sig2_bar_a_j[[iter]]))
    tau_list[[iter]] <- tau
    iter <- iter+1
    #print(tau)
  }
  
  if(!debias){
    sig2_b_j[[iter]][sig2_b_j[[iter]]< rho_N] <- 0
    out_sig2_b0 <- which(sig2_b_j[[iter]]!=0)
    in_sig2_b0 <- which(sig2_b_j[[iter]]==0)
    
    sig2_bar_a_j[[iter]][sig2_bar_a_j[[iter]]< rho_Na] <- 0
    out_sig2_a0 <- which(sig2_bar_a_j[[iter]]!=0)
    in_sig2_a0 <- which(sig2_bar_a_j[[iter]]==0)
  }
  
  sig2_b_j2 <- sig2_b_j[[iter]]
  which_temp <- which((sig2_b_j2 < rho_N2)&(sig2_b_j2 != 0))
  sig2_b_j2[which_temp] <- rho_N2
  
  sig2_bar_a_j2 <- sig2_bar_a_j[[iter]]
  which_temp <- which((sig2_bar_a_j2 < rho_Na2)&(sig2_bar_a_j2 != 0))
  sig2_bar_a_j2[which_temp] <- rho_Na2
  
  temp1_like <- list()
  temp2_like <- list()
  temp3_like <- list()
  
  for (s in 1:S) {
    part1 <- -log(1+exp(-sqrt(ksi_sij2[[iter]][[s]]))) # IJ
    part2 <- (Ysij[[s]]-1/2)*
      sweep(outer(q_mu_theta_is[[iter]][[s]], # IJ
                  q_mu_a_js[[iter]][,s],"*"),2, q_mu_b_js[[iter]][,s],"-") - #IJ-J
      (1/2)*sqrt(ksi_sij2[[iter]][[s]]) #IJ
    
    part31 <- q_mu_b_js[[iter]][,s]^2 + q_sig2_b_js[[iter]][,s] # J
    part32 <- -2*outer(q_mu_theta_is[[iter]][[s]], q_mu_a_js[[iter]][,s]*q_mu_b_js[[iter]][,s],'*') + #IJ
      outer(q_mu_theta_is[[iter]][[s]]^2 + (q_sig2_theta_is[[iter]][[s]]),#IJ
            (q_mu_a_js[[iter]][,s]^2 + q_sig2_a_js[[iter]][,s]),"*") - ksi_sij2[[iter]][[s]] #IJ
    part3 <- eta(ksi_sij2[[iter]][[s]])*sweep(part32,2,part31, "+")
    temp1_like[[s]] <- sum(part1 + part2 - part3)
    
    part4 <-  (Ns[s]+length(out_sig2_a0)+length(out_sig2_b0))*log(2*pi)
    part5 <- sum(log(sig2_theta[[iter]]) + 
                   (q_sig2_theta_is[[iter]][[s]] + 
                      (q_mu_theta_is[[iter]][[s]] - as.numeric(alpha_1v[[iter]]%*%group_matrix[s,]))^2)/sig2_theta[[iter]])
    part6 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0])
    part7 <- sum(log(sig2_b_j2[out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j2[out_sig2_b0])
    
    part8 <- sum(log(sig2_bar_a_j[[iter]][out_sig2_a0]) + 
                   (q_sig2_a_js[[iter]][,s][out_sig2_a0] + 
                      (q_mu_a_js[[iter]][,s][out_sig2_a0] - 
                         (gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0])^2)/sig2_bar_a_j[[iter]][out_sig2_a0] + 
                   2*log(pnorm((gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0]/sqrt(sig2_bar_a_j[[iter]][out_sig2_a0]))))
    
    part9 <- sum(log(sig2_bar_a_j2[out_sig2_a0]) + 
                   (q_sig2_a_js[[iter]][,s][out_sig2_a0] + 
                      (q_mu_a_js[[iter]][,s][out_sig2_a0] - 
                         (gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0])^2)/sig2_bar_a_j2[out_sig2_a0] + 
                   2*log(pnorm((gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0]/sqrt(sig2_bar_a_j2[out_sig2_a0]))))
    
    temp2_like[[s]] <- part4+part5+part6+part8
    temp3_like[[s]] <- part4+part5+part7+part9
  }
  
  likelihood <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)
  likelihood2 <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)
  
  BIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*log(sum(Ns))
  BIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*log(sum(Ns))
  GIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*c*log(sum(Ns))*log(log(sum(Ns)))
  GIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*c*log(sum(Ns))*log(log(sum(Ns)))
  # EBIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
  #   length(out_sig2_b0)*log(sum(Ns)) + 2*g*log(choose(J,length(out_sig2_b0)))
  
  all_return <- list(sig2_theta = sig2_theta[[iter]],
                     alpha_1v = alpha_1v[[iter]],
                     b_j = beta_jv[[iter]],
                     sig2_b_j = sig2_b_j[[iter]],
                     bar_a_j = gamma_jv[[iter]],
                     sig2_bar_a_j = sig2_bar_a_j[[iter]],
                     likelihood = likelihood,
                     likelihood2 = likelihood2,
                     BIC = BIC,
                     BIC2 = BIC2,
                     GIC = GIC,
                     GIC2 = GIC2,
                     iter = iter,
                     tau = tau,
                     tau_list = tau_list)
  return(all_return)
}

## ---- GVEM ab random impact ----
GVEM_ab_multitheta <- function(resp, group_matrix, all_start, sig_start, lambda, c, iter_criteria = 5e2, tau_criteria = 1e-3,rho_N, rho_N2, rho_Na, rho_Na2, debias = F,anchor){
  
  a_loss <- stan_model(model_code = a_loss_code)
  
  ## Ysij: a list with S elements, each with J columns, Ns rows
  Ysij <- split(resp%>%select(-starts_with("group")), f = resp$groupfull)
  Ysij <- lapply(Ysij, as.matrix)
  Ns <- as.vector(unlist(lapply(Ysij, nrow)))
  S <- length(unique(resp$groupfull))
  J <- ncol(resp%>%select(-starts_with("group")))
  d <- ncol(group_matrix) # 2*5*5： 1+4+4=9 | d <- 9
  group_matrix_tilde <- cbind(1,group_matrix)
  dim <- ncol(group_matrix)
  ksi_sij2 <- list()
  sig2_theta <- list()
  alpha_1v <- list()
  beta_jv <- list()
  sig2_b_j <- list()
  gamma_jv <- list()
  sig2_bar_a_j <- list()
  sig2_alpha_0 <- list()
  
  q_sig2_alpha_0s <- list()
  q_mu_alpha_0s <- list()
  q_sig2_theta_is <- list()
  q_mu_theta_is <- list()
  q_sig2_b_js <- list()
  q_mu_b_js <- list()
  q_sig2_bar_a_js <- list()
  q_mu_bar_a_js <- list()
  q_sig2_a_js <- list()
  q_mu_a_js <- list()
  
  # start value----
  ksi_sij2[[1]] <- lapply(Ns, function(s) {matrix(0.1, nrow = s, ncol = J)})
  sig2_alpha_0[[1]] <- all_start$sig2_alpha_0
  sig2_theta[[1]] <- all_start$sig2_theta
  alpha_1v[[1]] <- all_start$alpha_1v
  beta_jv[[1]] <- all_start$beta_jv
  sig2_b_j[[1]] <- sig_start$sig2_b_j
  gamma_jv[[1]] <- all_start$gamma_jv
  sig2_bar_a_j[[1]] <- sig_start$sig2_bar_a_j
  
  q_sig2_alpha_0s[[1]] <- rep(1,S)
  q_mu_alpha_0s[[1]] <- rep(0,S)
  q_sig2_theta_is[[1]] <- list()
  q_mu_theta_is[[1]] <- list()
  for (s in 1:S) {
    q_sig2_theta_is[[1]][[s]] <- rep(1,Ns[s])
    q_mu_theta_is[[1]][[s]] <- rep(1,Ns[s])
  }
  q_sig2_b_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_b_js[[1]] <- matrix(0,nrow = J,ncol = S)
  q_sig2_bar_a_js[[1]] <- matrix(NA,nrow = J,ncol = S)
  q_mu_bar_a_js[[1]] <- matrix(NA,nrow = J,ncol = S)
  q_sig2_a_js[[1]] <- matrix(1,nrow = J,ncol = S)
  q_mu_a_js[[1]] <- matrix(1,nrow = J,ncol = S)
  
  iter <- 1
  tau <- 1
  # iter_criteria = 6e2
  # tau_criteria = 1e-3
  
  rescale <- list()
  out_sig2_b0 <- which(sig2_b_j[[1]]!=0)
  in_sig2_b0 <- which(sig2_b_j[[1]]==0)
  out_sig2_a0 <- which(sig2_bar_a_j[[1]]!=0)
  in_sig2_a0 <- which(sig2_bar_a_j[[1]]==0)
  sig2_b_j[[1]][out_sig2_b0] <- 1
  sig2_bar_a_j[[1]][out_sig2_a0] <- 1
  
  q_sig2_a_js[[1]][in_sig2_a0,] <- 0
  q_mu_a_js[[1]][in_sig2_a0,] <- (gamma_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_a0,]
  q_sig2_b_js[[1]][in_sig2_b0,] <- 0
  q_mu_b_js[[1]][in_sig2_b0,] <- (beta_jv[[1]]%*%t(group_matrix_tilde))[in_sig2_b0,]
  
  while((tau >= tau_criteria)&&(iter<=iter_criteria)){
    ## initialize for E-step
    q_sig2_alpha_0s[[iter+1]] <- rep(NA,S)
    q_mu_alpha_0s[[iter+1]] <- rep(NA,S)
    q_sig2_theta_is[[iter+1]] <- list()
    q_mu_theta_is[[iter+1]] <- list()
    q_sig2_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_b_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_sig2_bar_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_bar_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_sig2_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    q_mu_a_js[[iter+1]] <- matrix(NA,nrow = J,ncol = S)
    
    # E-step----
    for (s in 1:S) {
      eta_temp_ij <- eta(ksi_sij2[[iter]][[s]]) #I*J: eta ksij[[s]]
      alpha_multi_temp <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
      beta_multi_temp <- beta_jv[[iter]]%*%t(group_matrix_tilde)
      gamma_multi_temp <- gamma_jv[[iter]]%*%t(group_matrix_tilde)
      
      temp1 <- Ns[s]*sig2_alpha_0[[iter]]+sig2_theta[[iter]]
      q_sig2_alpha_0s[[iter+1]][s] <- (sig2_alpha_0[[iter]]*sig2_theta[[iter]])/temp1
      q_sig2_alpha_0s[[iter]][s] <- q_sig2_alpha_0s[[iter+1]][s]
      
      q_mu_alpha_0s[[iter+1]][s] <- (sig2_alpha_0[[iter]]/temp1)*
        (sum(q_mu_theta_is[[iter]][[s]]) -  Ns[s]*alpha_multi_temp)
      q_mu_alpha_0s[[iter]][s] <- q_mu_alpha_0s[[iter+1]][s]
      
      temp2 <- (1+2*sig2_theta[[iter]]*
                  as.vector(eta_temp_ij%*%(q_mu_a_js[[iter]][,s]^2+q_sig2_a_js[[iter]][,s])))
      q_sig2_theta_is[[iter+1]][[s]] <- sig2_theta[[iter]]/temp2
      q_sig2_theta_is[[iter]][[s]] <- q_sig2_theta_is[[iter+1]][[s]]
      
      temp3 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*") # I*J matrix
      temp4 <- as.vector((Ysij[[s]]-1/2+2*temp3)%*%(q_mu_a_js[[iter]][,s]))*sig2_theta[[iter]]+
        alpha_multi_temp + q_mu_alpha_0s[[iter]][s]
      q_mu_theta_is[[iter+1]][[s]] <- temp4/temp2
      q_mu_theta_is[[iter]][[s]] <- q_mu_theta_is[[iter+1]][[s]]
      
      temp9 <- eta_temp_ij*(q_mu_theta_is[[iter]][[s]]^2 + q_sig2_theta_is[[iter]][[s]])
      temp9 <- temp9[,out_sig2_a0,drop=F]
      temp10 <- 1+2*sig2_bar_a_j[[iter]][out_sig2_a0]*colSums(temp9)
      q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0] <- sig2_bar_a_j[[iter]][out_sig2_a0]/temp10
      q_sig2_bar_a_js[[iter]][,s][out_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]
      
      temp11 <- sweep(eta_temp_ij,2,q_mu_b_js[[iter]][,s],"*")
      temp11 <- temp11[,out_sig2_a0,drop=F]
      temp12 <- colSums((Ysij[[s]][,out_sig2_a0,drop=F] - 1/2 + 2*temp11)*q_mu_theta_is[[iter]][[s]])
      temp13 <- gamma_multi_temp[,s][out_sig2_a0] + temp12*sig2_bar_a_j[[iter]][out_sig2_a0]
      q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0] <- temp13/temp10
      q_mu_bar_a_js[[iter]][,s][out_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_sig2_bar_a_js[[iter+1]][,s][in_sig2_a0] <- 0
      q_sig2_bar_a_js[[iter]][,s][in_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][in_sig2_a0]
      q_mu_bar_a_js[[iter+1]][,s][in_sig2_a0] <- gamma_multi_temp[,s][in_sig2_a0]
      q_mu_bar_a_js[[iter]][,s][in_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][in_sig2_a0]
      
      bar2_ratio <- (q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0]^2)/q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]
      temp14 <- (sqrt(bar2_ratio)*exp(-(bar2_ratio/2)))/
        (sqrt(2*pi)*pnorm(sqrt(bar2_ratio)))
      temp15 <- exp(-bar2_ratio)/(2*pi*(pnorm(sqrt(bar2_ratio)))^2)
      q_sig2_a_js[[iter+1]][,s][out_sig2_a0] <- q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0]*
        (1-temp14-temp15)
      q_sig2_a_js[[iter]][,s][out_sig2_a0] <- q_sig2_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_mu_a_js[[iter+1]][,s][out_sig2_a0] <- q_mu_bar_a_js[[iter+1]][,s][out_sig2_a0] +
        (sqrt(q_sig2_bar_a_js[[iter+1]][,s][out_sig2_a0])*exp(-(bar2_ratio/2)))/(sqrt(2*pi)*pnorm(sqrt(bar2_ratio)))
      q_mu_a_js[[iter]][,s][out_sig2_a0] <- q_mu_a_js[[iter+1]][,s][out_sig2_a0]
      
      q_sig2_a_js[[iter+1]][,s][in_sig2_a0] <- 0
      q_sig2_a_js[[iter]][,s][in_sig2_a0] <- q_sig2_a_js[[iter+1]][,s][in_sig2_a0]
      q_mu_a_js[[iter+1]][,s][in_sig2_a0] <- gamma_multi_temp[,s][in_sig2_a0]
      q_mu_a_js[[iter]][,s][in_sig2_a0] <- q_mu_a_js[[iter+1]][,s][in_sig2_a0]
      
      temp5 <- 1+2*sig2_b_j[[iter]][out_sig2_b0]*colSums(eta_temp_ij[,out_sig2_b0,drop=F])
      q_sig2_b_js[[iter+1]][,s][out_sig2_b0] <- sig2_b_j[[iter]][out_sig2_b0]/temp5
      q_sig2_b_js[[iter]][,s][out_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][out_sig2_b0]
      
      temp6 <- sweep(eta_temp_ij,2,q_mu_a_js[[iter]][,s],"*") # I*J matrix
      temp6 <- temp6[,out_sig2_b0,drop=F]
      temp7 <- colSums(Ysij[[s]][,out_sig2_b0,drop=F] - 1/2 - 2*temp6*q_mu_theta_is[[iter]][[s]])
      temp8 <- beta_multi_temp[,s][out_sig2_b0] - temp7*sig2_b_j[[iter]][out_sig2_b0]
      q_mu_b_js[[iter+1]][,s][out_sig2_b0] <- temp8/temp5
      q_mu_b_js[[iter]][,s][out_sig2_b0] <- q_mu_b_js[[iter+1]][,s][out_sig2_b0]
      
      q_sig2_b_js[[iter+1]][,s][in_sig2_b0] <- 0
      q_sig2_b_js[[iter]][,s][in_sig2_b0] <- q_sig2_b_js[[iter+1]][,s][in_sig2_b0]
      q_mu_b_js[[iter+1]][,s][in_sig2_b0] <- beta_multi_temp[,s][in_sig2_b0]
      q_mu_b_js[[iter]][,s][in_sig2_b0] <- q_mu_b_js[[iter+1]][,s][in_sig2_b0]
      
    }
    
    tau_M <- 1
    ## initialize for M-step
    sig2_theta[[iter+1]] <- sig2_theta[[iter]]
    alpha_1v[[iter+1]] <- alpha_1v[[iter]]
    beta_jv[[iter+1]] <- beta_jv[[iter]]
    sig2_b_j[[iter+1]] <- sig2_b_j[[iter]]
    gamma_jv[[iter+1]] <- gamma_jv[[iter]]
    sig2_bar_a_j[[iter+1]] <- sig2_bar_a_j[[iter]]
    sig2_alpha_0[[iter+1]] <- sig2_alpha_0[[iter]]
    ksi_sij2[[iter+1]] <- ksi_sij2[[iter]]
    
    while (tau_M > tau_criteria) {
      sig2_theta_old <- sig2_theta[[iter+1]]
      alpha_1v_old <- alpha_1v[[iter+1]]
      b_j_old <- beta_jv[[iter+1]]
      sig2_b_j_old<-sig2_b_j[[iter+1]]
      bar_a_j_old <- gamma_jv[[iter+1]]
      sig2_bar_a_j_old<-sig2_bar_a_j[[iter+1]]
      sig2_alpha_0_old <- sig2_alpha_0[[iter+1]]
      ksi_sij2_old <- ksi_sij2[[iter+1]]
      
      alpha_num_temp <- list()
      alpha_deno_temp <- list()
      beta_in_deno_temp <- list()
      beta_in_num_temp <- list()
      beta_deno_temp <- list()
      #
      gamma_in_num_temp <- list()
      gamma_in_deno_temp <- list()
      gamma_deno_temp <- list()
      gamma_num_temp <- list()
      sig2_theta_temp <- list()
      
      # M-step----
      for (s in 1:S) {
        #alpha_multi_temp_M <- as.numeric(alpha_1v[[iter]]%*%group_matrix[s,])
        #eta_temp_ij_M <- eta(ksi_sij2[[iter+1]][[s]])
        outer_temp1 <- outer(q_mu_theta_is[[iter+1]][[s]]^2 + q_sig2_theta_is[[iter+1]][[s]],
                             (q_mu_a_js[[iter+1]][,s]^2 + q_sig2_a_js[[iter+1]][,s]),"*")
        outer_temp2 <- outer(q_mu_theta_is[[iter+1]][[s]],
                             q_mu_a_js[[iter+1]][,s]*q_mu_b_js[[iter+1]][,s],'*')
        temp_1_M <- outer_temp1 - 2*outer_temp2
        temp_2_M <- q_mu_b_js[[iter+1]][,s]^2 + q_sig2_b_js[[iter+1]][,s]
        ksi_sij2[[iter+1]][[s]] <- sweep(temp_1_M,2,temp_2_M,"+")
        #PROBLEM SOLVED !!!!!-----
        eta_temp_ij_M <- eta(ksi_sij2[[iter+1]][[s]])
        
        alpha_deno_temp[[s]] <- Ns[s]*(as.matrix(group_matrix[s,])%*%(group_matrix[s,]))
        
        alpha_num_temp[[s]] <- as.matrix(group_matrix[s,])*(sum(q_mu_theta_is[[iter+1]][[s]]) -
                                                              Ns[s]*q_mu_alpha_0s[[iter+1]][s])
        
        beta_deno_temp[[s]] <- (as.matrix(group_matrix_tilde[s,])%*%(group_matrix_tilde[s,]))
        
        deno_temp <- array(beta_deno_temp[[s]], dim = c((1+dim),(1+dim),J))
        
        beta_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M), FUN = "*")
        #beta_in_deno_temp[[s]] <- colSums(2*eta_temp_ij_M)
        
        beta_in_num_temp[[s]] <- matrix(colSums(1/2 - Ysij[[s]] +
                                                  sweep((2*eta_temp_ij_M*q_mu_theta_is[[iter+1]][[s]]),
                                                        2,q_mu_a_js[[iter+1]][,s],"*")),ncol = 1)%*%group_matrix_tilde[s,]
        
        temp_3_M <- sweep(eta_temp_ij_M,2,q_mu_b_js[[iter+1]][,s],"*")
        temp_4_M <- Ysij[[s]]-1/2+2*temp_3_M
        temp_5_M <- temp_4_M*q_mu_theta_is[[iter+1]][[s]] #I*J * I
        gamma_in_num_temp[[s]] <- matrix(colSums(temp_5_M), ncol = 1)%*%group_matrix_tilde[s,]
        #(colSums(temp_5_M))
        
        temp_6_M <- (q_mu_theta_is[[iter+1]][[s]])^2 + q_sig2_theta_is[[iter+1]][[s]]
        gamma_in_deno_temp[[s]] <- sweep(deno_temp, MARGIN = 3, STATS = colSums(2*eta_temp_ij_M*temp_6_M), FUN = "*")
        #(colSums(2*eta_temp_ij_M*temp_6_M)) #I*J * I vector
      }
      
      a_optimal <- optimizing(a_loss,
                              dat = list(J = length(out_sig2_a0), S = S, M = dim+1, lambda = lambda,
                                         q_sig2_a_js = q_sig2_a_js[[iter+1]][out_sig2_a0,,drop=F],
                                         q_mu_a_js = q_mu_a_js[[iter+1]][out_sig2_a0,,drop=F],
                                         X_tilde_s = group_matrix_tilde),
                              init = list(sig2_bar_a_j = as.array(sig2_bar_a_j[[iter+1]][out_sig2_a0]),
                                          gamma_jv = as.matrix(gamma_jv[[iter+1]][out_sig2_a0,,drop =F])),
                              as_vector = FALSE)
      
      gamma_jv[[iter+1]][out_sig2_a0,] <- a_optimal$par$gamma_jv
      
      all_deno_g <- Reduce('+',gamma_in_deno_temp)
      if(length(in_sig2_a0)!=0){
        for (inj in 1:length(in_sig2_a0)) {
          gamma_jv[[iter+1]][in_sig2_a0[inj],] <- as.vector(solve(all_deno_g[,,in_sig2_a0[inj]])%*%matrix((Reduce('+',gamma_in_num_temp))[in_sig2_a0[inj],],ncol = 1))
        }
      }
      # bar_a_j[[iter+1]][out_sig2_a0] <- a_optimal$par[(length(out_sig2_a0)+1):(2*length(out_sig2_a0))]
      # bar_a_j[[iter+1]][in_sig2_a0] <- Reduce('+',gamma_in_num_temp)[in_sig2_a0]/Reduce('+',gamma_in_deno_temp)[in_sig2_a0]
      # bar_a_j[[iter+1]][in_sig2_a0][which(bar_a_j[[iter+1]][in_sig2_a0]<=0)] <- 1e-3
      
      sig2_bar_a_j[[iter+1]][out_sig2_a0] <- a_optimal$par$sig2_bar_a_j
      sig2_bar_a_j[[iter+1]][in_sig2_a0] <- 0
      
      gamma_jv[[iter+1]][anchor,-1] <- 0
      
      sig2_theta[[iter+1]] <- 1#Reduce("+",sig2_theta_temp)/sum(Ns)
      
      sig2_alpha_0[[iter+1]] <- mean(q_sig2_alpha_0s[[iter+1]]+q_mu_alpha_0s[[iter+1]]^2)
      
      alpha_1v[[iter+1]] <- as.vector(solve(Reduce('+',alpha_deno_temp))%*%Reduce('+',alpha_num_temp))
      
      beta_num_temp <- q_mu_b_js[[iter+1]]%*%group_matrix_tilde
      beta_jv[[iter+1]] <- t(solve(Reduce('+',beta_deno_temp))%*%t(beta_num_temp))#rowMeans(q_mu_b_js[[iter+1]])
      all_deno <- Reduce('+',beta_in_deno_temp)
      if(length(in_sig2_b0)!=0){
        for (inj in 1:length(in_sig2_b0)) {
          beta_jv[[iter+1]][in_sig2_b0[inj],] <- as.vector(solve(all_deno[,,in_sig2_b0[inj]])%*%matrix((Reduce('+',beta_in_num_temp))[in_sig2_b0[inj],],ncol = 1))
        }
      }
      # b_j[[iter+1]] <- rowMeans(q_mu_b_js[[iter+1]])
      # b_j[[iter+1]][in_sig2_b0] <- (Reduce('+',beta_in_num_temp))[in_sig2_b0]/(Reduce('+',beta_in_deno_temp))[in_sig2_b0]
      
      sig2_b_j[[iter+1]][out_sig2_b0] <- # sig2_b_j[[iter]]
        rowSums(q_sig2_b_js[[iter]][out_sig2_b0,,drop=F] +
                  ((beta_jv[[iter+1]]%*%t(group_matrix_tilde))[out_sig2_b0] -
                     q_mu_b_js[[iter]][out_sig2_b0,,drop=F])^2)/(S+2*lambda)
      
      sig2_b_j[[iter+1]][in_sig2_b0] <- 0
      
      beta_jv[[iter+1]][anchor,-1] <- 0
      
      tau_M <- max(#abs(ksi_sij2[[iter+1]] - ksi_sij2_old),
        abs(sig2_theta[[iter+1]] - sig2_theta_old),
        abs(alpha_1v[[iter+1]] - alpha_1v_old),
        abs(beta_jv[[iter+1]] - b_j_old),
        abs(sig2_b_j[[iter+1]]-sig2_b_j_old),
        abs(gamma_jv[[iter+1]] - bar_a_j_old),
        abs(sig2_bar_a_j[[iter+1]]-sig2_bar_a_j_old),
        abs(sig2_alpha_0[[iter+1]] - sig2_alpha_0_old))
      #print(tau_M)
    }
    ## update criteria
    tau <- max(abs(sig2_theta[[iter+1]] - sig2_theta[[iter]]),
               abs(alpha_1v[[iter+1]] - alpha_1v[[iter]]),
               abs(beta_jv[[iter+1]] - beta_jv[[iter]]),
               abs(sig2_b_j[[iter+1]]-sig2_b_j[[iter]]),
               abs(gamma_jv[[iter+1]] - gamma_jv[[iter]]),
               abs(sig2_bar_a_j[[iter+1]]-sig2_bar_a_j[[iter]]),
               abs(sig2_alpha_0[[iter+1]] - sig2_alpha_0[[iter]]))
    
    iter <- iter+1
    #print(tau)
  }
  
  if(!debias){
    sig2_b_j[[iter]][sig2_b_j[[iter]]< rho_N] <- 0
    out_sig2_b0 <- which(sig2_b_j[[iter]]!=0)
    in_sig2_b0 <- which(sig2_b_j[[iter]]==0)
    
    sig2_bar_a_j[[iter]][sig2_bar_a_j[[iter]]< rho_Na] <- 0
    out_sig2_a0 <- which(sig2_bar_a_j[[iter]]!=0)
    in_sig2_a0 <- which(sig2_bar_a_j[[iter]]==0)
  }
  
  sig2_b_j2 <- sig2_b_j[[iter]]
  which_temp <- which((sig2_b_j2 < rho_N2)&(sig2_b_j2 != 0))
  sig2_b_j2[which_temp] <- rho_N2
  
  sig2_bar_a_j2 <- sig2_bar_a_j[[iter]]
  which_temp <- which((sig2_bar_a_j2 < rho_Na2)&(sig2_bar_a_j2 != 0))
  sig2_bar_a_j2[which_temp] <- rho_Na2
  
  temp1_like <- list()
  temp2_like <- list()
  temp3_like <- list()
  
  for (s in 1:S) {
    part1 <- -log(1+exp(-sqrt(ksi_sij2[[iter]][[s]]))) # IJ
    part2 <- (Ysij[[s]]-1/2)*
      sweep(outer(q_mu_theta_is[[iter]][[s]], # IJ
                  q_mu_a_js[[iter]][,s],"*"),2, q_mu_b_js[[iter]][,s],"-") - #IJ-J
      (1/2)*sqrt(ksi_sij2[[iter]][[s]]) #IJ
    
    part31 <- q_mu_b_js[[iter]][,s]^2 + q_sig2_b_js[[iter]][,s] # J
    part32 <- -2*outer(q_mu_theta_is[[iter]][[s]], q_mu_a_js[[iter]][,s]*q_mu_b_js[[iter]][,s],'*') + #IJ
      outer(q_mu_theta_is[[iter]][[s]]^2 + (q_sig2_theta_is[[iter]][[s]]),#IJ
            (q_mu_a_js[[iter]][,s]^2 + q_sig2_a_js[[iter]][,s]),"*") - ksi_sij2[[iter]][[s]] #IJ
    part3 <- eta(ksi_sij2[[iter]][[s]])*sweep(part32,2,part31, "+")
    temp1_like[[s]] <- sum(part1 + part2 - part3)
    
    part4 <-  (Ns[s]+length(out_sig2_a0)+length(out_sig2_b0)+1)*log(2*pi) + log(sig2_alpha_0[[iter]]) +
      (q_sig2_alpha_0s[[iter]][[s]] + q_mu_alpha_0s[[iter]][[s]]^2)/sig2_alpha_0[[iter]]
    
    part5 <- sum(log(sig2_theta[[iter]]) + 
                   (q_sig2_theta_is[[iter]][[s]] + 
                      (q_mu_theta_is[[iter]][[s]] - as.numeric(alpha_1v[[iter]]%*%group_matrix[s,]))^2)/sig2_theta[[iter]])
    part6 <- sum(log(sig2_b_j[[iter]][out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j[[iter]][out_sig2_b0])
    part7 <- sum(log(sig2_b_j2[out_sig2_b0]) + 
                   (q_sig2_b_js[[iter]][,s][out_sig2_b0] + 
                      (q_mu_b_js[[iter]][,s][out_sig2_b0] - 
                         (beta_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_b0])^2)/sig2_b_j2[out_sig2_b0])
    
    part8 <- sum(log(sig2_bar_a_j[[iter]][out_sig2_a0]) + 
                   (q_sig2_a_js[[iter]][,s][out_sig2_a0] + 
                      (q_mu_a_js[[iter]][,s][out_sig2_a0] - 
                         (gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0])^2)/sig2_bar_a_j[[iter]][out_sig2_a0] + 
                   2*log(pnorm((gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0]/sqrt(sig2_bar_a_j[[iter]][out_sig2_a0]))))
    
    part9 <- sum(log(sig2_bar_a_j2[out_sig2_a0]) + 
                   (q_sig2_a_js[[iter]][,s][out_sig2_a0] + 
                      (q_mu_a_js[[iter]][,s][out_sig2_a0] - 
                         (gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0])^2)/sig2_bar_a_j2[out_sig2_a0] + 
                   2*log(pnorm((gamma_jv[[iter]]%*%t(group_matrix_tilde))[,s][out_sig2_a0]/sqrt(sig2_bar_a_j2[out_sig2_a0]))))
    
    temp2_like[[s]] <- part4+part5+part6+part8
    temp3_like[[s]] <- part4+part5+part7+part9
  }
  
  likelihood <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)
  likelihood2 <- Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)
  
  BIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*log(sum(Ns))
  BIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*log(sum(Ns))
  GIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*c*log(sum(Ns))*log(log(sum(Ns)))
  GIC2 <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp3_like)) +
    (length(out_sig2_a0)+length(out_sig2_b0))*c*log(sum(Ns))*log(log(sum(Ns)))
  # EBIC <- -2*(Reduce("+",temp1_like) - (1/2)*Reduce("+",temp2_like)) +
  #   length(out_sig2_b0)*log(sum(Ns)) + 2*g*log(choose(J,length(out_sig2_b0)))
  
  all_return <- list(sig2_theta = sig2_theta[[iter]],
                     alpha_1v = alpha_1v[[iter]],
                     b_j = beta_jv[[iter]],
                     sig2_b_j = sig2_b_j[[iter]],
                     bar_a_j = gamma_jv[[iter]],
                     sig2_bar_a_j = sig2_bar_a_j[[iter]],
                     sig2_alpha_0 = sig2_alpha_0[[iter]],
                     likelihood = likelihood,
                     likelihood2 = likelihood2,
                     BIC = BIC,
                     BIC2 = BIC2,
                     GIC = GIC,
                     GIC2 = GIC2,
                     iter = iter,
                     tau = tau)
  return(all_return)
}

# ---- MCMC code ----
## ---- MCMC bonly fixed impact ----
MCMC_b_code <- "
data {
  int<lower=0> n_student;             
  int<lower=0> n_item;               
  int<lower=0> n_group;               
  int<lower=1, upper=n_group> ind_g[n_student]; 
  int<lower=0,upper=1> Y[n_student, n_item];     
  int<lower=0> M;
  matrix[n_group, M] X_s;   
  int<lower=0> n_b_dif;             
  int<lower=1, upper=n_item> b_dif_ind[n_b_dif];
}

transformed data {
  matrix[n_group, M + 1] X_tilde_s = append_col(rep_vector(1, n_group), X_s);
  
  int<lower=0, upper=1> has_b_dif[n_item] = rep_array(0, n_item);
  int<lower=0, upper=n_b_dif> b_dif_pos[n_item] = rep_array(0, n_item);
  
  for (k in 1:n_b_dif) {
    int j = b_dif_ind[k];
    has_b_dif[j] = 1;
    b_dif_pos[j] = k; 
  }
}

parameters {
  vector[M] alpha1;                    
  vector[n_student] theta_raw;         
  
  vector[n_item] log_a;                
  
  matrix[n_item, M + 1] beta_j_coef;   
  vector<lower=0>[n_b_dif] sigma_b_j;  
  matrix[n_b_dif, n_group] b_js_raw;  
}

transformed parameters {
  vector[n_student] theta;
  for (i in 1:n_student) {
    int s = ind_g[i];
    theta[i] = X_s[s, ] * alpha1 + theta_raw[i];
  }
  
  matrix[n_item, n_group] a_js;
  for (j in 1:n_item) {
    for (s in 1:n_group) {
      a_js[j, s] = exp(log_a[j]);
    }
  }
  
  matrix[n_item, n_group] b_js;
  for (j in 1:n_item) {
    for (s in 1:n_group) {
      real mu_b = X_tilde_s[s, ] * beta_j_coef[j, ]';
      if (has_b_dif[j] == 1) {
        int k = b_dif_pos[j];
        b_js[j, s] = mu_b + b_js_raw[k, s] * sigma_b_j[k];
      } else {
        b_js[j, s] = mu_b;
      }
    }
  }
}

model {
  // Priors for abilities
  alpha1 ~ normal(0, 1);
  theta_raw ~ normal(0, 1);

  // Prior for item discrimination
  log_a ~ normal(0, 1);

  // Priors for difficulty
  to_vector(beta_j_coef) ~ normal(0, 1);
  sigma_b_j ~ gamma(1.5, 2);
  to_vector(b_js_raw) ~ normal(0, 1);

  // Likelihood
  for (i in 1:n_student) {
    int s = ind_g[i];
    vector[n_item] eta;
    for (j in 1:n_item) {
      eta[j] = a_js[j, s] * (theta[i] - b_js[j, s]);
    }
    Y[i, ] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  vector[n_student] log_lik;
  for (i in 1:n_student) {
    real lp_i = 0;
    int s = ind_g[i];
    for (j in 1:n_item) {
      real eta_ij = a_js[j, s] * (theta[i] - b_js[j, s]);
      lp_i += bernoulli_logit_lpmf(Y[i, j] | eta_ij);
    }
    log_lik[i] = lp_i;
  }
}
"

## ---- MCMC bonly random impact ----
MCMC_b_multilevel_code <- "
data {
  int<lower=0> n_student;              
  int<lower=0> n_item;       
  int<lower=0> n_group;        
  int<lower=1, upper=n_group> ind_g[n_student]; 
  int<lower=0,upper=1> Y[n_student, n_item];  
  int<lower=0> M;                    
  matrix[n_group, M] X_s;            
  
  int<lower=0> n_b_dif;               
  int<lower=1, upper=n_item> b_dif_ind[n_b_dif];
}

transformed data {
  matrix[n_group, M + 1] X_tilde_s = append_col(rep_vector(1, n_group), X_s);
  
  int<lower=0, upper=1> has_b_dif[n_item] = rep_array(0, n_item);
  int<lower=0, upper=n_b_dif> b_dif_pos[n_item] = rep_array(0, n_item);
  
  for (k in 1:n_b_dif) {
    int j = b_dif_ind[k];
    has_b_dif[j] = 1;
    b_dif_pos[j] = k; 
  }
}

parameters {
  vector[n_group] alpha0_s_raw;        
  real<lower=0> sigma_alpha0;          
  vector[M] alpha1;                   
  vector[n_student] theta_raw;        
  
  vector[n_item] log_a;               
  
  matrix[n_item, M + 1] beta_j_coef;   
  vector<lower=0>[n_b_dif] sigma_b_j;  
  matrix[n_b_dif, n_group] b_js_raw;
}

transformed parameters {
  vector[n_group] alpha0_s = alpha0_s_raw * sigma_alpha0;
  vector[n_student] theta;
  for (i in 1:n_student) {
    int s = ind_g[i];
    theta[i] = alpha0_s[s] + X_s[s, ] * alpha1 + theta_raw[i];
  }
  
  matrix[n_item, n_group] a_js;
  for (j in 1:n_item) {
    for (s in 1:n_group) {
      a_js[j, s] = exp(log_a[j]);
    }
  }
  
  matrix[n_item, n_group] b_js;
  for (j in 1:n_item) {
    for (s in 1:n_group) {
      real mu_b = X_tilde_s[s, ] * beta_j_coef[j, ]';
      if (has_b_dif[j] == 1) {
        int k = b_dif_pos[j];
        b_js[j, s] = mu_b + b_js_raw[k, s] * sigma_b_j[k];
      } else {
        b_js[j, s] = mu_b;
      }
    }
  }
}

model {
  alpha0_s_raw ~ normal(0, 1);
  sigma_alpha0 ~ gamma(1.5, 2);
  alpha1 ~ normal(0, 1);
  theta_raw ~ normal(0, 1);

  log_a ~ normal(0, 1);           

  to_vector(beta_j_coef) ~ normal(0, 1);
  sigma_b_j ~ gamma(1.5, 2);
  to_vector(b_js_raw) ~ normal(0, 1);

  for (i in 1:n_student) {
    int s = ind_g[i];
    vector[n_item] eta;
    for (j in 1:n_item) {
      eta[j] = a_js[j, s] * (theta[i] - b_js[j, s]);
    }
    Y[i, ] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  vector[n_student] log_lik;
  for (i in 1:n_student) {
    real lp_i = 0;
    int s = ind_g[i];
    for (j in 1:n_item) {
      real eta_ij = a_js[j, s] * (theta[i] - b_js[j, s]);
      lp_i += bernoulli_logit_lpmf(Y[i, j] | eta_ij);
    }
    log_lik[i] = lp_i;
  }
}
"

# ---- umbrella function ----
D2PL_rnd_gvem <- function(resp, group_matrix, all_lambda = seq(4,10,length.out = 20), cvalue = 0.04,
                          dif_type = "udif", impact_type = "fixed", anchor,
                          MCMC = FALSE, MCMC_iter = 1e4,
                          iter_criteria = 5e2, tau_criteria = 1e-3,
                          rho_N = 1e-3, rho_N2 = 0.1,
                          rho_Na = 1e-3, rho_Na2 = 0.1){
  
  dim <- ncol(group_matrix)
  group <- as.factor(resp$groupfull)
  J <- ncol(resp) - 1 
  invariance <- c('slopes','intercepts','free_means')
  S <- length(unique(resp$groupfull))
  
  c <- cvalue*S
  
  if(dif_type == "udif"){

    if(impact_type == "fixed"){
      
      all_start <- start_value(resp = resp, J = J, dim = dim, group = group, invariance = invariance,
                               dif_type = dif_type, impact_type = impact_type)
      
      GVEM_returns <- lapply(all_lambda, function(lambda) {
        GVEM_bonly_nomultitheta(resp = resp, group_matrix = group_matrix,
                                all_start = all_start,
                                sig_start = all_start,
                                lambda = lambda, c = c, iter_criteria = iter_criteria, 
                                tau_criteria = tau_criteria,
                                rho_N = rho_N, rho_N2 = rho_N2, anchor = anchor)})
      
      GVEM_returns_debias <- lapply(GVEM_returns, function(x){
        GVEM_bonly_nomultitheta(resp = resp, group_matrix = group_matrix,
                                all_start = all_start,
                                sig_start = x,
                                lambda = 0, c, iter_criteria = iter_criteria,
                                tau_criteria = tau_criteria,
                                rho_N = rho_N, rho_N2 = rho_N2, debias = T, anchor = anchor)})
      
      GIC_all <- lapply(GVEM_returns_debias, function(x) x[["GIC2"]])
      GVEM_min_GIC2 <- GVEM_returns_debias[[which.min(GIC_all)]]
      sig2_b <- GVEM_min_GIC2$sig2_b_j
      
      posterior <- NULL
      if(MCMC){
        library(rstan)
        options(mc.cores=parallel::detectCores()-6)
        rstan_options(auto_write = T)
        
        response <- resp[,-ncol(resp)]
        
        mod_rnd_b <- stan_model(model_code = MCMC_b_code)
        data_rnd_b  <- list(n_student = nrow(resp),
                            n_item = ncol(resp) - 1,
                            n_group = nrow(group_matrix),
                            ind_g = resp$groupfull,
                            Y = response,
                            M = ncol(group_matrix),
                            X_s = group_matrix,
                            n_b_dif = sum(sig2_b!=0),
                            b_dif_ind = which(sig2_b!=0))
        
        fit_rnd_b <- sampling(
          object = mod_rnd_b,
          data   = data_rnd_b,
          chains = 4,
          iter   = MCMC_iter,
          seed   = 123
        )
        
        posterior <- extract(fit_rnd_b)      # draws as lists/arrays
        
        # MCMC_sig2_b_j <- colMeans(posterior$sigma_b_j)
        # MCMC_beta_j <- colMeans(posterior$beta_j_coef)
      }
      
    }else if(impact_type == "random"){
      
      all_start <- start_value(resp = resp, J = J, dim = dim, group = group, invariance = invariance,
                               dif_type = dif_type, impact_type = impact_type)
      
      GVEM_returns <- lapply(all_lambda, function(lambda) {
        GVEM_bonly_multitheta(resp, group_matrix,
                              all_start = all_start,
                              sig_start = all_start,
                              lambda = lambda, c = c, rho_N = rho_N, rho_N2 = rho_N2, anchor = anchor)})
      
      GVEM_returns_debias <- lapply(GVEM_returns, function(x){
        GVEM_bonly_multitheta(resp = resp, group_matrix = group_matrix, 
                              all_start = all_start,
                              sig_start = x,
                              lambda = 0, rho_N = rho_N, rho_N2 = rho_N2, c = c, debias=T, anchor = anchor)})
      
      GIC_all <- lapply(GVEM_returns_debias, function(x) x[["GIC2"]])
      GVEM_min_GIC2 <- GVEM_returns_debias[[which.min(GIC_all)]]
      sig2_b <- GVEM_min_GIC2$sig2_b_j
      
      posterior <- NULL
      if(MCMC){
        library(rstan)
        options(mc.cores=parallel::detectCores()-6)
        rstan_options(auto_write = T)
        
        response <- resp[,-ncol(resp)]
        
        mod_rnd_b <- stan_model(model_code = MCMC_b_multilevel_code)
        data_rnd_b  <- list(n_student = nrow(resp),
                            n_item = ncol(resp) - 1,
                            n_group = nrow(group_matrix),
                            ind_g = resp$groupfull,
                            Y = response,
                            M = ncol(group_matrix),
                            X_s = group_matrix,
                            n_b_dif = sum(sig2_b!=0),
                            b_dif_ind = which(sig2_b!=0))
        
        fit_rnd_b <- sampling(
          object = mod_rnd_b,
          data   = data_rnd_b,
          chains = 4,
          iter   = MCMC_iter,
          seed   = 123
        )
        
        posterior <- extract(fit_rnd_b)      # draws as lists/arrays
      }
    }
  }else if(dif_type == "nudif"){

    if(impact_type == "fixed"){
      
      all_start <- start_value(resp = resp, J = J, dim = dim, group = group, invariance = invariance,
                               dif_type = dif_type, impact_type = impact_type)
      
      GVEM_returns <- lapply(all_lambda, function(lambda) {
        GVEM_ab_nomultitheta(resp = resp, group_matrix = group_matrix,
                             all_start = all_start,
                             sig_start = all_start,
                             lambda = lambda, c = c, iter_criteria = iter_criteria, 
                             tau_criteria = tau_criteria,
                             rho_N = rho_N, rho_N2 = rho_N2, rho_Na = rho_Na, rho_Na2 = rho_Na2, anchor = anchor)})
      
      GVEM_returns_debias <- lapply(GVEM_returns, function(x){
        GVEM_ab_nomultitheta(resp = resp, group_matrix = group_matrix,
                             all_start = all_start,
                             sig_start = x,
                             lambda = 0, c = c, iter_criteria = iter_criteria, 
                             tau_criteria = tau_criteria,
                             rho_N = rho_N, rho_N2 = rho_N2, rho_Na = rho_Na, rho_Na2 = rho_Na2, debias = T, anchor = anchor)})
      
      GIC_all <- lapply(GVEM_returns_debias, function(x) x[["GIC2"]])
      GVEM_min_GIC2 <- GVEM_returns_debias[[which.min(GIC_all)]]
      posterior <- NULL
      
    }else if(impact_type == "random"){
    
      all_start <- start_value(resp = resp, J = J, dim = dim, group = group, invariance = invariance,
                               dif_type = dif_type, impact_type = impact_type)
      
      GVEM_returns <- lapply(all_lambda, function(lambda) {
        GVEM_ab_multitheta(resp = resp, group_matrix = group_matrix,
                           all_start = all_start,
                           sig_start = all_start,
                           lambda = lambda, c = c, iter_criteria = iter_criteria, 
                           tau_criteria = tau_criteria,
                           rho_N = rho_N, rho_N2 = rho_N2, rho_Na = rho_Na, rho_Na2 = rho_Na2, anchor = anchor)})
      
      GVEM_returns_debias <- lapply(GVEM_returns, function(x){
        GVEM_ab_multitheta(resp = resp, group_matrix = group_matrix,
                           all_start = all_start,
                           sig_start = x,
                           lambda = 0, c = c, iter_criteria = iter_criteria, 
                           tau_criteria = tau_criteria,
                           rho_N = rho_N, rho_N2 = rho_N2, rho_Na = rho_Na, rho_Na2 = rho_Na2, debias=T, anchor = anchor)})
      GIC_all <- lapply(GVEM_returns_debias, function(x) x[["GIC2"]])
      GVEM_min_GIC2 <- GVEM_returns_debias[[which.min(GIC_all)]]
      posterior <- NULL
    }
  }
  return(list(gvem_all_lambda = GVEM_returns_debias,
              gvem_optimal_lambda = GVEM_min_GIC2,
              MCMC_est = posterior))
}

# ---- example use ----
# dif_type <- c("udif", "nudif")
# impact_type <- c("fixed","random")
# MCMC <- c(TRUE, FALSE)
# Users need to specify anchor items in the function using anchor = XX, XX is the number of items
# Note that the anchor is only on main effect (i.e., traditional DIF) but not on intersectional DIF

## ---- when dif_type == "udif" ----
# group_matrix <- readxl::read_excel("D2PL_rnd_groupmat.xlsx")
# write.csv(group_matrix,"D2PL_rnd_groupmat.csv",row.names = FALSE )
# group_matrix <- as.matrix(group_matrix)
# 
# resp <- readxl::read_excel("D2PL_rnd_resp_udif.xlsx")
# write.csv(resp,"D2PL_rnd_resp_udif.csv",row.names = FALSE )
# 
# fit_example_mcmc <- D2PL_rnd_gvem (resp, group_matrix, all_lambda = seq(1,10,length.out = 10), cvalue = 0.04,
#                               dif_type = "udif", impact_type = "fixed", anchor = c(1,2,3),
#                               MCMC = TRUE, MCMC_iter = 1e2,
#                               iter_criteria = 5e2, tau_criteria = 1e-3,
#                               rho_N = 1e-3, rho_N2 = 0.1,
#                               rho_Na = 1e-3, rho_Na2 = 0.1)
# ## the sig2_b_j is the item random effect on intercept
# ## not zero means this item has DIF, as indicated by 1 in the output below
# sign(fit_example$gvem_optimal_lambda$sig2_b_j)
# 
# ### ----we allow MCMC estimation under udif----
# ### with MCMC = TRUE we have more accurate estimations 
# ### NOTE the MCMC_iter = 1e2 in the example use is just for quick estimation
# ### fit_example$MCMC_est contains all the samplings
# ### USE THE FOLLOWING ITEM ESTIAMTIONS TO DRAW THE ICC
# ### all slope estimations, one estimation for each item, no variations across groups
# colMeans(fit_example$MCMC_est$a_js)[,1]
# ### all intercept estimations, each row is an item, each column is a group
# colMeans(fit_example$MCMC_est$b_js)
# 
# 
# ## ---- when dif_type == "nudif" ----
# ## the sig2_b_j is the item random effect on intercept
# ## the sig2_bar_a_j is the item random effect on intercept
# ## not zero means this item has DIF, as indicated by 1 in the output below
# group_matrix <- read.csv("D2PL_rnd_groupmat.csv")
# group_matrix <- as.matrix(group_matrix)
# resp <- read.csv("D2PL_rnd_resp_nudif.csv")
# write.csv(resp, "D2PL_rnd_resp_nudif.csv",row.names = FALSE )
# 
# fit_example_n <- D2PL_rnd_gvem (resp, group_matrix, all_lambda = seq(1,10,length.out = 10), cvalue = 0.04,
#                               dif_type = "nudif", impact_type = "fixed", anchor = 19:20,
#                               MCMC = FALSE, MCMC_iter = 1e2,
#                               iter_criteria = 5e2, tau_criteria = 1e-3,
#                               rho_N = 1e-3, rho_N2 = 0.1,
#                               rho_Na = 1e-3, rho_Na2 = 0.1)
# sign(fit_example$gvem_optimal_lambda$sig2_b_j)
# sign(fit_example$gvem_optimal_lambda$sig2_bar_a_j)

## we DO NOT allow MCMC estimation under nudif----
### SO NO ICC for nudif
### JUST flag which item has DIF on intercept and which on slope