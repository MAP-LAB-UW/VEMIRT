#' GVEM Algorithm for the Generalized Partial Credit Model
#'
#' @param data An \eqn{N\times J} matrix of item responses where 0 is the minimal partial credit score (missing responses should be coded as \code{NA}) 
#' @param model A \eqn{J\times K} matrix of loading indicators (K is the Number of latent dimension)(all items load on the only dimension by default)
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param SE Whether to calculate the standard errors
#' @param verbose Whether to show the progress
#' @param EFA Whether to rotate the output
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:
#'   \item{ ...$Sigma}{Group-level covariance matrices}#'   
#'   \item{ ...$MU}{Person-level posterior mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$ll}{Estimated lower bound of log-likelihood}
#'
#' @author Yijun Cheng <chengxb@uw.edu>
#' @export
#'
#' @examples
#' \dontrun{
#' with(MGPCM_data, MGPCM_gvem(data, model))}
#' 
MGPCM_gvem <- function(data, model = matrix(1, nrow = J, ncol = 4), group = rep(1, nrow(data)), iter = 2000, eps = 1e-5, SE = FALSE, verbose = TRUE, EFA = FALSE) {
  fn <- list(init.MGPCM_gvem, est.MGPCM_gvem)
  Y <- as.matrix(data)
  D <- as.matrix(model)
  N <- nrow(Y)
  group <- rep(1, nrow(data))
  X <- group
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  if (min(X) == 0)
    X <- X + 1
  output <- if (verbose) stderr() else nullfile()
  
  cat(file = output, 'Fitting the model for initial values...\n')
  init <- fn[[1]](Y, D, X, iter, eps, SE, EFA)
  cat(file = output, 'Preparing the output...\n')
  if(SE) cat(file = output, 'Calculating SE, please be patient...\n')
  result <- fn[[2]](init())
  
  if(verbose) {
    cat("\n")
    output <- capture.output(lst(a = result$a, b = result$b))
    cat(output, sep = "\n")
  }
  invisible(result)
}

library(torch)
library(tibble)
init.MGPCM_gvem <- function(Y, D, X, iter, eps, SE, EFA, ...) {
  
  N <- nrow(Y)
  n <- tapply(1:N, X, identity)
  J <- ncol(Y)
  K <- ncol(D)
  G <- max(X)
  Y.na <- is.na(Y)
  Y[Y.na] <- 0.5
  Y <- torch_tensor(Y)
  Y.mask <- torch_tensor(!Y.na)
  p <- groupmean(n, Y$masked_fill(!Y.mask, NaN)$float()$view(c(N, J)))
  
  P_list <- apply(Y, 2, max) #k-list
  P <- max(P_list) #
  P_range_max <- torch_arange(0, P, dtype = torch_int())  # Shape [P_max]
  P_range_matrix <- P_range_max$unsqueeze(1)$expand(c(J, P+1))$unsqueeze(1)[X]
  Y_matrix <- Y$unsqueeze(3)$expand(c(N,J,P+1))
  zero_mask <- torch_where(P_range_matrix == Y_matrix, 1, 0)
  # Initialize the matrix with zeros
  matrix <- torch_zeros(N, J, P+1)
  mask.b <- P_range_matrix <= torch_tensor(P_list)$unsqueeze(2)$unsqueeze(1)$expand(c(N,J,1))
  mask.bb <- mask.b[1,,]
  mask.k <- torch_where(zero_mask > 0, 0, mask.b)
  
  Mu <- torch_zeros(N, K)
  #xi <- torch_rand(N, J, P+1)
  xi <- torch_normal(mean = 0, std = 1, size = (N * J * (P+1)))$view(c(N, J, P+1))
  xi <- xi * mask.k
  #xi<-torch_tensor(xi_hat)# test
  
  #### initialize a
  total_score <- Y$sum(dim = 2)
  covariance <- ((Y - Y$mean(dim = 1)) * (total_score - total_score$mean())$unsqueeze(2))$sum(dim = 1) / (Y$size(1) - 1)
  a <- covariance / (Y$std(dim = 1) * total_score$std())
  size0 <- J*K
  random_matrx <- (torch_normal(mean=1, std=0.1, size= size0))$view(c(J, K))
  a <- a$unsqueeze(2)$expand_as(D)
  a <- a*random_matrx
  #a <- torch_rand(J, K)*(0.5)+0.5
  a.mask <- torch_tensor(D != 0)
  a <- a * a.mask
  #a <- torch_tensor(a_hat)# test
  
  
  # Generate random normal values for b_hat_full of shape (J, P)
  #b <- torch_rand(J, P)
  b <- torch_normal(mean = 0, std = 1, size = (J * P))$view(c(J, P))
  #b <- torch_tensor(b_hat) #test
  
  # Apply the mask to b_hat, keeping values where mask is TRUE and zero elsewhere
  zeros_column <- torch_zeros(b$size(1), 1)
  # Concatenate the column of zeros to the left of b_hat
  b <- torch_cat(list(zeros_column, b), dim = 2)*mask.bb
  
  n <- tapply(1:N, X, identity)
  
  k_diff <- (P_range_matrix-Y_matrix) * mask.b
  k_diff_square <- (k_diff^2)$unsqueeze(4)$unsqueeze(5)
  k_diff <- k_diff$unsqueeze(4)$unsqueeze(5)
  
  niter.init <- mstep.MGPCM_gvem()
 
  keep(Y, D, X, iter, eps, N, n, J, K, G, MU, SIGMA, Mu, Sigma, a, b, Y.mask, a.mask, mask.b,mask.bb, mask.k, eta, xi, niter.init, p,P, zero_mask, k_diff, k_diff_square,P_range_max, SE, EFA)
  init.new()
}

mstep.MGPCM_gvem <- function() {
  with(parent.frame(), {
    params.old <- NULL
    for (ite in 1:1000) {
      aG <- (a$unsqueeze(1))$unsqueeze(4)$unsqueeze(3)
      aG.t <- t(aG)
      b_ij <- (b*zero_mask)$sum(3)$unsqueeze(3)$tile(c(1, P+1))
      BB <- (b*mask.k-b_ij)*(!zero_mask)
      AG <- k_diff*aG[X]
      AG.t <- t(AG)
      
      eta <- torch_where(abs(xi) < 1e-3, 0.125, ((torch_exp(xi) - 1) / (4 * xi * (1 + torch_exp(xi)))))
      
      Sigma <- torch_stack(replicate(N, diag(K), simplify = F))
      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + (2 * (k_diff_square*eta$view(c(N, J, -1, 1, 1)) * (aG %*% aG.t)[X])$sum(3))$sum(2)
      SIGMA <- sym(SIGMA.inv$inverse())
      MU <- (SIGMA %*% ((((( 2 * eta *BB - 0.5))$view(c(N,J, -1, 1, 1)) ) * AG)$sum(3))$sum(2))$squeeze()
      #############################
      if(K == 1)  MU <- MU$unsqueeze(2)
      mu <- (a.mask * MU$unsqueeze(2))$unsqueeze(4)
      sigma.mu <- SIGMA$unsqueeze(2) * a.mask$unsqueeze(2) * a.mask$unsqueeze(3) + mu %*% t(mu)
      Sigma <- sym(groupmean(n, sigma.mu))
      #############################
      a <- 0.5 *((k_diff_square*eta$view(c(N, J, -1, 1, 1)) * (sigma.mu)$unsqueeze(3))$sum(3))$sum(1)$pinverse() %*% ((k_diff*( 2 * eta *BB - 0.5)$view(c(N,J, -1, 1, 1)) * mu$unsqueeze(3))$sum(3))$sum(1)
      a <- a$squeeze(length(dim(a)))* a.mask
    
      #update the difficulty parameters ----
      Y_matrix <- as.matrix(Y) + 1
      a_matrix <- as.matrix(a)
      MU_matrix <- as.matrix(MU)
      eta_numeric <- as.array(eta)
      b_hat_new_numeric <- as.matrix(b)
      
      for (j in 1:J) {
        for (k in 2:(P_list[j] + 1)) {
          bs1 <- 0
          bs2 <- 0
          for (i in 1:N) {
            kij <- Y_matrix[i, j]
            dw <- a_matrix[j,] %*% MU_matrix[i,] # Precompute dw once per i
            if (k != kij) {
              x <- eta_numeric[i, j, k]
              bs1 <- bs1 + x
              bs2 <- bs2 + 2 * x * (k - kij) * dw + 0.5 + 2 * x * b_hat_new_numeric[j, kij]
            }
            if(k==kij)
            {
              for (v in 1:(P_list[j] + 1))
              {
                if (v != k) {
                  x <- eta_numeric[i, j, v]
                  bs1 <- bs1 + x
                  bs2 <- bs2 - 2 * x * (v - k) * dw - 0.5 + 2 * x * b_hat_new_numeric[j, v]
                }
              }
            }
          }
          b_hat_new_numeric[j, k] <- 0.5 * bs2 / bs1
        }
      }
      
      # Update the original b_hat_new if needed
      b <- torch_tensor(b_hat_new_numeric)
      

      #############################
      DW2 <- (t((a$unsqueeze(1))$unsqueeze(4)[X])%*% mu)$squeeze()
      DW2 <- DW2$unsqueeze(3)$tile(c(1, 1, P+1))
      a_expanded = a$unsqueeze(1)[X] # (500, 20, 3)
      result = (a_expanded %*% SIGMA) %*%  t(a)$unsqueeze(1)
      Q <- result$diagonal(dim1=-2, dim2=-1)
      b_ij <- (b*zero_mask)$sum(3)$unsqueeze(3)$tile(c(1, P+1))*mask.bb
      
      VEC <- (k_diff$squeeze() * DW2)$squeeze() - (b-b_ij)
      XS <- VEC *VEC + k_diff_square$squeeze() * Q$unsqueeze(-1)
      xi <- torch_sqrt(torch_clamp(XS, min=0))
      xi <- torch_where(zero_mask==0, xi, 0)
      
      # params <- lapply(lst(SIGMA, MU, Sigma, Mu, a, b_hat_new), torch_clone)
      # if (!is.null(params.old) && all(distance(params, params.old) < eps))
      #   break
      params <- c(as.array(a),as.array(b))
      if(!is.null(params.old) && sum(((params.old-params)*(params.old-params))/( J*K + sum(P_list)))< eps)
        break
      params.old <- params
    }
    ite
  })
}

est.MGPCM_gvem <- function(e) {
  list2env(e, environment())
  niter <- c(0, 0)
  QXI <- torch_where(abs(xi*mask.bb) < 1e-3, 0.125, ((torch_exp(xi) - 1) / (4 * xi * (1 + torch_exp(xi)))))
  mu <- (a.mask * MU$unsqueeze(2))$unsqueeze(4)
  #p
  DW <- (t((a$unsqueeze(1))$unsqueeze(4)) %*% mu)$squeeze()
  DW <- DW$unsqueeze(3)$tile(c(1, 1, P+1))
  #q
  a_expanded = a$unsqueeze(1)[X]  
  result = (a_expanded %*% SIGMA) %*%  t(a)$unsqueeze(1)[X]
  Q <- result$diagonal(dim1=-2, dim2=-1)
  b_ij <- (b*zero_mask)$sum(3)$unsqueeze(3)$tile(c(1, P+1))
  BB <- (b*mask.b-b_ij)
  VEC <- (k_diff$squeeze() * DW)$squeeze() - BB
  XS <- VEC *VEC + k_diff_square$squeeze() * Q$unsqueeze(-1)
  #llk
  LLK <- (- QXI) * XS * mask.bb - 0.5 * VEC + QXI * (xi^2) + 0.5*xi - log(1 +exp(xi)* mask.bb  )
  LLK <- ((LLK*mask.bb)$sum(3)+log(2))$sum(2)
  ll <- (LLK  - 0.5*torch_sum(MU * MU, dim = 2) - 0.5*(SIGMA$diagonal(dim1 = -2,dim2 = -1))$sum(2))$sum(1)
  final.MGPCM_gvem()
}

final.MGPCM_gvem <- function(){
  with(parent.frame(), {
    #rotation
    
    sig_hat <- diag(K)
    a<- as.matrix(a)
    #promax rotation
    if(EFA == TRUE){
    if(K>1){
      rot <- (promax(a,m=4))
      G <- (rot$rotmat);
      a <- a %*% (G)
      sig_hat <- (solve(G)) %*% sig_hat %*% t(solve(G));
      MU <- MU %*% solve(G);
    }
    se <- sqrt(diag(sig_hat))
    for(d in 1:K)
    {
      sig_hat[,d] <- sig_hat[,d]/se;
      sig_hat[d,] <- sig_hat[d,]/se;
    }
    for(d in 1:K)
    {
      if(mean(a[,d])<0){
        a[,d] <- -a[,d];MU[,d] <- -MU[,d];
        sig_hat[,d] <- -sig_hat[,d];sig_hat[d,] <- -sig_hat[d,];}
    }
  }
    else{
      sigma <- as.matrix(Sigma[1,3,])
      lambda <- diag(1/(sqrt(diag(sigma))))
      sig_hat <- lambda %*% sigma %*% t(lambda)
      se <- sqrt(diag(diag(K)))
      a <- a * a.mask
    }
    #############################
    # permn=permutations(K,K)
    # error_last <- 10^6
    # permun <- 0
    # for(k in 1:K)
    # {
    #   error <- as.numeric((torch_tensor((a[,permn[k,]] - res1$a_true[1:J, ])^2))$sum())
    #   if(error_last > error)
    #   {
    #     permun <- k
    #   }
    # }
    # a <- a[,permn[permun,]]
    #############################
    vec <- function(x) {
      torch_cat(lapply(x, function(x) {
        x$view(-1)
      }))
    }
    if(SE){
    library(mvQuad)
    SE.level=4
    grid <- suppressMessages(createNIGrid(K, 'nHN', SE.level))
    rescale(grid, m = colMeans(as.array(MU)), C = sig_hat, dec.type = 1)

    z <- torch_tensor(getNodes(grid))
    w <- as.vector(getWeights(grid)) / as.array(torch_tensor(sig_hat)$det()$sqrt())
    b0 <- torch_tensor(b)[,2:(P+1)]$requires_grad_(T)
    zeros_column <- torch_zeros(J, 1)
    # Concatenate the column of zeros to the left of b_hat
    b <- torch_cat(list(zeros_column, b0), dim = 2)$contiguous()
    a <- torch_tensor(a)$requires_grad_(T)
    xi <- (z %*% t(a))$unsqueeze(3)$expand(c(dim(z)[1],J,P+1))*(P_range_max$unsqueeze(2)$expand(c(P+1, J))$t()$unsqueeze(1)$expand(c(dim(z)[1],J,P+1)))
    xi[,,2:(P+1)] <- xi[,,2:(P+1)] - b0$unsqueeze(1)
    xi <- torch_where(mask.bb, xi, -Inf)
    logP <- xi - torch_logsumexp(xi, dim = 3)$unsqueeze(-1)
    tmp <- ((zero_mask)$unsqueeze(2) * logP$unsqueeze(1))
    tmp <- torch_where(mask.bb, tmp, 0)$sum(3:4)
    l <- torch_logsumexp(tmp+log(torch_tensor(w))$unsqueeze(1), dim = 2)
    I <- Reduce(`+`, lapply(l$unbind(),function(l) {
      d <- vec(autograd_grad(l, list(a,b0), retain_graph = T))
      d$unsqueeze(2) %*% d$unsqueeze(1)
    }))
    SE.all <- I$pinverse()$diag()$sqrt()
    SE.a <- as.array(torch_zeros(a.mask$shape)$masked_scatter(a.mask, SE.all[1:(J*K)]))
    SE.b <- SE.all[(J*K+1):length(SE.all)]$view(c(J,-1))$masked_fill(!mask.bb[,2:(P+1)],0)
    b <- b[,2:(dim(b)[2])]
    c(lapply(lst(Sigma = sig_hat, MU, a, b, SE.a, SE.b,elbo = ll), as.array))
    }else {
      b <- b[,2:(dim(b)[2])]
      c(lapply(lst(Sigma = sig_hat, MU, a, b,elbo = ll), as.array))
    }
      
  })
}


#MGPCM_gvem(data, model = matrix(1, nrow = J, ncol = 4), group = rep(1, nrow(data)), iter = 20000, eps = 1e-5,  SE = FALSE, verbose = TRUE )



