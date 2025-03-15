
# ---- Required Packages -------------------------------------------------------
# Rcpp::sourceCpp("./_GvgemsgrmAlgorithm.cpp")
# Rcpp::sourceCpp("./_IwGvemGrmAlgorithm.cpp")

# ---- Call Function ----
GRM_gvem_method <- function(data, model = matrix(1, ncol(data)), iter = 200, tol = 1e-4, S = 10, M = 10, MinDim, MaxDim, verbose = FALSE, EFA = FALSE) {
  
  # -- load simulated data --
  y <- data
  y.na <- is.na(y)
  y[y.na] <- 0.5
  data[y.na] <- 0
  # missing data
  # compare those
  # 2 categories C2PL
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(model)
  
  R <- as.vector(apply(X = data, MARGIN = 2, FUN = function(x){length(unique(x))}))
  init_B <- matrix(NA, nrow=J, ncol=max(R)-1)
  for(j in 1:J){
    init_B[j,1:(R[j]-1)] <- qlogis(cumsum(table(y[,j]))/N)[1:(R[j]-1)]
  }
  init_ksi1 <- matrix(0.05, nrow=N, ncol=J)
  init_ksi2 <- matrix(0.05, nrow=N, ncol=J)
  output0 <- if (verbose) stderr() else nullfile()
  
  if (verbose) cat(file = output0, 'Fitting the model without regularization for initial values...\n')
  
  if(EFA == TRUE)
  { 
    estDim <- 0
    bic <- 1e100
    if(MinDim > MaxDim || MinDim == 0) 
    {
      warning("Please input valid MinDim and MaxDim!") 
      return(NULL)
    }
    for(dim in c(MinDim:MaxDim))
    {
      Mod <- matrix(1, nrow = J, ncol = dim)
      K <- ncol(Mod)
      a <- apply(X = y, MARGIN = 2, FUN = function(x){ cor(x, rowSums(y))})
      init_A <- matrix(a, nrow=J, ncol=K); init_A[Mod==0] <- 0
      init_sig <- diag(K)
      
      output0 <- vem_grm(y = y, R = R, old_A = init_A, old_B = init_B, old_sig = init_sig,
                        old_ksi1 = init_ksi1, old_ksi2 = init_ksi2, Mod = Mod, max_iter = iter,
                        tol_para = tol, stop_cri = 2, verbose = ifelse(verbose == TRUE, 1, 0))
      vbic <- output0$n2vlb + log(nrow(y))*(sum(output0$new_A!=0) + sum(output0$new_B!=0,na.rm=T) + (sum(output0$new_sig!=0) - ncol(output0$new_sig))/2) 
      if(vbic < bic)
      {
        estDim <- dim
        bic <- vbic
        output <- output0
      }
    }
    K <- estDim
  }
  if(EFA == FALSE)
  {
  Mod <-  model
  a <- apply(X = y, MARGIN = 2, FUN = function(x){ cor(x, rowSums(y))})
  init_A <- matrix(a, nrow=J, ncol=K); init_A[Mod==0] <- 0
  init_sig <- diag(K)

  output <- vem_grm(y = y, R = R, old_A = init_A, old_B = init_B, old_sig = init_sig,
                    old_ksi1 = init_ksi1, old_ksi2 = init_ksi2, Mod = Mod, max_iter = iter,
                    tol_para = tol, stop_cri = 2, verbose = ifelse(verbose == TRUE, 1, 0))
  vbic <- output$n2vlb + log(nrow(y))*(sum(output$new_A!=0) + sum(output$new_B!=0,na.rm=T) + (sum(output$new_sig!=0) - ncol(output$new_sig))/2) 
  }
  
  b <- output$new_B
  b.na <- is.na(b)
  b[b.na] <- 0
  
  sigma_n <- output$sigma_n
  
  #promax rotation
  if(EFA == TRUE){
    sig_hat <- diag(K)
    a <- output$new_A
    
    if(K>1){
      
      r <- try(
        rot <- (promax(a,m=4)), silent = TRUE
      )
      if (inherits(r, "try-error")) {
        warning("Please Check the Input Dimensions!")
        return(NULL)
      } 
      G <- (rot$rotmat);
      a <- a %*% (G)
      sig_hat <- (solve(G)) %*% sig_hat %*% t(solve(G));
      result <- apply(sigma_n, MARGIN = 3, function(x) {
        (solve(G)) %*% x %*% t(solve(G))
      })
      sigma_n <- array(result, dim = dim(sigma_n))
      MU <- output$mu_n %*% solve(G);
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
    sig_hat <- output$new_sig
    a <- output$new_A 
    MU <- output$mu_n
  }
  # if(K > 1){
  # result <- apply(sigma_n, MARGIN = 3, function(x) {
  #   diag(1 / sqrt(diag(x))) %*% x %*% t(diag(1 / sqrt(diag(x))))
  # })
  # }
  # Convert back to the original 3D shape
  # sigma_n <- array(result, dim = dim(sigma_n))
  
  result <- lst(
    a    = a,
    b    = b,
    SIGMA  = sig_hat,
    ksi1 = output$new_ksi1,
    ksi2 = output$new_ksi2,
    mu     = MU,
    dim  = K,
    sigma  = sigma_n,
    n2vlb    = output$n2vlb, # negative 2 variational lower bound
    iter         = output$it
  )

  # if(verbose) {
  #   cat("\n")
  #   output <- capture.output(lst(a = a, b = b, dim = K))
  #   cat(output, sep = "\n")
  # }
  return(result)
}

IwGRM_gvem_method <- function(data, model = matrix(1, ncol(data)), iter = 2000, tol = 1e-4, S = 10, M = 10, MinDim, MaxMim, verbose = TRUE, EFA = FALSE) {
  
  result <- GRM_gvem_method(data, model = model, iter = iter, tol = tol,MinDim, MaxMim, verbose = FALSE, EFA = FALSE )
  
  y <- data
  y.na <- is.na(y)
  y[y.na] <- 0.5
  data[y.na] <- 0
  
  Mod <- model
  R <- as.vector(apply(X = data, MARGIN = 2, FUN = function(x){length(unique(x))}))
  
  output0 <- if (verbose) stderr() else nullfile()
  
  if (verbose) cat(file = output0, 'Fitting the model without regularization for initial values...\n')

  output <- iwgvem_grm(y = y, R = R, old_A = result$a, old_B = result$b, old_sig = result$SIGMA,
                       mu_n = result$mu, sigma_n = result$sigma, Mod = Mod,
                       S = S, M = M,
                       beta1 = 0.9,    # exponential decay rates
                       beta2 = 0.999,  # exponential decay rates
                       eta_A   = 0.1,    # learning rate for Adaptive moment estimation
                       eta_gam = 0.001,  # learning rate for Adaptive moment estimation
                       eta_sig = 0.005,  # learning rate for Adaptive moment estimation
                       eps   = 0.001,
                       max_iter = iter,
                       tol = tol, 
                       verbose = ifelse(verbose == TRUE, 1, 0)
  )
  
  if (verbose) pb <- txtProgressBar(file = output0, 0, 1, style = 3)
  if (verbose) setTxtProgressBar(pb, pb$getVal() + 1)
  if (verbose) close(pb)
  b <- output$new_B
  b.na <- is.na(b)
  b[b.na] <- 0
  result <- lst(
    a    = output$new_A,
    b    = b,
    SIGMA  = output$new_sig,
    iter         = output$it
  )
  if(verbose) {
    cat("\n")
    output <- capture.output(lst(a = result$a, b = result$b))
    cat(output, sep = "\n")
  }
  invisible(result)
}
#' GVEM Algorithm for the Graded Response Model
#' @param data An \eqn{N\times J} matrix of item responses where 0 is the minimal partial credit score (missing responses should be coded as \code{NA}) 
#' @param model A \eqn{J\times K} matrix of loading indicators (K is the Number of latent dimension)(all items load on the only dimension by default)
#' @param criterion Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'BIC'}, or \code{'AIC'}
#' @param iter Maximum number of iterations
#' @param tol Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param S Sample size for approximating the expected lower bound (\code{'IWGVEM'} only)
#' @param M Sample size for approximating a tighter lower bound (\code{'IWGVEM'} only)
#' @param MinDim Minimum num of possible dimensions (\code{'EFA'} only)
#' @param MaxDim Maximum num of possible dimensions (\code{'EFA'} only)
#' @param verbose Whether to show the progress
#' @param EFA Whether to run EFA or CFA
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing the following elements:

#'   \item{ ...$SIGMA}{Person-level posterior covariance matrices}
#'   \item{ ...$MU}{Person-level posterior mean vectors}
#'   \item{ ...$Sigma}{Group-level covariance matrices}
#'   \item{ ...$Mu}{Group-level mean vectors}
#'   \item{ ...$ksi1}{Variational parameter 1}
#'   \item{ ...$ksi2}{Variational parameter 2}
#'   \item{ ...$dim}{Num of dimension between latent variables}
#'   \item{ ...$a}{Slopes}
#'   \item{ ...$b}{Intercepts}
#'   \item{ ...$n2vlb}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
#'   \item{iter}{Number(s) of iterations for initialization}
#'
#' @author Yijun Cheng <chengxb@uw.edu>
#' @export
#'
#' @examples
#' \dontrun{
#' with(MGRM_data, MGRM_gvem(data, method = "IWGVEM", model, EFA = FALSE))}

MGRM_gvem <- function(data, model = matrix(1, ncol(data)), method = "GVEM", iter = 200, tol = 1e-4,S = 10, M = 10, MinDim = 0, MaxDim = 0,  verbose = FALSE, EFA = FALSE) {
  fn <- switch(method,
               GVEM = GRM_gvem_method,
               IWGVEM = IwGRM_gvem_method,
               stop(paste0("Method '", method, "' not supported.")))
  fn(data, model = model, iter = iter, tol, S = 10, M = 10, MinDim, MaxDim, verbose = verbose, EFA = EFA)
}
