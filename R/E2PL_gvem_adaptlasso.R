#' Exploratory M2PL Analysis with Adaptive Lasso Penalty
#'
#'
#' @param u an \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param indic a \eqn{J \times K} \code{matrix} or a \code{data.frame} that
#' describes the factor loading structure of \eqn{J} items to \eqn{K} factors. It
#' consists of binary values where 0 refers to the item is irrelevant to this factor, and
#' 1 otherwise. For exploratory factor analysis with adaptive lasso penalty, \code{indic} should
#' include constraints on the a \eqn{K \times K} sub-matrix to ensure identifiability.
#' The remaining parts do not assume any pre-specified zero structure but instead, the
#' appropriate lasso penalty would recover the true zero structure. Also see \code{constrain}
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param constrain the constraint setting: \code{"C1"} or \code{"C2"}. To ensure
#' identifiability, \code{"C1"} sets a \eqn{K \times K} sub-matrix of \code{indic} to be an
#' identity matrix.This constraint anchor \eqn{K} factors by designating \eqn{K} items that load solely on each factor respectively.
#' Note that the \eqn{K \times K} matrix does not have to appear at the top of the \code{indic} matrix.
#' \code{"C2"} sets the \eqn{K \times K} sub-matrix to be a lower triangular matrix with the diagonal being ones. That is, there
#' are test items associated with each factor for sure and they may be associated
#' with other factors as well. Nonzero entries (in the lower triangular part) except for the diagonal entries of the
#' sub-matrix are penalized during the estimation procedure. For instance, assume \eqn{K=3}, then the \code{"C2"} constraint will
#' imply the following submatrix: \eqn{C2=\begin{bmatrix} 1 & 0 & 0\\ 1 & 1 & 0\\ 1 & 1 & 1\\\end{bmatrix}}. As shown, item 1 is allowed to only
#' load on the first factor, item 2 will for sure load on the second factor but it may also load on the first factor (hence a penalty is added
#' on the \eqn{(2,1)} element of \code{"C2"}, i.e., \eqn{C2_{2,1}} ). Item 3 will for sure load on the third factor but it may also load on the
#' first two factors. However, note that for all remaining items their loading vector will all be \eqn{(1, 1, 1)} hence indistinguishable from the
#' third anchor item. Therefore, we need to alert the algorithm that this third anchor item will for sure load on the third factor, and
#' whether or not it loads on the first two factors depends on the regularization results. Therefore, we need to specify
#' \code{"non_pen="} to identify the \eqn{K}th anchor item. Although, \code{"C2"} is much weaker than \code{"C1"}, it still ensures empirical identifiability. Default is \code{"C1"}.
#' During estimation, under both the \code{"C1"} and \code{"C2"} constraints, the population means and variances are constrained to be 0 and 1, respectively.
#' @param non_pen the index of an item that is associated with every factor under constraint \code{"C2"}.
#' For \code{C1}, the input can be \code{NULL}
#' @param gamma a numerical value of adaptive lasso parameter. Zou (2006) recommended three values, 0.5, 1, and 2.
#' The default value is 2.
#' @return a list containing the following objects:
#'   \item{ra}{item discrimination parameters, a \eqn{J \times K} \code{matrix}}
#'   \item{rb}{item difficulty parameters, vector of length \eqn{J}}
#'   \item{reta}{variational parameters \eqn{\eta(\xi)}, a \eqn{N \times J} matrix}
#'   \item{reps}{variational parameters \eqn{\xi}, a \eqn{N \times J} matrix}
#'   \item{rsigma}{population variance-covariance matrix, a \eqn{K \times K} matrix}
#'   \item{mu_i}{mean parameter for each person, a \eqn{K \times N} matrix}
#'   \item{sig_i}{covariance matrix for each person, a \eqn{K \times K \times N} array}
#'   \item{n}{the number of iterations for the EM cycle}
#'   \item{Q_mat}{factor loading structure, a \eqn{J \times K} matrix}
#'   \item{GIC}{model fit index}
#'   \item{AIC}{model fit index}
#'   \item{BIC}{model fit index}
#'   \item{lbd}{numerical value of lasso penalty parameter \eqn{\lambda}}
#' @references
#' Cho, A. E., Xiao, J., Wang, C., & Xu, G. (2022). Regularized Variational Estimation for Exploratory Item Factor Analysis. \emph{Psychometrika}. https://doi.org/10.1007/s11336-022-09874-6
#'
#' Zou, H. (2006). The adaptive LASSO and its oracle properties.  \emph{Journal of the American Statistical Association, 7}, 1011418–1429.
#'
#'
#' @author Jiaying Xiao <jxiao6@uw.edu>
#' @seealso \code{\link{E2PL_gvem_rot}}, \code{\link{E2PL_gvem_lasso}}, \code{\link{exampleIndic_efa2pl_c1}}, \code{\link{exampleIndic_efa2pl_c2}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(E2PL_data_C1, E2PL_gvem_adaptlasso(data, model, constrain = constrain, non_pen = non_pen, gamma=2))
#' with(E2PL_data_C2, E2PL_gvem_adaptlasso(data, model, constrain = constrain, non_pen = non_pen, gamma=2))}

#main function for gvem_2PLEFA_adaptlasso
E2PL_gvem_adaptlasso<-function(u,indic,max.iter=5000,constrain="C1",non_pen=NULL,gamma=2){
  #start=Sys.time()
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  if(constrain=="C1"){
    result=vem_2PLEFA_adaptive_const1_all(u,domain,indic,gamma,max.iter)
  }else{
    if(is.null(non_pen)){
      stop('non_pen argument is required for the C2 constraint',call.=FALSE)
    }else{
      result=vem_2PLEFA_adaptive_const2_all(u,domain,gamma,indic,non_pen,max.iter)
    }
  }
  if(result$lbd==0.1 || result$lbd==40){
    warning("The optimal penalty parameter may be out of range.",call. = FALSE)
  }
  if(result$n==max.iter){
    warning("The maximum number of EM cycles reached!",call.=FALSE)
  }
  #end=Sys.time()
  #duration=end-start
  #message(paste("Total Execution Time:", round(duration[[1]], 2),  units(duration)),"\n")
  new.vemirt_FA(result)
}

#adaptive lasso with constraint 1 function
vem_2PLEFA_adaptive_const1 <- function(u,new_a,new_b,eta,xi,Sigma, domain,lbd,indic,nopenalty_col,weights,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  #is_singular = 0
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  while(converge==1 && Matrix::rankMatrix(Sigma) == domain && n<max.iter){
    #update Sigma, MU, SIGMA
    par_Sigma = Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)


    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)

    par_a=new_a
    #update a
    new_a1=nalc12pl(u, indic, nopenalty_col, person, eta, new_b, SIGMA, MU)
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    #adaptive lasso penalty
    sdf=setdiff(1:item,nopenalty_col)
    new_a=paal2pl(u, domain, person, lbd, sdf, eta, new_a1, new_b, SIGMA, MU,weights)
    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")<0.0001){
      converge=0
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function for adaptive lasso l1: choose optimal lambda
vem_2PLEFA_adaptive_const1_all<-function(u,domain,indic,gamma,max.iter){
  lbd=seq(2,20,2)
  nopenalty_col=which(rowSums(indic)==1)
  person=dim(u)[1]
  item=dim(u)[2]
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  #create weights: use CFA
  wa=C2PL_gvem(u,indic,max.iter)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,max.iter)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
    Sigma=rl[[j]]$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }
  }
  id=which.min(gic)
  rs=C2PL_gvem(u,rl[[id]]$Q_mat,max.iter)
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)
}

#adaptive lasso with constraint 2 function
vem_2PLEFA_adaptive_const2 <- function(u, new_a,new_b,eta,xi,Sigma, domain,lbd,
                                       indic,nopenalty_col,weights,max.iter) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1


  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0

  par_Sigma = Sigma
  while(converge==1 && Matrix::rankMatrix(Sigma) == domain && n<max.iter){
    par_Sigma = Sigma
    #update MU, SIGMA, Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)


    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)

    par_a=new_a
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    new_a=nalc22pl(u, domain,new_a,nopenalty_col,lastone, person, eta, new_b, SIGMA, MU)
    #L1-penalty: off-diagnoal
    new_a=paalc22pl(u, new_a,indic,nopenalty_col,lastone, lbd, person, eta, new_b, SIGMA, MU,weights)
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    new_a=paalc22pl1(u, domain, item, person, lbd, eta, new_a, new_b, SIGMA, MU,pc,weights)

    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=0
    }
    n=n+1
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    #is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function: choose optimal lambda
vem_2PLEFA_adaptive_const2_all<-function(u,domain,gamma,indic,non_pen,max.iter){
  lbd=seq(2,20,2)
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=c(which(rowSums(indic)<domain),non_pen)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  #create weights: use CFA
  wa=C2PL_gvem(u,indic,max.iter)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],
                                        indic,nopenalty_col,weights,max.iter)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
    Sigma=rl[[j]]$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,max.iter)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }
  }
  id=which.min(gic)
  rs=C2PL_gvem(u, rl[[id]]$Q_mat,max.iter)
  rs$lbd=lbd[id]
  #rs$id=id
  return(rs)
}
