#' Likelihood Ratio Test for DIF Detection in M2PL Models
#'
#' @param data An \eqn{N\times J} binary matrix of item responses
#' @param model A \eqn{J\times K} binary matrix of loading indicators
#' @param group An \eqn{N} dimensional vector of group indicators (integers from \code{1} to \code{G})
#' @param unif Whether to detect uniform D2PL only
#'
#' @return A list:
#'   \item{Sigma}{Group-level posterior covariance matrices}
#'   \item{Mu}{Group-level posterior mean vectors}
#'   \item{a}{Slopes for group 1}
#'   \item{b}{Intercepts for group 1}
#'   \item{gamma}{D2PL parameters for the slopes}
#'   \item{beta}{D2PL parameters for the intercepts}
#'
#' @author Ruoyi Zhu <zhux0445@uw.edu>
#' @seealso \code{\link{D2PL_em}}, \code{\link{D2PL_pair_em}}, \code{\link{D2PL_gvem}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(D2PL_data, D2PL_lrt(data, model, group))}
D2PL_lrt <- function(data, model, group, unif = F){
  resp <- as.data.frame(data)
  indic <- t(model)
  Group <- as.factor(group)
  if (min(resp) == 0) {
    resp2 <- as.matrix(resp)
    resp <- resp + 1
  } else
    resp2 <- as.matrix(resp) - 1

  resp=resp2

  m=2 ##fixed, 2pl only
  N=nrow(resp)
  J=ncol(resp)
  domain=nrow(indic)
  y=length(unique(Group))
  y.allgroup=rbind(rep(0,y-1),diag(y-1))
  G=matrix(0,N,y-1)
  for (yy in 1:y){
    vec=which(Group==sort(unique(Group))[yy])
    for (i in 1:length(vec)){
      G[vec[i],]=y.allgroup[yy,]
    }
  }
  # defalt for no impact (when using mirt to estimate MLE, fix the mean and variance for all groups)
  COV <- matrix(TRUE,domain,domain); diag(COV)=FALSE
  model <- mirt.model(t(indic), COV=COV) ##

  #md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('slopes'))
  rownames(indic)=paste0("a",1:domain)
  if (unif){
    anchors0=c(1:J)
    diff=1
    while(diff>0){
      md.cons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors0]))
      d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
      anchors=which(d$adj_p>0.05)
      diff=length(anchors0)-length(anchors)
      anchors0=anchors
      #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
    }
    md.noncons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))
    dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])
    dif1.t=dif1[which(dif1$adj_p<0.05),]

    #refit
    if (length(rownames(dif1.t))==0){
      md.refit02 <-multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:J]))
    } else {
      md.refit02 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
    }
    #Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Mu=Mu,Sig=Sig,domain=domain,y=y
    est=cbind(coef(md.refit02,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,domain+1]-coef(md.refit02,simplify=T)[[1]]$items[,domain+1])
    }
    gra.est=as.matrix(est[,1:domain])
    grd.est=matrix(est[,(domain+1)],J,1)
    grgamma.est=array(0,dim=c((y-1),domain,J))
    grbeta.est=as.matrix(est[,(domain+1+1):(domain+1+1+y-1-1)])
    Sigma.est=matrix(0,domain*y,domain)
    Mu.est = numeric(domain*y)
    for (yy in 1:y){
      Sigma.est[((yy-1)*domain+1):(yy*domain),]=coef(md.refit02,simplify=T)[[yy]]$cov
      Mu.est[((yy-1)*domain+1):(yy*domain)]=coef(md.refit02,simplify=T)[[yy]]$means
    }
  } else {
    anchors0=c(1:J)
    diff=1
    while(diff>0){
      md.cons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors0]))
      d=DIF(md.cons0, which.par = c(rownames(indic),'d'), p.adjust = 'fdr',scheme = 'drop')
      anchors=which(d$adj_p>0.05)
      diff=length(anchors0)-length(anchors)
      anchors0=anchors
      #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
    }
    md.noncons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors]))
    Jt=c(1:J)[-anchors]
    dif1.t=NULL
    for (jj in Jt){
      an=rownames(indic)[which(indic[,jj]==1)]
      dif1=DIF(md.noncons0, which.par = c(an,'d'), p.adjust = 'fdr',scheme = 'add', items2test=jj)
      dif1.t=rbind(dif1.t,dif1[which(dif1$adj_p<0.05),])
    }

    #refit
    if (length(rownames(dif1.t))==0){
      md.refit02 <-multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:J]))
    } else {
      md.refit02 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
    }

    est=cbind(coef(md.refit02,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,1:domain]-coef(md.refit02,simplify=T)[[1]]$items[,1:domain])
    }
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,domain+1]-coef(md.refit02,simplify=T)[[1]]$items[,(domain+1)])
    }
    gra.est=as.matrix(est[,1:domain])
    grd.est=matrix(est[,domain+1],J,1)
    grgamma.est=array(double((y-1)*domain*J),dim = c((y-1),domain,J))
    for (yy in 1:(y-1)){
      grgamma.est[yy,,]=t(est[,(domain+1+(yy-1)*domain+1):(domain+1+(yy-1)*domain+domain)])
    }
    grbeta.est=as.matrix(est[,(domain+1+(y-1)*domain+1):(domain+1+(y-1)*domain+1+y-1-1)])
    Sigma.est=matrix(0,domain*y,domain)
    Mu.est = numeric(domain*y)
    for (yy in 1:y){
      Sigma.est[((yy-1)*domain+1):(yy*domain),]=coef(md.refit02,simplify=T)[[yy]]$cov
      Mu.est[((yy-1)*domain+1):(yy*domain)]=coef(md.refit02,simplify=T)[[yy]]$means
    }
  }
  #gra.est=matrix(0,J,domain)
  #grd.est=matrix(0,J,1)
  #grgamma.est=array(0,dim=c((y-1),domain,J))
  #grbeta.est=matrix(0,J,y-1)
  #Sigma.est=matrix(0,domain*y,domain)
  #Mu.est = numeric(domain*y)

  list(Sigma = aperm(array(Sigma.est, c(domain, y, domain)), c(2, 1, 3)), Mu = matrix(Mu.est, nrow = y, byrow = T),
       a = gra.est, b = -grd.est, gamma = abind(array(0, c(item, 1, domain)), aperm(grgamma.est, c(3, 1, 2)), along = 2),
       beta = cbind(0, grbeta.est))
}
