#' Confirmatory M2PL Analysis
#'
#' @param u an \eqn{N \times J} \code{matrix} or a \code{data.frame} that
#' consists of binary responses of \eqn{N} individuals to \eqn{J} items. The
#' missing values are coded as \code{NA}
#' @param indic a \eqn{J \times K} \code{matrix} or a \code{data.frame} that
#' describes the factor loading structure of \eqn{J} items to \eqn{K} factors. It
#' consists of binary values where 0 refers to the item is irrelevant with this factor,
#' 1 otherwise
#' @param max.iter the maximum number of iterations for the EM cycle; default is 5000
#' @param SE.est whether to estimate SE for item parameters using the updated
#' supplemented expectation maximization (USEM); default is FALSE
#'
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
#'   \item{SE}{Standard errors of item parameters, a \eqn{J \times (K+1)} matrix where the last column includes SE estimates for item difficulty parameters}
#'
#' @author Jiaying Xiao <jxiao6@uw.edu>
#' @seealso \code{\link{C3PL_sgvem}},\code{\link{importanceSampling}},\code{\link{C2PL_bs}}
#' @export
#'
#' @examples
#' \dontrun{
#' C2PL_gvem(exampleData_2pl, exampleIndic_cfa2pl)}
C2PL_gvem <- function(u,indic,max.iter=5000,SE.est=FALSE) {
  start=Sys.time()
  u=data.matrix(u)
  indic=data.matrix(indic)
  domain=dim(indic)[2]
  person=dim(u)[1]
  item=dim(u)[2]
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]] * indic
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  converge = 0
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  ntheta<-list()
  while(converge==0 && Matrix::rankMatrix(Sigma) == domain && n < max.iter){
    par_Sigma = Sigma
    #update MU, SIGMA, sigma, eta
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    xi[is.na(u)]=NA
    Spart=rs1$Spart
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)

    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part,na.rm = T)/colSums(2*eta,na.rm = T)
    par_a=new_a
    #update a
    new_a=acfa2(u, indic, person, item, domain, eta, new_b, SIGMA, MU,new_a)
    #2022/01/19:rescale a
    new_a=new_a%*%d_temp
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2") <0.0001){
      converge=1
    }
    n=n+1
    ntheta[[n]]<-list(new_a=new_a,new_b=new_b,eta = eta,Sigma = Sigma)
  }
  if (Matrix::rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
  }else{
    rsigma = Sigma
  }

  ntheta[[n]]<-list(new_a=new_a,new_b=new_b,eta = eta,Sigma = Sigma)
  se=matrix(NA,item,domain+1)
  if(SE.est==TRUE){
    bj2=2*colSums(eta)
    abj=-2*MU%*%eta*t(indic)
    a2=ea2(domain,person,item,eta,SIGMA,MU,indic)
    l=domain*item+item
    ic=matrix(0,l,l)
    #delta0=matrix(0,l,l)
    for (j in 1:item) {
      a21=a2[,,j]
      ic1=rbind(cbind(a21,abj[,j]),c(abj[,j],bj2[j]))
      ic[(1+(j-1)*(domain+1)):(j*(domain+1)),(1+(j-1)*(domain+1)):(j*(domain+1))]=ic1
    }
    d1=sum(indic)+item
    L=matrix(0,l,d1)
    L[cbind(which(diag(ic)!=0),1:d1)]=1
    newic=t(L)%*%ic%*%L
    par_mle=as.vector(t(cbind(new_a,new_b)))%*%L
    converge1 = 0
    n1=1

    delta_par=matrix(1,d1,d1)
    par1=par_mle
    iter=1:d1
    while(converge1==0 && n1 < n){
      delta=matrix(0,d1,d1)
      delta[-iter,]=delta_par[-iter,]
      par2=as.vector(t(cbind(ntheta[[n1]]$new_a,ntheta[[n1]]$new_b)))%*%L
      for (i in iter) {
        par1=par_mle
        par1[i]=par2[i]
        new_a1=t(matrix(par1%*%t(L),domain+1))[,1:domain]
        new_b1=t(matrix(par1%*%t(L),domain+1))[,-c(1:domain)]
        #update MU, SIGMA, sigma, eta
        par_Sigma1 = Sigma
        rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a1, new_b1)
        xi1=rs1$xi
        xi1[is.na(u)]=NA
        Spart1=rs1$Spart
        SIGMA1=rs1$SIGMA
        MU1=rs1$MU
        eta1=(exp(xi1)/(1+exp(xi1))-0.5)/(2*xi1)
        eta1=replace(eta1,abs(eta1)< 0.01,0.125)
        Sigma1=Spart1/person
        d_temp1=sqrt(diag(diag(Sigma1)))
        Sigma1=solve(d_temp1)%*%Sigma1%*%solve(d_temp1)

        #update b
        par_b1=new_b1
        b_part1 = t(new_a1%*%MU1)
        new_b1=colSums(0.5-u+2*eta1*b_part1,na.rm = T)/colSums(2*eta1,na.rm = T)
        #update a
        par_a1=new_a1
        new_a1=acfa2(u, indic, person, item, domain, eta1, new_b1, SIGMA1, MU1,new_a1)
        #2022/01/19:rescale a
        new_a1=new_a1%*%d_temp1
        par_new=as.vector(t(cbind(new_a1,new_b1)))%*%L
        #delta[,i]=(par_new-par_mle)/(par1[i]-par_mle[i])
        delta[i,]=(par_new-par_mle)/(par1[i]-par_mle[i])
        if (norm(as.vector(new_a1)-as.vector(par_a1),type="2")+norm(new_b1-par_b1,type="2")+
            norm(as.vector(Sigma1)-as.vector(par_Sigma1),type="2") < 0.01){
          iter=iter[!i]
        }
      }

      if( sum(iter)==0){
        converge1=1
      }
      delta_par=delta

      n1=n1+1

    }
    io1=solve(newic)%*%solve(diag(d1)-delta)
    isSymmetric(io1)
    vb=sqrt(diag(io1))
    vb[is.na(vb)]=0
    se=t(matrix(vb%*%t(L),domain+1))
  }
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat) - 2*lbound
  aic = 2*sum(Q_mat) -2*lbound
  if(n==max.iter){
    warning("The maximum number of EM cycles reached!",call.=FALSE)
  }
  end=Sys.time()
  duration=end-start
  cat(paste("Total Execution Time:", round(duration[[1]], 2),  units(duration)),"\n")
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic,SE=se))
}
