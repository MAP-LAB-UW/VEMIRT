# Written by Ruoyi Zhu
# Edited by Weicong Lyu

#Starting Values
DIF_init=function(resp,Group,indic,Unif){
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
  # default for no impact (when using mirt to estimate MLE, fix the mean and variance for all groups)
  COV <- matrix(TRUE,domain,domain); diag(COV)=FALSE
  model <- mirt.model(t(indic), COV=COV) ##
  if (Unif==T){
    md.noncons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE,invariance=c('slopes'))
    starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])

    }
    gra00=as.matrix(starting5new[,1:domain])
    rownames(gra00) <- c()
    colnames(gra00) <- c()
    grd00=matrix(starting5new[,domain+1],J,1)
    grgamma00=array(0,dim=c((y-1),domain,J))
    grbeta00=as.matrix(starting5new[,(domain+1+1):(domain+1+1+y-1-1)])
    rownames(grbeta00) <- c()
    colnames(grbeta00) <- c()
    #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
    Sigma0=matrix(0,domain*y,domain)
    #Mu0 = matrix(0,domain,y)
    Mu0 = numeric(domain*y)
    for (yy in 1:y){
      Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
      #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
    }
  } else {
    md.noncons0 <- multipleGroup(verbose = F, resp, model, group = Group,SE=TRUE)
    starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,1:domain]-coef(md.noncons0,simplify=T)[[1]]$items[,1:domain])
    }
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])

    }
    gra00=as.matrix(starting5new[,1:domain])
    rownames(gra00) <- c()
    colnames(gra00) <- c()
    grd00=matrix(starting5new[,domain+1],J,1)
    grgamma00=array(0,dim=c((y-1),domain,J))
    for (yy in 2:y){
      grgamma00[(yy-1),,]=t(starting5new[,(domain+2+(yy-2)*domain):(2*domain+1+(yy-2)*domain)])
    }
    grbeta00=as.matrix(starting5new[,(2*domain+2+(y-2)*domain):(2*domain+1+(y-2)*domain+(y-1))])
    rownames(grbeta00) <- c()
    colnames(grbeta00) <- c()
    #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
    Sigma0=matrix(0,domain*y,domain)
    Mu0 = matrix(0,domain,y)
    for (yy in 1:y){
      Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
      #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
    }
  }

  return(list(G=G,y=y,r=domain,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Sigma0=Sigma0,Mu0 =Mu0))
}

##################################################
##### Function: Reg_EM_DIF                      #####
##################################################
##### Inputs
##### resp: response vector of length L, L > 1
##### m: number of response categories
##### r: number of trait dimension
##### y: number of examinee groups y=3
##### N.vec: number of examinees in each group,
##### vector of y, e.g. N.vec=c(500,500,500)
##### Mu.list: prior mean vector for each group,
##### vector of r*y, e.g., the mean vector for
##### three groups are (0,0), (0.1,0.02), (-0.05,0.03),
##### then Mu.list=c(0,0,0.1,0.02,-0.05,0.03);Mu.list=c(mu100,mu200,mu300)
##### Sig.list: prior covariance matrix for each group,
##### matrix of r*y by r, each r by r matrix is the
##### covariance matrix of each group e.g. Sig100=matrix(c(1,0.8452613,0.8452613,1),2,2)
#####                                      Sig200=matrix(c(1.179328,1.065364,1.065364,1.179328),2,2);Sig300=matrix(c(0.9202015,0.8908855,0.8908855,0.9202015),2,2)
#####                                      Sig.list=rbind(Sig100,Sig200,Sig300)
##### gra00: starting values of a, should include the information of known loading structure
##### grd00: starting values of d
##### grbeta00: starting values of beta; Items with starting value zero are used as anchor. grbeta00=matrix(0,J,2)
##### grgamma00: starting values of gamma. Items with starting value zero are used as anchor.
##################################################
##### Outputs:
##### est: estimated parameter a and d.
##### Gamma: estimated parameter gamma
##### Beta: estimated parameter beta
##### iter: number of EM cycles
##### bic: BIC
##### means: estimated mean vector for each group
##### Covs: estimated covariance matrix for each group
##################################################
# 4

Reg_EM_DIF <- function(resp,Group,indic,Unif,eta,eps =1e-3,max.tol=1e-7,r,y,N.vec=N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list, M, c, ...)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  m=2 #fixed 2pl
  #m,r,y,N.vec,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL

  # Gauss-Hermite quadrature nodes
  X1=seq(-3,3,by=0.2)
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  ng <-  numeric(G)
  Xijk=array(double(N*J*m),dim = c(N,J,m))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:m){
        Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
      }
    }
  }
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) #y1, y2, y3
  # starting values
  gra=gra00
  grd=grd00
  grbeta=grbeta00 #matrix(0,J,2)
  grgamma=grgamma00
  #grgamma=grgamma00 #array(0,dim=c((y-1),r,J))
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)

  df.a <- df.d  <- df.gamma <- df.beta <- 1
  iter <- 0

  # regularied EM
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps | max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta

    # E STEP
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    #LiA=E_step0(resp=resp,N.vec=N.vec,X=X,y=y,G=G,y.allgroup=y.allgroup,Mu.list=c(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6]),Sig.list=rbind(Sig.est[1:2,],Sig.est[3:4,],Sig.est[5:6,]),gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    #scale
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }

    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    if (iter == M)
      break
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  iter0 <- iter
  # Re-estimation
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:r){
      for (nn in 1:(y-1)){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:(y-1)){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  }

  gra=gra00
  grd=grd00
  grgamma=grgamma00*sparsity1
  grbeta=grbeta00*sparsity2

  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta

    # E STEP
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    #scale
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }


    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
      P<- array(double(G*m*(y)),dim=c(G,m,y))
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    if (iter == M)
      break
    #print(c(2, iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2

  ll <- sum(lh)
  l0 <- l0norm
  AIC <- -2 * ll + l0 * 2
  BIC <- -2 * ll + l0 * log(N)
  GIC <- -2 * ll + c * l0 * log(N) * log(log(N))
  #BIC=-2*sum(lh)+l0norm*log(N)
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta, iter = c(iter0, iter), ll = ll, l0 = l0, aic = AIC, bic = BIC, gic = GIC, means=Mu.est,Covs=Sig.est))
}

# 5

Reg_EMM_DIF <- function(resp,Group,indic,eta,Unif=F,eps =1e-3,max.tol=1e-7,r,y,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list, M, c, ...)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  #m,r,y,N.vec,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL
  m=2 #fixed 2pl
  # Gauss-Hermite quadrature nodes
  X1=seq(-3,3,by=0.2)
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  ng <-  numeric(G)
  Xijk=array(double(N*J*m),dim = c(N,J,m))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:m){
        Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
      }
    }
  }
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) #y1, y2, y3
  # starting values
  gra=gra00
  grd=grd00
  grbeta=grbeta00 #matrix(0,J,2)
  grgamma=grgamma00 #array(0,dim=c((y-1),r,J))
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)

  df.a <- df.d  <- df.gamma <- df.beta <- 1
  iter <- 0

  # regularied EM
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps | max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta

    # E STEP
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    #LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    #scale
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }

    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    if (iter == M)
      break
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
    #print(gra)
    #print(grd)
    #print(grbeta)
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2

  ll <- sum(lh)
  l0 <- l0norm
  AIC <- -2 * ll + l0 * 2
  BIC <- -2 * ll + l0 * log(N)
  GIC <- -2 * ll + c * l0 * log(N) * log(log(N))
  #BIC=-2*sum(lh)+l0norm*log(N)
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter, ll = ll, l0 = l0, aic = AIC, bic = BIC, gic = GIC, means=Mu.est,Covs=Sig.est))
}


# 6
Reg_Adaptive_DIF <- function(resp,Group,indic,eta,lam=1,Unif=F,eps =1e-3,max.tol=1e-7,r,y,N.vec=N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list, M, c)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  m=2 #fixed 2pl
  #m,r,y,N.vec,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL

  # Gauss-Hermite quadrature nodes
  X1=seq(-3,3,by=0.2)
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  ng <-  numeric(G)
  Xijk=array(double(N*J*m),dim = c(N,J,m))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:m){
        Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
      }
    }
  }
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) #y1, y2, y3
  # starting values
  gra=gra00
  grd=grd00
  grbeta=grbeta00 #matrix(0,J,2)
  grgamma=grgamma00
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)

  df.a <- df.d  <- df.gamma <- df.beta <- 1
  iter <- 0

  # regularied EM
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps | max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta

    # E STEP
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    #scale
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }

    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstepadapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    if (iter == M)
      break
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  iter0 <- iter
  # Re-estimation
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:r){
      for (nn in 1:(y-1)){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:(y-1)){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  }

  gra=gra00
  grd=grd00
  grgamma=grgamma00*sparsity1
  grbeta=grbeta00*sparsity2
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta

    # E STEP
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((X-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    #scale
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)),((Xstar-matrix(Mu.est[((yy-1)*r+1):((yy-1)*r+r)], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }


    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
      P<- array(double(G*m*(y)),dim=c(G,m,y))
      #estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstepadapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    if (iter == M)
      break
    #print(c(2,iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2

  ll <- sum(lh)
  l0 <- l0norm
  AIC <- -2 * ll + l0 * 2
  BIC <- -2 * ll + l0 * log(N)
  GIC <- -2 * ll + c * l0 * log(N) * log(log(N))
  #BIC=-2*sum(lh)+l0norm*log(N)
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta, iter = c(iter0, iter), ll = ll, l0 = l0, aic = AIC, bic = BIC, gic = GIC, means=Mu.est,Covs=Sig.est))
}


#' EM Algorithms for DIF Detection in 2PL Models
#'
#' @param Y An \eqn{N\times J} binary matrix of item responses
#' @param D A \eqn{J\times K} binary matrix of loading indicators
#' @param X An \eqn{N} dimensional vector of group indicators (integers from \code{1} to \code{G})
#' @param method Estimation algorithm, one of \code{'EM'}, \code{'EMM'} and \code{'Adapt'}
#' @param unif Whether to detect uniform DIF only
#' @param Lambda0 A vector of \code{lambda0} values for \eqn{L_1} penalty (\code{lambda} is \code{sqrt(N) * lambda0})
#' @param criterion Information criterion for model selection, one of \code{'BIC'} (recommended), \code{'AIC'} and \code{'GIC'}
#' @param iter Maximum number of iterations
#' @param eps Termination criterion on numerical accuracy
#' @param c Constant for computing GIC
#' @param eta Tuning constant for adaptive lasso (\code{'Adapt'} only)
#'
#' @return An object of class \code{vemirt_DIF}, which is a list containing three elements:
#'   \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
#'   \item{best}{The location of \code{fit} in \code{all}}
#'   \item{all}{A list of models whose length is equal to \code{Lambda0}:}
#'   \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
#'   \item{ ...$lambda}{\code{sqrt(N) * lambda0}}
#'   \item{ ...$iter}{Number(s) of iterations}
#'   \item{ ...$Sigma}{Group-level posterior covariance matrices}
#'   \item{ ...$Mu}{Group-level posterior mean vectors}
#'   \item{ ...$a}{Slopes for group 1}
#'   \item{ ...$b}{Intercepts for group 1}
#'   \item{ ...$gamma}{DIF parameters for the slopes}
#'   \item{ ...$beta}{DIF parameters for the intercepts}
#'   \item{ ...$ll}{Log-likelihood}
#'   \item{ ...$l0}{Number of nonzero parameters in \code{gamma} and \code{beta}}
#'   \item{ ...$AIC}{Akaike Information Criterion}
#'   \item{ ...$BIC}{Bayesian Information Criterion}
#'   \item{ ...$GIC}{Generalized Information Criterion}
#'
#' @seealso \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(exampleDIF, em_DIF(Y, D, X))}
em_DIF <- function(Y, D, X, method = 'EMM', unif = F, Lambda0 = seq(0.1, 0.8, by = 0.1), criterion = 'BIC', iter = 1000, eps = 1e-3, c = 0.7, eta = 1) {
  if (!(method %in% c('EM', 'EMM', 'Adapt')))
    stop(paste0("Method '", method, "' not supported."))
  if (!(criterion %in% c('AIC', 'BIC', 'GIC')))
    stop(paste0("Information criterion '", criterion, "' not supported."))
  regdif <- switch(method, EM = Reg_EM_DIF, EMM = Reg_EMM_DIF, Adapt = Reg_Adaptive_DIF)
  X <- as.factor(X)
  resp <- as.data.frame(Y)
  if (min(resp) == 0) {
    resp2 <- as.matrix(resp)
    resp <- resp + 1
  } else
    resp2 <- as.matrix(resp) - 1
  indic <- t(D)
  Group <- X

  cat('Running multipleGroup() from mirt for initial values...\n')
  init=DIF_init(resp=resp2,Group=Group,indic=indic,Unif=unif)

  m=2 #fixed 2pl
  r=init$r
  y=init$y
  N.vec=as.vector(table(Group))
  gra00=init$gra00
  grd00=init$grd00
  grbeta00=init$grbeta00
  grgamma00=init$grgamma00
  grbetamle=init$grbeta00
  grgammamle=init$grgamma00
  Mu.list=init$Mu0
  Sig.list=init$Sigma0

  person=nrow(resp)
  item=ncol(resp)
  domain=nrow(indic)

  #y=length(unique(Group))
  #lbd.center=10+0.033*mean(N.vec)
  #lbd.vec=seq((lbd.center-10),(lbd.center+10),4)
  lbd.vec <- sqrt(person) * Lambda0
  bics=rep(0,length(lbd.vec))
  ADmat=array(double(item*(domain+1)*length(lbd.vec)),dim = c(item,(domain+1),length(lbd.vec)))
  Gammas=array(double((y-1)*domain*item*length(lbd.vec)),dim = c((y-1),domain,item,length(lbd.vec)))
  Betas=array(double(item*(y-1)*length(lbd.vec)),dim = c(item,(y-1),length(lbd.vec)))
  Mus=matrix(0,domain*y,length(lbd.vec))
  Sigs=array(double(domain*domain*y*length(lbd.vec)),dim = c(domain*y,domain,length(lbd.vec)))

  cat('Fitting the model using different lambdas...\n')
  pb <- txtProgressBar(0, length(lbd.vec), style = 3)
  result <- lapply(1:length(lbd.vec), function(k) {
    lambda <- lbd.vec[k]
    sim <- regdif(resp=resp,indic=indic,eta=lambda,lam=eta,Group=Group,Unif=unif,eps=eps,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list, M = iter, c = c)
    setTxtProgressBar(pb, pb$getVal() + 1)
    list(lambda0 = Lambda0[k], lambda = lambda, iter = sim$iter, Sigma = aperm(array(sim$Covs, c(domain, y, domain)), c(2, 1, 3)),
         Mu = matrix(sim$means, nrow = y, byrow = T), a = sim$est[, 1:domain], b = -sim$est[, domain + 1],
         gamma = abind(array(0, c(item, 1, domain)), aperm(sim$Gamma, c(3, 1, 2)), along = 2), beta = cbind(0, sim$Beta),
         ll = sim$ll, l0 = sim$l0, AIC = sim$aic, BIC = sim$bic, GIC = sim$gic)
  })
  close(pb)
  new_vemirt_DIF(0, result, person, criterion)
}


#' Likelihood Ratio Test for DIF Detection in 2PL Models
#'
#' @param Y An \eqn{N\times J} binary matrix of item responses
#' @param D A \eqn{J\times K} binary matrix of loading indicators
#' @param X An \eqn{N} dimensional vector of group indicators (integers from \code{1} to \code{G})
#' @param unif Whether to detect uniform DIF only
#'
#' @return A list:
#'   \item{Sigma}{Group-level posterior covariance matrices}
#'   \item{Mu}{Group-level posterior mean vectors}
#'   \item{a}{Slopes for group 1}
#'   \item{b}{Intercepts for group 1}
#'   \item{gamma}{DIF parameters for the slopes}
#'   \item{beta}{DIF parameters for the intercepts}
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
#'
#' @examples
#' \dontrun{
#' with(exampleDIF, lrt_DIF(Y, D, X))}
lrt_DIF <- function(Y, D, X, unif = F){
  X <- as.factor(X)
  resp <- as.data.frame(Y)
  if (min(resp) == 0) {
    resp2 <- as.matrix(resp)
    resp <- resp + 1
  } else
    resp2 <- as.matrix(resp) - 1
  indic <- t(D)
  Group <- X

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
