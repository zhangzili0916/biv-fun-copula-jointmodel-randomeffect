library(pracma)
library(cubature)
library(expm)
library("copula")
library(CVTuningCov)
library(JM)
library(joineR)
library(nlme)
library("abind")
library("magic")
library("matrixcalc")
library(Deriv)
library(mvtnorm)
library(emdbook)
library(latex2exp)
library(tikzDevice)
library(fda)
library(VineCopula)



GaussianHermite <- function(r)               #function for generating Gaussian hermite nodes and their corresponding weights
{
  esp <- 3.0e-14
  pim4 <- 1/sqrt(sqrt(pi)) 	
  maxit <- 1000
  x <- w <- rep(0,r) 
  m <- (r+1 - ((r+1)%%2))/2
  p <- 2*r
  q <- 2*r+1
  for(i in 1:m){
    if(i==1) z <- sqrt(q) -1.85575*q^(-.1667)
    if(i==2) z <- z-1.14^((.425)/z)
    if(i==3) z <- 1.86*z - 0.86*x[1]
    if(i==4) z <- 1.91*z - 0.91*x[2]
    if(i!=1 & i!=2 & i!=3 & i!=4)  z <- 2*z - x[i-2]
    its <- 1
    dif <- 1
    while(its<maxit & dif>esp){
      p1 <- pim4
      p2 <- 0
      for(j in 1:r){
        p3 <- p2
        p2 <- p1
        p1 <- z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3
      }
      pp <- sqrt(p)*p2
      z <- z - p1/pp
      dif <- abs(p1/pp)
      its <- its +1
    }
    x[i] <- z
    x[r+1-i] <- -z
    w[i] <- 2/(pp^2)
    w[r+1-i] <- w[i]
  }
  return(list(round(w,15),round(x,15)))
}


#A simulation study with unbalanced measurement time    t copula and Gaussian copula
#generate data with functional correlation rho(t)
Y.bivcopjoineR.rho.BS.unb=function(beta1,beta2,n,sigma,var.random,lambda,alpha,etapar,rate,trprob,error,
                                   genprob,ageprob,mi,tmax,obsprob,etaord,nu,copula,extra)
{
  schedule=seq(0,tmax,length=mi)
  Data=NULL
  G.eigen=eigen(var.random)
  G.sqrt=G.eigen$vectors%*%diag(sqrt(G.eigen$values))%*%t(G.eigen$vectors)
  for(subj in 1:n)
  {
    obsposi=c(1,rbinom((mi-1),1,obsprob))>0
    jit=runif((mi-1),min=-error,max=error)
    timepoint=c(schedule[1],schedule[2:mi]+jit)[obsposi]
    etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),norder=etaord,nbasis=length(etapar)))%*%etapar
    rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
    ni=length(timepoint)
    treat=rbinom(1,1,trprob)
    gender=rbinom(1,1,genprob)
    age=sample(c(0,1,2),1,prob=ageprob)
    Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
    Xi1=matrix(c(rep(1,ni),timepoint,rep(treat,ni),rep(gender,ni),rep(as.numeric(age==1),ni)
                 ,rep(as.numeric(age==2),ni)),ncol=6)
    Xi2=matrix(c(treat,gender,as.numeric(age==1),as.numeric(age==2)),ncol=4)
    ceni=rexp(1,rate)
    bi=G.sqrt%*%rnorm(2)
    Yimean=Xi1%*%beta1+Zi1%*%bi
    Ci=lambda*exp(beta2[1]*treat+beta2[2]*gender
                  +beta2[3]*Xi2[3]+beta2[4]*Xi2[4]+alpha[1]*bi[1])
    Ti=0
    Yisur=0
    Yiobs=0
    timepoint[ni+1]=Inf
    j=1
    while((j<(ni+1))&(Ti<=timepoint[j]))
    {
      paramar=list(list(min=0,max=1))
      paramar[[2]]=list(mean=Yimean[j],sd=sigma)
      if(copula=="tcop")
      {
        Joint=rMvdc(1,mvdc(copula = ellipCopula(family = "t", dim=2,df=nu,dispstr = "un",
                                                param =rhoty[j]),margins = c("unif","norm"),paramMargins=paramar))
      }
      if(copula=="Gaucop")
      {
        Joint=rMvdc(1,mvdc(copula = ellipCopula(family = "normal", dim=2,df=nu,dispstr = "un",
                                                param =rhoty[j]),margins = c("unif","norm"),paramMargins=paramar))
      }
      Uij=Joint[1]
      Yisur[j]=Joint[2]
      Fij=1-exp(-Ci/(alpha[2]*bi[2])*(exp(alpha[2]*bi[2]*timepoint[j+1])-exp(alpha[2]*bi[2]*timepoint[j])))
      if(Uij>Fij)
      {
        j=j+1
        Ti=timepoint[j]
      }
      if(Uij<=Fij)
      {
        Tij=1/(alpha[2]*bi[2])*log(exp(alpha[2]*bi[2]*timepoint[j])-alpha[2]*bi[2]*log(1-Uij)/Ci)
        Ti=Tij
      }
    }
    Tiobs=min(Ti,ceni,(tmax+extra))
    indi=as.numeric(Ti<=ceni&&Ti<=(tmax+extra))
    dimyi=sum(timepoint<=Tiobs)
    Yiobs=Yisur[1:dimyi]
    Datai=matrix(0,nrow=dimyi,ncol=13)
    Datai[,1]=Yiobs
    Datai[,2]=timepoint[1:dimyi]
    Datai[,3]=subj
    Datai[,4]=Ti 
    Datai[,5]=ceni
    Datai[,6]=Tiobs
    Datai[,7]=indi
    Datai[,8]=treat
    Datai[,9]=gender
    Datai[,10]=Xi2[,3]
    Datai[,11]=Xi2[,4] 
    Datai[,12]=bi[1]
    Datai[,13]=bi[2]
    Data=rbind(Data,Datai)
  }
  Data=data.frame(Data)
  colnames(Data)=c('resp','ti','subject','surti','centi','obsti','indicator'
                   ,'treat','gender','middle','young','bi0','bi1')
  return(Data)
}



#log-likelihood function for the bivariate functional Gaussian copula joint model 
#Estimate etaty by B-spline first then transfer to rhoty in Gaussian copula   
appGQbivGaucoploglik.funrhoty=function(theta,data,m,tmax,neta,etaord)         
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  etapar=theta[14:(13+neta)]
  lambda=theta[14+neta]
  sigma=theta[15+neta]
  alpha=theta[16+neta]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      timepoint=data$ti[data$subject==sub]
      etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=neta,norder=etaord))%*%etapar
      rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
      ni=length(timepoint)
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],data$gender[data$subject==sub],
                   data$middle[data$subject==sub],data$young[data$subject==sub]),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(data$middle[data$subject==sub]),unique(data$young[data$subject==sub]))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      V11=Zi1%*%D%*%t(Zi1)+sigma^2*diag(ni)
      V22=D
      V12=Zi1%*%D
      V21=t(V12)
      mubiconyi=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)
      Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      Zyij=(matrix(rep(data$resp[data$subject==sub],m^2),ncol=m^2)-
              matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)/sigma
      postijcondimean=1
      if(ni>1)
      {
        for(k in 1:(ni-1))
        {
          Sijk=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*timepoint[k+1])-exp(alpha*bij[2,]*timepoint[k])))
          Zyijk=Zyij[k,]
          mutijkcon=rhoty[k]*Zyijk
          sigmatijkcon=sqrt(1-rhoty[k]^2)
          ftijkcon=pnorm((qnorm(Sijk)+mutijkcon)/sigmatijkcon)
          postijcondimean=ftijkcon*postijcondimean  
        }
      }
      if(obsti!=timepoint[ni])
      {
        hijni=Cij*exp(alpha*bij[2,]*obsti)
        Sijni=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-exp(alpha*bij[2,]*timepoint[ni])))
        fijni=hijni*Sijni
        Zyijni=Zyij[ni,]
        mutijnicon=rhoty[ni]*Zyijni
        sigmatijnicon=sqrt(1-rhoty[ni]^2)
        ftijnicon=indi*dnorm(-qnorm(Sijni),mutijnicon,sigmatijnicon)*fijni/dnorm(-qnorm(Sijni))+
          (1-indi)*(pnorm((qnorm(Sijni)+mutijnicon)/sigmatijnicon))
      } else {ftijnicon=1}
      posticondimean=sum(wij*postijcondimean*ftijnicon)
      ll=ll+(-ni/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1))%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}

#Estimate by etaty=0, which means conditional independence between the two sub-models
appGQbivGaucoploglik.0rhoty=function(theta,data,m)        
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  lambda=theta[14]
  sigma=theta[15]
  alpha=theta[16]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      timepoint=data$ti[data$subject==sub]
      ni=length(timepoint)
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],data$gender[data$subject==sub],
                   data$middle[data$subject==sub],data$young[data$subject==sub]),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(data$middle[data$subject==sub]),unique(data$young[data$subject==sub]))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      V11=Zi1%*%D%*%t(Zi1)+sigma^2*diag(ni)
      V22=D
      V12=Zi1%*%D
      V21=t(V12)
      mubiconyi=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)
      Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      postijcondimean=1
      if(ni>1)
      {
        for(k in 1:(ni-1))
        {
          Sijk=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*timepoint[k+1])-exp(alpha*bij[2,]*timepoint[k])))
          ftijkcon=Sijk
          postijcondimean=ftijkcon*postijcondimean  
        }
      }
      if(obsti!=timepoint[ni])
      {
        hijni=Cij*exp(alpha*bij[2,]*obsti)
        Sijni=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-exp(alpha*bij[2,]*timepoint[ni])))
        fijni=hijni*Sijni
        ftijnicon=indi*fijni+(1-indi)*Sijni
      } else {ftijnicon=1}
      posticondimean=sum(wij*postijcondimean*ftijnicon)
      ll=ll+(-ni/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1))%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}


#loglikelihood equivalent to etaty=0 as in "appGQbivGaucoploglik.0rhoty"
appGQbivGaucoploglik.uncor=function(theta,data,m)        
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  lambda=theta[14]
  sigma=theta[15]
  alpha=theta[16]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      timepoint=data$ti[data$subject==sub]
      ni=length(timepoint)
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],data$gender[data$subject==sub],
                   data$middle[data$subject==sub],data$young[data$subject==sub]),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(data$middle[data$subject==sub]),unique(data$young[data$subject==sub]))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      V11=Zi1%*%D%*%t(Zi1)+sigma^2*diag(ni)
      V22=D
      V12=Zi1%*%D
      V21=t(V12)
      mubiconyi=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)
      Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      hij=Cij*exp(alpha*bij[2,]*obsti)
      Sij=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-1))
      fij=hij*Sij
      ftijcon=indi*fij+(1-indi)*Sij
      posticondimean=sum(wij*ftijcon)
      ll=ll+(-ni/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1))%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}



#log-likelihood function for the bivariate functional t copula joint model 
#Estimate etaty by B-spline first then transfer to rhoty in t copula  
appGQbivtcoploglik.funrhoty=function(theta,data,m,tmax,neta,etaord,nu)      
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  etapar=theta[14:(13+neta)]
  lambda=theta[14+neta]
  sigma=theta[15+neta]
  alpha=theta[16+neta]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      timepoint=data$ti[data$subject==sub]
      etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=neta,norder=etaord))%*%etapar
      rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
      ni=length(timepoint)
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],data$gender[data$subject==sub],
                   data$middle[data$subject==sub],data$young[data$subject==sub]),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(data$middle[data$subject==sub]),unique(data$young[data$subject==sub]))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      V11=Zi1%*%D%*%t(Zi1)+sigma^2*diag(ni)
      V22=D
      V12=Zi1%*%D
      V21=t(V12)
      mubiconyi=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)
      Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      Wyij=qt(pnorm((matrix(rep(data$resp[data$subject==sub],m^2),ncol=m^2)-
                       matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)/sigma),df=nu)
      postijcondimean=1
      if(ni>1)
      {
        for(k in 1:(ni-1))
        {
          Sijk=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*timepoint[k+1])-exp(alpha*bij[2,]*timepoint[k])))
          Wyijk=Wyij[k,]
          mutijkcon=rhoty[k]*Wyijk
          sigmatijkcon=sqrt((nu+Wyijk^2)*(1-rhoty[k]^2)/(nu+1))
          ftijkcon=pt((qt(Sijk,df=nu)+mutijkcon)/sigmatijkcon,df=nu+1)
          postijcondimean=ftijkcon*postijcondimean  
        }
      }
      if(obsti!=timepoint[ni])
      {
        hijni=Cij*exp(alpha*bij[2,]*obsti)
        Sijni=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-exp(alpha*bij[2,]*timepoint[ni])))
        fijni=hijni*Sijni
        Wyijni=Wyij[ni,]
        mutijnicon=rhoty[ni]*Wyijni
        sigmatijnicon=sqrt((nu+Wyijni^2)*(1-rhoty[ni]^2)/(nu+1))
        ftijnicon=indi*(dt((-qt(Sijni,df=nu)-mutijnicon)/sigmatijnicon,df=nu+1)/sigmatijnicon)*fijni/dt(-qt(Sijni,df=nu),df=nu)+
          (1-indi)*pt((qt(Sijni,df=nu)+mutijnicon)/sigmatijnicon,df=nu+1)
      } else {ftijnicon=1}
      posticondimean=sum(wij*postijcondimean*ftijnicon)
      ll=ll+(-ni/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1))%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}


#Fitting by nlm estimation of observed data by integrate out bi in ni measurements, functional rhoty, linear population mean
EstY.bivcop.funrhoty.nlm=function(initpara,truepara,N,Jointdata,tmax,neta,m,etaord,nu,copula)
{
  paranum=length(truepara)
  gradmatrix=matrix(0,nrow=N,ncol=paranum)
  estmatrix=matrix(0,nrow=N,ncol=paranum)
  sdmatrix=matrix(0,nrow=N,ncol=paranum)
  samplestmatrix=matrix(0,nrow=1,ncol=paranum)
  biasmatrix=matrix(0,nrow=N,ncol=paranum)
  RMSEmatrix=matrix(0,nrow=1,ncol=paranum)
  CPmatrix=matrix(0,nrow=N,ncol=paranum)
  ECPmatrix=matrix(0,nrow=N,ncol=paranum)
  i=0;k=1;logliki=0;datasetid=NULL;iter=0
  colnames(gradmatrix)=c('gbeta11','gbeta21','gbeta31','gbeta41','gbeta51','gbeta61','gbeta12','gbeta22','gbeta32',
                        'gbeta42','gD11','gD22','gD12',rep('getapar',neta),'glambda','gsigma','galpha')
  colnames(estmatrix)=c('beta11','beta21','beta31','beta41','beta51','beta61','beta12','beta22','beta32','beta42',
                        'D11','D22','D12',rep('etapar',neta),'lambda','sigma','alpha')
  colnames(sdmatrix)=c('sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51','sdbeta61','sdbeta12','sdbeta22',
                       'sdbeta32','sdbeta42','sdD11','sdD22','sdD12',rep('sdetapar',neta),'sdlambda','sdsigma','sdalpha')
  colnames(samplestmatrix)=c('ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51','ssdbeta61','ssdbeta12',
                        'ssdbeta22','ssdbeta32','ssdbeta42','ssdD11','ssdD22','ssdD12',rep('ssdetapar',neta),
                        'ssdlambda','ssdsigma','ssdalpha')
  colnames(biasmatrix)=c('bbeta11','bbeta21','bbeta31','bbeta41','bbeta51','bbeta61','bbeta12','bbeta22','bbeta32',
                        'bbeta42','bD11','bD22','bD12',rep('betapar',neta),'blambda','bsigma','balpha')
  colnames(RMSEmatrix)=c('rmbeta11','rmbeta21','rmbeta31','rmbeta41','rmbeta51','rmbeta61','rmbeta12','rmbeta22','rmbeta32',
                        'rmbeta42','rmD11','rmD22','rmD12',rep('rmetapar',neta),'rmlambda','rmsigma','rmalpha')
  colnames(CPmatrix)=c('cpbeta11','cpbeta21','cpbeta31','cpbeta41','cpbeta51','cpbeta61','cpbeta12','cpbeta22',
                       'cpbeta32','cpbeta42','cpD11','cpD22','cpD12',rep('cpetapar',neta),'cplambda','cpsigma','cpalpha')
  colnames(ECPmatrix)=c('ecpbeta11','ecpbeta21','ecpbeta31','ecpbeta41','ecpbeta51','ecpbeta61',
                        'ecpbeta12','ecpbeta22','ecpbeta32','ecpbeta42','ecpD11','ecpD22','ecpD12',
                        rep('ecpetapar',neta),'ecplambda','ecpsigma','ecpalpha')
  while(k<=N)
  {
    i=i+1
    if(copula=="tcop")
    {
      estGHOb=nlm(f=appGQbivtcoploglik.funrhoty,p=initpara,data=Jointdata[[i]],m=m,tmax=tmax,
                  hessian=T,neta=neta,etaord=etaord,nu=nu,iterlim=1000)
    }
    if(copula=="Gaucop")
    {
      estGHOb=nlm(f=appGQbivGaucoploglik.funrhoty,p=initpara,data=Jointdata[[i]],m=m,tmax=tmax,
                  hessian=T,neta=neta,etaord=etaord,iterlim=1000)
    }
    if(copula=="uncor")
    {
      estGHOb=nlm(f=appGQbivGaucoploglik.0rhoty,p=initpara,data=Jointdata[[i]],m=m,hessian=T,iterlim=1000)
    }
    if(copula=="uncor2")
    {
      estGHOb=nlm(f=appGQbivGaucoploglik.uncor,p=initpara,data=Jointdata[[i]],m=m,hessian=T,iterlim=1000)
    }
    if(estGHOb$code==1)
    {
      estmatrix[k,]=estGHOb$estimate
      sdmatrix[k,]=sqrt(diag(solve(estGHOb$hessian)))
      gradmatrix[k,]=estGHOb$gradient
      biasmatrix[k,]=estmatrix[k,]-truepara
      logliki[k]=-estGHOb$minimum
      iter[k]=estGHOb$iter
      print(c(i,k,iter[k]))
      datasetid=c(datasetid,i)
      k=k+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  RMSEmatrix[1,]=colMeans((estmatrix-matrix(rep(truepara,N),nrow=N,byrow=T))^2)^{0.5}
  eupmatrix=estmatrix+matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  elowmatrix=estmatrix-matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  upmatrix=estmatrix+qnorm(0.975)*sdmatrix
  lowmatrix=estmatrix-qnorm(0.975)*sdmatrix
  for(j in 1:N)
  {
    ECPmatrix[j,]=as.numeric(elowmatrix[j,]<truepara&truepara<eupmatrix[j,])
    CPmatrix[j,]=as.numeric(lowmatrix[j,]<truepara&truepara<upmatrix[j,])
  }
  results=list(gradmatrix,estmatrix,sdmatrix,colMeans(estmatrix),colMeans(sdmatrix),samplestmatrix,
               colMeans(biasmatrix),RMSEmatrix,colMeans(ECPmatrix),colMeans(CPmatrix),logliki,i,datasetid,iter)
  return(results)
}



#dynamic prediction for survival probabilities
#posterior distribution of random effect bi, i.e., f(bi|ti,yi), for bivaraite functional Gaussian copula joint model
fbiposGau=function(beta1,beta2,D11,D22,D12,etapar,lambda,sigma,alpha,etaord,data,dynati,bi,m,tmax)
{
  timepoint=data$ti
  ni=length(timepoint)
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  Xi1=matrix(c(rep(1,ni),timepoint,data$treat,data$gender,data$middle,data$young),ncol=6)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  Xi2=c(unique(data$treat),unique(data$gender),unique(data$middle),unique(data$young))
  etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=length(etapar),norder=etaord))%*%etapar
  rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
  obsti=unique(data$obsti)
  if(dynati<obsti) {indi=0
  ti=dynati} else {indi=unique(data$indicator)
  ti=obsti}
  dimyi=sum(timepoint<=ti)
  V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+sigma^2*diag(dimyi)
  V22=D
  V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
  V21=t(V12)
  mubiconyi=V21%*%solve(V11)%*%(data$resp-Xi1%*%beta1)[1:dimyi,]
  Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
  Ci=lambda*exp(c(Xi2%*%beta2)+alpha*bi[1])
  ftikconyibi=1
  if(dimyi>1)
  {
    for(k in 1:(dimyi-1))
    {
      Sik=exp(-Ci/(alpha*bi[2])*(exp(alpha*bi[2]*timepoint[k+1])-exp(alpha*bi[2]*timepoint[k])))   
      Zyik=(data$resp-Xi1%*%beta1-Zi1%*%bi)[k,]/sigma
      mutikconyibi=rhoty[k]*Zyik
      sigmatikconyibi=sqrt(1-rhoty[k]^2)
      Stikconyibi=pnorm((qnorm(Sik)+mutikconyibi)/sigmatikconyibi)
      ftikconyibi=ftikconyibi*Stikconyibi 
    }
  }
  if(ti!=timepoint[dimyi])
  {
    hini=Ci*exp(alpha*bi[2]*ti)
    Sini=exp(-Ci/(alpha*bi[2])*(exp(alpha*bi[2]*ti)-exp(alpha*bi[2]*timepoint[dimyi])))   
    fini=hini*Sini
    Zyini=(data$resp-Xi1%*%beta1-Zi1%*%bi)[dimyi,]/sigma
    mutiniconyibi=rhoty[dimyi]*Zyini
    sigmatiniconyibi=sqrt(1-rhoty[dimyi]^2)
    ftiniconyibi=indi*dnorm(-qnorm(Sini),mutiniconyibi,sigmatiniconyibi)*fini/dnorm(-qnorm(Sini))+
      (1-indi)*(pnorm((qnorm(Sini)+mutiniconyibi)/sigmatiniconyibi))
  } else {ftiniconyibi=1}
  fticonyibi=ftikconyibi*ftiniconyibi
  ftijkconyibi=1
  roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
  roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
  bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
  wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
  Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
  if(dimyi>1)
  {
    for(k in 1:(dimyi-1))
    {
      Sijk=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*timepoint[k+1])-exp(alpha*bij[2,]*timepoint[k])))   
      Zyijk=(matrix(rep(data$resp,m^2),ncol=m^2)-matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[k,]/sigma
      mutijkconyibi=rhoty[k]*Zyijk
      sigmatijkconyibi=sqrt(1-rhoty[k]^2)
      Stijkconyibi=pnorm((qnorm(Sijk)+mutijkconyibi)/sigmatijkconyibi)
      ftijkconyibi=ftijkconyibi*Stijkconyibi 
    }
  }
  if(ti!=timepoint[dimyi])
  {
    hijni=Cij*exp(alpha*bij[2,]*ti)
    Sijni=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*ti)-exp(alpha*bij[2,]*timepoint[dimyi])))   
    fijni=hijni*Sijni
    Zyijni=(matrix(rep(data$resp,m^2),ncol=m^2)-matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[dimyi,]/sigma
    mutijniconyibi=rhoty[dimyi]*Zyijni
    sigmatijniconyibi=sqrt(1-rhoty[dimyi]^2)
    ftijniconyibi=indi*dnorm(-qnorm(Sijni),mutijniconyibi,sigmatijniconyibi)*fijni/dnorm(-qnorm(Sijni))+
      (1-indi)*(pnorm((qnorm(Sijni)+mutijniconyibi)/sigmatijniconyibi))
  } else {ftijniconyibi=1}
  fticonyi=sum(ftijkconyibi*ftijniconyibi*wij)
  fbiconyiti=fticonyibi*dmvnorm(bi,c(mubiconyi),Sigmabiconyi)/fticonyi
  return(-fbiconyiti)
}

#posterior distribution of random effect bi, i.e., f(bi|ti,yi), for bivaraite functional t copula joint model
fbipostcop=function(beta1,beta2,D11,D22,D12,etapar,lambda,sigma,alpha,etaord,data,dynati,bi,m,tmax,nu)
{
  timepoint=data$ti
  ni=length(timepoint)
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  Xi1=matrix(c(rep(1,ni),timepoint,data$treat,data$gender,data$middle,data$young),ncol=6)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  Xi2=c(unique(data$treat),unique(data$gender),unique(data$middle),unique(data$young))
  etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=length(etapar),norder=etaord))%*%etapar
  rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
  obsti=unique(data$obsti)
  if(dynati<obsti) {indi=0
  ti=dynati} else {indi=unique(data$indicator)
  ti=obsti}
  dimyi=sum(timepoint<=ti)
  V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+sigma^2*diag(dimyi)
  V22=D
  V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
  V21=t(V12)
  mubiconyi=V21%*%solve(V11)%*%(data$resp-Xi1%*%beta1)[1:dimyi,]
  Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
  Ci=lambda*exp(c(Xi2%*%beta2)+alpha*bi[1])
  ftikconyibi=1
  if(dimyi>1)
  {
    for(k in 1:(dimyi-1))
    {
      Sik=exp(-Ci/(alpha*bi[2])*(exp(alpha*bi[2]*timepoint[k+1])-exp(alpha*bi[2]*timepoint[k])))   
      Wyik=qt(pnorm((data$resp-Xi1%*%beta1-Zi1%*%bi)[k,]/sigma),df=nu)
      mutikconyibi=rhoty[k]*Wyik
      sigmatikconyibi=sqrt((nu+Wyik^2)*(1-rhoty[k]^2)/(nu+1))
      Stikconyibi=pt((qt(Sik,df=nu)+mutikconyibi)/sigmatikconyibi,df=nu+1)
      ftikconyibi=ftikconyibi*Stikconyibi 
    }
  }
  if(ti!=timepoint[dimyi])
  {
    hini=Ci*exp(alpha*bi[2]*ti)
    Sini=exp(-Ci/(alpha*bi[2])*(exp(alpha*bi[2]*ti)-exp(alpha*bi[2]*timepoint[dimyi])))   
    fini=hini*Sini
    Wyini=qt(pnorm((data$resp-Xi1%*%beta1-Zi1%*%bi)[dimyi,]/sigma),df=nu)
    mutiniconyibi=rhoty[dimyi]*Wyini
    sigmatiniconyibi=sqrt((nu+Wyini^2)*(1-rhoty[dimyi]^2)/(nu+1))
    ftiniconyibi=(1-indi)*pt((qt(Sini,df=nu)+mutiniconyibi)/sigmatiniconyibi,df=nu+1)+
      indi*(dt((-qt(Sini,df=nu)-mutiniconyibi)/sigmatiniconyibi,df=nu+1)/sigmatiniconyibi)*fini/dt(-qt(Sini,df=nu),df=nu)
  } else {ftiniconyibi=1}
  fticonyibi=ftikconyibi*ftiniconyibi
  ftijkconyibi=1
  roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
  roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
  bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
  wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
  Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
  if(dimyi>1)
  {
    for(k in 1:(dimyi-1))
    {
      Sijk=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*timepoint[k+1])-exp(alpha*bij[2,]*timepoint[k])))   
      Wyijk=qt(pnorm((matrix(rep(data$resp,m^2),ncol=m^2)-matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[k,]/sigma),df=nu)
      mutijkconyibi=rhoty[k]*Wyijk
      sigmatijkconyibi=sqrt((nu+Wyijk^2)*(1-rhoty[k]^2)/(nu+1))
      Stijkconyibi=pt((qt(Sijk,df=nu)+mutijkconyibi)/sigmatijkconyibi,df=nu+1)
      ftijkconyibi=ftijkconyibi*Stijkconyibi 
    }
  }
  if(ti!=timepoint[dimyi])
  {
    hijni=Cij*exp(alpha*bij[2,]*ti)
    Sijni=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*ti)-exp(alpha*bij[2,]*timepoint[dimyi])))   
    fijni=hijni*Sijni
    Wyijni=qt(pnorm((matrix(rep(data$resp,m^2),ncol=m^2)-matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[dimyi,]/sigma),df=nu)
    mutijniconyibi=rhoty[dimyi]*Wyijni
    sigmatijniconyibi=sqrt((nu+Wyijni^2)*(1-rhoty[dimyi]^2)/(nu+1))
    ftijniconyibi= (1-indi)*pt((qt(Sijni,df=nu)+mutijniconyibi)/sigmatijniconyibi,df=nu+1)+
      indi*(dt((-qt(Sijni,df=nu)-mutijniconyibi)/sigmatijniconyibi,df=nu+1)/sigmatijniconyibi)*fijni/dt(-qt(Sijni,df=nu),df=nu)
  } else {ftijniconyibi=1}
  fticonyi=sum(ftijkconyibi*ftijniconyibi*wij)
  fbiconyiti=fticonyibi*dmvnorm(bi,c(mubiconyi),Sigmabiconyi)/fticonyi
  return(-fbiconyiti)
}



#posterior distribution of random effect bi, i.e., f(bi|ti,yi), for regular joint model (conditional independence)
#This is essential rho(t)=0 under the bivariate functional Gaussain copula joint model
fbiposuncor=function(beta1,beta2,D11,D22,D12,lambda,sigma,alpha,data,dynati,bi,m)
{
  timepoint=data$ti
  ni=length(timepoint)
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  Xi1=matrix(c(rep(1,ni),timepoint,data$treat,data$gender,data$middle,data$young),ncol=6)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  Xi2=c(unique(data$treat),unique(data$gender),unique(data$middle),unique(data$young))
  obsti=unique(data$obsti)
  if(dynati<obsti) {indi=0
  ti=dynati} else {indi=unique(data$indicator)
  ti=obsti}
  dimyi=sum(timepoint<=ti)
  V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+sigma^2*diag(dimyi)
  V22=D
  V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
  V21=t(V12)
  mubiconyi=V21%*%solve(V11)%*%(data$resp-Xi1%*%beta1)[1:dimyi,]
  Sigmabiconyi=V22-V21%*%solve(V11)%*%V12
  Ci=lambda*exp(c(Xi2%*%beta2)+alpha*bi[1])
  hi=Ci*exp(alpha*bi[2]*ti)
  Si=exp(-Ci/(alpha*bi[2])*(exp(alpha*bi[2]*ti)-1))
  fi=hi*Si
  fticon=indi*fi+(1-indi)*Si
  postijcondimean=0
  roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
  roti=eigen(Sigmabiconyi)$vector%*%diag(sqrt(eigen(Sigmabiconyi)$value))%*%roti1
  bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubiconyi)
  wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
  Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
  hij=Cij*exp(alpha*bij[2,]*ti)
  Sij=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*ti)-1))
  fij=hij*Sij
  ftijcon=indi*fij+(1-indi)*Sij
  postijcondimean=sum(ftijcon*wij)
  fbicon=fticon*dmvnorm(bi,c(mubiconyi),Sigmabiconyi)/postijcondimean
  return(-fbicon)
}


#maximise f(bi|ti,yi) to obtain hat^{bi} for bivaraite copula joint model
dynabi=function(beta1,beta2,D11,D22,D12,etapar,lambda,sigma,alpha,etaord,data,dynati,m,tmax,copula,nu,mtool)
{
  timepoint=data$ti
  ni=length(timepoint)
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  Xi1=matrix(c(rep(1,ni),timepoint,data$treat,data$gender,data$middle,data$young),ncol=6)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  Xi2=c(unique(data$treat),unique(data$gender),unique(data$middle),unique(data$young))
  ti=min(dynati,unique(data$obsti))
  dimyi=sum(timepoint<=ti)
  V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+sigma^2*diag(dimyi)
  V22=D
  V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
  V21=t(V12)
  mubiconyi=V21%*%solve(V11)%*%(data$resp-Xi1%*%beta1)[1:dimyi,]
  if(copula=="Gaucop")
  {
    if(mtool=="nlm")
    {
    bihat=nlm(fbiposGau,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,
              sigma=sigma,alpha=alpha,etaord=etaord,data=data,dynati=dynati,m=m,tmax=tmax,hessian=T,iterlim=10000)$est
    } else{
    bihat=optim(fbiposGau,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,
              sigma=sigma,alpha=alpha,etaord=etaord,data=data,dynati=dynati,m=m,tmax=tmax,hessian=T,control=list(maxit=20000))$par
    }
  }
  if(copula=="tcop")
  {
    if(mtool=="nlm")
    {
    bihat=nlm(fbipostcop,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,
              sigma=sigma,alpha=alpha,etaord=etaord,data=data,dynati=dynati,m=m,tmax=tmax,nu=nu,hessian=T,iterlim=10000)$est
    } else{
    bihat=optim(fbipostcop,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,
        sigma=sigma,alpha=alpha,etaord=etaord,data=data,dynati=dynati,m=m,tmax=tmax,nu=nu,hessian=T,control=list(maxit=20000))$par
    }
  }
  if(copula=="uncor") #This is essential rho(t)=0 under the bivariate functional Gaussain copula joint model
  {
    if(mtool=="nlm")
    {
    bihat=nlm(fbiposuncor,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,lambda=lambda,
              sigma=sigma,alpha=alpha,data=data,dynati=dynati,m=m,hessian=T,iterlim=10000)$est
    } else{
      optim(fbiposuncor,p=c(mubiconyi),beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,lambda=lambda,
          sigma=sigma,alpha=alpha,data=data,dynati=dynati,m=m,hessian=T,control=list(maxit=20000))$par
    }
  }
  results=list(bihat,c(mubiconyi))
  return(results)
}


#predict survival probabilities
predSur=function(beta1,beta2,D11,D22,D12,etapar,lambda,sigma,alpha,etaord,data,dynati,predinterv,m,tmax,acc,copula,nu,mtool)
{
  timepoint=data$ti
  ni=length(timepoint)
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  Xi1=matrix(c(rep(1,ni),timepoint,data$treat,data$gender,data$middle,data$young),ncol=6)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  Xi2=c(unique(data$treat),unique(data$gender),unique(data$middle),unique(data$young))
  ti=seq(min(dynati,unique(data$obsti)),predinterv,by=acc)
  dimyi=sum(timepoint<=min(dynati,unique(data$obsti)))
  outbi=dynabi(beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,sigma=sigma,
               alpha=alpha,etaord=etaord,data=data,dynati=dynati,m=m,tmax=tmax,copula=copula,nu=nu,mtool=mtool)
  hatbi=outbi[[1]]
  hatbiconyi=outbi[[2]]
  Ci=lambda*exp(c(Xi2%*%beta2)+alpha*hatbi[1])
  Sicon=exp(-Ci/(alpha*hatbi[2])*(exp(alpha*hatbi[2]*ti)-exp(alpha*hatbi[2]*timepoint[dimyi]))) 
  if(copula=="Gaucop")
  {
    etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=length(etapar),norder=etaord))%*%etapar
    rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
    Zyi=(data$resp-Xi1%*%beta1-Zi1%*%hatbi)[dimyi,]/sigma
    muticonyibi=rhoty[dimyi]*Zyi
    sigmaticonyibi=sqrt(1-rhoty[dimyi]^2)
    predi=pnorm((qnorm(Sicon)+muticonyibi)/sigmaticonyibi)/pnorm((qnorm(Sicon[1])+muticonyibi)/sigmaticonyibi)
  }
  if (copula=="tcop")
  {
    etaty=eval.basis(timepoint,create.bspline.basis(c(0,tmax+0.2),nbasis=length(etapar),norder=etaord))%*%etapar
    rhoty=(exp(2*etaty)-1)/(exp(2*etaty)+1)
    Wyi=qt(pnorm((data$resp-Xi1%*%beta1-Zi1%*%hatbi)[dimyi,]/sigma),df=nu)
    muticonyibi=rhoty[dimyi]*Wyi
    sigmaticonyibi=sqrt((nu+Wyi^2)*(1-rhoty[dimyi]^2)/(nu+1))
    predi=pt((qt(Sicon,df=nu)+muticonyibi)/sigmaticonyibi,df=nu+1)/pt((qt(Sicon[1],df=nu)+muticonyibi)/sigmaticonyibi,df=nu+1)
  }
  if(copula=="uncor")
  {
    predi=Sicon/Sicon[1]
  }
  results=list(ti,Sicon/Sicon[1],predi,dimyi,hatbi,hatbiconyi)
  return(results)
}

#calculate AUC and PE for censored data

#function for providing some imputs for "dynaAUC.cen" and "dynaPE.cen" functions later on.
dynasurgroup.cen=function(Data,beta1,beta2,D11,D22,D12,etapar,lambda,sigma,alpha,etaord,dynati,predinterv,m,
                          tmax,copula,nu,mtool)
{
  allprob=0
  j=0
  i=1
  ind=0
  obst=0
  for(i in unique(Data$subj))
  {
    Datai=Data[Data$subj==i,]
    obsti=unique(Datai$obsti)
    indi=unique(Datai$indicator)
    if(obsti>dynati)
    {
      j=j+1
      allprob[j]=predSur(beta1=beta1,beta2=beta2,D11=D11,D22=D22,D12=D12,etapar=etapar,lambda=lambda,
                            sigma=sigma,alpha=alpha,etaord=etaord,data=Datai,dynati=dynati,
                            predinterv=predinterv,m=m,tmax=tmax,acc=predinterv-dynati,copula=copula,nu=nu,mtool=mtool)[[3]][-1]
      ind[j]=indi;obst[j]=obsti
    }
  }
  results=list(j,allprob,ind,obst)
  return(results)
}

#function for calculate AUC
dynaAUC.cen=function(probs,ind,obst,dynati,predinterv)
{
  upcount=0;lowcount=0
  for(i in 1:length(probs))
  {
    if(obst[i]<=predinterv&obst[i]>dynati&ind[i]==1)
    {
      upcount=sum((obst>predinterv)*(probs[i]<probs))+sum((obst>obst[i]&obst<=predinterv&ind==0)*(probs[i]<probs)*probs)+
        upcount
      lowcount=sum(obst>predinterv)+sum((obst>obst[i]&obst<=predinterv&ind==0)*probs)+lowcount
    }
    if(obst[i]<=predinterv&obst[i]>dynati&ind[i]==0)
    {
      upcount=sum((obst>predinterv)*(probs[i]<probs)*(1-probs[i]))+
        sum((obst>obst[i]&obst<=predinterv&ind==0)*(probs[i]<probs)*(1-probs[i])*probs)+upcount
      lowcount=sum((obst>predinterv)*(1-probs[i]))+sum((obst>obst[i]&obst<=predinterv&ind==0)*(1-probs[i])*probs)+
        lowcount
    }
  }
  return(upcount/lowcount)
}


#function for calculate PE
dynaPE.cen=function(probs,ind,obst,dynati,predinterv)
{
     upcount=sum((obst>predinterv)*(1-probs)^2+ind*(obst<predinterv)*(0-probs)^2+(1-ind)*(obst<predinterv)*
                   (probs*(1-probs)^2+(1-probs)*(0-probs)^2))
    lowcount=sum(obst>dynati)
  return(upcount/lowcount)
}
