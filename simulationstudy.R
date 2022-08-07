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

