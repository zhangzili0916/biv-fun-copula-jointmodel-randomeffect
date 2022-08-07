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



GaussianHermite <- function(r){
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
