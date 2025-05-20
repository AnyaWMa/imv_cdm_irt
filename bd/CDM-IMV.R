# wenchao's code
imv<-function (resp, pv1, pv2, eps = 1e-06) 
{
  pv1 <- ifelse(pv1 < eps, eps, pv1)
  pv2 <- ifelse(pv2 < eps, eps, pv2)
  pv1 <- ifelse(pv1 > 1 - eps, 1 - eps, pv1)
  pv2 <- ifelse(pv2 > 1 - eps, 1 - eps, pv2)
  ll <- function(x, p) {
    z <- log(p) * resp + log(1 - p) * (1 - resp)
    z <- sum(z)/length(x)
    exp(z)
  }
  loglik1 <- ll(resp, pv1)
  loglik2 <- ll(resp, pv2)
  getcoins <- function(a) {
    f <- function(p, a) abs(p * log(p) +
                              (1 - p) * log(1 - p) - log(a))
    nlminb(0.5, f, lower = 0.001, upper = 0.999, a = a)$par
  }
  c1 <- getcoins(loglik1)
  c2 <- getcoins(loglik2)
  ew <- function(p1, p0) (p1 - p0)/p0
  imv <- ew(c2, c1)
  imv
}
set.seed(123445)

N <- 500
Q <- sim30GDINA$simQ
J <- nrow(Q)
K <- ncol(Q)
library(mirt)
library(GDINA)
rr <- c(0,.25,.5,.75)
imv.vec <- matrix(NA,3,4)
for(i in 1:4){
  gs <- matrix(runif(J*2,0,0.3),J,2)
  cutoffs <- qnorm(c(1:K)/(K+1))
  m <- rep(0,K)
  
  vcov <- matrix(rr[i],K,K)
  diag(vcov) <- 1
  sim <- simGDINA(N,Q,gs.parm = gs, att.dist = "mvnorm",
                    mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))
  # test data - same p for each cell
  sim.test <- simGDINA(N,Q,catprob.parm = sim$catprob.parm,attribute = sim$attribute)
  cdmfit <- GDINA(sim$dat,Q,mono.constraint = TRUE)
  
  irtfit<-mirt(data.frame(sim$dat),1,'2PL')
  
  # predictions from irt
  
  th.est<-fscores(irtfit)
  est<-coef(irtfit,simplify=TRUE,IRTpars=TRUE)$items
  k<-outer(th.est[,1],est[,2],'-')
  k<-matrix(est[,1],nrow=N,ncol=J,byrow=TRUE)*k
  irtp<-1/(1+exp(-k))
  
  #prediction from cdm based on MAP
  gr <- GDINA:::matchMatrix(as.matrix(attributepattern(K)),as.matrix(personparm(cdmfit,"MAP")[,1:K]))
  cdmp <- t(cdmfit$LC.prob[,gr])
  
  
  #prediction from true cdm
  truep <- sim$LCprob.parm[sim$att.group,]
  
  # prediction from cdm based on mp
  lg <- cdmfit$technicals$reduced.LG
  d <- cdmfit$delta.parm
  whichjk <- apply(Q,1,function(x)which(x==1))
  mp <- personparm(cdmfit,"mp")
  pmp <- matrix(NA,N,J)
  for(ii in 1:nrow(mp)){
    x <- mp[ii,]
    x[x<1e-10] <- 1e-10
    p <- vector(length=J)
    for(jj in 1:J){
      if(length(whichjk[[jj]])==1){
        aj <- c(1,x[whichjk[[jj]]])
      }else{
        aj <- exp(lg[[jj]]%*%log(x[whichjk[[jj]]]))
      } 
      aj[1] <- 1
      p[jj] <- sum(aj * d[[jj]])
    }
    pmp[ii,] <- p
  }
  
  # calculating IMV - CDM prediction vs IRT prediction
  imv.vec[1,i] <- imv(as.numeric(sim.test$dat),
      pv1=as.numeric(irtp),
      pv2=as.numeric(cdmp))
  
  # calculating IMV - CDM MAP prediction vs CDM mp prediction
  
  imv.vec[2,i] <- imv(as.numeric(sim.test$dat),
                      pv1=as.numeric(cdmp),
                      pv2=as.numeric(pmp))
  
  # calculating IMV - CDM MAP prediction vs true prediction
  
  imv.vec[3,i] <- imv(as.numeric(sim.test$dat),
                      pv1=as.numeric(cdmp),
                      pv2=as.numeric(truep))
}
imv.vec



