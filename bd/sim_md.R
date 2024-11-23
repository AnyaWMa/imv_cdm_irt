##equal share of items load on two dimensions
source("00funs.R")
simfun<-function(r) {
    N<-5000
    library(MASS)
    ##theta
    th<-mvrnorm(N,mu=rep(0,2),Sigma=matrix(c(1,r,r,1),2,2))
    ##skills
    logit<-function(x,a=1.7) 1/(1+exp(-a*x))
    s1<-rbinom(N,1,logit(th[,1]+1))
    s2<-rbinom(N,1,logit(th[,1]))
    s3<-rbinom(N,1,logit(th[,1]-1))
    s4<-rbinom(N,1,logit(th[,2]+1))
    s5<-rbinom(N,1,logit(th[,2]))
    s6<-rbinom(N,1,logit(th[,2]-1))
    sk<-cbind(s1,s2,s3,s4,s5,s6)
    
    ##qmatrix
    ##between-item MD
    ## qm<-matrix(
    ##     c(1,0,0,0,0,0,
    ##       0,1,0,0,0,0,
    ##       0,0,1,0,0,0,
    ##       1,1,0,0,0,0,
    ##       0,1,1,0,0,0,
    ##       1,0,1,0,0,0,
    ##       1,1,1,0,0,0,
    ##       0,0,0,1,0,0,
    ##       0,0,0,0,1,0,
    ##       0,0,0,0,0,1,
    ##       0,0,0,1,1,0,
    ##       0,0,0,0,1,1,
    ##       0,0,0,1,0,1,
    ##       0,0,0,1,1,1
    ##       ),14,6,byrow=TRUE)
    ##within-item MD
    S<-TRUE
    while (S) {
        qm<-matrix(rbinom(30*6,1,.5),30,6)
        S<-any(rowMeans(qm)==0)
    }        
    qm0<-qm

    ##response probabilities for with g=s=0.1
    pL<-respL<-list()
    g<-runif(nrow(qm0),min=0,max=.35)
    s<-runif(nrow(qm0),min=0,max=.35)
    for (i in 1:nrow(qm0)) {
        ii<-which(qm[i,]==1)
        z<-sk[,ii,drop=FALSE]
        rm<-rowMeans(z)
        p.tmp<-ifelse(rm==1,1-s[i],g[i])
        respL[[i]]<-rbinom(N,1,p.tmp)
        pL[[i]]<-p.tmp
    }
    resp<-do.call("cbind",respL)
    p.true<-do.call("cbind",pL)
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    
    p.irt<-irt.pr(resp)
    p.cdm<-cdm.pr(resp,qm0)
    
    L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    true<-cor(z)[2,3:4]
    
    ##cv imv values
    om<-oos.compare(resp,qm0,nfolds=5)
    list(t=true,om=om)

}

rs<-sort(runif(250,-1,1))
library(parallel)
L<-mclapply(rs,simfun,mc.cores=4)

tr<-lapply(L,function(x) x$t)
tr<-do.call("rbind",tr)
om<-sapply(L,function(x) x$om)

par(mgp=c(2,1,0),mfrow=c(1,2))
##
irt<-tr[,1]
cdm<-tr[,2]
plot(NULL,xlim=c(-1,1),ylim=0:1,xlab='correlation between dimensions',ylab='correlation between estimates and true probabilities')
pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...,lwd=2)
}
lines(pf(rs,irt))
lines(pf(rs,cdm,col='red'))
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"))
##
plot(NULL,xlim=c(-1,1),ylim=c(-.05,.15),xlab='correlation between dimensions',ylab='IMV')
lines(pf(rs,om))

