source("00funs.R")


################################3
##simulate data with irt

irt.sim<-function(N) {
    th<-rnorm(N)
    b<-sort(rnorm(30))
    th.mat<-matrix(th,length(th),length(b),byrow=FALSE)
    b.mat<-matrix(b,length(th),length(b),byrow=TRUE)
    k<-th.mat-b.mat
    p<-1/(1+exp(-k))
    resp<-p
    for (i in 1:ncol(resp)) resp[,i]<-rbinom(nrow(resp),1,p[,i])
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    p.true<-p
    ##make up q matrix
    sk<-rnorm(5)
    p<-outer(b,sk,'-')
    p<-apply(p,2,function(x) 1/(1+exp(-x)))
    qm0<-p
    for (i in 1:nrow(qm0)) {
        S<-0
        while (S==0) {
            qm0[i,]<-rbinom(ncol(qm0),1,p[i,])
            S<-sum(qm0[i,])
        }            
    }         
    ##
    ## p.irt<-irt.pr(resp)
    ## p.cdm<-cdm.pr(resp,qm0)
    ## coors<-list()
    ## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
    ## ##note that third row >> fourth row!
    ## sapply(coors,function(x) x[1,])
    ##cv imv values
    oos.compare(resp,qm0,nfolds=4)
}

N<-runif(5,3,4.5)
N<-round(10^N)
om<-sapply(N,irt.sim)

plot(N,om,pch=19); abline(h=0)
