source("00funs.R")


################################3
##simulate data with irt

irt.sim<-function(a,N=1000) {
    th<-rnorm(N)
    b<-sort(rnorm(10,sd=.3))
    th.mat<-matrix(th,length(th),length(b),byrow=FALSE)
    b.mat<-matrix(b,length(th),length(b),byrow=TRUE)
    k<-th.mat-b.mat
    p<-1/(1+exp(-k))
    resp<-p
    for (i in 1:ncol(resp)) resp[,i]<-rbinom(nrow(resp),1,p[,i])
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    p.true<-p
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    ##make up q matrix
    sk<-rnorm(3)
    p<-outer(b,sk,'-')
    p<-apply(p,2,function(x) 1/(1+exp(-a*x)))
    qm<-p
    S<-TRUE
    while (S) {
        for (i in 1:nrow(qm)) qm[i,]<-rbinom(ncol(qm),1,p[i,])
        S<-any(c(colMeans(qm),rowMeans(qm))==0)
    }            
    
    p.irt<-irt.pr(resp)
    p.cdm<-cdm.pr.marg(resp,qm)
    
    L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    true<-cor(z)[2,3:4]
    
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true)
    list(t=true,om=om)
}

a<-sort(runif(100,min=.1,max=2))
library(parallel)
L<-mclapply(a,irt.sim,mc.cores=10)


tr<-lapply(L,function(x) x$t)
tr<-do.call("rbind",tr)
om<-sapply(L,function(x) x$om)



#pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario2.pdf",width=6,height=3)
par(mgp=c(2,1,0),mfrow=c(1,2),mar=c(3,3,1,1),oma=rep(.5,4))
##
irt<-tr[,1]
cdm<-tr[,2]
plot(NULL,xlim=c(0,2),ylim=0:1,xlab='a',ylab='r(true,est)')
pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...,lwd=2)
}
lines(pf(a,irt))
lines(pf(a,cdm,col='red'))
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"),title="est")
##
plot(NULL,xlim=c(0,2.5),ylim=c(-.25,0),xlab='a',ylab='IMV(IRT,CDM)')
lines(pf(a,om))
#dev.off()


irt.sim<-function(n,N=1000) {
    th<-rnorm(N)
    b<-sort(rnorm(50,sd=1))
    th.mat<-matrix(th,length(th),length(b),byrow=FALSE)
    b.mat<-matrix(b,length(th),length(b),byrow=TRUE)
    k<-th.mat-b.mat
    p<-1/(1+exp(-k))
    resp<-p
    for (i in 1:ncol(resp)) resp[,i]<-rbinom(nrow(resp),1,p[,i])
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    p.true<-p
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')

    ##make up q matrix
    out<-list()
    a<-c(0,runif(1,min=0,max=1.5))
    qmL<-list()
    for (ii in 1:length(a)) {
        sk<-rnorm(n)
        p<-outer(b,sk,'-')
        p<-apply(p,2,function(x) 1/(1+exp(-a[ii]*x)))
        qm<-p
        S<-TRUE
        while (S) {
            for (i in 1:n) qm[,i]<-rbinom(nrow(qm),1,p[,i])
            S<-any(c(colMeans(qm),rowMeans(qm))==0)
        }            
        qmL[[ii]]<-qm
    }

    om<-oos.compare.newresp.2q(resp,qmL[[1]],qmL[[2]],truep=p.true)
    c(n,a[2],om)
}
library(parallel)
nn<-c(rep(3,10),rep(6,10))
z<-mclapply(nn,irt.sim,mc.cores=10,N=2500)



