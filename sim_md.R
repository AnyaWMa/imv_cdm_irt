##equal share of items load on two dimensions
simfun<-function(r) {
    N<-1000
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
    qm<-matrix(
        c(1,0,0,0,0,0,
          0,1,0,0,0,0,
          0,0,1,0,0,0,
          1,1,0,0,0,0,
          0,1,1,0,0,0,
          1,0,1,0,0,0,
          1,1,1,0,0,0,
          0,0,0,1,0,0,
          0,0,0,0,1,0,
          0,0,0,0,0,1,
          0,0,0,1,1,0,
          0,0,0,0,1,1,
          0,0,0,1,0,1,
          0,0,0,1,1,1
          ),14,6,byrow=TRUE)
    qm0<-rbind(qm,qm,qm)
    
    ##response probabilities for with g=s=0.1
    pL<-respL<-list()
    for (j in 1:3) {
        p<-resp<-list()
        g<-s<-0.1*j
        for (i in 1:nrow(qm)) {
            ii<-which(qm[i,]==1)
            z<-sk[,ii,drop=FALSE]
            rm<-rowMeans(z)
            p.tmp<-ifelse(rm==1,1-s,g)
            resp[[i]]<-rbinom(N,1,p.tmp)
            p[[i]]<-p.tmp
        }
        respL[[j]]<-do.call("cbind",resp)
        pL[[j]]<-do.call("cbind",p)
    }
    resp<-do.call("cbind",respL)
    p.true<-do.call("cbind",pL)
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp))
    
    source("00funs.R")
    p.irt<-irt.pr(resp)
    p.cdm<-cdm.pr(resp,qm0)
    
    L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    cor(z)
}

rs<-sort(runif(100,-1,1))
L<-lapply(rs,simfun)

irt<-sapply(L,function(x) x[2,3])
cdm<-sapply(L,function(x) x[2,4])
par(mgp=c(2,1,0))
plot(NULL,xlim=c(-1,1),ylim=0:1,xlab='corrrelation between dimensions',ylab='correlation between estimates and true probabilities')
pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...,lwd=2)
}
lines(pf(rs,irt))
lines(pf(rs,cdm,col='red'))
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"))

