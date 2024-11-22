##as opposed to sim_md.R, here most of the items load on the first dimension.
simfun<-function(r) {
    N<-5000
    library(MASS)
    ##theta
    th<-mvrnorm(N,mu=rep(0,2),Sigma=matrix(c(1,r,r,1),2,2))
    ##skills
    logit<-function(x,a=1.7) 1/(1+exp(-a*x))
    s1<-rbinom(N,1,logit(th[,1]+.5))
    s2<-rbinom(N,1,logit(th[,1]))
    s3<-rbinom(N,1,logit(th[,1]-.5))
    s4<-rbinom(N,1,logit(th[,1]-1))
    s5<-rbinom(N,1,logit(th[,2]))
    s6<-rbinom(N,1,logit(th[,2]-1))
    sk<-cbind(s1,s2,s3,s4,s5,s6)
    
    ##qmatrix
    qm<-matrix(
        c(1,0,0,0,0,0,
          0,1,0,0,0,0,
          0,0,1,0,0,0,
          0,0,0,1,0,0,
          1,1,0,0,0,0,
          0,1,1,0,0,0,
          0,0,1,1,0,0,
          1,0,1,0,0,0,
          1,0,0,1,0,0,
          1,1,1,0,0,0,
          0,1,1,1,0,0,
          1,1,1,1,0,0,
          0,0,0,1,0,0,
          0,0,0,0,1,0,
          0,0,0,0,0,1,
          0,0,0,0,1,1
          ),16,6,byrow=TRUE)
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

    th2.index<-which(qm0[,5]>0 | qm0[,6]>0)
    th1.index<-which(!(qm0[,5]>0 | qm0[,6]>0))

    L<-list(as.matrix(resp[,th1.index]),p.true[,th1.index],p.irt[,th1.index],p.cdm[,th1.index])
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    c1<-cor(z)

    L<-list(as.matrix(resp[,th2.index]),p.true[,th2.index],p.irt[,th2.index],p.cdm[,th2.index])
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    c2<-cor(z)
    list(c1,c2)
}

rs<-sort(runif(50,-1,1))
L<-lapply(rs,simfun)

th1<-lapply(L,function(x) x[[1]])
th2<-lapply(L,function(x) x[[2]])

pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...,lwd=2)
}
##
f<-function(L) {
    irt<-sapply(L,function(x) x[2,3])
    cdm<-sapply(L,function(x) x[2,4])
    lines(pf(rs,irt))
    lines(pf(rs,cdm,col='red'))
}
par(mgp=c(2,1,0),mfrow=c(1,2))
plot(NULL,xlim=c(-1,1),ylim=0:1,xlab='corrrelation between dimensions',ylab='correlation between estimates and true probabilities')
f(th1)
plot(NULL,xlim=c(-1,1),ylim=0:1,xlab='corrrelation between dimensions',ylab='correlation between estimates and true probabilities')
f(th2)
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"))

