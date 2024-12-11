##equal share of items load on two dimensions
source("00funs.R")
simfun<-function(r,N=1000) {
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
    
    S<-TRUE
    while (S) {
        qm<-matrix(rbinom(50*6,1,.5),50,6)
        ##block out loadings across thetas
        for (i in 1:25) qm[i,4:6]<-0
        for (i in 26:50) qm[i,1:3]<-0
        S<-any(c(colMeans(qm),rowMeans(qm))==0)
    }        

    ##response probabilities for with g=s=0.1
    pL<-respL<-list()
    g<-runif(nrow(qm),min=0,max=.35)
    s<-runif(nrow(qm),min=0,max=.35)
    for (i in 1:nrow(qm)) {
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
    p.cdm<-cdm.pr.marg(resp,qm)
    
    ## L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    ## L<-lapply(L,as.numeric)
    ## z<-do.call("cbind",L)
    ## true<-cor(z)[2,3:4]
    ## coors<-list()
    ## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
    ## ss<-sapply(coors,function(x) x[1,])
    ## true<-rowMeans(ss)[3:4]
    
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true)
                                        #list(t=true,om=om)
    om
}

rs<-sort(runif(500,-1,1))
library(parallel)
L<-mclapply(rs,simfun,mc.cores=10)

##save.image(file="scenario3.Rdata")


pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario3.pdf",width=6,height=3)
par(mgp=c(2,1,0),mfrow=c(1,2),mar=c(3,3,1,1),oma=rep(.5,4))
plot(NULL,xlim=c(-1,1),ylim=c(0,.6),xlab=expression(rho),ylab='r(true,est)')
pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...,lwd=2)
}
f<-function(out,...) {
    z<-do.call("rbind",out)
    irt<-z[,4]
    cdm<-z[,5]
    irt.p<-z[,6]
    cdm.p<-z[,7]
    pf(rs,irt,col='blue',...)
    pf(rs,cdm,col='red',...)
    pf(rs,irt.p,col='blue',...,lty=2)
    pf(rs,cdm.p,col='red',...,lty=2)
}
#f(out1)
f(L)
legend("topright",bty='n',lty=c(1,1,2,2),col=c("blue","red","blue","red"),cex=.7,ncol=2,
       c("(resp,IRT)","(resp,CDM)","(True,IRT)","(True,CDM)"))
#####
plot(NULL,xlim=c(c(-1,1)),ylab="IMV",xlab='a',ylim=c(0,.2))
f<-function(out,...) {
    om<-do.call("rbind",out)
    abline(h=0)
    pf(rs,om[,1],...)
    pf(rs,om[,2],col='blue',...,lty=2)
    pf(rs,om[,3],col='red',...,lty=2)
}
#f(out1,lty=1)
f(L)
legend("topright",bty='n',lty=c(1,2,2),col=c("black","blue","red"),cex=.7,
       c("(IRT,CDM)","(IRT,True)","(CDM,True)")
       )
##
dev.off()
