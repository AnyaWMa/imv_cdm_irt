##equal share of items load on two dimensions
remotes::install_github("hansorlee/irwpkg")
source("00funs.R")
simfun<-function(r,N=500,bound=NULL) {
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
    ##
    S<-TRUE
    while (S) {
        I_k <- diag(6)
        
        random_matrix<-matrix(rbinom(44*6,1,.5),44,6)
        ##block out loadings across thetas
        for (i in 1:22) random_matrix[i,4:6]<-0
        for (i in 23:44) random_matrix[i,1:3]<-0
        combined_matrix <- rbind(I_k, random_matrix)
        
        qm <- combined_matrix[sample(nrow(combined_matrix)), ]
        S <-any(c(colMeans(qm),rowMeans(qm))==0)
    }        
    ##response probabilities 
    pL<-respL<-list()
    if (is.null(bound)) {
        g<-runif(nrow(qm),min=0,max=.35)
        s<-runif(nrow(qm),min=0,max=.35)
    } else {
        g<-rep(bound,nrow(qm))
        s<-rep(bound,nrow(qm))
    }
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
    ##
    p.irt<-irt.pr(resp)
    p.cdm<-cdm.pr.marg(resp,qm)
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true)
    ##
    om
}

rs<-sort(runif(100,-1,1))
library(parallel)
#L<-mclapply(rs,simfun,mc.cores=10)

L<-list()
for (bound in c(.1,.2,.3)) L[[as.character(bound)]]<-mclapply(rs,simfun,mc.cores=10,bound=bound)

#L[["rs"]] <- rs
#save(L , file = "../simulation_data/scenario2_dina.RData")

pdf("../plots/scenario2_update.pdf",width=6,height=3)
#par(mgp=c(2,1,0),mfrow=c(1,2),mar=c(3,3,1,1),oma=rep(.5,4))
par(mgp=c(2,1,0),mfrow=c(3,2),mar=c(3,3,1,1),oma=rep(.5,4))
for (i in 1:length(L)) {
    plot(NULL,xlim=c(-1,1),ylim=c(0,.6),xlab=expression(rho),ylab='RMSE(true,est)')
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
    f(L[[i]])
    legend("topright",bty='n',lty=c(1,1,2,2),col=c("blue","red","blue","red"),cex=.7,ncol=2,
           c("(resp,IRT)","(resp,CDM)","(True,IRT)","(True,CDM)"))
#####
    plot(NULL,xlim=c(c(-1,1)),ylab="IMV",xlab=expression(rho),ylim=c(0,.2))
    f<-function(out,...) {
        om<-do.call("rbind",out)
        abline(h=0)
        pf(rs,om[,1],...)
        pf(rs,om[,2],col='blue',...,lty=2)
        pf(rs,om[,3],col='red',...,lty=2)
    }
                                        #f(out1,lty=1)
    f(L[[i]])
    legend("topright",bty='n',lty=c(1,2,2),col=c("black","blue","red"),cex=.7,
           c("(IRT,CDM)","(IRT,True)","(CDM,True)")
           )
    ##
}
dev.off()
