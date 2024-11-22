source("00funs.R")


################################3
##simulate data with naughty cdm
cdm.sim<-function(N) {
    ##skills
    sk<-rbinom(5*N,1,.5)
    sk<-matrix(sk,nrow=N,ncol=5)
    ##qmatrix
    #qm<-diag(5)
    qm<-matrix(rbinom(10*5,1,.5),10,5)
    qm0<-rbind(qm,qm,qm)
    ##response probabilities 
    pL<-respL<-list()
    g<-c(.1,.2,.3)
    for (j in 1:3) {
        p<-resp<-list()
        for (i in 1:nrow(qm)) {
            ii<-which(qm[i,]==1)
            z<-sk[,ii,drop=FALSE]
            rm<-rowMeans(z)
            p.tmp<-ifelse(rm==1,1-g[j],g[j])
            resp[[i]]<-rbinom(N,1,p.tmp)
            p[[i]]<-p.tmp
        }
        respL[[j]]<-do.call("cbind",resp)
        pL[[j]]<-do.call("cbind",p)
    }
    resp<-do.call("cbind",respL)
    p.true<-do.call("cbind",pL)
    resp<-as.data.frame(resp)
    names(resp)<-paste("i",1:ncol(resp),sep='')
    
    ## p.irt<-irt.pr(resp)
    ## p.cdm<-cdm.pr(resp,qm0)
    
    ## ##using true values
    ## coors<-list()
    ## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
    ## ##note that third row << fourth row
    ## sapply(coors,function(x) x[1,])
    
    ##cv imv values
    oos.compare(resp,qm0,nfolds=4)
}

N<-runif(5,3,4.5)
N<-round(10^N)
om<-sapply(N,cdm.sim)

plot(N,om,pch=19); abline(h=0)
