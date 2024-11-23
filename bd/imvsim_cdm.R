##idea: you are taking within-item correlations below. might be better to take items across all correlations


source("00funs.R")


################################3
##simulate data with naughty cdm
cdm.sim<-function(a,N=5000) {
    ##skills
                                        #sk<-rbinom(5*N,1,.5)
                                        #sk<-matrix(sk,nrow=N,ncol=5)
    th<-rnorm(N)
    sk<-runif(5)
    p<-a*(outer(th,sk-mean(sk),'-'))
    p<-apply(p,2,function(x) 1/(1+exp(-x)))
    sk<-apply(p,2,function(x) rbinom(nrow(p),1,p))
          
    ##qmatrix
    #qm<-diag(5)
    S<-TRUE
    while (S) {
        qm<-matrix(rbinom(5*30,1,.5),30,5)
        S<-any(c(rowSums(qm),colSums(qm))==0)
    }            
    qm0<-qm
    ##response probabilities 
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
    
    ## ##using true values
    coors<-list()
    for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
    ##note that third row << fourth row
    s<-sapply(coors,function(x) x[1,])
    true<-rowMeans(s)[3:4]
    
    ##cv imv values
    om<-oos.compare(resp,qm0,nfolds=5)
    list(t=true,om=om)
}

## N<-runif(50,3,4.5)
## N<-round(10^N)
## om<-sapply(N,cdm.sim)

a<-sort(runif(100,min=0,max=3))
library(parallel)
out<-mclapply(a,cdm.sim,mc.cores=4)


#pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario1.pdf",width=4,height=2.5)
par(mgp=c(2,1,0),mfrow=c(1,2))

####
tr<-lapply(out,function(x) x$t)

z<-do.call("rbind",tr)
plot(a,z[,1],type='l',col='black',ylim=0:1)
lines(a,z[,2],type='l',col='red')


#####
om<-lapply(out,function(x) x$om)
df<-data.frame(a=a,om=unlist(om))

plot(df$a,df$om,pch=19,col='lightgray',cex=.5,ylab="IMV",xlab='a')
abline(h=0)

m<-loess(om~a,df)
pr<-predict(m,df$a,se=TRUE)
lines(df$a,pr$fit,col='red',lwd=TRUE)
cc<-col2rgb("red")
polygon(c(df$a,rev(df$a)),c(pr$fit+1.96*pr$se.fit,rev(pr$fit-1.96*pr$se.fit)),border=NA,col=rgb(cc[1],cc[2],cc[3],max=255,alpha=30))
#dev.off()
