source("00funs.R")


################################3
##simulate data with naughty cdm
cdm.sim<-function(a,N=1000,sk.offset=0,nsk=6) {
    ##skills
                                        #sk<-rbinom(5*N,1,.5)
                                        #sk<-matrix(sk,nrow=N,ncol=5)
    th<-rnorm(N)
    sk<-runif(nsk)
    p<-a*(outer(th,sk-mean(sk),'-'))
    p<-p+sk.offset ##controlling prevalence of skills
    p<-apply(p,2,function(x) 1/(1+exp(-(x))))
    sk<-apply(p,2,function(x) rbinom(nrow(p),1,p))
          
    ##qmatrix
    S<-TRUE
    while (S) {
        qm<-matrix(rbinom(50*nsk,1,.65),50,nsk)
        S<-any(c(colMeans(qm),rowMeans(qm))==0)
    }        

    ##response probabilities 
    pL<-respL<-list()
    g<-runif(nrow(qm),min=0,max=.35)
    s<-runif(nrow(qm),min=0,max=.35)
    #s<-rep(.05,nrow(qm))
    #g<-rep(.05,nrow(qm))
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
    
    ## p.irt<-irt.pr(resp)
    ## p.cdm<-cdm.pr.marg(resp,qm)
    
    ## ## L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    ## ## L<-lapply(L,as.numeric)
    ## ## z<-do.call("cbind",L)
    ## ## true<-cor(z)[2,3:4]
    ## coors<-list()
    ## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
    ## ss<-sapply(coors,function(x) x[1,])
    ## true<-rowMeans(ss)[3:4]
    
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true)
    #list(t=true,om=om)
    om
}
a<-sort(runif(50,min=0,max=3))
library(parallel)
out1<-mclapply(a,cdm.sim,mc.cores=10)
out2<-mclapply(a,cdm.sim,mc.cores=10,sk.offset=1.5)
out<-list(out1=out1,out2=out2)
save(out,file="scenario1.Rdata")



#pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario1.pdf",width=6,height=3)
par(mgp=c(2,1,0),mfrow=c(1,2),mar=c(3,3,1,1),oma=rep(.5,4))
####
plot(NULL,xlim=c(0,3),ylim=c(0,.5),xlab='a',ylab='rmse(oos resp,est)')
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
    pf(a,irt,col='blue',...)
    pf(a,cdm,col='red',...)
    pf(a,irt.p,col='blue',...,lty=2)
    pf(a,cdm.p,col='red',...,lty=2)
}
#f(out1,lty=1)
f(out2)
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"),title="est")
#####
plot(NULL,xlim=c(0,3),ylab="IMV",xlab='a',ylim=c(0,.15))
f<-function(out,...) {
    om<-do.call("rbind",out)
    abline(h=0)
    pf(a,om[,1],...)
    pf(a,om[,2],col='blue',...)
    pf(a,om[,3],col='red',...)
}
#f(out1,lty=1)
f(out2)
legend("topright",bty='n',lwd=1,lty=c(1,2),title=expression(delta),legend=c(0,1.5))
##
#dev.off()



## ##anya work:
## cdm.sim.fast <- function(a, N = 500, sk.offset = 0) {
##   th <- rnorm(N)
##   sk <- runif(10)
##   p <- a * (outer(th, sk - mean(sk), "-")) + sk.offset
##   p <- 1 / (1 + exp(-p))
##   sk <- matrix(rbinom(length(p), 1, p), nrow = N)
  
##   repeat {
##     qm <- matrix(rbinom(50 * 10, 1, .65), 50, 10)
##     if (all(colMeans(qm) > 0) && all(rowMeans(qm) > 0)) break
##   }
  
##   g <- runif(nrow(qm), min = 0, max = 0.35)
##   s <- runif(nrow(qm), min = 0, max = 0.35)
##   respL <- vector("list", nrow(qm))
##   pL <- vector("list", nrow(qm))
  
##   for (i in seq_len(nrow(qm))) {
##     ii <- which(qm[i, ] == 1)
##     z <- sk[, ii, drop = FALSE]
##     rm <- rowMeans(z)
##     p.tmp <- ifelse(rm == 1, 1 - s[i], g[i])
##     respL[[i]] <- rbinom(N, 1, p.tmp)
##     pL[[i]] <- p.tmp
##   }
  
##   resp <- do.call(cbind, respL)
##   p.true <- do.call(cbind, pL)
  
##   resp <- as.data.frame(resp)
##   names(resp) <- paste("i", 1:ncol(resp), sep = "")
  
##   p.irt <- irt.pr(resp)
##   p.cdm <- cdm.pr.marg(resp, qm)
  
##   z <- cbind(as.numeric(as.matrix(resp)),
##              as.numeric(as.matrix(p.true)),
##              as.numeric(p.irt),
##              as.numeric(p.cdm))
##   true <- cor(z)[2, 3:4]
  
##   om <- oos.compare.newresp(resp, qm, truep = p.true)
##   list(t = true, om = om)
## }

## a<-sort(runif(100,min=0,max=3))
## library(parallel)
## out1<-mclapply(a,cdm.sim.fast,mc.cores=10)
## out2<-mclapply(a,cdm.sim.fast,mc.cores=10,sk.offset=1.5)
