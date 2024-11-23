source("00funs.R")


################################3
##simulate data with naughty cdm
cdm.sim<-function(a,N=1000,sk.offset=0) {
    ##skills
                                        #sk<-rbinom(5*N,1,.5)
                                        #sk<-matrix(sk,nrow=N,ncol=5)
    th<-rnorm(N)
    sk<-runif(10)
    p<-a*(outer(th,sk-mean(sk),'-'))
    p<-p+sk.offset ##controlling prevalence of skills
    p<-apply(p,2,function(x) 1/(1+exp(-(x))))
    sk<-apply(p,2,function(x) rbinom(nrow(p),1,p))
          
    ##qmatrix
    S<-TRUE
    while (S) {
        qm<-matrix(rbinom(50*10,1,.65),50,10)
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
    
    p.irt<-irt.pr(resp)
    p.cdm<-cdm.pr(resp,qm)
    
    L<-list(as.matrix(resp),p.true,p.irt,p.cdm)
    L<-lapply(L,as.numeric)
    z<-do.call("cbind",L)
    true<-cor(z)[2,3:4]
    
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true)
    list(t=true,om=om)
}

a<-sort(runif(100,min=0,max=3))
library(parallel)
out1<-mclapply(a,cdm.sim,mc.cores=10)
out2<-mclapply(a,cdm.sim,mc.cores=10,sk.offset=1.5)
#out<-list(out1=out1,out2=out2)
#save(out,"scenario1.Rdata")

pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario1.pdf",width=6,height=3)
par(mgp=c(2,1,0),mfrow=c(1,2),mar=c(3,3,1,1),oma=rep(.5,4))
####
plot(NULL,xlim=c(0,3),ylim=c(0,1),xlab='a',ylab='cor(true,est)')
f<-function(out,...) {
    tr<-lapply(out,function(x) x$t)
    z<-do.call("rbind",tr)
    irt<-z[,1]
    cdm<-z[,2]
    pf<-function(x,y,...) {
        m<-loess(y~x)
        lines(x,predict(m),...,lwd=2)
    }
    lines(pf(a,irt,...))
    lines(pf(a,cdm,col='red',...))
}
f(out1,lty=1)
f(out2,lty=2)
legend("bottomright",bty='n',fill=c("black","red"),c("irt","cdm"),title="est")

#####
plot(NULL,xlim=c(0,3),ylab="IMV",xlab='a',ylim=c(-.075,.1))
f<-function(out,...) {
    om<-lapply(out,function(x) x$om)
    df<-data.frame(a=a,om=unlist(om))
    abline(h=0)
    m<-loess(om~a,df)
    pr<-predict(m,df$a,se=TRUE)
    lines(df$a,pr$fit,lwd=TRUE,...)
    cc<-col2rgb("red")
    #polygon(c(df$a,rev(df$a)),c(pr$fit+1.96*pr$se.fit,rev(pr$fit-1.96*pr$se.fit)),border=NA,col=rgb(cc[1],cc[2],cc[3],max=255,alpha=30))
}
f(out1,lty=1)
f(out2,lty=2)
legend("topright",bty='n',lwd=1,lty=c(1,2),title=expression(delta),legend=c(0,1.5))
##
dev.off()
