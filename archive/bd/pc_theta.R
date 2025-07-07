a<-c(2,1,0,0,0,0)
N<-5000
nsk<-6
sk.offset<-0

##skills
                                        #sk<-rbinom(5*N,1,.5)
                                        #sk<-matrix(sk,nrow=N,ncol=5)
th<-rnorm(N)
sk<-runif(nsk)
p<-outer(th,sk-mean(sk),'-')
for (i in 1:ncol(p)) p[,i]<-a[i]*p[,i]
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

library(GDINA)
m <- GDINA(resp,qm,'DINA')
map <- personparm(m, what = "mp")[,1:ncol(qm)]

fa<-factanal(map,factors=1,scores='regression')
cor(fa$scores,th)
