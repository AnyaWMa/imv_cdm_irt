source("00funs.R")
library(GDINA)
resp<-frac20$dat
qm<-frac20$Q

set.seed(1010101)
oos.compare(resp,qm,5,modeltype="2PL")
##imv=0.017 which is way different than what anya got.

##random qmatrix
getq<-function(qm) {
    nr<-nrow(qm)
    nc<-ncol(qm)
    S<-TRUE
    while (S) {
                                        #qm<-matrix(rbinom(nr*nc,1,.65),nr,nc)
        cm<-colMeans(qm)
        for (i in 1:nr) qm[i,]<-rbinom(nc,1,cm)
        S<-any(c(colMeans(qm),rowMeans(qm))==0)
    }        
    qm
}
om<-numeric()
for (i in 1:3) {
    qm2<-getq(qm)
    om[i]<-oos.compare.2q(resp,qm2,qm,nfolds=5)
}
om

