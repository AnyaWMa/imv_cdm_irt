## source("00funs.R")

## ################################3
## ##simulate data with naughty cdm
## N<-1000
## ##skills
## sk<-rbinom(5*N,1,.5)
## sk<-matrix(sk,nrow=N,ncol=5)
## ##qmatrix
## qm<-diag(5)
## qm0<-rbind(qm,qm,qm)
## ##response probabilities 
## pL<-respL<-list()
## for (j in 1:3) {
##     p<-resp<-list()
##     g<-s<-0.1*j
##     for (i in 1:nrow(qm)) {
##         ii<-which(qm[i,]==1)
##         z<-sk[,ii,drop=FALSE]
##         rm<-rowMeans(z)
##         p.tmp<-ifelse(rm==1,1-s,g)
##         resp[[i]]<-rbinom(N,1,p.tmp)
##         p[[i]]<-p.tmp
##     }
##     respL[[j]]<-do.call("cbind",resp)
##     pL[[j]]<-do.call("cbind",p)
## }
## resp<-do.call("cbind",respL)
## p.true<-do.call("cbind",pL)
## resp<-as.data.frame(resp)
## names(resp)<-paste("i",1:ncol(resp))

## p.irt<-irt.pr(resp)
## p.cdm<-cdm.pr(resp,qm0)

## coors<-list()
## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
## ##note that third row << fourth row
## sapply(coors,function(x) x[1,])

## ################################3
## ##simulate data with rasch model
## th<-rnorm(1000)
## b<-sort(rnorm(15))
## th.mat<-matrix(th,length(th),length(b),byrow=FALSE)
## b.mat<-matrix(b,length(th),length(b),byrow=TRUE)
## k<-th.mat-b.mat
## p<-1/(1+exp(-k))
## resp<-p
## for (i in 1:ncol(resp)) resp[,i]<-rbinom(nrow(resp),1,p[,i])
## resp<-as.data.frame(resp)
## names(resp)<-paste("i",1:ncol(resp))
## p.true<-p

## qm<-matrix(0,ncol=3,nrow=5)
## L<-list(qm,qm,qm)
## for (i in 1:3) L[[i]][,i]<-1
## qm0<-do.call("rbind",L)

## p.irt<-irt.pr(resp)
## p.cdm<-cdm.pr(resp,qm0)

## coors<-list()
## for (i in 1:ncol(resp)) coors[[i]]<-cor(cbind(resp[,i],p.true[,i],p.irt[,i],p.cdm[,i]))
## ##note that third row >> fourth row!
## sapply(coors,function(x) x[1,])
