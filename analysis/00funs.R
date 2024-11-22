irt.pr<-function(resp,modeltype='Rasch',th.type='EAP') {
    ##get irt pvalues
    library(mirt)
    m<-mirt(resp,1,modeltype)
    th<-fscores(m)
    co<-coef(m,IRTpars=TRUE,simplify=TRUE)$items
    x<-list()
    for (i in 1:nrow(co)) {
        z<-co[i,1]*(th[,1]-co[i,2])
        x[[i]]<-1/(1+exp(-1*z))
    }
    do.call("cbind",x)
}

cdm.pr<-function(resp,qm0,modeltype="DINA") {
    ##cdm pvalues
    library(GDINA)
    m <- GDINA(resp,qm0,model=modeltype)
    map <- personparm(m, what = "EAP")[,1:ncol(qm0)]
    gs<-do.call("rbind",coef(m))
    p<-list()
    for (i in 1:nrow(qm0)) {
        ii<-which(qm0[i,]==1)
        z<-map[,ii,drop=FALSE]
        rm<-rowMeans(z)
        p[[i]]<-ifelse(rm==1,gs[i,2],gs[i,1])
    }
    p.cdm<-do.call("cbind",p)
}
