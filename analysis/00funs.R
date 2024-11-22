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
    gs<-coef(m,what='gs')
    p<-list()
    for (i in 1:nrow(qm0)) {
        ii<-which(qm0[i,]==1)
        z<-map[,ii,drop=FALSE]
        rm<-rowMeans(z)
        p[[i]]<-ifelse(rm==1,1-gs[i,2],gs[i,1])
    }
    p.cdm<-do.call("cbind",p)
}

oos.compare<-function(resp,qm0,nfolds) {
    id<-1:nrow(resp)
    item<-names(resp)
    L<-list()
    for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],resp=resp[,i])
    df<-do.call("rbind",L)
    df$gr<-sample(1:nfolds,nrow(df),replace=TRUE)
    
    om<-numeric()
    for (gr in unique(df$gr)) {
        oos<-df[df$gr==gr,]
        ins<-df[df$gr!=gr,]
        x<-irw::long2resp(ins)
        id<-x$id
        x<-x[,names(resp)]
        ##
        p.cdm<-cdm.pr(x,qm0)
        L<-list()
        for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.cdm=p.cdm[,i])
        p<-do.call("rbind",L)
        oos<-merge(oos,p)
        ##
        p.irt<-irt.pr(x)
        L<-list()
        for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.irt=p.irt[,i])
        p<-do.call("rbind",L)
        oos<-merge(oos,p)
        ##
        om[gr]<-imv::imv.binary(oos$resp,oos$p.irt,oos$p.cdm)
    }
    mean(om)
}
