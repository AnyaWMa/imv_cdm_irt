irt.pr<-function(resp,modeltype='Rasch',th.type='EAP') {
    ##get irt pvalues
    library(mirt)
    m<-mirt(resp,1,modeltype)
    th<-fscores(m,method=th.type)
    co<-coef(m,IRTpars=TRUE,simplify=TRUE)$items
    x<-list()
    for (i in 1:nrow(co)) {
        z<-co[i,1]*(th[,1]-co[i,2])
        x[[i]]<-1/(1+exp(-1*z))
    }
    do.call("cbind",x)
}

oos.compare<-function(resp,qm,nfolds=5,modeltype,estmethod = "mp") {
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
        x<-irw::long2resp(df)
        id<-x$id
        x<-x[,names(resp)]
        ##
        p.cdm<-cdm.pr.marg.update(x,qm, modeltype, estmethod)
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

oos.compare.update <-function(resp,qm,truep, nfolds=5,modeltype,estmethod = "mp") {
  id<-1:nrow(resp)
  item<-names(resp)
  L<-list()
  for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],resp=resp[,i],truep=truep[,i])
  df<-do.call("rbind",L)
  df$gr<-sample(1:nfolds,nrow(df),replace=TRUE)
  om<-numeric()
  om.cdm<-numeric()
  om.irt<-numeric()
  for (gr in unique(df$gr)) {
    oos<-df[df$gr==gr,]
    ins<-df[df$gr!=gr,]
    x<-irw::long2resp(df)
    id<-x$id
    x<-x[,names(resp)]
    ##
    p.cdm<-cdm.pr.marg.update(x,qm, modeltype, estmethod)
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
    om.cdm[gr]<-imv::imv.binary(oos$resp,oos$p.cdm,oos$truep)
    om.irt[gr]<-imv::imv.binary(oos$resp,oos$p.irt,oos$truep)
  }
  c(mean(om), mean(om.cdm), mean(om.irt))
  
}

oos.compare.newresp<-function(resp,qm,truep, modeltype = "DINA", irtmodel = "Rasch", estmethod = "mp") {
    id<-1:nrow(resp)
    item<-names(resp)
    L<-list()
    for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],resp=resp[,i],truep=truep[,i])
    df<-do.call("rbind",L)
    om<-numeric()
    x<-irw::long2resp(df)
    id<-x$id
    print(x)
    x<-x[,names(resp)]
    ##
    p.cdm<-cdm.pr.marg.update(x,qm, modeltype, estmethod)
    L<-list()
    for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.cdm=p.cdm[,i])
    p<-do.call("rbind",L)
    df<-merge(df,p)
    ##
    p.irt<-irt.pr(x, modeltype = irtmodel)
    L<-list()
    for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.irt=p.irt[,i])
    p<-do.call("rbind",L)
    df<-merge(df,p)
    ##
    df$resp<-rbinom(nrow(df),1,df$truep) ##obliterate old response
    om<-imv::imv.binary(df$resp,df$p.irt,df$p.cdm)
    oracle.irt<-imv::imv.binary(df$resp,df$p.irt,df$truep)
    oracle.cdm<-imv::imv.binary(df$resp,df$p.cdm,df$truep)
    ##
    rms<-function(x) sqrt(mean(x^2))
    rms.irt<-rms(df$resp-df$p.irt)
    rms.cdm<-rms(df$resp-df$p.cdm)
    rms.irt.p<-rms(df$truep-df$p.irt)
    rms.cdm.p<-rms(df$truep-df$p.cdm)
    ##
    c(om=om,oracle.irt=oracle.irt,oracle.cdm=oracle.cdm,
      rms.irt=rms.irt,rms.cdm=rms.cdm,rms.irt.p=rms.irt.p,rms.cdm.p=rms.cdm.p)
}


oos.compare.2q<-function(resp,qm1,qm2,nfolds=5) {
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
        x<-irwpkg::irw_long2resp(ins)
        id<-x$id
        x<-x[,names(resp)]
        ##
        p.cdm<-cdm.pr.marg(x,qm1)
        L<-list()
        for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.cdm1=p.cdm[,i])
        p<-do.call("rbind",L)
        oos<-merge(oos,p)
        ##
        p.cdm<-cdm.pr.marg(x,qm2)
        L<-list()
        for (i in 1:ncol(resp)) L[[i]]<-data.frame(id=id,item=names(resp)[i],p.cdm2=p.cdm[,i])
        p<-do.call("rbind",L)
        oos<-merge(oos,p)
        ##
        om[gr]<-imv::imv.binary(oos$resp,oos$p.cdm1,oos$p.cdm2)
    }
    mean(om)
}

cdm.pr.marg<-function(resp,qm,modeltype="DINA") {
  ##cdm pvalues
  library(GDINA)
  m <- GDINA(resp,qm,modeltype)
  map <- personparm(m, what ="mp")[,1:ncol(qm)]
  gs<-coef(m,what='gs')
  p<-list()
  
  for (i in 1:nrow(qm)) {
    ii<-which(qm[i,]==1)
    z<-map[,ii,drop=FALSE]
    rm<-apply(z,1,prod)
    
    if (modeltype == "DINO"){
      rm_dino <- 1 - apply(1 - z, 1, prod)
      p[[i]] <- (1 - gs[i,2]) * rm_dino + (1 - rm_dino) * gs[i,1]
    } else {
      p[[i]]<-(1-gs[i,2])*rm + (1-rm)*gs[i,1]
    }
  }
  p.cdm<-do.call("cbind",p)
}

cdm.pr.marg2<-function(resp,qm,modeltype="GDINA") { #really only works with GDINA
   # if (modeltype!="GDINA") stop("must be GDINA")
    ##
    library(GDINA)
    m <- GDINA(resp,qm,modeltype)
    map <- personparm(m, what = "mp")[,1:ncol(qm)]
    co<-coef(m) #gs<-coef(m,what='gs')
    p<-list()
    for (i in 1:nrow(qm)) {
        ii<-which(qm[i,]==1)
        z<-map[,ii,drop=FALSE]
        nms<-names(co[[i]])
        cats<-gsub(")","",gsub("P(","",nms,fixed=TRUE),fixed=TRUE)
        cats<-strsplit(cats,"")
        cats<-lapply(cats,as.numeric)
        cats<-do.call("rbind",cats)
        pr<-list()
        for (j in 1:nrow(cats)) {
            y<-cats[j,]
            z2<-z
        for (k in 1:ncol(z2)) if (y[k]==1) z2[,k]<-z[,k] else z2[,k]<-1-z[,k]
            pr[[j]]<-apply(z2,1,prod)
        }
        pr<-do.call("cbind",pr)
        p[[i]]<-pr %*% matrix(co[[i]],ncol=1)
    }
    p.cdm<-do.call("cbind",p)
}


cdm.pr.marg.update<-function(resp,qm,modeltype="DINA",estmethod = "mp") { 
  library(GDINA)
  m <- GDINA(resp,qm,modeltype, mono.constr = TRUE)
  if (!GDINA::extract(m, "convergence")) {
    print("CDM non-converge")
    return(NA)
  } 
  map <- personparm(m, what = estmethod)[,1:ncol(qm), drop = FALSE]
  co<-coef(m) #gs<-coef(m,what='gs')
  p<-list()
  for (i in 1:nrow(qm)) {
    ii<-which(qm[i,]==1)
    z<-map[,ii,drop=FALSE]
    nms<-names(co[[i]])
    cats<-gsub(")","",gsub("P(","",nms,fixed=TRUE),fixed=TRUE)
    cats<-strsplit(cats,"")
    cats<-lapply(cats,as.numeric)
    cats<-do.call("rbind",cats)
    pr<-list()
    for (j in 1:nrow(cats)) {
      y<-cats[j,]
      z2<-z
      for (k in 1:ncol(z2)) if (y[k]==1) z2[,k]<-z[,k] else z2[,k]<-1-z[,k]
      pr[[j]]<-apply(z2,1,prod)
    }
    pr<-do.call("cbind",pr)
    p[[i]]<-pr %*% matrix(co[[i]],ncol=1)
  }
  p.cdm<-do.call("cbind",p)
}


check.qm.complete <- function(qm) {
  K <- ncol(qm)  # Number of skills (attributes)
  
  # Find rows that have exactly one '1' (single-attribute items)
  single_attribute_rows <- apply(qm, 1, function(row) ifelse(sum(row) == 1, which(row == 1), NA))
  
  # Extract unique skills (columns) covered by single-attribute items
  unique_single_skills <- unique(na.omit(single_attribute_rows))
  
  # Check if all K skills are covered
  return(length(unique_single_skills) == K)
    
}

