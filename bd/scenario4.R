remotes::install_github("hansorlee/irwpkg")
source("00funs.R") ##https://github.com/AnyaWMa/IRW-Qmatrix/blob/main/bd/00funs.R
library(GDINA)
################################3
##simulate data with naughty cdm
cdm.sim<-function(a,N=500,sk.offset=0,nsk=6,bound = 0.1,prevalence=0.35, modeltype = "DINA",estmethod="MAP") {
    ##skills
    th<-rnorm(N)
    sk<-runif(nsk)
    p<-a*(outer(th,sk-mean(sk),'-'))
    p<-p+sk.offset ##controlling prevalence of skills
    p<-apply(p,2,function(x) 1/(1+exp(-(x))))
    sk<-apply(p,2,function(x) rbinom(nrow(p),1,p))
    ##qmatrix
    S<-TRUE
    
    n_valid <- 44     # desired number of valid rows
    
    # Initialize an empty matrix to store valid rows
    valid_rows <- matrix(nrow = 0, ncol = nsk)
    
    # Loop until we have 50 valid rows
    while (nrow(valid_rows) < n_valid) {
      #prevalence of skills in items
      #ps <-runif(1,min=prevalence-0.15,max=prevalence + 0.15)
      ps <- rnorm(1, mean = prevalence, sd = 0.15)
      candidate <- rbinom(nsk, 1, ps)  # Generate one candidate row
      if (any(is.na(candidate)) || all(candidate == 0)) next
      valid_rows <- rbind(valid_rows, candidate)
    }
    
    # Create the identity matrix
    I_k <- diag(nsk)
    # Combine the identity matrix with the 50 valid rows to form the Q-matrix
    qm <- rbind(I_k, valid_rows)
    # Display the resulting Q-matrix
    rownames(qm) <- NULL
  
    ##response probabilities 
    pL<-respL<-list()
    if (is.null(bound)) {
        g<-runif(nrow(qm),min=0,max=.35)
        s<-runif(nrow(qm),min=0,max=.35)
    } else {
        g<-rep(bound,nrow(qm))
        s<-rep(bound,nrow(qm))
    }
    
    gs <- cbind(unlist(g), unlist(s))
    if (modeltype == "GDINA"){
      simD <- simGDINA(N,qm,gs.parm = gs, model = "GDINA",attribute = sk)
      dat <- extract(simD,"dat")
      resp <- as.data.frame(dat)
      names(resp)<-paste("i",1:ncol(resp),sep='')
      
      J <- nrow(qm)  # number of items
      co <- simD$catprob.parm
      p<-list()
      
      for (i in 1:nrow(qm)) {
        ii<-which(qm[i,]==1)
        z<-sk[,ii,drop=FALSE]
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
          #print(z2)
          pr[[j]]<-apply(z2,1,prod)
          
        }
        
        pr<-do.call("cbind",pr)
        p[[i]]<-pr %*% matrix(co[[i]],ncol=1)
        
      }
      p.true <- do.call("cbind", p)
      
    }
    else {
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
    }
    ##cv imv values
    om<-oos.compare.newresp(resp,qm,truep=p.true, modeltype = modeltype, estmethod = estmethod)
    ##
    om
}
a<-sort(runif(n = 100,min=0,max=3))
library(parallel)

out2dina<-list()
for (i in c(.2, 0.4, 0.6)) out2dina[[as.character(i)]]<-mclapply(a,cdm.sim,mc.cores=10,sk.offset=1.5,bound = 0.1, prevalence=i,  modeltype = "GDINA",estmethod = "mp")
#save(out2 , file = "../simulation_data/out_dina_offset15_ps035.RData")
out2 <- out2dina
#pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario1.pdf",width=6,height=3)
par(mgp = c(2,1,0), mfrow = c(3,2), mar = c(3,3,2,1), oma = rep(0.5,4))

panel_titles <- c("RMSE ps = 0.2", "IMV ps = 0.2", "RMSE ps = 0.4", "IMV ps = 0.4", "RMSE ps = 0.6", "IMV ps = 0.6")

for (i in 1:length(out2)) {
  # RMSE Plot
  plot(NULL, xlim = c(0,3), ylim = c(0,0.6), xlab = 'a', ylab = 'RMSE(oos resp, est)')
  title(panel_titles[2*i - 1], line = 0.5, cex.main = 1)
  mtext(paste0("(", letters[2*i - 1], ")"), side = 3, adj = 0, line = -1.2, cex = 0.9)
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x)
    lines(x, predict(m), ..., lwd = 2)
  }
  
  f <- function(out, ...) {
    z <- do.call("rbind", out)
    irt <- z[,4]
    cdm <- z[,5]
    irt.p <- z[,6]
    cdm.p <- z[,7]
    pf(a, irt, col = 'blue', ...)
    pf(a, cdm, col = 'red', ...)
    pf(a, irt.p, col = 'blue', ..., lty = 2)
    pf(a, cdm.p, col = 'red', ..., lty = 2)
  }
  
  f(out2[[i]])
  legend("topright", bty = 'n', lty = c(1,1,2,2),
         col = c("blue","red","blue","red"), cex = .7,
         legend = c("(resp,IRT)","(resp,CDM)","(True,IRT)","(True,CDM)"))
  
  # IMV Plot
  plot(NULL, xlim = c(0,3), ylab = "IMV", xlab = 'a', ylim = c(-0.05, 0.2))
  title(panel_titles[2*i], line = 0.5, cex.main = 1)
  mtext(paste0("(", letters[2*i], ")"), side = 3, adj = 0, line = -1.2, cex = 0.9)
  
  f <- function(out, ...) {
    om <- do.call("rbind", out)
    abline(h = 0)
    pf(a, om[,1], ...)
    pf(a, om[,2], col = 'blue', ..., lty = 2)
    pf(a, om[,3], col = 'red', ..., lty = 2)
  }
  
  f(out2[[i]])
  legend("topright", bty = 'n', lty = c(1,2,2),
         col = c("black","blue","red"), cex = .7,
         legend = c("(IRT,CDM)","(IRT,True)","(CDM,True)"))
}
##legend("topright",bty='n',lwd=1,lty=c(1,2),title=expression(delta),legend=c(0,1.5))
##dev.off()
