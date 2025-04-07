remotes::install_github("hansorlee/irwpkg")
source("00funs.R") ##https://github.com/AnyaWMa/IRW-Qmatrix/blob/main/bd/00funs.R
library(GDINA)
library(parallel)
################################3
##simulate data with naughty cdm
cdm.sim<-function(a,N=500,sk.offset=0,nsk=6,bound=NULL, modeltype = "DINA",estmethod="MAP") {
    ##skills
    th<-rnorm(N)
    sk<-runif(nsk)
    p<-a*(outer(th,sk-mean(sk),'-'))
    p<-p+sk.offset ##controlling prevalence of skills
    p<-apply(p,2,function(x) 1/(1+exp(-(x))))
    sk<-apply(p,2,function(x) rbinom(nrow(p),1,p))
    ## qmatrix
    S<-TRUE
    
    n_valid <- 44     # desired number of valid rows
    
    # Initialize an empty matrix to store valid rows
    valid_rows <- matrix(nrow = 0, ncol = nsk)
    
    # Loop until we have 50 valid rows
    while (nrow(valid_rows) < n_valid) {
      #prevalence of skills in items
      ps <- rnorm(1, mean = 0.35, sd = 0.15)
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

simulate_scenario_1 <- function(modeltype = "GDINA", estmethod = "mp") {
  # Generate sorted vector a
  a <- sort(runif(n = 100, min = 0, max = 3))
  
  # Initialize output list
  out2 <- list()
  
  # Loop over bounds and apply cdm.sim in parallel
  for (i in c(0.1, 0.2, 0.3)) {
    out2[[as.character(i)]] <- mclapply(
      a,
      cdm.sim,
      mc.cores = 10,
      sk.offset = 1.5,
      bound = i,
      modeltype = modeltype,
      estmethod = estmethod
    )
  }
  
  # Return both a and out2
  return(list(a = a, out2 = out2))
}

## Simulate GDINA
result_gdina_1 <- simulate_scenario_1(modeltype = "GDINA", estmethod = "mp")
save(result_gdina_1 , file = "../simulation_data/scenario1_out_gdina_offset15_ps035_mp.RData")

## Simulate GDINA
result_dina_1 <- simulate_scenario_1(modeltype = "DINA", estmethod = "mp")
save(result_dina_1 , file = "../simulation_data/scenario1_out_dina_offset15_ps035_mp.RData")
#save(out2 , file = "../simulation_data/out_dina_offset15_ps035.RData")

simulate_scenario_2 <- function() {
  # Generate sorted vector a
  a <- sort(runif(n = 100, min = 0, max = 3))
  
  # Initialize output list
  out2 <- list()
  
  # Loop over model types and estimation methods
  for (i in c("DINA", "GDINA")) {
    print(i)
    for (j in c("mp", "MAP")) {
      print(j)
      key <- paste(i, j, sep = "_")
      out2[[key]] <- mclapply(
        a,
        cdm.sim,
        mc.cores = 10,
        sk.offset = 1.5,
        bound = 0.1,
        modeltype = i,
        estmethod = j
      )
    }
  }
  
  # Return both a and out2
  return(list(a = a, out2 = out2))
}

result_compare_estimator <- simulate_scenario_2()
save(result_compare_estimator , file = "../simulation_data/scenario1_result_compare_estimator.RData")
#pdf("/home/bdomingu/Dropbox/Apps/Overleaf/CDM_predictions/scenario1.pdf",width=6,height=3)

plot_rmse_imv_panels <- function(a, out2) {
  par(mgp = c(2, 1, 0), mfrow = c(3, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4))
  
  panel_titles <- c("RMSE (High item quality)", "IMV (High item quality)", 
                    "RMSE (Medium item quality)", "IMV (Medium item quality)", 
                    "RMSE (Low item quality)", "IMV (Low item quality)")
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x)
    lines(x, predict(m), ..., lwd = 2)
  }
  
  for (i in seq_along(out2)) {
    
    ## RMSE Plot
    plot(NULL, xlim = c(0, 3), ylim = c(0, 0.6), 
         xlab = 'a', ylab = 'RMSE(oos resp, est)')
    title(panel_titles[2 * i - 1], line = 0.5, cex.main = 1)
    mtext(paste0("(", letters[2 * i - 1], ")"), side = 3, adj = 0, cex = 1)
    
    f_rmse <- function(out, ...) {
      z <- do.call("rbind", out)
      irt <- z[, 4]
      cdm <- z[, 5]
      irt.p <- z[, 6]
      cdm.p <- z[, 7]
      pf(a, irt, col = 'blue', ...)
      pf(a, cdm, col = 'red', ...)
      pf(a, irt.p, col = 'blue', lty = 2, ...)
      pf(a, cdm.p, col = 'red', lty = 2, ...)
    }
    
    f_rmse(out2[[i]])
    
    legend("topright", bty = 'n', lty = c(1, 1, 2, 2),
           col = c("blue", "red", "blue", "red"), cex = 1,
           legend = c("(resp,IRT)", "(resp,CDM)", "(True,IRT)", "(True,CDM)"))
    
    ## IMV Plot
    plot(NULL, xlim = c(0, 3), ylim = c(-0.05, 0.2),
         xlab = 'a', ylab = 'IMV')
    title(panel_titles[2 * i], line = 0.5, cex.main = 1)
    mtext(paste0("(", letters[2 * i], ")"), side = 3, adj = 0, cex = 1)
    
    f_imv <- function(out, ...) {
      om <- do.call("rbind", out)
      abline(h = 0)
      pf(a, om[, 1], ...)
      pf(a, om[, 2], col = 'blue', lty = 2, ...)
      pf(a, om[, 3], col = 'red', lty = 2, ...)
    }
    
    f_imv(out2[[i]])
    
    legend("topright", bty = 'n', lty = c(1, 2, 2),
           col = c("black", "blue", "red"), cex = 1,
           legend = c("(IRT,CDM)", "(IRT,True)", "(CDM,True)"))
  }
}
pdf("../plots/simulation1_gdina_mp.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_gdina_1$a, result_gdina_1$out2)
if (!is.null(file)) dev.off()

pdf("../plots/simulation1_dina_mp.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_dina_1$a, result_dina_1$out2)
if (!is.null(file)) dev.off()

plot_imv_only_panels <- function(a, out2) {
  # Open new graphics window to reset layout (works in RStudio)
  if (dev.cur() != 1) dev.new(width = 10, height = 8)
  
  par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4), mgp = c(2, 1, 0))
  
  panel_titles <- c(
    "DINA with mp", 
    "DINA with MAP", 
    "G-DINA with mp", 
    "G-DINA with MAP"
  )
  
  keys <- c("DINA_mp", "DINA_MAP", "GDINA_mp", "GDINA_MAP")
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x)
    lines(x, predict(m), ..., lwd = 2)
  }
  
  for (i in seq_along(keys)) {
    plot(NULL, xlim = c(0, 3), ylim = c(-0.05, 0.2),
         xlab = 'a', ylab = 'IMV', main = panel_titles[i], cex.main = 1)
    mtext(paste0("(", letters[i], ")"), side = 3, adj = 0, cex = 1)
    abline(h = 0)
    
    out <- out2[[keys[i]]]
    if (!is.null(out)) {
      om <- do.call("rbind", out)
      pf(a, om[, 1], col = 'black')
      pf(a, om[, 2], col = 'blue', lty = 2)
      pf(a, om[, 3], col = 'red', lty = 2)
      
      legend("topright", bty = 'n', lty = c(1, 2, 2),
             col = c("black", "blue", "red"), cex = 0.9,
             legend = c("(IRT,CDM)", "(IRT,True)", "(CDM,True)"))
    } else {
      warning(paste("No data for", keys[i]))
    }
  }
}

pdf("../plots/simulation1_dina_gdina_mp_MAP.pdf", width = 6, height = 6)
plot_imv_only_panels(result_compare_estimator$a, result_compare_estimator$out2)
if (!is.null(file)) dev.off()


##legend("topright",bty='n',lwd=1,lty=c(1,2),title=expression(delta),legend=c(0,1.5))
##dev.off()

##legend("topright",bty='n',lwd=1,lty=c(1,2),title=expression(delta),legend=c(0,1.5))
##dev.off()
