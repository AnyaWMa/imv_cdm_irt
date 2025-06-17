remotes::install_github("hansorlee/irwpkg")
o
source("00funs.R") ##https://github.com/AnyaWMa/IRW-Qmatrix/blob/main/bd/00funs.R
library(GDINA)
library(parallel)
################################3
##simulate data with naughty cdm
set.seed(123)
Q <- sim30GDINA$simQ
J <- nrow(Q)
K <- ncol(Q)

cdm.sim<-function(a, N=500, modeltype = "DINA") {
  Q <- sim30GDINA$simQ
  J <- nrow(Q)
  K <- ncol(Q)
  
  # set item difficulty
  gs <- matrix(runif(J*2,0,0.3),J,2)
  cutoffs <- qnorm(c(1:K)/(K+1))
  m <- rep(0,K)
  
  vcov <- matrix(a,K,K)
  diag(vcov) <- 1
  
  sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "mvnorm",
                  mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))
  
  # test data - same p for each cell
  sim.test <- simGDINA(N,Q, model = modeltype, catprob.parm = sim$catprob.parm,attribute = sim$attribute)
  
  #prediction from true cdm
  truep <- sim$LCprob.parm[sim$att.group,]
  
  irtfit<-mirt(data.frame(sim$dat),1,'2PL')
  
  # predictions from irt
  th.est<-fscores(irtfit)
  est<-coef(irtfit,simplify=TRUE,IRTpars=TRUE)$items
  k<-outer(th.est[,1],est[,2],'-')
  k<-matrix(est[,1],nrow=N,ncol=J,byrow=TRUE)*k
  irtp<-1/(1+exp(-k))
  
  cdmp <- cdm.pr.marg.update(sim$dat,Q,modeltype,estmethod = "MAP")
  
  pmp <- cdm.pr.marg.update(sim$dat,Q,modeltype,estmethod = "mp")
  
  om<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(cdmp))
  om.mp <- imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(pmp))
  oracle.irt<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(truep))
  oracle.cdm<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(cdmp),as.numeric(truep))
  oracle.cdm.mp<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(pmp),as.numeric(truep))
  
  rms<-function(x) sqrt(mean(x^2))
  rms.irt<-rms(as.numeric(sim.test$dat)-as.numeric(irtp))
  rms.cdm<-rms(as.numeric(sim.test$dat)-as.numeric(cdmp))
  rms.cdm.mp <- rms(as.numeric(sim.test$dat)-as.numeric(pmp))
  rms.irt.p<-rms(as.numeric(truep)-as.numeric(irtp))
  rms.cdm.p<-rms(as.numeric(truep)-as.numeric(cdmp))
  rms.cdm.mp.p<-rms(as.numeric(truep)-as.numeric(pmp))
  
  c(om=om,
    oracle.irt=oracle.irt,
    oracle.cdm=oracle.cdm,
    rms.irt=rms.irt,
    rms.cdm=rms.cdm,
    rms.irt.p=rms.irt.p,
    rms.cdm.p=rms.cdm.p, 
    om.mp = om.mp, 
    oracle.cdm.mp = oracle.cdm.mp, 
    rms.cdm.mp = rms.cdm.mp,
    rms.cdm.mp.p = rms.cdm.mp.p
    )
}

simulate_scenario_1 <- function(modeltype = "GDINA") {
  # Generate sorted vector a
  a <- sort(runif(n = 100, min = 0, max = 0.8))
  
  # Initialize output list
  out2 <- list()
  
  # Loop over bounds and apply cdm.sim in parallel
  for (i in c(200, 500, 1000)) {
    print(i)
    out2[[as.character(i)]] <- mclapply(
      a,
      cdm.sim,
      N = i,
      modeltype = modeltype
    )
  }
  
  # Return both a and out2
  return(list(a = a, out2 = out2))
}


## Simulate DINA
result_dina_1 <- simulate_scenario_1(modeltype = "DINA")
save(result_dina_1 , file = "simulation_data_update/scenario1_out_dina.RData")

## Simulate GDINA
result_gdina_1 <- simulate_scenario_1(modeltype = "GDINA")
save(result_gdina_1 , file = "simulation_data_update/scenario1_out_gdina.RData")


## Plotting
plot_rmse_imv_panels <- function(a, out2) {
  par(mgp = c(2, 1, 0), mfrow = c(3, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4))
  
  panel_titles <- c("RMSE (N = 200)", "IMV (N = 200)", 
                    "RMSE (N = 500)", "IMV (N = 500)", 
                    "RMSE (N = 1000)", "IMV (N = 1000)")
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x, family = "symmetric")
    lines(x, predict(m), ..., lwd = 2)
  }
  
  for (i in seq_along(out2)) {
    
    ## RMSE Plot
    plot(NULL, xlim = c(0, 0.8), ylim = c(0, 0.6), 
         xlab = expression(rho), ylab = 'RMSE(oos resp, est)')
    title(panel_titles[2 * i - 1], line = 0.5, cex.main = 1)
    #mtext(paste0("(", letters[2 * i - 1], ")"), side = 3, adj = 0, cex = 1)
    
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
    plot(NULL, xlim = c(0, 0.8), ylim = c(-0.05, 0.2),
         xlab = expression(rho), ylab = 'IMV')
    title(panel_titles[2 * i], line = 0.5, cex.main = 1)
    #mtext(paste0("(", letters[2 * i], ")"), side = 3, adj = 0, cex = 1)
    
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
pdf("plots_update/simulation1_gdina_map.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_gdina_1$a, result_gdina_1$out2)
if (!is.null(file)) dev.off()

pdf("plots_update/simulation1_dina_map.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_dina_1$a, result_dina_1$out2)
if (!is.null(file)) dev.off()

# for mp ploting
plot_rmse_imv_panels_mp <- function(a, out2) {
  par(mgp = c(2, 1, 0), mfrow = c(3, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4))
  
  panel_titles <- c("RMSE (N = 200)", "IMV (N = 200)", 
                    "RMSE (N = 500)", "IMV (N = 500)", 
                    "RMSE (N = 1000)", "IMV (N = 1000)")
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x,  family = "symmetric")
    lines(x, predict(m), ..., lwd = 2)
  }
  
  for (i in seq_along(out2)) {
    
    ## RMSE Plot
    plot(NULL, xlim = c(0, 0.8), ylim = c(0, 0.6), 
         xlab = expression(rho), ylab = 'RMSE(oos resp, est)')
    title(panel_titles[2 * i - 1], line = 0.5, cex.main = 1)
    #mtext(paste0("(", letters[2 * i - 1], ")"), side = 3, adj = 0, cex = 1)
    
    f_rmse <- function(out, ...) {
      z <- do.call("rbind", out)
      irt <- z[, 4]
      cdm <- z[, 10]
      irt.p <- z[, 6]
      cdm.p <- z[, 11]
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
    plot(NULL, xlim = c(0, 0.8), ylim = c(-0.05, 0.2),
         xlab = expression(rho), ylab = 'IMV')
    title(panel_titles[2 * i], line = 0.5, cex.main = 1)
    #mtext(paste0("(", letters[2 * i], ")"), side = 3, adj = 0, cex = 1)
    
    f_imv <- function(out, ...) {
      om <- do.call("rbind", out)
      abline(h = 0)
      pf(a, om[, 8], ...)
      pf(a, om[, 2], col = 'blue', lty = 2, ...)
      pf(a, om[, 9], col = 'red', lty = 2, ...)
    }
    
    f_imv(out2[[i]])
    
    legend("topright", bty = 'n', lty = c(1, 2, 2),
           col = c("black", "blue", "red"), cex = 1,
           legend = c("(IRT,CDM)", "(IRT,True)", "(CDM,True)"))
  }
}

pdf("plots_update/simulation1_gdina_mp.pdf", width = 6, height = 8)
plot_rmse_imv_panels_mp(result_gdina_1$a, result_gdina_1$out2)
if (!is.null(file)) dev.off()

pdf("plots_update/simulation1_dina_mp.pdf", width = 6, height = 8)
plot_rmse_imv_panels_mp(result_dina_1$a, result_dina_1$out2)
if (!is.null(file)) dev.off()

plot_imv_only_panels <- function(result_dina, result_gdina) {
  result_dina_a <- result_dina$a
  result_gdina_a <- result_gdina$a
  
  result_dina_om <- do.call("rbind", result_dina$out2[['500']])
  result_gdina_om <- do.call("rbind", result_gdina$out2[['500']])
  
  par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4), mgp = c(2, 1, 0))
  
  panel_titles <- c("DINA with MAP", "G-DINA with MAP", "DINA with mp", "G-DINA with mp")
  
  pf <- function(x, y, ...) {
    m <- loess(y ~ x, family = "symmetric")
    lines(x, predict(m), ..., lwd = 2)
  }
  
  # Define source data for each panel
  plot_configs <- list(
    list(a = result_dina_a, om = result_dina_om, cols = c(1, 2, 3)),
    list(a = result_gdina_a, om = result_gdina_om, cols = c(1, 2, 3)),
    list(a = result_dina_a, om = result_dina_om, cols = c(8, 2, 9)),
    list(a = result_gdina_a, om = result_gdina_om, cols = c(8, 2, 9))
  )
  
  for (i in seq_along(plot_configs)) {
    config <- plot_configs[[i]]
    a <- config$a
    om <- config$om
    cols <- config$cols
    
    plot(NULL, xlim = c(0, 0.8), ylim = c(-0.05, 0.2),
         xlab = expression(rho), ylab = "IMV", main = panel_titles[i], cex.main = 1)
  #  mtext(paste0("(", letters[i], ")"), side = 3, adj = 0, cex = 1)
    abline(h = 0)
    
    pf(a, om[, cols[1]], col = 'black')
    pf(a, om[, cols[2]], col = 'blue', lty = 2)
    pf(a, om[, cols[3]], col = 'red', lty = 2)
    
    legend("topright", bty = 'n', lty = c(1, 2, 2),
           col = c("black", "blue", "red"), cex = 0.9,
           legend = c("(IRT,CDM)", "(IRT,True)", "(CDM,True)"))
  }
}

pdf("plots_update/simulation1_dina_gdina_mp_MAP.pdf", width = 6, height = 6)
plot_imv_only_panels(result_dina_1, result_gdina_1)
if (!is.null(file)) dev.off()

