#remotes::install_github("hansorlee/irwpkg")

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
  
  warning1 <- c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)
  warning2 <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,-Inf, -Inf)
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
  
  # --- build "free" confirmatory model + priors (all items load on all factors) ---
  cov_terms <- character(0)
  if (K >= 2) {
    for (k1 in 1:(K-1)) for (k2 in (k1+1):K) cov_terms <- c(cov_terms, paste0("F", k1, "*F", k2))
  }
  
  prior_terms <- paste0("(1-", J, ", a", 1:K, ", lnorm, 0.0, 1.0)", collapse = ",")
  
  s_prior <- paste0(
    paste0("F", 1:K, " = 1-", J, collapse = ",\n"),
    if (K >= 2) paste0("\nCOV = ", paste(cov_terms, collapse = ",")) else "",
    "\nPRIOR = ", prior_terms
  )
  
  mod_free_prior <- mirt.model(s_prior)
  
  # --- fit M2PL-free with prior ---
  m2pl_free_prior_fit <- try(
    mirt(data.frame(sim$dat), mod_free_prior, itemtype = "2PL",
         method = "MHRM", verbose = FALSE),
    silent = TRUE
  )
  
  ## --- guard ---
  if (inherits(m2pl_free_prior_fit, "try-error") || !extract.mirt(m2pl_free_prior_fit, "converged")) {
    print("M2PL-free-with-prior not converged")
    return(warning1)
  }
  
  ## --- extract predicted probabilities ---
  th_free_prior <- fscores(m2pl_free_prior_fit, full.scores = TRUE)
  p_all <- probtrace(m2pl_free_prior_fit, th_free_prior)
  
  pcols <- grep("\\.P\\.1$", colnames(p_all))
  if (length(pcols) != J) pcols <- grep("\\.P\\.2$", colnames(p_all))
  
  p_m2pl_free_prior <- p_all[, pcols, drop = FALSE]
  colnames(p_m2pl_free_prior) <- colnames(data.frame(sim$dat))
  
  irtp = p_m2pl_free_prior
  
  cdmp <- cdm.pr.marg.update(sim$dat,Q,modeltype,estmethod = "MAP")
  
  if (length(cdmp)==1) {
    return (warning2)
  }
  
  pmp <- cdm.pr.marg.update(sim$dat,Q,modeltype,estmethod = "mp")
  if (length(pmp)==1) {
    return (warning2)
  }
  
  om<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(cdmp))
  om.mp <- imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(pmp))
  oracle.irt<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(irtp),as.numeric(truep))
  oracle.cdm<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(cdmp),as.numeric(truep))
  oracle.cdm.mp<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(pmp),as.numeric(truep))
  oracle.map.mp<-imv::imv.binary(as.numeric(sim.test$dat),as.numeric(cdmp),as.numeric(pmp))
  
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
    rms.cdm.mp.p = rms.cdm.mp.p, 
    oracle.map.mp = oracle.map.mp
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
result_dina_1_multi_irt <- simulate_scenario_1(modeltype = "DINA")
save(result_dina_1_multi_irt , file = "simulation_data_update/scenario1_out_dina_multi_irt.RData")

sim_sizes <- c("200", "500", "1000")

## Plotting
plot_rmse_imv_panels <- function(a, out2) {
  par(mgp = c(2, 1, 0), mfrow = c(3, 2), mar = c(3, 3, 2, 1), oma = rep(0.5, 4))
  
  panel_titles <- c("RMSE (N = 200)", "IMV (N = 200)", 
                    "RMSE (N = 500)", "IMV (N = 500)", 
                    "RMSE (N = 1000)", "IMV (N = 1000)")
  
  pf <- function(x, y, ...) {
    valid <- is.finite(y)
    x <- x[valid]
    y <- y[valid]
    m <- loess(y ~ x)
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
      pf(a, irt, col = '#33BFD5', ...)
      pf(a, cdm, col = '#D22730', ...)
      pf(a, irt.p, col = '#33BFD5', lty = 2, ...)
      pf(a, cdm.p, col = '#D22730', lty = 2, ...)
    }
    
    f_rmse(out2[[i]])
    
    legend("topright", bty = 'n', lty = c(1, 1, 2, 2),
           col = c("#33BFD5", "#D22730", "#33BFD5", "#D22730"), cex = 1,
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
      pf(a, om[, 2], col = '#33BFD5', lty = 2, ...)
      pf(a, om[, 3], col = '#D22730', lty = 2, ...)
    }
    
    f_imv(out2[[i]])
    
    legend("topright", bty = 'n', lty = c(1, 2, 2),
           col = c("black", "#33BFD5", "#D22730"), cex = 1,
           legend = c("(IRT,CDM)", "(IRT,True)", "(CDM,True)"))
  }
}


pdf("plots_update/simulation1_dina_map_multi_irt.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_dina_1_multi_irt$a, result_dina_1_multi_irt$out2)
if (!is.null(file)) dev.off()
