remotes::install_github("hansorlee/irwpkg")
o
source("bd/00funs.R") ##https://github.com/AnyaWMa/IRW-Qmatrix/blob/main/bd/00funs.R
library(GDINA)
library(parallel)
################################
##simulate data with naughty cdm
set.seed(123)
Q <- sim30GDINA$simQ
J <- nrow(Q)
K <- ncol(Q)

# Simulation 2
cdm.sim.2<-function(a, N=500, modeltype = "DINA", attribute = "mvnorm") {
  Q <- sim30GDINA$simQ
  J <- nrow(Q)
  K <- ncol(Q)
  
  # set item difficulty
  gs <- matrix(runif(J*2,0,0.3),J,2)
  cutoffs <- qnorm(c(1:K)/(K+1))
  m <- rep(0,K)
  
  vcov <- matrix(a,K,K)
  diag(vcov) <- 1
  
  if (attribute == "independent"){
    sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "mvnorm", mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))
  }
  else if (attribute == "higher.order"){
      theta <- rnorm(N)
      lambda <- data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))
      sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "higher.order",  higher.order.parm = list(theta = theta,lambda = lambda))
    }
  else if (attribute == "linear") {
    linear <- list(c(1,2),
                   c(2,3),
                   c(3,4),
                   c(4,5))
    struc <- att.structure(linear,K)
    sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "categorical",att.prior = struc$att.prob)
  } else if (attribute == "converg"){
    converg <- list(c(1,2),
                   c(2,3),
                   c(2,4),
                   c(4,5), 
                   c(3,5))
    struc <- att.structure(converg,K)
    sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "categorical",att.prior = struc$att.prob)
  } else if (attribute == "diverg") {
    diverg <- list(c(1,2),
                   c(2,3),
                   c(1,4),
                   c(4,5))
    struc <- att.structure(diverg,K)
    sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "categorical",att.prior = struc$att.prob)
  } else {
    unstructured <- list(c(1,2),
                   c(1,3),
                   c(1,4),
                   c(1,5))
    struc <- att.structure(unstructured,K)
    sim <- simGDINA(N,Q,gs.parm = gs, model = modeltype, att.dist = "categorical",att.prior = struc$att.prob)
  }
  
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
    rms.cdm.mp.p = rms.cdm.mp.p,
    attribute = attribute
  )
}

simulate_scenario_2_attribute <- function(modeltype = "GDINA") {
  # Generate sorted vector a
  a <- sort(runif(n = 10, min = 0, max = 0))
  
  # Initialize output list
  out2 <- list()
  
  # Loop over bounds and apply cdm.sim in parallel
  types <- c("independent", "higher.order", "linear", "converg", "diverg", "unstructured")
  for (i in c(1, 2, 3, 4, 5, 6)) {
    print(i)
    out2[[as.character(i)]] <- mclapply(
      a,
      cdm.sim.2,
      N = 500,
      modeltype = modeltype,
      attribute = types[[i]]
    )
  }
  
  # Return both a and out2
  return(out2 = out2)
}

df_all <- map2_dfr(result_gdina_2 , names(result_gdina_2), function(attr_list, attr_name) {
  map_dfr(seq_along(attr_list), function(i) {
    as_tibble_row(attr_list[[i]]) %>%
      mutate(attribute = attr_name, rep = i)
  })
})

df_all_clean <- df_all %>%
  mutate(across(
    where(is.character) & !any_of(c("attribute", "rep")),
    as.numeric
  ))

df_summary <- df_all_clean %>%
  group_by(attribute) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
    .names = "{.col}_{.fn}"
  ), .groups = "drop")
## Simulate DINA
result_gdina_2 <- simulate_scenario_2_attribute(modeltype = "GDINA")
#save(result_gdina_2 , file = "simulation_data_update/scenario1_out_dina.RData")


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
         xlab = 'a', ylab = 'RMSE(oos resp, est)')
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
         xlab = 'a', ylab = 'IMV')
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
pdf("plots_update/simulation2_gdina_map.pdf", width = 6, height = 8)
plot_rmse_imv_panels(result_gdina_2$a, result_gdina_2$out2)
if (!is.null(file)) dev.off()

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
         xlab = 'a', ylab = 'RMSE(oos resp, est)')
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
         xlab = 'a', ylab = 'IMV')
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


pdf("plots_update/simulation2_gdina_mp.pdf", width = 6, height = 8)
plot_rmse_imv_panels_mp(result_gdina_2$a, result_gdina_2$out2)
if (!is.null(file)) dev.off()
