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
  
  #warning1 <- c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)
  #warning2 <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,-Inf, -Inf)
  
  ## =========================================================================
  ## [Added by Yiqing] Change the length of warning match the length of output
  ## =========================================================================
  warning1 <- rep(Inf, 24)
  warning2 <- rep(-Inf, 24)
  ############################################################################
  
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
  
  if (!extract.mirt(irtfit , 'converged')) {
    print("irt not converged")
    return (warning1)
  }
  
  # predictions from irt
  th.est<-fscores(irtfit)
  est<-coef(irtfit,simplify=TRUE,IRTpars=TRUE)$items
  k<-outer(th.est[,1],est[,2],'-')
  k<-matrix(est[,1],nrow=N,ncol=J,byrow=TRUE)*k
  irtp<-1/(1+exp(-k))
  
  ## ============================================================
  ## [Added by Yiqing] Unrestricted K-dimensional M2PL (no Q constraints)
  ## Goal: fit an exploratory multidimensional 2PL on the *same training data*
  ##       used in Anya's simulation (sim$dat), then extract predicted P(Y=1)
  ##       for downstream IMV comparisons on the OOS test set (sim.test$dat).
  ## Note: For K>1, we do NOT use the 1D hand-coded logistic formula; instead,
  ##       we extract probabilities via probtrace(), which is parameterization-safe.
  ## Output: p_m2pl_free (N x J) predicted probabilities on the training persons.
  ## ============================================================
  
  # K is already defined: K <- ncol(Q), I used MHRM here for better convergence
  # in higher dimensions; can switch to others if preferred.
  m2pl_free_fit <- try(
    mirt(data.frame(sim$dat), K, itemtype = "2PL",
         method = "MHRM", verbose = FALSE),
    silent = TRUE
  )
  
  ## --- guard: M2PL-free convergence ---
  if (inherits(m2pl_free_fit, "try-error") || !extract.mirt(m2pl_free_fit, "converged")) {
    print("M2PL-free not converged")
    return(warning1)
  }
  
  ## --- extract predicted probabilities (align with Anya's probtrace usage) ---
  th_free <- fscores(m2pl_free_fit, full.scores = TRUE)  # N x K
  p_all   <- probtrace(m2pl_free_fit, th_free)           # N x (2J)
  
  pcols <- grep("\\.P\\.1$", colnames(p_all))
  if (length(pcols) != J) pcols <- grep("\\.P\\.2$", colnames(p_all))  # fallback
  
  p_m2pl_free <- p_all[, pcols, drop = FALSE]            # N x J
  colnames(p_m2pl_free) <- colnames(data.frame(sim$dat))
  
  ## sanity check: predicted probabilities in [0,1]
  cat("[M2PL-free] p range:", range(p_m2pl_free), "\n")
  ## End of added block
  ####################################################################################
  
  ## ============================================================
  ## [Added by Yiqing] Unrestricted K-dimensional M2PL with PRIORS
  ## Goal: same as M2PL-free, but add weakly-informative priors on slopes
  ##       to stabilize high-dimensional estimation (K=5) without Q constraints.
  ## Note: PRIOR is applied to a1..aK across all items. Here we use lognormal
  ##       to encourage positive slopes (as in our earlier 2D setup).
  ## Output: p_m2pl_free_prior (N x J)
  ## ============================================================
  
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
  
  # optional sanity
  # cat("[M2PL-free-prior] p range:", range(p_m2pl_free_prior), "\n")
  ## End of added block
  ####################################################################################
  
  ## ============================================================
  ## [Added by Yiqing] Q-aligned confirmatory K-dimensional M2PL
  ## Goal: fit a multidimensional 2PL whose loading pattern is aligned with
  ##       the CDM Q-matrix (Q=0 => loading fixed 0; Q=1 => freely estimated),
  ##       then extract predicted P(Y=1) via probtrace().
  ## Output: p_m2pl_q (N x J) predicted probabilities on the training persons.
  ## ============================================================
  
  # --- build confirmatory mirt model syntax from Q ---
  lines <- character(0)
  for (k in 1:K) {
    items_k <- which(Q[, k] == 1L)
    if (length(items_k) == 0) {
      print(paste0("Q-aligned M2PL: dimension ", k, " has no items (all zeros)."))
      return(warning1)
    }
    lines <- c(lines, paste0("F", k, " = ", paste(items_k, collapse = ",")))
  }
  
  # allow factor covariances (to match generating correlation structure)
  if (K >= 2) {
    cov_terms <- c()
    for (k1 in 1:(K - 1)) for (k2 in (k1 + 1):K) cov_terms <- c(cov_terms, paste0("F", k1, "*F", k2))
    lines <- c(lines, paste0("COV = ", paste(cov_terms, collapse = ",")))
  }
  
  mod_q <- mirt.model(paste(lines, collapse = "\n"))
  
  # --- fit Q-aligned M2PL ---
  m2pl_q_fit <- try(
    mirt(data.frame(sim$dat), mod_q, itemtype = "2PL",
         method = "MHRM", verbose = FALSE),
    silent = TRUE
  )
  
  ## --- guard: M2PL-Q convergence ---
  if (inherits(m2pl_q_fit, "try-error") || !extract.mirt(m2pl_q_fit, "converged")) {
    print("M2PL-Q (Q-aligned) not converged")
    return(warning1)
  }
  
  ## --- extract predicted probabilities ---
  th_q  <- fscores(m2pl_q_fit, full.scores = TRUE)  # N x K
  p_all <- probtrace(m2pl_q_fit, th_q)              # N x (2J)
  
  pcols <- grep("\\.P\\.1$", colnames(p_all))
  if (length(pcols) != J) pcols <- grep("\\.P\\.2$", colnames(p_all))  # fallback
  
  p_m2pl_q <- p_all[, pcols, drop = FALSE]          # N x J
  colnames(p_m2pl_q) <- colnames(data.frame(sim$dat))
  
  # optional quick sanity
  cat("[M2PL-Q] p range:", range(p_m2pl_q), "\n")
  ## End of added block
  ####################################################################################
  
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
  
  ## ============================================================
  ## [Added by Yiqing] Compare new M2PL predictions with CDM / oracle
  ## ============================================================
  
  # IMV: new M2PL vs CDM (MAP)
  om.m2free <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_free), as.numeric(cdmp))
  om.m2q    <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_q),    as.numeric(cdmp))
  
  # Oracle IMV: new M2PL vs true probability
  oracle.m2free <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_free), as.numeric(truep))
  oracle.m2q    <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_q),    as.numeric(truep))
  
  # m2free with prior
  om.m2free.prior <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_free_prior), as.numeric(cdmp))
  oracle.m2free.prior <- imv::imv.binary(as.numeric(sim.test$dat), as.numeric(p_m2pl_free_prior), as.numeric(truep))
  rms.m2free.prior <- rms(as.numeric(sim.test$dat) - as.numeric(p_m2pl_free_prior))
  rms.m2free.prior.p <- rms(as.numeric(truep) - as.numeric(p_m2pl_free_prior))
  
  # RMSE vs test responses
  rms.m2free <- rms(sim.test$dat - p_m2pl_free)
  rms.m2q    <- rms(sim.test$dat - p_m2pl_q)
  
  # RMSE vs oracle probability
  rms.m2free.p <- rms(truep - p_m2pl_free)
  rms.m2q.p    <- rms(truep - p_m2pl_q)
  # End of added block
  ####################################################################################
  
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
    oracle.map.mp = oracle.map.mp,
    
    ## ---- [Added by Yiqing] new M2PL outputs ----
    om.m2free = om.m2free,
    oracle.m2free = oracle.m2free,
    rms.m2free = rms.m2free,
    rms.m2free.p = rms.m2free.p,
    om.m2q = om.m2q,
    oracle.m2q = oracle.m2q,
    rms.m2q = rms.m2q,
    rms.m2q.p = rms.m2q.p,
    
    om.m2free.prior = om.m2free.prior,
    oracle.m2free.prior = oracle.m2free.prior,
    rms.m2free.prior = rms.m2free.prior,
    rms.m2free.prior.p = rms.m2free.prior.p
    ## -------------------------------------------
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
      modeltype = modeltype, mc.cores = 1
    )
  }
  
  # Return both a and out2
  return(list(a = a, out2 = out2))
}


## Simulate DINA
##### Rest of your original code #####