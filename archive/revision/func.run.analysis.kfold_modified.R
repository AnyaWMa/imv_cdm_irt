## OOS prediction
func.run.analysis.kfold <- function(df, cdm_model = "DINA", irt_model = "2PL", if_oos = TRUE, kfold = 5, cdm_pr_estimate = cdm.pr.marg.update) {
  n_items <- n_distinct(df$item)
  item_list <- sort(unique(df$item))
  
  df.cdm.structure <- func.convert.Q.data.structure.with.k.fold(df, kfold = kfold)
  
  func.get.prediction.result <- function(df.tmp, dat, qm){
    
    df.ids <- df.tmp %>% group_by(id) %>%
      tally()
    
    df.dat <- dat %>%
      dplyr ::select(-gr) %>%
      dplyr :: arrange(item) %>%
      pivot_wider(names_from = item, values_from = resp) %>%
      dplyr ::arrange(id)
    
    df.dat.no.id <- df.dat %>% dplyr :: select(-id)
    
    # IRT models
    # ---- add by Yiqing：support m2pl_prior and m2pl_q_constrained，output prediction.irt ----
    if (irt_model %in% c("m2pl_prior")) {
      
      library(mirt)
      
      resp_mat <- df.dat.no.id
      Q <- qm
      
      ############# fit model with priors (m2pl_prior) #############
      J <- ncol(resp_mat)
      K <- ncol(Q)
      
      # COV terms
      cov_terms <- character(0)
      if (K >= 2) {
        for (k1 in 1:(K-1)) for (k2 in (k1+1):K) cov_terms <- c(cov_terms, paste0("F", k1, "*F", k2))
      }
      
      # PRIOR terms: (1-J, a1, lnorm, 0, 1), ... (1-J, aK, lnorm, 0, 1)
      prior_terms <- paste0("(1-", J, ", a", 1:K, ", lnorm, 0.0, 1.0)", collapse = ",")
      
      s_prior <- paste0(
        paste0("F", 1:K, " = 1-", J, collapse = "\n"),
        if (K >= 2) paste0("\nCOV = ", paste(cov_terms, collapse = ",")) else "",
        "\nPRIOR = ", prior_terms
      )
      
      modspec_prior <- mirt.model(s_prior)
      
      mod_m2pl_prior <- mirt(
        data     = resp_mat,
        model    = modspec_prior,
        itemtype = "2PL",
        method   = "MHRM",
        verbose  = TRUE,
        technical = list(NCYCLES = 2000)
      )
      
      if (!extract.mirt(mod_m2pl_prior, "converged")) {
        print("irt not converged")
        return(NULL)
      }
      
      # ---- predicted P(Y=1), N x J ----
      theta_prior <- fscores(mod_m2pl_prior, full.scores = TRUE)
      p_all_prior <- probtrace(mod_m2pl_prior, Theta = theta_prior)
      P_irt_prior <- p_all_prior[, grep("\\.P\\.1$", colnames(p_all_prior)), drop = FALSE]
      colnames(P_irt_prior) <- colnames(resp_mat)
      
      result.irt <- as.data.frame(P_irt_prior) %>%
        mutate(id = df.dat$id) %>%
        pivot_longer(cols = colnames(resp_mat), names_to = "item", values_to = "prediction.irt") %>%
        left_join(df.tmp, by = c("id","item")) %>%
        dplyr ::select(id, item, resp, gr, prediction.irt)
      
    } else if (irt_model %in% c("m2pl_q_constrained")) {
      
      library(mirt)
      
      resp_mat <- df.dat.no.id
      Q <- qm
      
      ############# fit model with Q-constraints (m2pl_q_constrained) #############
      J <- ncol(resp_mat)
      K <- ncol(Q)
      
      lines <- character(0)
      for (k in 1:K) {
        items_k <- which(Q[, k] != 0)
        if (length(items_k) == 0) stop(paste0("Q-constrained: Factor F", k, " has 0 items (empty column in Q)."))
        lines <- c(lines, paste0("F", k, " = ", paste(items_k, collapse = ",")))
      }
      
      # allow factor correlations
      if (K >= 2) {
        cov_terms <- character(0)
        for (k1 in 1:(K-1)) for (k2 in (k1+1):K) cov_terms <- c(cov_terms, paste0("F", k1, "*F", k2))
        lines <- c(lines, paste0("COV = ", paste(cov_terms, collapse = ",")))
      }
      
      modspec_q <- mirt.model(paste(lines, collapse = "\n"))
      
      mod_m2pl_q <- mirt(
        data     = resp_mat,
        model    = modspec_q,
        itemtype = "2PL",
        method   = "MHRM",
        verbose  = TRUE,
        technical = list(NCYCLES = 2000)
      )
      
      if (!extract.mirt(mod_m2pl_q, "converged")) {
        print("irt not converged")
        return(NULL)
      }
      
      # ---- predicted P(Y=1), N x J ----
      theta_q <- fscores(mod_m2pl_q, full.scores = TRUE)
      p_all_q <- probtrace(mod_m2pl_q, Theta = theta_q)
      P_irt_q <- p_all_q[, grep("\\.P\\.1$", colnames(p_all_q)), drop = FALSE]
      colnames(P_irt_q) <- colnames(resp_mat)
      
      result.irt <- as.data.frame(P_irt_q) %>%
        mutate(id = df.dat$id) %>%
        pivot_longer(cols = colnames(resp_mat), names_to = "item", values_to = "prediction.irt") %>%
        left_join(df.tmp, by = c("id","item")) %>%
        dplyr ::select(id, item, resp, gr, prediction.irt)
      
    } else {
      # ---- original 1d 2pl  ----
      mod.irt <- mirt(df.dat.no.id, 1, itemtype = irt_model)
      
      if (!extract.mirt(mod.irt, "converged")) {
        print("irt not converged")
        return(NULL)
      }
      
      thetas <- as.data.frame(fscores(mod.irt))
      
      params <- data.frame(coef(mod.irt, IRTpars = TRUE, simplify = TRUE))
      
      result.irt <- tibble(th = thetas$F1, id = df.dat$id) %>%
        left_join(df.tmp) %>%
        left_join(rownames_to_column(params, var = "item")) %>%
        mutate(prediction.irt = prob.irt(items.a, items.b, th)) %>%
        dplyr ::select(id, item, resp, gr, prediction.irt)
    }
    
    # CDM model
    predict.matrix.cdm.map <- cdm_pr_estimate(df.dat.no.id, qm, model = cdm_model, "MAP")
    if (is.null(predict.matrix.cdm.map)) {
      print("cdm map not converged")
      return(NULL)
    }
    
    predict.matrix.cdm.mp <- cdm_pr_estimate(df.dat.no.id, qm, model = cdm_model, "mp")
    if (is.null(predict.matrix.cdm.mp)) {
      print("cdm mp not converged")
      return(NULL)
    }
    
    
    result.cdm.map <- as.data.frame(predict.matrix.cdm.map) %>%
      set_names(item_list) %>%
      mutate(id = df.ids$id) %>%
      pivot_longer(cols = item_list, names_to = "item", values_to = "prediction.cdm.map")
    
    result.cdm.mp <- as.data.frame(predict.matrix.cdm.mp) %>%
      set_names(item_list) %>%
      mutate(id = df.ids$id) %>%
      pivot_longer(cols = item_list, names_to = "item", values_to = "prediction.cdm.mp")
    
    result.compared <- result.cdm.map %>%
      left_join(result.cdm.mp) %>%
      left_join(df.tmp) %>%
      left_join(result.irt)
  }
  
  list.imv.1 <- list()
  list.imv.2 <- list()
  list.imv.3 <- list()
  
  if (!if_oos) {
    result.compared <- func.get.prediction.result(
      df.cdm.structure$dat_kfold, 
      df.cdm.structure$dat_kfold, 
      df.cdm.structure$Q
    )
    
    if (is.null(result.compared)) return(NA)
    
    df.compare <- result.compared %>% 
      filter(!is.na(prediction.irt), 
             !is.na(prediction.cdm.map), 
             !is.na(prediction.cdm.mp))
    
    list.imv.1 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.irt, df.compare$prediction.cdm.map))
    list.imv.2 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.irt, df.compare$prediction.cdm.mp))
    list.imv.3 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.cdm.map, df.compare$prediction.cdm.mp))
    
  } else {
    for (i in seq(1, kfold, 1)) {
      df.tmp <- df.cdm.structure$dat_kfold
      
      df.train <- df.tmp %>% 
        filter(gr != i) 
      
      result.compared <- func.get.prediction.result(df.tmp, df.train, df.cdm.structure$Q)
      
      if (is.null(result.compared)) return(next)
      
      df.compare <- result.compared %>% 
        filter(gr == i) %>% 
        filter(!is.na(prediction.irt), !is.na(prediction.cdm.map), !is.na(prediction.cdm.mp))
      
      imv.1 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.irt, df.compare$prediction.cdm.map))
      imv.2 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.irt, df.compare$prediction.cdm.mp))
      imv.3 <- c(imv::imv.binary(df.compare$resp, df.compare$prediction.cdm.map, df.compare$prediction.cdm.mp))
      
      list.imv.1 <- c(list.imv.1, imv.1) 
      list.imv.2 <- c(list.imv.2, imv.2) 
      list.imv.3 <- c(list.imv.3, imv.3) 
      
    }
  }
  
  return(c(
    mean(unlist(list.imv.1)), sd(unlist(list.imv.1)),
    mean(unlist(list.imv.2)), sd(unlist(list.imv.2)),
    mean(unlist(list.imv.3)), sd(unlist(list.imv.3))
  ))
}
