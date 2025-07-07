cdm.pr.marg.response <- function(resp, qm, 
                                 modeltype = c("DINA", "DINO", "GDINA"), 
                                 method = "mp", 
                                 cv = FALSE, 
                                 nfolds = 5, 
                                 seed = 123) {
  # Ensure modeltype is one of the allowed values
  modeltype <- match.arg(modeltype)
  
  # Load required packages
  library(GDINA)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(mirt)
  
  #### Branch 1: Overfitting (Full Data) Mode
  if (!cv) {
    # 1. Fit CDM on full data and compute predicted probabilities.
    m <- GDINA(resp, qm, modeltype)
    map <- personparm(m, what = method)[, 1:ncol(qm)]
    p <- list()
    
    if(modeltype %in% c("DINA", "DINO")){
      gs <- coef(m, what = "gs")
      for(i in 1:nrow(qm)){
        ii <- which(qm[i, ] == 1)
        z <- map[, ii, drop = FALSE]
        if(modeltype == "DINO"){
          rm_dino <- 1 - apply(1 - z, 1, prod)
          p[[i]] <- (1 - gs[i, 2]) * rm_dino + (1 - rm_dino) * gs[i, 1]
        } else {
          rm <- apply(z, 1, prod)
          p[[i]] <- (1 - gs[i, 2]) * rm + (1 - rm) * gs[i, 1]
        }
      }
    } else if(modeltype == "GDINA"){
      co <- coef(m)
      for(i in 1:nrow(qm)){
        ii <- which(qm[i, ] == 1)
        z <- map[, ii, drop = FALSE]
        nms <- names(co[[i]])
        cats <- gsub(")", "", gsub("P(", "", nms, fixed = TRUE), fixed = TRUE)
        cats <- strsplit(cats, "")
        cats <- lapply(cats, as.numeric)
        cats <- do.call("rbind", cats)
        pr <- list()
        for(j in 1:nrow(cats)){
          y <- cats[j, ]
          z2 <- z
          for(k in 1:ncol(z2)){
            z2[, k] <- ifelse(y[k] == 1, z[, k], 1 - z[, k])
          }
          pr[[j]] <- apply(z2, 1, prod)
        }
        pr <- do.call("cbind", pr)
        p[[i]] <- pr %*% matrix(co[[i]], ncol = 1)
      }
    } else {
      stop("Invalid model type. Please choose 'DINA', 'DINO' or 'GDINA'.")
    }
    
    # Combine predicted vectors into a matrix.
    p.cdm <- do.call("cbind", p)
    colnames(p.cdm) <- colnames(resp)
    
    # 2. Convert CDM predictions to long format and merge with true responses.
    pred_cdm_long <- as.data.frame(p.cdm) %>%
      tibble::rownames_to_column(var = "id") %>%
      pivot_longer(cols = -id, names_to = "item", values_to = "predicted_cdm")
    
    true_resp_long <- as.data.frame(resp) %>%
      tibble::rownames_to_column(var = "id") %>%
      pivot_longer(cols = -id, names_to = "item", values_to = "true_resp")
    
    pred_cdm_long <- left_join(pred_cdm_long, true_resp_long, by = c("id", "item"))
    
    #### 3. Fit a 2PL IRT model on the full data and compute IRT predicted probabilities.
    mod.irt <- mirt(resp, 1, itemtype = "2PL")
    thetas <- as.data.frame(fscores(mod.irt))
    params <- data.frame(coef(mod.irt, IRTpars = TRUE, simplify = TRUE))
    
    # Use your tested code block (do not modify)
    pred_irt_long <- tibble(th = thetas[, 1], id = rownames(resp)) %>%
      tidyr::crossing(tibble::rownames_to_column(params, var = "item")) %>%
      mutate(prediction_irt = prob.irt(items.a, items.b, th)) %>%
      dplyr::select(id, item, prediction_irt)
    
    #### 4. Merge the CDM and IRT predictions.
    combined_long <- left_join(pred_cdm_long, pred_irt_long, by = c("id", "item"))
    
    #### 5. Compute IMV statistic over the full data.
    imv_stats <- imv.binary(as.numeric(combined_long$true_resp),
                            as.numeric(combined_long$prediction_irt),
                            as.numeric(combined_long$predicted_cdm))
    
    return(list(cv = FALSE,
                p.cdm = p.cdm,
                pred_cdm_long = pred_cdm_long,
                pred_irt_long = pred_irt_long,
                combined_long = combined_long,
                imv_stats = imv_stats))
    
  } else {
    #### Branch 2: Cross-Validation Mode
    set.seed(seed)
    n <- nrow(resp)
    J <- ncol(resp)
    folds_mat <- matrix(NA, nrow = n, ncol = J)
    for (i in 1:n) {
      folds_mat[i, ] <- sample(rep(1:nfolds, length.out = J))
    }
    
    predictions_by_fold <- vector("list", nfolds)
    fold_predictions_df <- data.frame()
    
    for (f in 1:nfolds) {
      resp_train <- resp
      mask_indices <- (folds_mat == f)
      resp_train[mask_indices] <- NA
      
      m <- GDINA(resp_train, qm, modeltype)
      full_map <- personparm(m, what = method)[, 1:ncol(qm)]
      p_pred <- matrix(NA, nrow = n, ncol = J)
      
      if (modeltype %in% c("DINA", "DINO")) {
        gs <- coef(m, what = "gs")
        for (i in 1:nrow(qm)) {
          ii <- which(qm[i, ] == 1)
          z <- full_map[, ii, drop = FALSE]
          if (modeltype == "DINO") {
            rm_dino <- 1 - apply(1 - z, 1, prod)
            p_item <- (1 - gs[i, 2]) * rm_dino + (1 - rm_dino) * gs[i, 1]
          } else {
            rm <- apply(z, 1, prod)
            p_item <- (1 - gs[i, 2]) * rm + (1 - rm) * gs[i, 1]
          }
          p_pred[, i] <- p_item
        }
      } else if (modeltype == "GDINA") {
        co <- coef(m)
        for (i in 1:nrow(qm)) {
          ii <- which(qm[i, ] == 1)
          z <- full_map[, ii, drop = FALSE]
          nms <- names(co[[i]])
          cats <- gsub(")", "", gsub("P(", "", nms, fixed = TRUE), fixed = TRUE)
          cats <- strsplit(cats, "")
          cats <- lapply(cats, as.numeric)
          cats <- do.call("rbind", cats)
          pr <- list()
          for (j in 1:nrow(cats)) {
            y <- cats[j, ]
            z2 <- z
            for (k in 1:ncol(z2)) {
              z2[, k] <- ifelse(y[k] == 1, z[, k], 1 - z[, k])
            }
            pr[[j]] <- apply(z2, 1, prod)
          }
          pr <- do.call("cbind", pr)
          p_item <- pr %*% matrix(co[[i]], ncol = 1)
          p_pred[, i] <- p_item
        }
      } else {
        stop("Invalid model type. Please choose 'DINA', 'DINO' or 'GDINA'.")
      }
      
      heldout_idx <- which(mask_indices, arr.ind = TRUE)
      heldout_predictions <- p_pred[mask_indices]
      heldout_true <- resp[mask_indices]
      
      if (!is.null(rownames(resp))) {
        ids <- rownames(resp)[heldout_idx[, "row"]]
      } else {
        ids <- heldout_idx[, "row"]
      }
      if (!is.null(colnames(resp))) {
        items <- colnames(resp)[heldout_idx[, "col"]]
      } else {
        items <- heldout_idx[, "col"]
      }
      
      df_fold <- data.frame(
        fold = f,
        row = heldout_idx[, "row"],
        col = heldout_idx[, "col"],
        id = as.character(ids),
        item = as.character(items),
        true_resp = heldout_true,
        predicted_cdm = heldout_predictions,
        stringsAsFactors = FALSE
      )
      
      predictions_by_fold[[f]] <- df_fold
      fold_predictions_df <- rbind(fold_predictions_df, df_fold)
    }
    
    #### IRT predictions (same as above)
    mod.irt <- mirt(resp, 1, itemtype = "2PL")
    thetas <- as.data.frame(fscores(mod.irt))
    params <- data.frame(coef(mod.irt, IRTpars = TRUE, simplify = TRUE))
    # Use the tested code block for IRT predictions:
    pred_irt_long <- tibble(th = thetas[, 1], id = rownames(resp)) %>%
      tidyr::crossing(tibble::rownames_to_column(params, var = "item")) %>%
      mutate(prediction_irt = prob.irt(items.a, items.b, th)) %>%
      dplyr::select(id, item, prediction_irt)
    
    for (f in seq_along(predictions_by_fold)) {
      predictions_by_fold[[f]] <- predictions_by_fold[[f]] %>%
        mutate(id = as.character(id),
               item = as.character(item)) %>%
        left_join(pred_irt_long, by = c("id", "item"))
    }
    
    fold_predictions_df <- fold_predictions_df %>%
      mutate(id = as.character(id),
             item = as.character(item)) %>%
      left_join(pred_irt_long, by = c("id", "item"))
    
    # Compute IMV for each fold, then average.
    imv_folds <- sapply(predictions_by_fold, function(df) {
      imv.binary(as.numeric(as.character(df$true_resp)),
                 as.numeric(as.character(df$prediction_irt)),
                 as.numeric(as.character(df$predicted_cdm)))
    })
    imv_stats <- mean(imv_folds, na.rm = TRUE)
    
    return(list(cv = TRUE,
                predictions_by_fold = predictions_by_fold,
                fold_predictions_df = fold_predictions_df,
                folds_mat = folds_mat,
                pred_irt_long = pred_irt_long,
                imv_stats = imv_stats))
  }
}

#-------------------------------
# Helper function for IRT probability (must be defined before use)
prob.irt <- function(a, b, theta) {
  1 / (1 + exp(-a * (theta - b)))
}

# load the dataset
setwd("~/desktop")
load("CDM_AlcoholApplication_Data.Rdata")
df <- data0
resp <- df[,41:80]
qm <- Q0

# imv results function
get_imv_stats <- function(resp, qm, 
                          methods = c("mp", "MAP"), 
                          model_types = c("DINA", "DINO", "GDINA"),
                          cv = FALSE, nfolds = 5, seed = 123) {
  # Create an empty list to store results
  imv_results <- list()
  
  for (meth in methods) {
    imv_results[[meth]] <- list()
    for (mod in model_types) {
      cat("Processing method =", meth, "and model =", mod, "\n")
      # Call the merged function for this method/model combination.
      res <- cdm.pr.marg.response(resp, qm, modeltype = mod, method = meth, cv = cv, nfolds = nfolds, seed = seed)
      # Store the computed imv_stats in the list
      imv_results[[meth]][[mod]] <- res$imv_stats
    }
  }
  return(imv_results)
}

# For overfitting (full-data) mode:
imv_overfit <- get_imv_stats(resp, qm, methods = c("mp", "MAP"), model_types = c("DINA", "DINO", "GDINA"), cv = FALSE)
print(imv_overfit)

# For cross-validation mode:
imv_cv <- get_imv_stats(resp, qm, methods = c("mp", "MAP"), model_types = c("DINA", "DINO", "GDINA"), cv = TRUE, nfolds = 5, seed = 123)
print(imv_cv)

