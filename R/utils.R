# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

metrics_regression <- function(y_true, y_pred) {
  ok <- !is.na(y_true) & !is.na(y_pred)
  scored <- sum(ok)
  total  <- length(y_true)
  coverage <- if (total > 0) scored / total else NA_real_
  
  if (scored == 0) {
    return(c(
      MSE      = NA_real_,
      RMSE     = NA_real_,
      MAE      = NA_real_,
      R2       = NA_real_,
      Scored   = scored,
      Total    = total,
      Coverage = coverage
    ))
  }
  
  mse  <- mean((y_true[ok] - y_pred[ok])^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(y_true[ok] - y_pred[ok]))
  
  ybar <- mean(y_true[ok])
  
  sse <- sum((y_true[ok] - y_pred[ok])^2)
  sst <- sum((y_true[ok] - ybar)^2)
  r2  <- if (sst > 0) 1 - sse / sst else NA_real_
  
  c(
    MSE      = mse,
    RMSE     = rmse,
    MAE      = mae,
    R2       = r2,
    Scored   = scored,
    Total    = total,
    Coverage = coverage
  )
}

compute_auc <- function(y_true, scores) {
  y_true <- as.integer(y_true)
  pos <- which(y_true == 1L)
  neg <- which(y_true == 0L)
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  
  r <- rank(scores, ties.method = "average")
  auc <- (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
  as.numeric(auc)
}

classification_report <- function(
    y_true,
    y_score,
    threshold = 0.5,
    header = NULL,
    quiet = FALSE
) {
  y_true <- as.integer(y_true)
  
  # If we scored nothing in this fold, return NA metrics safely
  if (length(y_true) == 0L || length(y_score) == 0L) {
    out <- list(
      accuracy  = NA_real_,
      precision = NA_real_,
      recall    = NA_real_,
      f1        = NA_real_,
      auc       = NA_real_
    )
    return(invisible(out))
  }
  
  y_pred <- as.integer(y_score >= threshold)
  
  tp <- sum(y_pred == 1 & y_true == 1)
  tn <- sum(y_pred == 0 & y_true == 0)
  fp <- sum(y_pred == 1 & y_true == 0)
  fn <- sum(y_pred == 0 & y_true == 1)
  
  acc  <- mean(y_pred == y_true)
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
    2 * prec * rec / (prec + rec)
  } else {
    NA_real_
  }
  auc  <- compute_auc(y_true, y_score)
  
  if (!quiet) {
    if (!is.null(header)) cat(header, "\n", sep = "")
    cat(sprintf("  Accuracy : %.6f\n", acc))
    cat(sprintf("  Precision: %.6f\n", prec))
    cat(sprintf("  Recall   : %.6f\n", rec))
    cat(sprintf("  F1 Score : %.6f\n", f1))
    cat(sprintf("  AUC      : %.6f\n", auc))
  }
  
  invisible(list(
    accuracy  = acc,
    precision = prec,
    recall    = rec,
    f1        = f1,
    auc       = auc
  ))
}

# ------------------------------------------------------------
# Helper: collapse rare factor levels in training data
# - Replaces low-frequency levels with "Other".
# - Applied to TRAIN only to avoid data leakage.
# ------------------------------------------------------------
collapse_rare_levels_train <- function(
    train_df,
    cols,
    threshold = 10L,
    other_label = "Other"
) {
  if (length(cols) == 0L) return(train_df)
  
  for (col in cols) {
    x_chr <- as.character(train_df[[col]])
    freq <- table(x_chr, useNA = "no")
    rare_levels <- names(freq[freq < threshold])
    
    x_chr[x_chr %in% rare_levels] <- other_label
    
    x_fac <- factor(x_chr)
    
    # Ensure "Other" exists as a level even if not used
    if (!(other_label %in% levels(x_fac))) {
      levels(x_fac) <- c(levels(x_fac), other_label)
    }
    
    train_df[[col]] <- x_fac
  }
  
  train_df
}

# ------------------------------------------------------------
# Helper: align test factor levels to training factor levels
# - If train has a factor column, force test to that factor's levels.
# - Any unseen test level becomes NA, so model.frame(..., na.omit)
#   can drop those rows before prediction/scoring.
# ------------------------------------------------------------
align_test_levels_to_train <- function(
    test_df,
    train_df,
    cols,
    other_label = "Other"
) {
  if (length(cols) == 0L) return(test_df)
  
  for (col in cols) {
    if (!is.factor(train_df[[col]])) next
    
    tr_levels <- levels(train_df[[col]])
    if (!(other_label %in% tr_levels)) tr_levels <- c(tr_levels, other_label)
    
    te_chr <- as.character(test_df[[col]])
    te_chr[!(te_chr %in% tr_levels) & !is.na(te_chr)] <- other_label
    
    test_df[[col]] <- factor(te_chr, levels = tr_levels)
  }
  
  test_df
}

# ------------------------------------------------------------
# Helper: detect categorical predictors
# - Converts character/logical predictors to factors.
# - Excludes the response variable (yName).
# - Returns updated data and names of categorical predictors.
# ------------------------------------------------------------
get_categorical_predictors <- function(df, yName) {
  preds <- setdiff(names(df), yName)
  
  # Convert character/logical predictors to factor so we detect them
  for (nm in preds) {
    if (is.character(df[[nm]]) || is.logical(df[[nm]])) {
      df[[nm]] <- as.factor(df[[nm]])
    }
  }
  
  cat_cols <- preds[vapply(df[preds], is.factor, logical(1))]
  list(df = df, cat_cols = cat_cols)
}


# ------------------------------------------------------------
# Unified bootstrap(): k-fold CV for regression OR classification
# - task chooses metrics: "regression" or "classification"
# - method chooses missingness strategy: "AC" or "CC" or "PREFILL"
# ------------------------------------------------------------
bootstrap <- function(
    data,
    yName,
    k = 5,
    
    # Task type
    task = c("regression", "classification"),
    
    # Missing-data / fitting strategy
    method = c("AC", "CC", "PREFILL", "TOWER"),
    
    # Classification-specific
    threshold = 0.5,
    
    # Reproducibility
    seed = 42,
    
    # PREFILL-specific arguments
    impute_method = c("mice", "amelia", "missforest", "complete"),
    mice_method = NULL,     # mice-specific
    m = 5,
    use_dummies = FALSE,
    
    # TOWER-specific arguments
    tower_regFtnName = "lm",
    tower_opts       = list(),
    tower_scaling    = NULL,
    tower_yesYVal    = NULL,

    ...
) {
  task   <- match.arg(task)
  method <- match.arg(method)
  set.seed(seed)
  
  if (k < 2) stop("k must be at least 2.")
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!(yName %in% names(data))) stop(paste("Column", yName, "not found in data."))
  
  # -----------------------------
  # Missingness handling
  # -----------------------------
  data_used <- if (method == "CC") na.omit(data) else data
  
  n <- nrow(data_used)
  if (n < k) stop("Not enough rows for the requested number of folds.")
  
  folds <- sample(rep(1:k, length.out = n))
  
  # -----------------------------
  # Metric storage
  # -----------------------------
  if (task == "classification") {
    acc_vec  <- numeric(k)
    prec_vec <- numeric(k)
    rec_vec  <- numeric(k)
    f1_vec   <- numeric(k)
    auc_vec  <- numeric(k)
  } else {
    mse_vec  <- numeric(k)
    rmse_vec <- numeric(k)
    mae_vec  <- numeric(k)
    r2_vec   <- numeric(k)
  }
  
  # -----------------------------
  # Fold loop
  # -----------------------------
  for (i in seq_len(k)) {
    test_idx  <- which(folds == i)
    train_idx <- which(folds != i)
    
    train_dat <- data_used[train_idx, , drop = FALSE]
    test_dat  <- data_used[test_idx,  , drop = FALSE]
    
    # ------------------------------------------------------------
    # Fold-wise categorical preprocessing (NO leakage)
    # - Convert character/logical -> factor (predictors only)
    # - Identify categorical predictors
    # - Collapse rare levels in TRAIN
    # - Align TEST levels to TRAIN levels (unseen -> "Other")
    # ------------------------------------------------------------
    tmp_tr <- get_categorical_predictors(train_dat, yName)
    train_dat <- tmp_tr$df
    cat_cols  <- tmp_tr$cat_cols
    
    tmp_te <- get_categorical_predictors(test_dat, yName)
    test_dat <- tmp_te$df
    
    train_dat <- collapse_rare_levels_train(
      train_df = train_dat,
      cols = cat_cols,
      threshold = 10L,
      other_label = "Other"
    )
    test_dat <- align_test_levels_to_train(
      test_df = test_dat,
      train_df = train_dat,
      cols = cat_cols,
      other_label = "Other"
    )
    
    # -------------------------
    # CC path
    # -------------------------
    if (method == "CC") {
      form <- reformulate(setdiff(names(train_dat), yName), response = yName)
      
      # Build model.frame on TEST and drop rows with any NA 
      mf_te <- model.frame(form, test_dat, na.action = na.omit)
      
      if (nrow(mf_te) == 0L) {
        if (task == "classification") {
          acc_vec[i]  <- NA_real_
          prec_vec[i] <- NA_real_
          rec_vec[i]  <- NA_real_
          f1_vec[i]   <- NA_real_
          auc_vec[i]  <- NA_real_
        } else {
          mse_vec[i]  <- NA_real_
          rmse_vec[i] <- NA_real_
          mae_vec[i]  <- NA_real_
          r2_vec[i]   <- NA_real_
        }
        next
      }
      
      if (task == "classification") {
        fit <- glm(form, data = train_dat, family = binomial())
        
        y_true  <- model.response(mf_te)
        y_score <- as.numeric(predict(fit, newdata = mf_te, type = "response"))
      } else {
        fit <- lm(form, data = train_dat)
        
        y_true <- model.response(mf_te)
        y_pred <- as.numeric(predict(fit, newdata = mf_te))
      }
      
    } else if (method == "AC") {
      # -------------------------
      # AC path
      # -------------------------
      obj <- lm_ac(data_used, yName = yName, holdout = test_idx, ...)
      
      te_raw <- obj$testing_data
      mf_te  <- model.frame(obj$terms, te_raw, na.action = na.omit)
      
      if (nrow(mf_te) == 0L) {
        if (task == "classification") {
          acc_vec[i]  <- NA_real_
          prec_vec[i] <- NA_real_
          rec_vec[i]  <- NA_real_
          f1_vec[i]   <- NA_real_
          auc_vec[i]  <- NA_real_
        } else {
          mse_vec[i]  <- NA_real_
          rmse_vec[i] <- NA_real_
          mae_vec[i]  <- NA_real_
          r2_vec[i]   <- NA_real_
        }
        next
      }
      
      y_true <- model.response(mf_te)
      
      if (task == "classification") {
        y_score <- as.numeric(predict(obj, mf_te))
      } else {
        y_pred <- as.numeric(predict(obj, mf_te))
      }
      
    } else if (method == "TOWER") {
      # -------------------------
      # TOWER path (lm_tower)
      # -------------------------
      
      # y may have NAs in test; drop those rows for scoring only
      y_all <- test_dat[[yName]]
      ok_y  <- !is.na(y_all)
      
      if (!any(ok_y)) {
        if (task == "classification") {
          acc_vec[i]  <- NA_real_
          prec_vec[i] <- NA_real_
          rec_vec[i]  <- NA_real_
          f1_vec[i]   <- NA_real_
          auc_vec[i]  <- NA_real_
        } else {
          mse_vec[i]  <- NA_real_
          rmse_vec[i] <- NA_real_
          mae_vec[i]  <- NA_real_
          r2_vec[i]   <- NA_real_
        }
        next
      }
      
      # Fit tower on TRAIN
      fit_tw <- lm_tower(
        train_dat,
        yName      = yName,
        regFtnName = tower_regFtnName,
        opts       = tower_opts,
        scaling    = tower_scaling,
        yesYVal    = tower_yesYVal
      )
      
      # Prepare test predictors 
      x_te <- test_dat[ok_y, , drop = FALSE]
      x_te[[yName]] <- NULL
      
      # Predict
      if (task == "classification") {
        y_true  <- y_all[ok_y]
        y_score <- as.numeric(predict(fit_tw, newdata = x_te))
      } else {
        y_true <- y_all[ok_y]
        y_pred <- as.numeric(predict(fit_tw, newdata = x_te))
      } 
    } else {
      # -------------------------
      # PREFILL path (lm_prefill)
      # -------------------------
      impute_method <- tolower(impute_method)
      impute_method <- match.arg(impute_method)
      
      m_int <- suppressWarnings(as.integer(m)[1L])
      if (is.na(m_int) || m_int < 1L) stop("m must be a positive integer scalar (e.g., m = 5).")
      
      # Capture user-specified args (e.g., method="pmm") and forward them
      dots <- list(...)
      
      # Forward mice_method as mice::mice(method=...)
      if (impute_method == "mice" && !is.null(mice_method)) {
        dots$method <- mice_method
      }
      
      # Split fold
      train_fold <- data_used[train_idx, , drop = FALSE]
      test_fold  <- data_used[test_idx,  , drop = FALSE]
      
      # Convert predictors to factor where needed
      preds_all <- setdiff(names(train_fold), yName)
      for (nm in preds_all) {
        if (is.character(train_fold[[nm]]) || is.logical(train_fold[[nm]])) {
          train_fold[[nm]] <- as.factor(train_fold[[nm]])
        }
        if (is.character(test_fold[[nm]]) || is.logical(test_fold[[nm]])) {
          test_fold[[nm]] <- as.factor(test_fold[[nm]])
        }
      }
      
      # Identify categorical predictors based on TRAIN
      tmp_tr <- get_categorical_predictors(train_fold, yName)
      train_fold <- tmp_tr$df
      cat_cols   <- tmp_tr$cat_cols
      
      # Collapse rare levels in TRAIN, align TEST to TRAIN
      train_fold <- collapse_rare_levels_train(
        train_df = train_fold,
        cols = cat_cols,
        threshold = 10L,
        other_label = "Other"
      )
      test_fold <- align_test_levels_to_train(
        test_df = test_fold,
        train_df = train_fold,
        cols = cat_cols,
        other_label = "Other"
      )
      
      # Rebuild data_fold with factor levels shrunk to TRAIN levels
      data_fold <- data_used
      
      # Put back non-categorical columns normally
      data_fold[train_idx, ] <- train_fold
      data_fold[test_idx,  ] <- test_fold
      
      # Now forcibly rebuild each categorical column with TRAIN's levels only
      for (col in cat_cols) {
        lv <- levels(train_fold[[col]])  # includes "Other" from helper
        
        col_chr <- rep(NA_character_, nrow(data_fold))
        col_chr[train_idx] <- as.character(train_fold[[col]])
        col_chr[test_idx]  <- as.character(test_fold[[col]])
        
        data_fold[[col]] <- factor(col_chr, levels = lv)
      }
      
      if (impute_method == "mice" && is.null(dots$nnet.MaxNWts)) {
        dots$nnet.MaxNWts <- 50000
      }
      
      args_pf <- c(list(
        data_fold,
        yName = yName,
        impute_method = impute_method,
        m = m_int,
        use_dummies = use_dummies,
        holdout = test_idx
      ), dots)
      
      obj_pf <- do.call(lm_prefill, args_pf)
      
      te_raw <- obj_pf$testing_data
      mf_te  <- model.frame(obj_pf$formula, te_raw, na.action = na.omit)
      
      if (nrow(mf_te) == 0L) {
        if (task == "classification") {
          acc_vec[i]  <- NA_real_
          prec_vec[i] <- NA_real_
          rec_vec[i]  <- NA_real_
          f1_vec[i]   <- NA_real_
          auc_vec[i]  <- NA_real_
        } else {
          mse_vec[i]  <- NA_real_
          rmse_vec[i] <- NA_real_
          mae_vec[i]  <- NA_real_
          r2_vec[i]   <- NA_real_
        }
        next
      }
      
      y_true <- model.response(mf_te)
      
      if (task == "classification") {
        y_score <- as.numeric(predict(obj_pf, newdata = mf_te))
      } else {
        y_pred  <- as.numeric(predict(obj_pf, newdata = mf_te))
      }
    }
    
    
    # -------------------------
    # Metrics
    # -------------------------
    if (task == "classification") {
      if (is.factor(y_true)) {
        if (nlevels(y_true) != 2L) stop("Classification requires a binary response.")
        y_bin <- as.integer(y_true == levels(y_true)[2L])
      } else {
        y_bin <- as.integer(y_true)
      }
      
      met <- classification_report(
        y_true    = y_bin,
        y_score   = y_score,
        threshold = threshold,
        header    = NULL,
        quiet     = TRUE
      )
      
      acc_vec[i]  <- met$accuracy
      prec_vec[i] <- met$precision
      rec_vec[i]  <- met$recall
      f1_vec[i]   <- met$f1
      auc_vec[i]  <- met$auc
      
    } else {
      met <- metrics_regression(y_true, y_pred)
      mse_vec[i]  <- met["MSE"]
      rmse_vec[i] <- met["RMSE"]
      mae_vec[i]  <- met["MAE"]
      r2_vec[i]   <- met["R2"]
    }
  }
  
  # -----------------------------
  # Output
  # -----------------------------
  if (task == "classification") {
    list(
      task = task,
      method = method,
      k = k,
      threshold = threshold,
      
      accuracy_mean  = mean(acc_vec,  na.rm = TRUE),
      precision_mean = mean(prec_vec, na.rm = TRUE),
      recall_mean    = mean(rec_vec,  na.rm = TRUE),
      f1_mean        = mean(f1_vec,   na.rm = TRUE),
      auc_mean       = mean(auc_vec,  na.rm = TRUE),
      
      accuracy  = acc_vec,
      precision = prec_vec,
      recall    = rec_vec,
      f1        = f1_vec,
      auc       = auc_vec
    )
  } else {
    list(
      task = task,
      method = method,
      k = k,
      
      MSE_mean  = mean(mse_vec,  na.rm = TRUE),
      RMSE_mean = mean(rmse_vec, na.rm = TRUE),
      MAE_mean  = mean(mae_vec,  na.rm = TRUE),
      R2_mean   = mean(r2_vec,   na.rm = TRUE),
      
      MSE  = mse_vec,
      RMSE = rmse_vec,
      MAE  = mae_vec,
      R2   = r2_vec
    )
  }
}
