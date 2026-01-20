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
# Helper: align test factor levels to training factor levels
# - If train has a factor column, force test to that factor's levels.
# - Any unseen test level becomes NA, so model.frame(..., na.omit)
#   can drop those rows before prediction/scoring.
# ------------------------------------------------------------
align_test_levels_to_train <- function(train_dat, test_dat, yName) {
  pred_names <- setdiff(names(train_dat), yName)
  
  for (nm in pred_names) {
    if (is.factor(train_dat[[nm]])) {
      
      # Ensure test is factor-like
      if (is.character(test_dat[[nm]]) || is.logical(test_dat[[nm]])) {
        test_dat[[nm]] <- as.factor(test_dat[[nm]])
      } else if (!is.factor(test_dat[[nm]])) {
        test_dat[[nm]] <- as.factor(test_dat[[nm]])
      }
      
      train_lvls <- levels(train_dat[[nm]])
      test_chr   <- as.character(test_dat[[nm]])
      
      # Unseen levels -> NA
      unseen <- !(test_chr %in% train_lvls)
      if (any(unseen, na.rm = TRUE)) {
        test_chr[unseen] <- NA_character_
      }
      
      test_dat[[nm]] <- factor(test_chr, levels = train_lvls)
    }
  }
  
  test_dat
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
    method = c("AC", "CC", "PREFILL"),
    
    # Classification-specific
    threshold = 0.5,
    
    # Reproducibility
    seed = 42,
    
    # PREFILL-specific arguments
    impute_method = c("mice", "amelia", "missforest", "complete"),
    m = 5,
    use_dummies = FALSE,
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
    
    # -------------------------
    # CC path
    # -------------------------
    if (method == "CC") {
      form <- reformulate(setdiff(names(train_dat), yName), response = yName)
      
      # Convert character/logical predictors to factor in train
      for (nm in setdiff(names(train_dat), yName)) {
        if (is.character(train_dat[[nm]]) || is.logical(train_dat[[nm]])) {
          train_dat[[nm]] <- as.factor(train_dat[[nm]])
        }
      }
      
      # Align TEST factor levels to TRAIN factor levels; unseen -> NA
      test_dat2 <- align_test_levels_to_train(train_dat, test_dat, yName)
      
      # Build model.frame on TEST and drop rows that became NA (unseen levels, etc.)
      mf_te <- model.frame(form, test_dat2, na.action = na.omit)
      
      # If no scorable rows remain in this fold, store NAs and continue
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
        # Logistic regression for CC classification
        fit <- glm(form, data = train_dat, family = binomial())
        
        y_true  <- model.response(mf_te)
        y_score <- as.numeric(predict(fit, newdata = mf_te, type = "response"))
      } else {
        # Linear regression for CC regression
        fit <- lm(form, data = train_dat)
        
        y_true <- model.response(mf_te)
        y_pred <- as.numeric(predict(fit, newdata = mf_te))
      }
      
      # -------------------------
      # AC path
      # -------------------------
    } else if (method == "AC") {
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
      
    } else {
      # -------------------------
      # PREFILL path (lm_prefill)
      # -------------------------
      
      # 1) Normalize impute_method (case-insensitive)
      impute_method <- tolower(impute_method)
      impute_method <- match.arg(impute_method)
      
      # 2) Force m to be a single positive integer
      m_int <- suppressWarnings(as.integer(m)[1L])
      if (is.na(m_int) || m_int < 1L) stop("m must be a positive integer scalar (e.g., m = 5).")
      
      preds_all <- setdiff(names(data_used), yName)
      for (nm in preds_all) {
        if (is.character(data_used[[nm]]) || is.logical(data_used[[nm]])) {
          data_used[[nm]] <- as.factor(data_used[[nm]])
        }
      }
      
      # Run lm_prefill
      obj_pf <- lm_prefill(
        data_used,
        yName = yName,
        impute_method = impute_method,
        m = m_int,
        use_dummies = use_dummies,
        holdout = test_idx,
        ...
      )
      
      # Training split used by lm_prefill (use for factor level alignment)
      train_pf <- obj_pf$data
      
      # Ensure train_pf has factor types 
      preds_pf <- setdiff(names(train_pf), yName)
      for (nm in preds_pf) {
        if (is.character(train_pf[[nm]]) || is.logical(train_pf[[nm]])) {
          train_pf[[nm]] <- as.factor(train_pf[[nm]])
        }
      }
      
      # Raw test split from lm_prefill
      te_raw <- obj_pf$testing_data
      
      # Coerce test predictor types to factor when applcable
      for (nm in preds_pf) {
        if (is.factor(train_pf[[nm]])) {
          if (is.character(te_raw[[nm]]) || is.logical(te_raw[[nm]]) || !is.factor(te_raw[[nm]])) {
            te_raw[[nm]] <- as.factor(te_raw[[nm]])
          }
        }
      }
      
      # Align test factor levels to training factor levels (unseen → NA)
      te_aligned <- align_test_levels_to_train(train_pf, te_raw, yName)
      
      # Build test model.frame using lm_prefill's formula; drop rows made NA by alignment
      mf_te <- model.frame(obj_pf$formula, te_aligned, na.action = na.omit)
      
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
        # lm_prefill currently fits lm(); treat output as score until glm support is added
        y_score <- as.numeric(predict(obj_pf, newdata = mf_te))
      } else {
        y_pred  <- as.numeric(predict(obj_pf, newdata = mf_te))
      }
    }
    
    # -------------------------
    # Metrics 
    # -------------------------
    if (task == "classification") {
      # Convert y_true to {0,1}
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
