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
    test_split = 0,
    
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
  dots <- list(...)
  set.seed(seed)
  
  if (k < 2) stop("k must be at least 2.")
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!(yName %in% names(data))) stop(paste("Column", yName, "not found in data."))
  if (!is.numeric(test_split) || length(test_split) != 1L || is.na(test_split) ||
      test_split < 0 || test_split >= 1) {
    stop("test_split must be a single number in [0, 1).")
  }
  
  # -----------------------------
  # Missingness handling
  # -----------------------------
  data_all <- data
  n_total <- nrow(data_all)
  if (n_total < k) stop("Not enough rows for the requested number of folds.")
  
  if (test_split > 0) {
    n_test <- floor(n_total * test_split)
    if (n_test < 1L) n_test <- 1L
    if (n_test >= n_total) stop("test_split leaves no rows for training.")
    outer_test_idx <- sort(sample.int(n_total, size = n_test))
    outer_train_idx <- setdiff(seq_len(n_total), outer_test_idx)
  } else {
    outer_train_idx <- seq_len(n_total)
    outer_test_idx <- integer(0)
  }
  
  data_used <- data_all[outer_train_idx, , drop = FALSE]
  
  n <- nrow(data_used)
  if (n < k) stop("Not enough rows for the requested number of folds.")
  
  folds <- sample(rep(1:k, length.out = n))
  
  fit_and_predict_split <- function(train_dat, test_dat, data_reference, holdout_idx) {
    empty_res <- list(model = NULL, y_true = NULL, y_pred = NULL, y_score = NULL)
    
    if (method == "CC") {
      tmp_tr <- get_categorical_predictors(train_dat, yName)
      train_proc <- tmp_tr$df
      cat_cols <- tmp_tr$cat_cols
      tmp_te <- get_categorical_predictors(test_dat, yName)
      test_proc <- tmp_te$df
      train_proc <- collapse_rare_levels_train(
        train_df = train_proc,
        cols = cat_cols,
        threshold = 10L,
        other_label = "Other"
      )
      test_proc <- align_test_levels_to_train(
        test_df = test_proc,
        train_df = train_proc,
        cols = cat_cols,
        other_label = "Other"
      )
      
      form <- reformulate(setdiff(names(train_proc), yName), response = yName)
      train_cc <- stats::na.omit(train_proc)
      if (nrow(train_cc) == 0L) return(empty_res)
      
      model <- if (task == "classification") {
        glm(form, data = train_cc, family = binomial())
      } else {
        lm(form, data = train_cc)
      }
      
      mf_te <- model.frame(form, test_proc, na.action = na.omit)
      if (nrow(mf_te) == 0L) {
        return(list(model = model, y_true = NULL, y_pred = NULL, y_score = NULL))
      }
      
      if (task == "classification") {
        return(list(
          model = model,
          y_true = model.response(mf_te),
          y_pred = NULL,
          y_score = as.numeric(predict(model, newdata = mf_te, type = "response"))
        ))
      }
      
      return(list(
        model = model,
        y_true = model.response(mf_te),
        y_pred = as.numeric(predict(model, newdata = mf_te)),
        y_score = NULL
      ))
    }
    
    if (method == "AC") {
      model <- do.call(lm_ac, c(list(data = data_reference, yName = yName, holdout = holdout_idx), dots))
      te_raw <- model$testing_data
      mf_te <- model.frame(model$terms, te_raw, na.action = na.omit)
      if (nrow(mf_te) == 0L) {
        return(list(model = model, y_true = NULL, y_pred = NULL, y_score = NULL))
      }
      
      if (task == "classification") {
        return(list(
          model = model,
          y_true = model.response(mf_te),
          y_pred = NULL,
          y_score = as.numeric(predict(model, mf_te))
        ))
      }
      
      return(list(
        model = model,
        y_true = model.response(mf_te),
        y_pred = as.numeric(predict(model, mf_te)),
        y_score = NULL
      ))
    }
    
    if (method == "TOWER") {
      tmp_tr <- get_categorical_predictors(train_dat, yName)
      train_proc <- tmp_tr$df
      cat_cols <- tmp_tr$cat_cols
      tmp_te <- get_categorical_predictors(test_dat, yName)
      test_proc <- tmp_te$df
      train_proc <- collapse_rare_levels_train(
        train_df = train_proc,
        cols = cat_cols,
        threshold = 10L,
        other_label = "Other"
      )
      test_proc <- align_test_levels_to_train(
        test_df = test_proc,
        train_df = train_proc,
        cols = cat_cols,
        other_label = "Other"
      )
      
      form <- reformulate(setdiff(names(train_proc), yName), response = yName)
      model <- lm_tower(
        train_proc,
        yName      = yName,
        regFtnName = tower_regFtnName,
        opts       = tower_opts,
        scaling    = tower_scaling,
        yesYVal    = tower_yesYVal
      )
      
      mf_te <- model.frame(form, test_proc, na.action = na.omit)
      if (nrow(mf_te) == 0L) {
        return(list(model = model, y_true = NULL, y_pred = NULL, y_score = NULL))
      }
      
      x_te <- mf_te
      x_te[[yName]] <- NULL
      
      if (task == "classification") {
        return(list(
          model = model,
          y_true = model.response(mf_te),
          y_pred = NULL,
          y_score = as.numeric(predict(model, newdata = x_te))
        ))
      }
      
      return(list(
        model = model,
        y_true = model.response(mf_te),
        y_pred = as.numeric(predict(model, newdata = x_te)),
        y_score = NULL
      ))
    }
    
    impute_method_local <- tolower(impute_method)
    impute_method_local <- match.arg(
      impute_method_local,
      choices = c("mice", "amelia", "missforest", "complete")
    )
    m_int <- suppressWarnings(as.integer(m)[1L])
    if (is.na(m_int) || m_int < 1L) stop("m must be a positive integer scalar (e.g., m = 5).")
    
    train_fold <- train_dat
    test_fold <- test_dat
    preds_all <- setdiff(names(train_fold), yName)
    for (nm in preds_all) {
      if (is.character(train_fold[[nm]]) || is.logical(train_fold[[nm]])) {
        train_fold[[nm]] <- as.factor(train_fold[[nm]])
      }
      if (is.character(test_fold[[nm]]) || is.logical(test_fold[[nm]])) {
        test_fold[[nm]] <- as.factor(test_fold[[nm]])
      }
    }
    
    tmp_tr <- get_categorical_predictors(train_fold, yName)
    train_fold <- tmp_tr$df
    cat_cols <- tmp_tr$cat_cols
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
    
    data_fold <- data_reference
    data_fold[setdiff(seq_len(nrow(data_reference)), holdout_idx), ] <- train_fold
    data_fold[holdout_idx, ] <- test_fold
    
    for (col in cat_cols) {
      lv <- levels(train_fold[[col]])
      col_chr <- rep(NA_character_, nrow(data_fold))
      col_chr[setdiff(seq_len(nrow(data_reference)), holdout_idx)] <- as.character(train_fold[[col]])
      col_chr[holdout_idx] <- as.character(test_fold[[col]])
      data_fold[[col]] <- factor(col_chr, levels = lv)
    }
    
    prefill_dots <- dots
    if (impute_method_local == "mice" && !is.null(mice_method)) {
      prefill_dots$method <- mice_method
    }
    if (impute_method_local == "mice" && is.null(prefill_dots$nnet.MaxNWts)) {
      prefill_dots$nnet.MaxNWts <- 50000
    }
    
    model <- do.call(lm_prefill, c(list(
      data_fold,
      yName = yName,
      impute_method = impute_method_local,
      m = m_int,
      use_dummies = use_dummies,
      holdout = holdout_idx
    ), prefill_dots))
    
    te_raw <- model$testing_data
    mf_te <- model.frame(model$formula, te_raw, na.action = na.omit)
    if (nrow(mf_te) == 0L) {
      return(list(model = model, y_true = NULL, y_pred = NULL, y_score = NULL))
    }
    
    if (task == "classification") {
      return(list(
        model = model,
        y_true = model.response(mf_te),
        y_pred = NULL,
        y_score = as.numeric(predict(model, newdata = mf_te))
      ))
    }
    
    list(
      model = model,
      y_true = model.response(mf_te),
      y_pred = as.numeric(predict(model, newdata = mf_te)),
      y_score = NULL
    )
  }

  # -----------------------------
  # Missingness storage for training folds
  # -----------------------------
  train_missing_tables <- vector("list", k)
  test_missing_tables <- vector("list", k)
  fold_models <- vector("list", k)
  
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
    # Training-fold missingness table
    # - Includes predictors and target
    # - Stores percent missingness for this fold's training data
    # -------------------------
    train_missing_tables[[i]] <- data.frame(
      fold = i,
      variable = names(train_dat),
      role = ifelse(names(train_dat) == yName, "target", "predictor"),
      missing_pct = colMeans(is.na(train_dat)) * 100,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    test_missing_tables[[i]] <- data.frame(
      fold = i,
      variable = names(test_dat),
      role = ifelse(names(test_dat) == yName, "target", "predictor"),
      missing_pct = colMeans(is.na(test_dat)) * 100,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("\nTraining-data missingness for fold %d:\n", i))
    print(train_missing_tables[[i]], row.names = FALSE)
    cat(sprintf("\nTesting-data missingness for fold %d:\n", i))
    print(test_missing_tables[[i]], row.names = FALSE)
    
    split_res <- fit_and_predict_split(
      train_dat = train_dat,
      test_dat = test_dat,
      data_reference = data_used,
      holdout_idx = test_idx
    )
    fold_models[[i]] <- split_res$model
    
    
    # -------------------------
    # Metrics
    # -------------------------
    if (task == "classification") {
      if (is.null(split_res$y_true) || length(split_res$y_true) == 0L) {
        acc_vec[i]  <- NA_real_
        prec_vec[i] <- NA_real_
        rec_vec[i]  <- NA_real_
        f1_vec[i]   <- NA_real_
        auc_vec[i]  <- NA_real_
        next
      }
      
      y_true <- split_res$y_true
      y_score <- split_res$y_score
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
      if (is.null(split_res$y_true) || length(split_res$y_true) == 0L) {
        mse_vec[i]  <- NA_real_
        rmse_vec[i] <- NA_real_
        mae_vec[i]  <- NA_real_
        r2_vec[i]   <- NA_real_
        next
      }
      
      y_true <- split_res$y_true
      y_pred <- split_res$y_pred
      met <- metrics_regression(y_true, y_pred)
      mse_vec[i]  <- met["MSE"]
      rmse_vec[i] <- met["RMSE"]
      mae_vec[i]  <- met["MAE"]
      r2_vec[i]   <- met["R2"]
    }
  }
  
  # -----------------------------
  # Average missingness across training folds
  # -----------------------------
  train_missing_all <- do.call(rbind, train_missing_tables)
  
  train_missing_average <- aggregate(
    missing_pct ~ variable + role,
    data = train_missing_all,
    FUN = mean
  )
  
  train_missing_average <- train_missing_average[
    order(train_missing_average$role, train_missing_average$variable),
  ]
  
  cat("\nAverage training-data missingness across folds:\n")
  print(train_missing_average, row.names = FALSE)
  
  test_missing_all <- do.call(rbind, test_missing_tables)
  test_missing_average <- aggregate(
    missing_pct ~ variable + role,
    data = test_missing_all,
    FUN = mean
  )
  test_missing_average <- test_missing_average[
    order(test_missing_average$role, test_missing_average$variable),
  ]
  
  cat("\nAverage testing-data missingness across folds:\n")
  print(test_missing_average, row.names = FALSE)
  
  final_model <- NULL
  test_results <- NULL
  if (test_split > 0) {
    final_res <- fit_and_predict_split(
      train_dat = data_used,
      test_dat = data_all[outer_test_idx, , drop = FALSE],
      data_reference = data_all,
      holdout_idx = outer_test_idx
    )
    final_model <- final_res$model
    n_test_scored <- if (is.null(final_res$y_true)) 0L else length(final_res$y_true)
    
    if (task == "classification") {
      if (is.null(final_res$y_true) || n_test_scored == 0L) {
        test_results <- list(
          accuracy = NA_real_,
          precision = NA_real_,
          recall = NA_real_,
          f1 = NA_real_,
          auc = NA_real_,
          n_test_scored = 0L,
          n_test_dropped = length(outer_test_idx)
        )
      } else {
        y_true_test <- final_res$y_true
        if (is.factor(y_true_test)) {
          if (nlevels(y_true_test) != 2L) stop("Classification requires a binary response.")
          y_bin_test <- as.integer(y_true_test == levels(y_true_test)[2L])
        } else {
          y_bin_test <- as.integer(y_true_test)
        }
        
        met_test <- classification_report(
          y_true = y_bin_test,
          y_score = final_res$y_score,
          threshold = threshold,
          header = NULL,
          quiet = TRUE
        )
        test_results <- c(met_test, list(
          n_test_scored = n_test_scored,
          n_test_dropped = length(outer_test_idx) - n_test_scored
        ))
      }
    } else {
      if (is.null(final_res$y_true) || n_test_scored == 0L) {
        test_results <- list(
          MSE = NA_real_,
          RMSE = NA_real_,
          MAE = NA_real_,
          R2 = NA_real_,
          n_test_scored = 0L,
          n_test_dropped = length(outer_test_idx)
        )
      } else {
        met_test <- metrics_regression(final_res$y_true, final_res$y_pred)
        test_results <- list(
          MSE = as.numeric(met_test["MSE"]),
          RMSE = as.numeric(met_test["RMSE"]),
          MAE = as.numeric(met_test["MAE"]),
          R2 = as.numeric(met_test["R2"]),
          n_test_scored = n_test_scored,
          n_test_dropped = length(outer_test_idx) - n_test_scored
        )
      }
    }
  }
  
  split_info <- list(
    n_total = n_total,
    n_train = length(outer_train_idx),
    n_test = length(outer_test_idx),
    train_indices = outer_train_idx,
    test_indices = outer_test_idx
  )

  # -----------------------------
  # Output
  # -----------------------------
  if (task == "classification") {
    list(
      task = task,
      method = method,
      k = k,
      threshold = threshold,
      test_split = test_split,
      split_info = split_info,
      
      accuracy_mean  = mean(acc_vec,  na.rm = TRUE),
      precision_mean = mean(prec_vec, na.rm = TRUE),
      recall_mean    = mean(rec_vec,  na.rm = TRUE),
      f1_mean        = mean(f1_vec,   na.rm = TRUE),
      auc_mean       = mean(auc_vec,  na.rm = TRUE),
      
      accuracy  = acc_vec,
      precision = prec_vec,
      recall    = rec_vec,
      f1        = f1_vec,
      auc       = auc_vec,
      fold_models = fold_models,
      final_model = final_model,
      test_results = test_results,
      training_missingness_by_fold = train_missing_tables,
      training_missingness_average = train_missing_average
    )
  } else {
    list(
      task = task,
      method = method,
      k = k,
      test_split = test_split,
      split_info = split_info,
      
      MSE_mean  = mean(mse_vec,  na.rm = TRUE),
      RMSE_mean = mean(rmse_vec, na.rm = TRUE),
      MAE_mean  = mean(mae_vec,  na.rm = TRUE),
      R2_mean   = mean(r2_vec,   na.rm = TRUE),
      
      MSE  = mse_vec,
      RMSE = rmse_vec,
      MAE  = mae_vec,
      R2   = r2_vec,
      fold_models = fold_models,
      final_model = final_model,
      test_results = test_results,
      training_missingness_by_fold = train_missing_tables,
      training_missingness_average = train_missing_average
    )
  }
}
