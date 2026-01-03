# ============================================================
# lm_ac: Linear Model using Available Cases (pairwise deletion)
# ------------------------------------------------------------
# Implements regression when X or y contain NAs by computing
# each element of X'X and X'y as the *average* of intact pairs,
# then solving the normal equations directly on those averages
# ============================================================

# -----------------------------------------------------------------------
# Constructor: Initializes an lm_ac object with formula and raw data.
# -----------------------------------------------------------------------
lm_ac <- function(data, yName, holdout = NULL, ...) {
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame.")
  if (!(yName %in% names(data))) stop(paste("Column", yName, "not found in data."))
  if (!is.numeric(data[[yName]])) stop("Response variable must be numeric.")
  if (nrow(data) == 0) stop("Input data is empty.")
  
  # Holdout set
  if (!is.null(holdout)) {
    if (is.logical(holdout) && length(holdout) == nrow(data)) {
      training_data <- data[!holdout, , drop = FALSE]
      testing_data  <- data[ holdout, , drop = FALSE]
    } else if (is.numeric(holdout) && all(holdout %in% seq_len(nrow(data)))) {
      training_data <- data[-holdout, , drop = FALSE]
      testing_data  <- data[ holdout, , drop = FALSE]
    } else {
      stop("Holdout is invalid")
    }
  } else {
    training_data <- data
    testing_data  <- NULL
  }
  
  if (nrow(training_data) == 0) stop("No training rows available.")
  
  # Convert character/logical predictors to factors
  for (col in names(training_data)) {
    if (is.character(training_data[[col]]) || is.logical(training_data[[col]])) {
      training_data[[col]] <- as.factor(training_data[[col]])
    }
  }
  
  # Regression formula
  formula <- reformulate(setdiff(names(training_data), yName), response = yName)
  
  # model.frame/matrix with na.pass
  mf <- model.frame(formula, training_data, na.action = na.pass, ...)
  y  <- model.response(mf)
  X  <- model.matrix(attr(mf, "terms"), mf, na.action = na.pass)
  
  trm  <- terms(mf)
  coln <- colnames(X)
  
  XtX_avg <- ac_mul(X)
  Xty_avg <- ac_vec(X, y)
  
  if (anyNA(XtX_avg) || anyNA(Xty_avg)) {
    stop("Missing values in system matrices: need at least one intact pair per term (and intercept).")
  }
  
  beta_hat <- tryCatch(
    solve(XtX_avg, Xty_avg),
    error = function(e) stop("XtX is singular; regression cannot proceed.")
  )
  
  obj <- list(
    data         = training_data,
    testing_data = testing_data,
    yName        = yName,
    formula      = formula,
    terms        = trm,
    fit_obj      = list(
      coef     = as.vector(beta_hat),
      colnames = coln
    )
  )
  class(obj) <- "lm_ac"
  obj
}

ac_mul <- function(A) {
  if (!is.matrix(A)) A <- as.matrix(A)
  p <- ncol(A)
  result <- matrix(NA_real_, p, p)
  for (i in 1:p) {
    xi <- A[, i]
    for (j in i:p) {
      xj <- A[, j]
      valid_idx <- !is.na(xi) & !is.na(xj)
      if (any(valid_idx)) {
        val <- mean(xi[valid_idx] * xj[valid_idx])
        result[i, j] <- val
        result[j, i] <- val
      }
    }
  }
  colnames(result) <- colnames(A)
  rownames(result) <- colnames(A)
  result
}

ac_vec <- function(X, y) {
  if (!is.matrix(X)) X <- as.matrix(X)
  p <- ncol(X)
  result <- numeric(p)
  for (i in 1:p) {
    xi <- X[, i]
    valid_idx <- !is.na(xi) & !is.na(y)
    result[i] <- if (any(valid_idx)) mean(xi[valid_idx] * y[valid_idx]) else NA_real_
  }
  names(result) <- colnames(X)
  result
}

summary.lm_ac <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  coefs <- object$fit_obj$coef
  names(coefs) <- object$fit_obj$colnames
  print(coefs)
  invisible(coefs)
}

predict.lm_ac <- function(object, newdata, ...) {
  stopifnot(inherits(object, "lm_ac"))
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  if (missing(newdata)) stop("'newdata' is required for predict.lm_ac().")
  
  for (col in names(newdata)) {
    if (is.character(newdata[[col]]) || is.logical(newdata[[col]])) {
      newdata[[col]] <- as.factor(newdata[[col]])
    }
  }
  
  mf_new <- model.frame(object$terms, newdata, na.action = na.pass, ...)
  X_new  <- model.matrix(object$terms, mf_new, na.action = na.pass)
  
  as.vector(X_new %*% object$fit_obj$coef)
}

# ============================================================
# bootstrap = k-fold CV for CLASSIFICATION (CC vs AC)
# ============================================================

bootstrap <- function(
    data,
    yName,
    k = 5,
    threshold = 0.5,
    method = c("AC", "CC"),
    seed = 42,
    ...
) {
  method <- match.arg(method)
  set.seed(seed)
  
  if (k < 2) stop("k must be at least 2.")
  
  if (method == "CC") {
    data_used <- na.omit(data)
  } else {
    data_used <- data
  }
  
  n <- nrow(data_used)
  if (n < k) stop("Not enough rows for the requested number of folds.")
  
  folds <- sample(rep(1:k, length.out = n))
  
  acc_vec  <- numeric(k)
  prec_vec <- numeric(k)
  rec_vec  <- numeric(k)
  f1_vec   <- numeric(k)
  auc_vec  <- numeric(k)
  
  for (i in seq_len(k)) {
    test_idx  <- which(folds == i)
    train_idx <- which(folds != i)
    
    train_dat <- data_used[train_idx, , drop = FALSE]
    test_dat  <- data_used[test_idx,  , drop = FALSE]
    
    if (method == "CC") {
      form <- reformulate(setdiff(names(train_dat), yName), response = yName)
      fit  <- lm(form, data = train_dat)
      
      y_true  <- test_dat[[yName]]
      y_score <- as.numeric(predict(fit, newdata = test_dat))
    } else {
      holdout_set <- test_idx
      obj <- lm_ac(data_used, yName = yName, holdout = holdout_set, ...)
      
      te_raw <- obj$testing_data
      mf_te  <- model.frame(obj$terms, te_raw, na.action = na.omit)
      if (nrow(mf_te) == 0L) {
        acc_vec[i]  <- NA_real_
        prec_vec[i] <- NA_real_
        rec_vec[i]  <- NA_real_
        f1_vec[i]   <- NA_real_
        auc_vec[i]  <- NA_real_
        next
      }
      
      y_true <- model.response(mf_te)
      X_te   <- model.matrix(obj$terms, mf_te, na.action = na.pass)
      X_te   <- X_te[, obj$fit_obj$colnames, drop = FALSE]
      y_score <- as.numeric(X_te %*% obj$fit_obj$coef)
    }
    
    if (is.factor(y_true)) {
      if (nlevels(y_true) != 2L) stop("For k-fold CV classification, y must be binary.")
      y_bin <- as.integer(y_true == levels(y_true)[2L])
    } else {
      y_bin <- as.integer(y_true)
    }
    
    m <- classification_report(
      y_true = y_bin,
      y_score = y_score,
      threshold = threshold,
      header = NULL,
      quiet = TRUE
    )
    
    acc_vec[i]  <- m$accuracy
    prec_vec[i] <- m$precision
    rec_vec[i]  <- m$recall
    f1_vec[i]   <- m$f1
    auc_vec[i]  <- m$auc
  }
  
  list(
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
}

kfold_regression <- function(
    data,
    yName,
    k      = 5,
    method = c("AC", "CC"),
    seed   = 42,
    ...
) {
  method <- match.arg(method)
  set.seed(seed)
  
  if (k < 2) stop("k must be at least 2.")
  
  if (method == "CC") {
    data_used <- na.omit(data)
  } else {
    data_used <- data
  }
  
  n <- nrow(data_used)
  if (n < k) stop("Not enough rows for the requested number of folds.")
  
  folds <- sample(rep(1:k, length.out = n))
  
  mse_vec  <- numeric(k)
  rmse_vec <- numeric(k)
  mae_vec  <- numeric(k)
  r2_vec   <- numeric(k)
  
  for (i in seq_len(k)) {
    test_idx  <- which(folds == i)
    train_idx <- which(folds != i)
    
    train_dat <- data_used[train_idx, , drop = FALSE]
    test_dat  <- data_used[test_idx,  , drop = FALSE]
    
    if (method == "CC") {
      form <- reformulate(setdiff(names(train_dat), yName), response = yName)
      fit  <- lm(form, data = train_dat)
      
      y_true <- test_dat[[yName]]
      y_pred <- as.numeric(predict(fit, newdata = test_dat))
    } else {
      holdout_set <- test_idx
      obj <- lm_ac(data_used, yName = yName, holdout = holdout_set, ...)
      
      te_raw <- obj$testing_data
      mf_te  <- model.frame(obj$terms, te_raw, na.action = na.omit)
      if (nrow(mf_te) == 0L) {
        mse_vec[i]  <- NA_real_
        rmse_vec[i] <- NA_real_
        mae_vec[i]  <- NA_real_
        r2_vec[i]   <- NA_real_
        next
      }
      
      y_true <- model.response(mf_te)
      te_clean <- mf_te
      y_pred <- predict(obj, te_clean)
    }
    
    m <- metrics_regression(y_true, y_pred)
    mse_vec[i]  <- m["MSE"]
    rmse_vec[i] <- m["RMSE"]
    mae_vec[i]  <- m["MAE"]
    r2_vec[i]   <- m["R2"]
  }
  
  list(
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