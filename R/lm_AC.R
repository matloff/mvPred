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
  # ----------------------------
  # Input validation
  # ----------------------------
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame.")
  if (!(yName %in% names(data))) stop(paste("Column", yName, "not found in data."))
  if (!is.numeric(data[[yName]])) stop("Response variable must be numeric.")
  if (nrow(data) == 0) stop("Input data is empty.")
  
  # ----------------------------
  # Holdout set for bootstrap or validation
  # ----------------------------
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
  
  # ----------------------------
  # Convert character/logical predictors to factors
  # ----------------------------
  for (col in names(training_data)) {
    if (is.character(training_data[[col]]) || is.logical(training_data[[col]])) {
      training_data[[col]] <- as.factor(training_data[[col]])
    }
  }
  
  # ----------------------------
  # Construct regression formula: y ~ all other predictors
  # ----------------------------
  formula <- reformulate(setdiff(names(training_data), yName), response = yName)
  
  # ----------------------------
  # Build model.frame and model.matrix with na.pass
  #   (this ensures rows with NAs are retained for AC math)
  # ----------------------------
  mf <- model.frame(formula, training_data, na.action = na.pass, ...)
  y  <- model.response(mf)
  X  <- model.matrix(attr(mf, "terms"), mf, na.action = na.pass)
  
  # store the terms for consistent newdata handling
  trm  <- terms(mf)
  coln <- colnames(X)
  
  # ----------------------------
  # Compute average-based XtX and Xty (pairwise deletion)
  # ----------------------------
  XtX_avg <- ac_mul(X)
  Xty_avg <- ac_vec(X, y)
  
  if (anyNA(XtX_avg) || anyNA(Xty_avg)) {
    stop("
      Missing values in system matrices: 
      need at least one intact pair per term (and intercept).
    ")
  }
  
  # ----------------------------
  # Solve system and catch singular matrix errors
  # ----------------------------
  beta_hat <- tryCatch(
    solve(XtX_avg, Xty_avg),
    error = function(e) {
      # try a sequence of larger ridge values
      lam_seq <- 10 ^ seq(-6, -1, by = 1)  # 1e-6, 1e-5, ..., 1e-1
      for (lam in lam_seq) {
        M <- XtX_avg + diag(lam, ncol(XtX_avg))
        ans <- try(solve(M, Xty_avg), silent = TRUE)
        if (!inherits(ans, "try-error")) return(ans)
      }
      stop("Even ridge-regularized system is too ill-conditioned.")
    }
  )
  
  # ----------------------------
  # Return lm_ac object
  # ----------------------------
  obj <- list(
    data         = training_data,
    testing_data = testing_data,  # raw holdout (predict() will handle transform)
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

# -----------------------------------------------------------------------
# ac_mul: Compute (1/n_ij) * sum[X_i * X_j] over available pairs
# For each pair of columns (i,j), only rows where both are non-NA are used.
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# ac_vec: Compute (1/n_i) * sum[X_i * y] over available pairs
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# summary.lm_ac: Print coefficients (same shape as lm())
# -----------------------------------------------------------------------
summary.lm_ac <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  coefs <- object$fit_obj$coef
  names(coefs) <- object$fit_obj$colnames
  print(coefs)
  invisible(coefs)
}

# -----------------------------------------------------------------------
# predict.lm_ac: Given newdata (with possible NAs), rebuild design matrix
# using stored terms, and compute y_hat = X_new %*% β_hat (na.pass here).
# NOTE: In evaluation, we clean test rows to ensure no NAs (see test func).
# -----------------------------------------------------------------------
predict.lm_ac <- function(object, newdata, ...) {
  stopifnot(inherits(object, "lm_ac"))
  if (missing(newdata)) stop("'newdata' is required")
  
  newdata <- as.data.frame(newdata)
  for (nm in names(newdata)) {
    if (is.character(newdata[[nm]]) || is.logical(newdata[[nm]])) {
      newdata[[nm]] <- as.factor(newdata[[nm]])
    }
  }
  
  mf_new <- model.frame(object$terms, newdata, na.action = na.pass, ...)
  if (nrow(mf_new) == 0L) return(numeric(0))
  X_new  <- model.matrix(object$terms, mf_new, na.action = na.pass)
  
  # align to training columns
  train_cols <- object$fit_obj$colnames
  miss <- setdiff(train_cols, colnames(X_new))
  if (length(miss)) {
    X_new <- cbind(X_new, matrix(0, nrow(X_new), length(miss),
                                 dimnames = list(NULL, miss)))
  }
  X_new <- X_new[, train_cols, drop = FALSE]
  
  as.vector(X_new %*% object$fit_obj$coef)
}


# -----------------------------------------------------------------------
# k-fold cross-validation for lm_ac
#   - Splits data into k folds
#   - For each fold: train on k-1 folds, test on the held-out fold
#   - Aggregates coef, MSPE, and MAPE across folds
# -----------------------------------------------------------------------
bootstrap <- function(
    data, yName, coef_name, k = 5, ...
) {
  n <- nrow(data)
  if (n < 2) stop("Not enough rows for cross-validation.")
  if (k < 2 || k > n) stop("'k' must be between 2 and nrow(data).")
  
  # Randomly assign each row to one of k folds
  fold_ids <- sample(rep(seq_len(k), length.out = n))
  
  coefs <- numeric(k)
  mspe_values <- numeric(k)
  mape_values <- numeric(k)
  
  for (fold in seq_len(k)) {
    holdout_set <- which(fold_ids == fold)
    
    # Fit lm_ac with this fold as the holdout
    obj <- lm_ac(data, yName, holdout = holdout_set, ...)
    
    beta <- obj$fit_obj$coef
    names(beta) <- obj$fit_obj$colnames
    coefs[fold] <- if (coef_name %in% names(beta)) beta[[coef_name]] else NA_real_
    
    # Clean test (no NAs w.r.t. model terms)
    te_raw <- obj$testing_data
    mf_te  <- model.frame(obj$terms, te_raw, na.action = na.omit)
    if (nrow(mf_te) == 0L) {
      mspe_values[fold] <- NA_real_
      mape_values[fold] <- NA_real_
      next
    }
    
    y_true <- model.response(mf_te)
    X_te   <- model.matrix(obj$terms, mf_te, na.action = na.pass)
    
    ## 🔑 HERE is the snippet you asked about:
    # Make test design matrix compatible with training cols
    train_cols <- obj$fit_obj$colnames
    miss <- setdiff(train_cols, colnames(X_te))
    if (length(miss)) {
      X_te <- cbind(
        X_te,
        matrix(0, nrow(X_te), length(miss),
               dimnames = list(NULL, miss))
      )
    }
    # Now align column order
    X_te <- X_te[, train_cols, drop = FALSE]
    ## 🔑 snippet ends
    
    y_pred <- as.numeric(X_te %*% obj$fit_obj$coef)
    mspe_values[fold] <- mean((y_true - y_pred)^2)
    mape_values[fold] <- mean(abs((y_true - y_pred) / y_true))
  }
  
  avg_coef <- mean(coefs, na.rm = TRUE)
  se <- sd(coefs, na.rm = TRUE) / sqrt(sum(!is.na(coefs)))
  CI <- c(avg_coef - (1.96 * se), avg_coef + (1.96 * se))
  
  list(
    coef_mean     = avg_coef,
    coef_se       = se,
    coef_CI       = CI,
    coefs         = coefs,
    mspe          = mean(mspe_values, na.rm = TRUE),
    mape          = mean(mape_values, na.rm = TRUE),
    mspe_per_fold = mspe_values,
    mape_per_fold = mape_values,
    folds         = fold_ids
  )
}