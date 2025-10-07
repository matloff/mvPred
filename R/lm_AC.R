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
    error = function(e) stop("XtX is singular; regression cannot proceed.")
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
# For each pair of columns (i,j), only rows where both are non-NA
# are used, and we take the *mean* of the products. That yields
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
# For each predictor column i, only rows with non-NA X_i and y
# are used; we take the mean of the products to estimate E[X_i y].
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
# predict.lm_ac: Given newdata (with possible NAs in predictors),
# rebuild design matrix using stored terms, and compute y_hat = X_new %*% β_hat.
# We use na.pass to retain rows with missing predictors.
# -----------------------------------------------------------------------
predict.lm_ac <- function(object, newdata, ...) {
  stopifnot(inherits(object, "lm_ac"))
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  if (missing(newdata)) stop("'newdata' is required for predict.lm_ac().")
  
  # Convert character/logical predictors to factors
  for (col in names(newdata)) {
    if (is.character(newdata[[col]]) || is.logical(newdata[[col]])) {
      newdata[[col]] <- as.factor(newdata[[col]])
    }
  }
  
  # Build model.frame and model.matrix with na.pass
  mf_new <- model.frame(object$terms, newdata, na.action = na.pass, ...)
  X_new  <- model.matrix(object$terms, mf_new, na.action = na.pass)
  
  as.vector(X_new %*% object$fit_obj$coef)
}

# -----------------------------------------------------------------------
# bootstrap: Performs bootstrap estimate of standard error for lm_ac
# Samples holdout sets repeatedly, fits model, collects coefs & errors.
# -----------------------------------------------------------------------
bootstrap <- function(
    data, yName, coef_name, reps = 100, holdout_size = 0.2, ...
) {
  n <- nrow(data)
  if (n < 2) stop("Not enough rows for bootstrap.")
  if (holdout_size <= 0 || holdout_size >= 1) stop("'holdout_size' should be in (0,1).")
  
  coefs <- numeric(reps)
  mspe_values <- numeric(reps)
  mape_values <- numeric(reps)
  size <- max(1L, floor(holdout_size * n))
  
  for (i in seq_len(reps)) {
    holdout_set <- sample.int(n, size = size)
    obj <- lm_ac(data, yName, holdout = holdout_set, ...)
    
    beta <- obj$fit_obj$coef
    names(beta) <- obj$fit_obj$colnames
    coefs[i] <- if (coef_name %in% names(beta)) beta[[coef_name]] else NA_real_
    
    te <- obj$testing_data
    te <- te[!is.na(te[[yName]]), , drop = FALSE]
    if (nrow(te) == 0) {
      mspe_values[i] <- NA
      mape_values[i] <- NA
      next
    }
    y_true <- te[[yName]]
    y_pred <- predict(obj, te)
    ok <- !is.na(y_pred)
    if (any(ok)) {
      mspe_values[i] <- mean((y_true[ok] - y_pred[ok])^2)
      mape_values[i] <- mean(abs((y_true[ok] - y_pred[ok]) / y_true[ok]))
    } else {
      mspe_values[i] <- NA
      mape_values[i] <- NA
    }
  }
  
  avg_coef <- mean(coefs, na.rm = TRUE)
  se <- sd(coefs, na.rm = TRUE) / sqrt(sum(!is.na(coefs)))
  CI <- c(avg_coef - (1.96 * se), avg_coef + (1.96 * se))
  
  list(
    coef_mean = avg_coef,
    coef_se   = se,
    coef_CI   = CI,
    coefs     = coefs,
    mspe      = mean(mspe_values, na.rm = TRUE),
    mape      = mean(mape_values, na.rm = TRUE)
  )
}