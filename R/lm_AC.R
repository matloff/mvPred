library(regtools)

# lm_AC: Linear Model using Available Cases (pairwise deletion)
# ------------------------------------------------------------
# Implements regression when X or y contain NAs by computing
# each element of X'X and X'y as the *average* of intact pairs,
# then solving the normal equations directly on those averages
# (per Prof. Matloff’s Option 2).

# -----------------------------------------------------------------------
# Constructor: Initializes an lm_AC object with formula and raw data.
# -----------------------------------------------------------------------
lm_AC <- function(data, yName, ...) {
  # ----------------------------
  # Input validation
  # ----------------------------
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame.")
  if (!(yName %in% names(data))) stop(paste("Column", yName, "not found in data."))
  if (!is.numeric(data[[yName]])) stop("Response variable must be numeric.")
  if (nrow(data) == 0) stop("Input data is empty.")
  
  # ----------------------------
  # Convert character/logical to factor, then to dummy variables
  # ----------------------------
  for (col in names(data)) {
    if (is.character(data[[col]]) || is.logical(data[[col]])) {
      data[[col]] <- as.factor(data[[col]])
    }
  }
  data <- regtools::factorsToDummies(data)
  
  # # ----------------------------
  # # 3. Remove columns with all NA
  # # ----------------------------
  # keep_cols <- sapply(data, function(col) sum(!is.na(col)) > 0)
  # data <- data[, keep_cols, drop = FALSE]
  # if (!(yName %in% names(data))) stop("Response variable became empty after dropping NA columns.")
  
  # ----------------------------
  # Clean column names (in case of dummies)
  # ----------------------------
  names(data) <- make.names(names(data), unique = TRUE)
  
  # ----------------------------
  # Construct formula
  # ----------------------------
  formula <- reformulate(setdiff(names(data), yName), response = yName)

  # ----------------------------
  # Build model.frame and model.matrix with do.call to allow user-specified options
  # ----------------------------
  mf_args <- list(formula = formula, data = data, na.action = NULL, ...)
  mf <- do.call(model.frame, mf_args)
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)

  # ----------------------------
  # Drop constant columns (0 variance)
  # ----------------------------
  # const_idx <- apply(X, 2, function(col) var(col, na.rm = TRUE) == 0)
  # if (any(const_idx)) {
  #   X <- X[, !const_idx, drop = FALSE]
  # }

  # ----------------------------
  # 8. Compute average-based XtX and Xty
  # ----------------------------
  XtX_avg <- ac_mul(X)
  Xty_avg <- ac_vec(X, y)

  # ----------------------------
  # 9. Solve system and catch singular matrix errors
  # ----------------------------
  if (anyNA(XtX_avg) || anyNA(Xty_avg)) {
    stop("Missing values in system matrices: need at least one intact pair per term.")
  }

  beta_hat <- tryCatch(
    solve(XtX_avg, Xty_avg),
    error = function(e) stop("XtX is singular; regression cannot proceed.")
  )

  # ----------------------------
  # 10. Return model
  # ----------------------------
  obj <- list(
    data     = data,
    yName    = yName,
    formula  = formula,
    fit_obj  = list(
      coef     = beta_hat,
      colnames = colnames(X)
    )
  )
  class(obj) <- "lm_AC"
  obj
}


# -----------------------------------------------------------------------
# ac_mul: Compute (1/n_i_j) * sum[X_i * X_j] over available pairs
# For each pair of columns (i,j), only rows where both are non-NA
# are used, and we take the *mean* of the products.  That yields
# an estimate of E[X_i X_j] (prof’s Option 2), not the raw sum.
# -----------------------------------------------------------------------
ac_mul <- function(A) {
  if (!is.matrix(A)) A <- as.matrix(A)
  p <- ncol(A)
  result <- matrix(NA, p, p)
  
  for (i in 1:p) {
    for (j in i:p) {
      # find rows where both X[,i] and X[,j] are observed
      valid_idx <- which(!is.na(A[,i]) & !is.na(A[,j]))
      if (length(valid_idx) > 0) {
        # average of X_i * X_j over intact rows:
        result[i,j] <- mean(A[valid_idx, i] * A[valid_idx, j])
        result[j,i] <- result[i,j]  # symmetry
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
  p <- ncol(X)
  result <- numeric(p)
  
  for (i in 1:p) {
    valid_idx <- which(!is.na(X[,i]) & !is.na(y))
    if (length(valid_idx) > 0) {
      # average of X_i * y over intact rows:
      result[i] <- mean(X[valid_idx, i] * y[valid_idx])
    } else {
      result[i] <- NA
    }
  }
  
  result
}

# -----------------------------------------------------------------------
# summary.lm_AC: Print coefficients (same shape as lm())
# -----------------------------------------------------------------------
summary.lm_AC <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  coefs <- object$fit_obj$coef
  names(coefs) <- object$fit_obj$colnames
  print(coefs)
  invisible(coefs)
}

# -----------------------------------------------------------------------
# predict.lm_AC: Given newdata (with possible NAs in predictors),
# build design matrix and compute y_hat = X_new %*% β_hat
# -----------------------------------------------------------------------
predict.lm_AC <- function(object, newdata, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  X_new <- model.matrix(object$formula, newdata)
  as.vector(X_new %*% object$fit_obj$coef)
}