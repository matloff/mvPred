# lm_AC: Linear Model using Available Cases (pairwise deletion)
# ------------------------------------------------------------
# Implements regression when X or y contain NAs by computing
# each element of X'X and X'y as the "average" of intact pairs,
# then solving the normal equations directly on those averages

# -----------------------------------------------------------------------
# Constructor: Initializes an lm_AC object with formula and raw data.
# -----------------------------------------------------------------------
lm_AC <- function(yName, data) {
  formula <- reformulate(".", response = yName)
  obj <- list(
    data    = data,        # original data.frame (may contain NAs)
    yName   = yName,       # name of the response variable
    formula = formula,     # will be used to build X and y
    fit_obj = NULL         # placeholder for fitted model results
  )
  class(obj) <- "lm_AC"
  obj
}

# -----------------------------------------------------------------------
# ac_mul: Compute (1/n_i_j) * sum[X_i * X_j] over available pairs
# For each pair of columns (i,j), only rows where both are non-NA
# are used, and we take the mean of the products. That yields
# an estimate of E[X_i X_j].
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

# dispatch generic
fit <- function(object, ...) UseMethod("fit")

# -----------------------------------------------------------------------
# fit.lm_AC: Fit the model using Option 2 (average over intact pairs)
#
# 1. Build model frame & design matrix X, response y (allow NAs).
# 2. Compute pairwise averages XtX_avg = (1/n_ij) Σ X_i X_j,
#    and Xty_avg = (1/n_i) Σ X_i y.
# 3. Solve β = (XtX_avg)^{-1} Xty_avg, which is equivalent to
#    solving on sums, since the n-scaling cancels out.
# -----------------------------------------------------------------------
fit.lm_AC <- function(object, ...) {
  data    <- object$data
  formula <- object$formula
  
  # pull out X and y, preserving NA rows
  mf <- model.frame(formula, data, na.action = NULL)
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  
  # compute pairwise-average matrices per Option 2
  XtX_avg <- ac_mul(X)
  Xty_avg <- ac_vec(X, y)
  
  # stop if any element is still NA (i.e. no intact pairs for some term)
  if (anyNA(XtX_avg) || anyNA(Xty_avg)) {
    stop("Missing values in system matrices: need at least one intact pair per term.")
  }
  
  # solve normal equations on averages—multiplier cancels out
  beta_hat <- solve(XtX_avg, Xty_avg)
  
  # store results
  object$fit_obj <- list(
    coef    = beta_hat,
    colnames = colnames(X),
    formula  = formula
  )
  object
}

# -----------------------------------------------------------------------
# summary.lm_AC: Print coefficients (same shape as lm())
# -----------------------------------------------------------------------
summary.lm_AC <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted; call fit() first.")
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
  if (is.null(object$fit_obj)) stop("Model not fitted; call fit() first.")
  X_new <- model.matrix(object$formula, newdata)
  as.vector(X_new %*% object$fit_obj$coef)
}
