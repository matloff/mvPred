# lm_AC: Linear Model using Available Cases (pairwise deletion)

# -----------------------------------------------------------------------
# Constructor: Initializes an lm_AC object with formula and raw data.
# -----------------------------------------------------------------------
lm_AC <- function(yName, data) {
  formula <- reformulate(".", response = yName)
  obj <- list(
    data = data,        # original dataset
    yName = yName,      # name of the target variable
    formula = formula,  # model formula y ~ .
    fit_obj = NULL      # will store model fit results
  )
  class(obj) <- "lm_AC"
  obj
}


# -----------------------------------------------------------------------
# ac_mul: Computes X^T X using pairwise complete observations (Available Cases)
# For each (i, j), uses only rows where both X[,i] and X[,j] are not NA
# -----------------------------------------------------------------------
ac_mul <- function(A) {
  if (!is.matrix(A)) A <- as.matrix(A)
  
  p <- ncol(A)
  result <- matrix(NA, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in i:p) {
      # Find rows where both columns i and j are observed
      valid_idx <- which(!is.na(A[, i]) & !is.na(A[, j]))
      if (length(valid_idx) > 0) {
        value <- sum(A[valid_idx, i] * A[valid_idx, j])
        result[i, j] <- value
        result[j, i] <- value  # enforce symmetry
      } else {
        result[i, j] <- result[j, i] <- NA
      }
    }
  }
  
  colnames(result) <- colnames(A)
  rownames(result) <- colnames(A)
  return(result)
}

fit <- function(object, ...) {
  UseMethod("fit")
}

# -----------------------------------------------------------------------
# fit.lm_AC: Fits the linear model using pairwise deletion (Available Cases)
# Computes XtX and Xty manually using available cases for each term
# Solves: beta = (X^T X)^(-1) X^T y
# -----------------------------------------------------------------------
fit.lm_AC <- function(object, ...) {
  data <- object$data
  formula <- object$formula
  
  # Extract design matrix X and response vector y
  mf <- model.frame(formula, data, na.action = NULL)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  
  # Compute X^T X using available (non-NA) data per entry
  XtX <- ac_mul(X)
  
  # Compute X^T y using available data per feature
  Xty <- numeric(ncol(X))
  for (i in 1:ncol(X)) {
    valid_idx <- which(!is.na(X[, i]) & !is.na(y))
    Xty[i] <- sum(X[valid_idx, i] * y[valid_idx])
  }
  
  # Fail early if NA entries prevent solution
  if (anyNA(XtX) || anyNA(Xty)) stop("Missing values in system matrices.")
  
  # Solve normal equations for beta
  beta_hat <- solve(XtX, Xty)
  
  # Store coefficients and model metadata
  object$fit_obj <- list(
    coef = beta_hat,
    colnames = colnames(X),
    formula = formula
  )
  object
}

# -----------------------------------------------------------------------
# summary.lm_AC: Prints model coefficients
# -----------------------------------------------------------------------
summary.lm_AC <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  coefs <- object$fit_obj$coef
  names(coefs) <- object$fit_obj$colnames
  print(coefs)
  invisible(coefs)
}

# -----------------------------------------------------------------------
# predict.lm_AC: Makes predictions using the fitted coefficients
# Simply performs: y_hat = X_new %*% beta_hat
# -----------------------------------------------------------------------
predict.lm_AC <- function(object, newdata, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  X_new <- model.matrix(object$formula, newdata)
  as.vector(X_new %*% object$fit_obj$coef)
}