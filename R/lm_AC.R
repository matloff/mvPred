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