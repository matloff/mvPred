library(regtools)

# lm_ac: Linear Model using Available Cases (pairwise deletion)
# ------------------------------------------------------------
# Implements regression when X or y contain NAs by computing
# each element of X'X and X'y as the *average* of intact pairs,
# then solving the normal equations directly on those averages
# (per Prof. Matloff’s Option 2).

# -----------------------------------------------------------------------
# Constructor: Initializes an lm_ac object with formula and raw data.
# -----------------------------------------------------------------------
lm_ac <- function(data, yName, holdout = NULL, ...) {
  # ----------------------------
  # Input validation
  # ----------------------------
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame.")
  
  if (!(yName %in% names(data))){
    stop(paste("Column", yName, "not found in data.")) 
  }
  
  if (!is.numeric(data[[yName]])) stop("Response variable must be numeric.")
  if (nrow(data) == 0) stop("Input data is empty.")

  # ----------------------------
  # Holdout set for bootstrap
  # ----------------------------

  if (!is.null(holdout)) {
    if (!is.logical(holdout) || length(holdout) != nrow(data)) {
      stop("Holdout is invalid")
    }
    training_data <- data[!holdout, 1:ncol(data), drop = FALSE]
    testing_data <- data[holdout, 1:ncol(data), drop = FALSE]
  } else {
    training_data <- data
    testing_data <- NULL
  }
  
  # ----------------------------
  # Convert character/logical to factor, then to dummy variables
  # ----------------------------
  for (col in names(data)) {
    if (is.character(data[[col]]) || is.logical(data[[col]])) {
      data[[col]] <- as.factor(data[[col]])
    }
  }
  
  data <-  as.data.frame(regtools::factorsToDummies(data))
  
  # ----------------------------
  # Clean column names (in case of dummies)
  # ----------------------------
  names(data) <- make.names(names(data), unique = TRUE)
  
  # ----------------------------
  # Construct formula
  # ----------------------------
  formula <- reformulate(setdiff(names(data), yName), response = yName)

  # ----------------------------
  # Build model.frame and model.matrix with do.call to allow user-specified
  # options 
  # ----------------------------
  mf_args <- list(formula = formula, data = data, na.action = NULL, ...)
  mf <- do.call(model.frame, mf_args)
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)

  # ----------------------------
  # Compute average-based XtX and Xty
  # ----------------------------
  XtX_avg <- ac_mul(X)
  Xty_avg <- ac_vec(X, y)

  # ----------------------------
  # Solve system and catch singular matrix errors
  # ----------------------------
  if (anyNA(XtX_avg) || anyNA(Xty_avg)) {
    stop("
      Missing values in system matrices: 
      need at least one intact pair per term.
    ") 
  }

  beta_hat <- tryCatch(
    solve(XtX_avg, Xty_avg),
    error = function(e) stop("XtX is singular; regression cannot proceed.")
  )

  # ----------------------------
  # Return model
  # ----------------------------
  obj <- list(
    data     = training_data,
    testing_data = testing_data,
    yName    = yName,
    formula  = formula,
    fit_obj  = list(
      coef     = beta_hat,
      colnames = colnames(X)
    )
  )
  class(obj) <- "lm_ac"
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
# build design matrix and compute y_hat = X_new %*% β_hat
# -----------------------------------------------------------------------
predict.lm_ac <- function(object, newdata, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted.")
  X_new <- model.matrix(object$formula, newdata)
  as.vector(X_new %*% object$fit_obj$coef)
}

# -----------------------------------------------------------------------
# bootstrap: Performs bootstrap estimate of standard error for lm_ac
# -----------------------------------------------------------------------

bootstrap <- function(
  data, yName, coef_name, reps = 100, holdout_size = 0.2, ...
) {
  coefs <- numeric(reps)
  mspe_values <- numeric(reps)
  mape_values <- numeric(reps)
  size <- holdout_size * nrow(data)

  for(i in 1:reps) {
    holdout_set <- sample(nrow(data), size = size)

    obj <- lm_ac(data, yName, holdout = holdout_set, ...)

    if (is.null(obj)) {
      coefs[i] <- NA
      mspe_values[i] <- NA
      mape_values[i] <- NA
      next
    }

    coef <- coef(obj)
    if (coef_name %in% names(coef)) {
        coefs[i] <- coef[coef_name]
    } else {
        coefs[i] <- NA
    }

    # Calculate prediction errors
    y_actual <- obj$testing_data[[yName]]
    y_pred <- predict.lm_ac(obj, newdata = obj$testing_data)

    mspe_values[i] <- mean((y_actual - y_pred)^2)
    mape_values[i] <- mean(abs((y_actual - y_pred) / y_actual))
  }

  # Calculate confidence interval
  avg_coef <- mean(coefs)
  
  se <- sd(coefs) / sqrt(reps)
  CI <- c(avg_coef - (1.96 * se), avg_coef + (1.96 * se))

  # Summarize results
  results <- list(
    coef_mean = avg_coef,
    coef_se = se,
    coef_CI = CI,
    coefs = coefs,
    mspe = mean(mspe_values),
    mape = mean(mape_values)
  )
  return(results)
}
