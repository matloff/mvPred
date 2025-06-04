# source("lm_AC_2.R")       # Contains the lm_ac() function using available-case analysis
# source("lm_prefill.R")    # Contains the lm_prefill() function for imputing missing data before modeling

# Define a fallback operator: returns 'a' if it's not NULL, otherwise returns 'b'.
# This helps us avoid writing repetitive conditional checks for missing objects.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Function to compare different modeling approaches for handling missing data:
# 1. Available-case analysis via lm_ac()
# 2. Multiple/Single imputation via lm_prefill()
# 3. Standard complete-case regression via lm()
lm_compare <- function(data, yName, impute_method = "mice", m = 5, use_dummies = FALSE, metrics = c("rmse", "r2", "beta")) {
  
  # Ensure the specified response variable exists in the dataset
  if (!(yName %in% names(data))) stop("yName not found in data.")
  
  # Remove all rows with any missing values to create a clean dataset for training/testing
  # This ensures a fair and consistent test set across all methods
  complete_data <- na.omit(data)
  if (nrow(complete_data) < 10) warning("Few complete cases. Results may not generalize.")
  
  # Randomly split the complete data into training (80%) and testing (20%) sets
  set.seed(42)
  idx <- sample(seq_len(nrow(complete_data)), size = floor(0.8 * nrow(complete_data)))
  train_data <- complete_data[idx, ]
  test_data  <- complete_data[-idx, ]
  y_true     <- test_data[[yName]]  # Ground truth for model evaluation
  
  # --------------------------
  # Train and evaluate lm_ac (Available Case analysis)
  # --------------------------
  # Fit a model that handles missing values in predictors/response via pairwise deletion
  model_ac <- lm_ac(train_data, yName)
  pred_ac <- predict(model_ac, test_data)
  
  # Extract estimated coefficients and assign readable names
  beta_ac <- model_ac$fit_obj$coef
  names(beta_ac) <- model_ac$fit_obj$colnames
  
  # --------------------------
  # Train and evaluate lm_prefill (Imputation-based regression)
  # --------------------------
  # Impute missing values in the training set and fit a linear model to the imputed data
  model_pre <- lm_prefill(train_data, yName, impute_method = impute_method, m = m, use_dummies = use_dummies)
  pred_pre <- predict(model_pre, newdata = test_data)
  
  # Aggregate coefficients from multiple imputations if necessary
  if (impute_method %in% c("mice", "amelia")) {
    fits <- model_pre$fit_obj$analyses %||% model_pre$fit_obj
    beta_pre <- Reduce("+", lapply(fits, coef)) / m  # Average coefficients across imputations
  } else {
    beta_pre <- coef(model_pre$fit_obj)
  }
  
  # --------------------------
  # Train and evaluate standard lm (Complete-case analysis)
  # --------------------------
  # Fit a conventional linear model using only complete cases (same as train_data)
  formula <- reformulate(setdiff(names(train_data), yName), response = yName)
  model_lm <- lm(formula, data = train_data)
  pred_lm <- predict(model_lm, newdata = test_data)
  beta_lm <- coef(model_lm)
  
  # --------------------------
  # Metric definitions for evaluating predictions
  # --------------------------
  # RMSE: measures average prediction error magnitude
  rmse <- function(y, y_hat) sqrt(mean((y - y_hat)^2))
  
  # R²: measures proportion of variance explained by the model
  r2 <- function(y, y_hat) 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
  
  # --------------------------
  # Construct and return results
  # --------------------------
  result <- list()
  
  if ("rmse" %in% metrics) {
    result$RMSE <- c(
      ac = rmse(y_true, pred_ac),
      prefill = rmse(y_true, pred_pre),
      lm = rmse(y_true, pred_lm)
    )
  }
  
  if ("r2" %in% metrics) {
    result$R2 <- c(
      ac = r2(y_true, pred_ac),
      prefill = r2(y_true, pred_pre),
      lm = r2(y_true, pred_lm)
    )
  }
  
  if ("beta" %in% metrics) {
    # Union of all coefficient names across models ensures consistent column comparison
    all_names <- union(names(beta_ac), union(names(beta_pre), names(beta_lm)))
    result$Beta <- data.frame(
      term = all_names,
      ac = beta_ac[all_names],
      prefill = beta_pre[all_names],
      lm = beta_lm[all_names]
    )
  }
  
  return(result)
}

# Example usage:
# results <- lm_compare(airquality, yName = "Ozone", impute_method = "mice", m = 5)
# results$RMSE   # View RMSE for each method
# results$R2     # View R² for each method
# results$Beta   # View estimated coefficients from all three approaches
