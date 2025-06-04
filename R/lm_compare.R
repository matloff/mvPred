# Load implementations for available-case regression and imputation-based regression
# source("lm_AC_2.R")       # Contains the lm_ac() function using available-case analysis
# source("lm_prefill.R")    # Contains the lm_prefill() function for imputing missing data before modeling

# Fallback operator: if `a` is not NULL, return `a`; otherwise, return `b`
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Function to compare the performance of various fitted regression models
lm_compare <- function(models, data, yName, test_frac = 0.2, metrics = c("rmse", "r2", "beta")) {
  if (!(yName %in% names(data))) stop("yName not found in input data.")
  
  # Define evaluation metrics for model performance
  rmse <- function(y, y_hat) sqrt(mean((y - y_hat)^2))
  r2   <- function(y, y_hat) 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
  
  # Generate a test set using only complete cases to ensure fair evaluation across models
  complete_data <- na.omit(data)
  if (nrow(complete_data) < 10) warning("Few complete cases. Results may not generalize.")
  set.seed(42)  # For reproducibility of the train-test split
  idx <- sample(seq_len(nrow(complete_data)), size = floor((1 - test_frac) * nrow(complete_data)))
  train_data_lm <- complete_data[idx, ]
  test_data     <- complete_data[-idx, ]
  y_true        <- test_data[[yName]]
  
  # Initialize containers for predictions and beta coefficients from all models
  preds <- list()
  betas <- list()
  
  # Iterate through the user-provided models and extract predictions and coefficients
  for (name in names(models)) {
    model <- models[[name]]
    
    # Handle lm_ac models
    if (inherits(model, "lm_ac")) {
      preds[[name]] <- predict(model, test_data)
      beta <- model$fit_obj$coef
      names(beta) <- model$fit_obj$colnames
      betas[[name]] <- beta
      
      # Handle lm_prefill models
    } else if (inherits(model, "lm_prefill")) {
      preds[[name]] <- predict(model, newdata = test_data)
      
      # If multiple imputations were used, average coefficients across imputations
      if (!is.null(model$fit_obj$analyses)) {
        fits <- model$fit_obj$analyses
        betas[[name]] <- Reduce("+", lapply(fits, coef)) / length(fits)
      } else if (is.list(model$fit_obj)) {
        fits <- model$fit_obj
        betas[[name]] <- Reduce("+", lapply(fits, coef)) / length(fits)
      } else {
        betas[[name]] <- coef(model$fit_obj)
      }
      
      # Skip any model that doesn’t match known types
    } else {
      warning(paste("Skipping unknown model type:", name))
    }
  }
  
  # Internally train and evaluate a standard complete-case lm model
  formula <- reformulate(setdiff(names(train_data_lm), yName), response = yName)
  model_lm <- lm(formula, data = train_data_lm)
  preds$lm <- predict(model_lm, newdata = test_data)
  betas$lm <- coef(model_lm)
  
  # Collect evaluation results for each requested metric
  result <- list()
  
  if ("rmse" %in% metrics) {
    result$RMSE <- sapply(preds, function(p) rmse(y_true, p))
  }
  
  if ("r2" %in% metrics) {
    result$R2 <- sapply(preds, function(p) r2(y_true, p))
  }
  
  if ("beta" %in% metrics) {
    all_terms <- Reduce(union, lapply(betas, names))
    beta_df <- data.frame(term = all_terms)
    for (name in names(betas)) {
      beta_df[[name]] <- betas[[name]][all_terms]
    }
    result$Beta <- beta_df
  }
  
  return(result)
}


# Example usage

# Load the dataset
# data(airquality)
# 
# # Fit models externally using available-case and imputation methods
# model_ac <- lm_ac(airquality, "Ozone")
# model_pre <- lm_prefill(airquality, "Ozone", impute_method = "mice", m = 5)
# 
# # Compare performance using the compare function
# results <- lm_compare(
#   models = list(ac = model_ac, prefill = model_pre),
#   data = airquality,
#   yName = "Ozone"
# )
# 
# # View comparison metrics
# results$RMSE
# results$R2
# results$Beta

