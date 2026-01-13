#library(mice)
#library(Amelia)
#library(missForest)

# Helper: Detect categorical columns and optionally convert to dummies
detect_and_prepare_noms <- function(data, use_dummies = FALSE) {
  noms <- c()
  for (col in names(data)) {
    if (is.factor(data[[col]]) || is.character(data[[col]]) || is.logical(data[[col]])) {
      noms <- c(noms, col)
      if (is.character(data[[col]]) || is.logical(data[[col]])) {
        data[[col]] <- as.factor(data[[col]])
      }
    }
  }
  if (use_dummies && length(noms) > 0) {
    if (!requireNamespace("regtools", quietly = TRUE)) {
      stop("Package 'regtools' needed for dummy variable creation. Please install it.")
    }
    data <- regtools::factorsToDummies(data)
    noms <- NULL
  }
  list(data = data, noms = noms)
}
# ==============================================================================
# PREPROCESSING FUNCTION
# ==============================================================================

preprocess_data <- function(
    train,
    validation,
    target,
    task = c("classification", "regression"),
    categorical_cols = NULL
) {
  task <- match.arg(task)
  
  # Classification-specific preprocessing
  if (task == "classification" && !is.null(categorical_cols)) {
    for (col in categorical_cols) {
      # Convert to character for manipulation
      train[[col]] <- as.character(train[[col]])
      validation[[col]] <- as.character(validation[[col]])
      
      # Identify rare levels in training (frequency < 10)
      freq <- table(train[[col]])
      rare_levels <- names(freq[freq < 10])
      
      # Replace rare levels with "Other" in training
      train[[col]][train[[col]] %in% rare_levels] <- "Other"
      
      # Replace rare levels with "Other" in validation
      validation[[col]][validation[[col]] %in% rare_levels] <- "Other"
      
      # Align levels: both train and validation should have same factor levels
      lvls <- union(unique(train[[col]]), "Other")
      train[[col]] <- factor(train[[col]], levels = lvls)
      validation[[col]] <- factor(validation[[col]], levels = lvls)
    }
  }
  
  # Regression: ensure numeric columns are numeric
  if (task == "regression") {
    numeric_cols <- setdiff(names(train), c(target, categorical_cols))
    for (col in numeric_cols) {
      if (!is.numeric(train[[col]])) {
        train[[col]] <- as.numeric(as.character(train[[col]]))
      }
      if (!is.numeric(validation[[col]])) {
        validation[[col]] <- as.numeric(as.character(validation[[col]]))
      }
    }
  }
  
  return(list(train = train, validation = validation))
}

# ==============================================================================
# METRICS FUNCTIONS
# ==============================================================================

compute_metrics_classification <- function(true, pred_class, pred_prob) {
  true_pos <- sum(pred_class == 1 & true == 1)
  false_pos <- sum(pred_class == 1 & true == 0)
  true_neg <- sum(pred_class == 0 & true == 0)
  false_neg <- sum(pred_class == 0 & true == 1)
  
  div <- function(a, b) ifelse(b == 0, NA, a / b)
  
  accuracy  <- div(true_pos + true_neg, true_pos + false_pos + true_neg + false_neg)
  precision <- div(true_pos, true_pos + false_pos)
  recall    <- div(true_pos, true_pos + false_neg)
  f1        <- div(2 * precision * recall, precision + recall)
  
  return(list(
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1 = f1
  ))
}

compute_metrics_regression <- function(true, pred) {
  mse <- mean((true - pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(true - pred))
  r_squared <- 1 - sum((true - pred)^2) / sum((true - mean(true))^2)
  
  return(list(
    mse = mse,
    rmse = rmse,
    mae = mae,
    r_squared = r_squared
  ))
}

# Helper function to calculate missingness percentage
missing_percent <- function(df) {
  sapply(df, function(col) {
    sum(is.na(col)) / length(col) * 100
  })
}

# ==============================================================================
# BOOTSTRAP FUNCTION
# ==============================================================================

bootstrap <- function(
    data,
    yName,
    task = c("classification", "regression"),
    categorical_cols = NULL,
    impute_method = "mice",
    m = 5,
    k = 5,
    use_dummies = FALSE,
    ...
) {
  
  task <- match.arg(task)
  
  # Validate inputs
  if (!yName %in% names(data)) {
    stop("Target variable '", yName, "' not found in data. Available columns: ", 
         paste(names(data), collapse = ", "))
  }
  
  n <- nrow(data)
  set.seed(123)
  
  # Create k folds
  fold_ids <- sample(rep(1:k, length.out = n))
  
  # Initialize storage
  metric_results <- vector("list", k)
  train_missing_list <- vector("list", k)
  validation_missing_list <- vector("list", k)
  
  # Loop through each fold
  for (i in 1:k) {
    cat("\n=== Processing fold", i, "of", k, "===\n")
    
    # Split data
    train <- data[fold_ids != i, ]
    validation <- data[fold_ids == i, ]
    
    cat("Train size:", nrow(train), "| Validation size:", nrow(validation), "\n")
    
    # Preprocess data
    preprocessed <- preprocess_data(
      train = train,
      validation = validation,
      target = yName,
      task = task,
      categorical_cols = categorical_cols
    )
    train <- preprocessed$train
    validation <- preprocessed$validation
    
    # Store missingness before removing NAs
    train_missing_list[[i]] <- missing_percent(train)
    validation_missing_list[[i]] <- missing_percent(validation)
    
    # Fit model using lm_prefill
    
    cat("Fitting model...\n")
    lm_model <- lm_prefill(
        train,
        yName,  
        impute_method = impute_method,
        m = m,
        use_dummies = use_dummies,
        ...
      )
    
    if (is.null(lm_model)) {
      cat("Skipping fold", i, "due to model fitting error\n")
      next
    }
    
    cat("Model fitted successfully\n")
    
    # Filter validation data before prediction
    # Remove incomplete rows
    validation_complete <- validation[complete.cases(validation), ]
    
    cat("Complete validation cases:", nrow(validation_complete), "\n")
    
    if (nrow(validation_complete) == 0) {
      cat("No complete cases in validation fold", i, ", skipping\n")
      next
    }
    
    # Prepare data for prediction (remove target)
    validation_pred <- validation_complete
    validation_pred[[yName]] <- NULL  
    
    # Ensure validation_pred has the same structure as training data (minus target)
    train_cols <- setdiff(names(train), yName)
    validation_pred <- validation_pred[, train_cols, drop = FALSE]
    
    cat("Predicting...\n")
    # Predict
    # Predict
    predictions <- predict(lm_prefill_model, newdata = validation_pred)
    
    if (is.null(predictions)) {
      cat("Skipping fold", i, "due to prediction error\n")
      next
    }
    
    cat("Predictions made successfully. Length:", length(predictions), "\n")
    
    # Get true values (already filtered to complete cases)
    true_values <- validation_complete[[yName]]
    
    # Compute metrics based on task
    if (task == "classification") {
      predicted_class <- ifelse(predictions > 0.5, 1, 0)
      metric_results[[i]] <- compute_metrics_classification(
        true_values,
        predicted_class,
        predictions
      )
    } else {
      metric_results[[i]] <- compute_metrics_regression(true_values, predictions)
    }
  }
  
  # Compute average metrics
  # Filter out NULL results from failed folds
  valid_results <- metric_results[!sapply(metric_results, is.null)]
  
  if (length(valid_results) == 0) {
    stop("All folds failed. Cannot compute metrics.")
  }
  
  if (length(valid_results) < k) {
    cat("Warning:", k - length(valid_results), "folds failed\n")
  }
  
  metrics_df <- do.call(rbind, lapply(valid_results, as.data.frame))
  avg_metrics <- colMeans(metrics_df, na.rm = TRUE)
  
  # Average missingness across folds
  avg_train_missing <- colMeans(do.call(rbind, train_missing_list), na.rm = TRUE)
  avg_validation_missing <- colMeans(do.call(rbind, validation_missing_list), na.rm = TRUE)
  
  # Return results
  return(list(
    metrics = avg_metrics,
    metrics_per_fold = metrics_df,
    train_missing = avg_train_missing,
    validation_missing = avg_validation_missing
  ))
}
# Internal constructor
lm_prefill_create <- function(data, yName, holdout = NULL) {
  formula <- reformulate(".", response = yName)
  
  if (!is.null(holdout)) {
    if (is.logical(holdout) && length(holdout) == nrow(data)) {
      training_data <- data[!holdout, 1:ncol(data), drop = FALSE]
      testing_data <- data[holdout, 1:ncol(data), drop = FALSE]
    } else if (is.numeric(holdout) && all(holdout %in% seq_len(nrow(data)))) {
      training_data <- data[-holdout, 1:ncol(data), drop = FALSE]
      testing_data  <- data[holdout, 1:ncol(data), drop = FALSE]
    } else {
      stop("Holdout is invalid")
    }
    
  } else {
    training_data <- data
    testing_data <- NULL
  }
  
  # drop NAs in testing data
  testing_data <- na.omit(testing_data)
  obj <- list(
    data = training_data,
    testing_data = testing_data,
    yName = yName,
    formula = formula,
    fit_obj = NULL,
    imputed_data = NULL,
    impute_method = NULL
  )
  class(obj) <- "lm_prefill"
  obj
}

# Main user-facing function
lm_prefill <- function(data, yName, impute_method = "mice", m = 5, use_dummies = FALSE, holdout = NULL, ...) {
  dots <- list(...)
  
  # Argument names for each function
  mice_formals <- names(formals(mice::mice))
  amelia_formals <- names(formals(Amelia::amelia))
  missforest_formals <- names(formals(missForest::missForest))
  lm_formals <- names(formals(stats::lm))
  
  # Remove arguments handled by our function
  handled <- c("impute_method", "m", "use_dummies", "data", "yName")
  
  # Special handling: if user passed method="qr" or "model.frame", it's for lm(), not mice
  method_for_lm <- FALSE
  if ("method" %in% names(dots)) {
    if (dots$method %in% c("qr", "model.frame")) {
      method_for_lm <- TRUE
    }
  }
  
  # Imputation arguments (for the chosen method)
  if (impute_method == "mice") {
    impute_args <- dots[names(dots) %in% mice_formals & !(names(dots) %in% handled)]
    if (method_for_lm) impute_args$method <- NULL
  } else if (impute_method == "amelia") {
    impute_args <- dots[names(dots) %in% amelia_formals & !(names(dots) %in% handled)]
  } else if (impute_method == "missforest") {
    impute_args <- dots[names(dots) %in% missforest_formals & !(names(dots) %in% handled)]
  } else {
    impute_args <- list()
  }
  
  # Modeling arguments (for lm)
  fit_args <- dots[names(dots) %in% lm_formals & !(names(dots) %in% handled)]
  if (!method_for_lm) fit_args$method <- NULL
  
  obj <- lm_prefill_create(data, yName, holdout = holdout)
  obj <- impute.lm_prefill(obj, impute_method = impute_method, m = m, use_dummies = use_dummies, impute_args = impute_args)
  obj <- fit.lm_prefill(obj, fit_args = fit_args)
  obj
}

# Impute Data
impute.lm_prefill <- function(object, impute_method = "mice", m = 5, use_dummies = FALSE, impute_args = list()) {
  if (impute_method == "complete") {
    object$imputed_data <- stats::na.omit(object$data)
    object$impute_method <- "complete"
  } else if (impute_method == "mice") {
    data_for_impute <- as.data.frame(object$data)
    mice_args <- c(list(data = data_for_impute, m = m), impute_args)
    object$imputed_data <- do.call(mice::mice, mice_args)
    object$impute_method <- "mice"
  } else if (impute_method == "amelia") {
    prep <- detect_and_prepare_noms(object$data, use_dummies = use_dummies)
    object$data <- prep$data
    noms <- prep$noms
    amelia_args <- c(list(x = object$data, m = m), impute_args)
    if (!is.null(noms) && length(noms) > 0) {
      amelia_args$noms <- noms
    }
    object$imputed_data <- do.call(Amelia::amelia, amelia_args)
    object$impute_method <- "amelia"
  } else if (impute_method == "missforest") {
    data_for_impute <- as.data.frame(object$data)
    mf_args <- c(list(xmis = data_for_impute), impute_args)
    mf_result <- do.call(missForest::missForest, mf_args)
    object$imputed_data <- mf_result$ximp
    object$impute_method <- "missforest"
    object$missforest_OOBerror <- mf_result$OOBerror
  } else {
    stop("Unknown imputation method.")
  }
  object
}

# Fit Model
fit.lm_prefill <- function(object, fit_args = list()) {
  if (is.null(object$imputed_data)) stop("No imputed data found. Call impute() first.")
  if (object$impute_method %in% c("complete", "missforest")) {
    lm_args <- c(list(formula = object$formula, data = object$imputed_data), fit_args)
    object$fit_obj <- do.call(stats::lm, lm_args)
  } else if (object$impute_method == "mice") {
    fits <- lapply(1:object$imputed_data$m, function(i) {
      dat <- mice::complete(object$imputed_data, action = i)
      lm_args <- c(list(formula = object$formula, data = dat), fit_args)
      do.call(stats::lm, lm_args)
    })
    object$fit_obj <- list(analyses = fits)
  } else if (object$impute_method == "amelia") {
    fits <- lapply(object$imputed_data$imputations, function(dat) {
      lm_args <- c(list(formula = object$formula, data = dat), fit_args)
      do.call(stats::lm, lm_args)
    })
    object$fit_obj <- fits
  } else {
    stop("Unknown imputation method.")
  }
  object
}

# Summary method returning lm-style summaries for each imputed dataset
summary.lm_prefill <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  if (object$impute_method %in% c("complete", "missforest")) {
    return(summary(object$fit_obj))
  } else if (object$impute_method == "mice") {
    summaries <- lapply(object$fit_obj$analyses, summary)
    class(summaries) <- "summary.lm_prefill_mice"
    return(summaries)
  } else if (object$impute_method == "amelia") {
    summaries <- lapply(object$fit_obj, summary)
    class(summaries) <- "summary.lm_prefill_amelia"
    return(summaries)
  } else {
    stop("Unknown imputation method.")
  }
}

# Predict method averaging predictions across imputations
predict.lm_prefill <- function(object, newdata, type = "response", ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  missing_rows <- !stats::complete.cases(newdata)
  if (any(missing_rows)) {
    warning("Some rows in newdata have missing values. These rows will be omitted from prediction.")
    newdata <- newdata[!missing_rows, , drop = FALSE]
  }
  if (nrow(newdata) == 0) return(numeric(0))
  if (object$impute_method %in% c("complete", "missforest")) {
    stats::predict(object$fit_obj, newdata = newdata, type = type, ...)
  } else if (object$impute_method == "mice") {
    fits <- object$fit_obj$analyses
    preds <- do.call(cbind, lapply(fits, function(fit) stats::predict(fit, newdata = newdata, type = type, ...)))
    if (is.null(dim(preds))) mean(preds) else rowMeans(preds)
  } else if (object$impute_method == "amelia") {
    preds <- do.call(cbind, lapply(object$fit_obj, function(fit) stats::predict(fit, newdata = newdata, type = type, ...)))
    if (is.null(dim(preds))) mean(preds) else rowMeans(preds)
  } else {
    stop("Unknown imputation method in object.")
  }
}

#Example Usage
#lm_obj <- lm_prefill(data,"mpg",method = "mice", m = 5)
#summary(lm_obj)      # calls summary.lm_prefill()
#predict(lm_obj, newdata)  # calls predict.lm_prefill()

