#Constructor for Complete Cases and MI
lm_prefill <- function(yName, data) {
  # yName: character string, the name of the target variable
  # data: data.frame or tibble

  # Build the formula: yName ~.
  formula <- reformulate(".", response = yName)

  obj <- list(
    data = data,
    yName = yName,
    formula = formula,
    fit_obj = NULL,
    imputed_data = NULL,
    impute_method = NULL
  )
  class(obj) <- "lm_prefill"
  obj
}

#Impute Data
impute.lm_prefill <- function(object, method = "mice", ...) {
  if (method == "complete") {
    object$imputed_data <- na.omit(object$data)
    object$impute_method <- "complete"
  } else if (method == "mice") {
    library(mice)
    object$imputed_data <- mice(object$data, m = 5, ...)
    object$impute_method <- "mice"
  } else if (method == "amelia") {
    library(Amelia)
    object$imputed_data <- amelia(object$data, m = 5, ...)
    object$impute_method <- "amelia"
  } else if (method == "missforest") {
    library(missForest)
    # missForest requires a data.frame, not a tibble
    data_for_impute <- as.data.frame(object$data)
    mf_result <- missForest(data_for_impute, ...)
    object$imputed_data <- mf_result$ximp  # imputed data frame
    object$impute_method <- "missforest"
    object$missforest_OOBerror <- mf_result$OOBerror  # optional: store OOB error
  } else {
    stop("Unknown imputation method.")
  }
  object
}

#fit data

fit.lm_prefill <- function(object, ...) {
  if (is.null(object$imputed_data)) stop("No imputed data found. Call impute() first.")
  if (object$impute_method == "complete" || object$impute_method == "missforest") {
    object$fit_obj <- lm(object$formula, data = object$imputed_data)
  } else if (object$impute_method == "mice") {
    object$fit_obj <- with(object$imputed_data, lm(formula = object$formula))
  } else if (object$impute_method == "amelia") {
    fits <- lapply(object$imputed_data$imputations, function(dat) {
      lm(object$formula, data = dat)
    })
    object$fit_obj <- fits
  } else {
    stop("Unknown imputation method.")
  }
  object
}

#Summary function
summary.lm_prefill <- function(object, ...) {
    if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
    
    if (object$impute_method %in% c("complete", "missforest")) {
        # Single imputation: just show the summary
        return(summary(object$fit_obj))
    } else if (object$impute_method == "mice") {
        # MICE: fit_obj is a mira object
        return(summary(object$fit_obj))
    } else if (object$impute_method == "amelia") {
        # Amelia: fit_obj is a list of lm objects
        cat("Summary for first imputed dataset (Amelia):\n")
        return(summary(object$fit_obj[[1]]))
    } else {
        stop("Unknown imputation method.")
    }
}



#Predict data

predict.lm_prefill <- function(object, newdata, type = "response", ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  
  # Optionally handle missing values in newdata
  missing_rows <- apply(newdata, 1, function(row) any(is.na(row)))
  if (any(missing_rows)) {
    warning("Some rows in newdata have missing values. These rows will be omitted from prediction.")
    newdata <- newdata[!missing_rows, , drop = FALSE]
  }
  if (nrow(newdata) == 0) return(numeric(0))
  
  # Prediction logic for each imputation method
  if (object$impute_method %in% c("complete", "missforest")) {
    # Both use a single fitted lm object
    predict(object$fit_obj, newdata = newdata, type = type, ...)
    
  } else if (object$impute_method == "mice") {
    fits <- object$fit_obj$analyses
    preds <- do.call(cbind, lapply(fits, function(fit) predict(fit, newdata = newdata, type = type, ...)))
    if (is.null(dim(preds))) mean(preds) else rowMeans(preds)
    
  } else if (object$impute_method == "amelia") {
    preds <- do.call(cbind, lapply(object$fit_obj, function(fit) predict(fit, newdata = newdata, type = type, ...)))
    if (is.null(dim(preds))) mean(preds) else rowMeans(preds)
    
  } else {
    stop("Unknown imputation method in object.")
  }
}

#Example Usage
# Step 1: Create lm_prefill object with target variable "mpg"
lm_obj <- lm_prefill(yName = "mpg", data = mtcars_missing)

# Step 2: Impute missing data using one of the methods, e.g., "mice"
lm_obj <- impute(lm_obj, method = "mice", m = 5, maxit = 5, seed = 123)

# Step 3: Fit the linear model on the imputed data
lm_obj <- fit(lm_obj)

# Step 4: Summarize the fitted model
summary(lm_obj)
                                   
