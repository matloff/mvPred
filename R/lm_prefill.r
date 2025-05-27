library(mice)
library(Amelia)
library(missForest)

# Internal constructor (hidden from user)
lm_prefill_create <- function(yName, data) {
  formula <- reformulate(".", response = yName)
  environment(formula) <- environment()
  
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

# Impute Data
impute.lm_prefill <- function(object, method = "mice", m = 5, ...) {
  if (method == "complete") {
    object$imputed_data <- na.omit(object$data)
    object$impute_method <- "complete"
    
  } else if (method == "mice") {
    data_for_impute <- as.data.frame(object$data)
    object$imputed_data <- mice(data_for_impute, m = m, ...)
    object$impute_method <- "mice"
    
  } else if (method == "amelia") {
    object$imputed_data <- amelia(object$data, m = m, ...)
    object$impute_method <- "amelia"
    
  } else if (method == "missforest") {
    data_for_impute <- as.data.frame(object$data)
    mf_result <- missForest(data_for_impute, ...)
    object$imputed_data <- mf_result$ximp
    object$impute_method <- "missforest"
    object$missforest_OOBerror <- mf_result$OOBerror
    
  } else {
    stop("Unknown imputation method.")
  }
  object
}

# Fit Model
fit.lm_prefill <- function(object, ...) {
  if (is.null(object$imputed_data)) stop("No imputed data found. Call impute() first.")
  
  if (object$impute_method %in% c("complete", "missforest")) {
    object$fit_obj <- lm(object$formula, data = object$imputed_data)
    
  } else if (object$impute_method == "mice") {
    fits <- lapply(1:object$imputed_data$m, function(i) {
      dat <- complete(object$imputed_data, action = i)
      lm(object$formula, data = dat)
    })
    object$fit_obj <- as.mira(fits)
    
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

# Summary method returning lm-style summaries for each imputed dataset
summary.lm_prefill <- function(object, ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  
  if (object$impute_method %in% c("complete", "missforest")) {
    return(summary(object$fit_obj))
    
  } else if (object$impute_method == "mice") {
    # mira object: list of lm fits in $analyses
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


# User-facing function: creates object, imputes, fits, returns summary
lm_prefill <- function(yName, data, method = "mice", m = 5, ...) {
  obj <- lm_prefill_create(yName, data)
  obj <- impute.lm_prefill(obj, method = method, m = m, ...)
  obj <- fit.lm_prefill(obj)
  # Return the full object, not summary
  obj
}

# Predict method averaging predictions across imputations
predict.lm_prefill <- function(object, newdata, type = "response", ...) {
  if (is.null(object$fit_obj)) stop("Model not fitted. Call fit() first.")
  
  # Remove rows with missing values in newdata
  missing_rows <- apply(newdata, 1, function(row) any(is.na(row)))
  if (any(missing_rows)) {
    warning("Some rows in newdata have missing values. These rows will be omitted from prediction.")
    newdata <- newdata[!missing_rows, , drop = FALSE]
  }
  if (nrow(newdata) == 0) return(numeric(0))
  
  if (object$impute_method %in% c("complete", "missforest")) {
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
#lm_obj <- lm_prefill("mpg", data, method = "mice", m = 5)
#summary(lm_obj)      # calls summary.lm_prefill()
#predict(lm_obj, newdata)  # calls predict.lm_prefill()
