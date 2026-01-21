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


