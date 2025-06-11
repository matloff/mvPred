# Load implementations for available-case regression and imputation-based
# regression source("lm_AC_2.R")       # Contains the lm_ac() function using
# available-case analysis source("lm_prefill.R")    # Contains the lm_prefill()
# function for imputing missing data before modeling 

# Fallback operator: if `a` is not NULL, return `a`; otherwise, return `b` 
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Function to compare the performance of various fitted regression models
lm_compare <- function(models, yName,
                       metrics = c("rmse", "r2", "beta")) {
  
  rmse <- function(y, y_hat) sqrt(mean((y - y_hat)^2))
  r2   <- function(y, y_hat) 1 - sum((y - y_hat)^2) /
    sum((y - mean(y))^2)
  
  preds <- list(); betas <- list()
  
  for (nm in names(models)) {
    mdl <- models[[nm]]
    test <- mdl$testing_data
    if (is.null(test)) {
      warning(paste("Model", nm, "has no testing_data. Skipping."))
      next
    }
    
    ## Predictions
    preds[[nm]] <- predict(mdl, newdata = test)
    
    ## Coefficients (works for both lm_ac and pooled‑MI objects)
    if (inherits(mdl, "lm_ac")) {
      beta <- mdl$fit_obj$coef
      names(beta) <- mdl$fit_obj$colnames
      betas[[nm]] <- beta
    } else if (inherits(mdl, "lm_prefill")) {
      if (!is.null(mdl$fit_obj$analyses)) {          # “mice” list
        fits <- mdl$fit_obj$analyses
        betas[[nm]] <- Reduce("+", lapply(fits, coef)) / length(fits)
      } else if (is.list(mdl$fit_obj)) {             # list of lm objects
        fits <- mdl$fit_obj
        betas[[nm]] <- Reduce("+", lapply(fits, coef)) / length(fits)
      } else {                                       # single lm
        betas[[nm]] <- coef(mdl$fit_obj)
      }
    }
  }
  
  out <- list()
  
  if ("rmse" %in% metrics) {
    out$RMSE <- mapply(function(p, nm)
      rmse(models[[nm]]$testing_data[[yName]], p),
      preds, names(preds))
  }
  
  if ("r2" %in% metrics) {
    out$R2 <- mapply(function(p, nm)
      r2(models[[nm]]$testing_data[[yName]], p),
      preds, names(preds))
  }
  
  if ("beta" %in% metrics) {
    all_terms <- Reduce(union, lapply(betas, names))
    beta_df   <- data.frame(term = all_terms)
    for (nm in names(betas))
      beta_df[[nm]] <- betas[[nm]][all_terms]
    out$Beta <- beta_df
  }
  
  out
}


# data(airquality)
# set.seed(42)

# n   <- nrow(airquality)
# idx <- sample(n, size = floor(0.20 * n))    # 20 % hold‑out row numbers

# holdout_vec <- rep(FALSE, n)
# holdout_vec[idx] <- TRUE

# ## Fit the models
# model_ac  <- lm_ac(airquality, "Ozone", holdout = holdout_vec)
# model_pre <- lm_prefill(airquality, "Ozone",
#                         impute_method = "mice", m = 5,
#                         holdout = holdout_vec)

# ## Compare
# results <- lm_compare(list(ac = model_ac,
#                            prefill = model_pre),
#                       yName = "Ozone")

# results$RMSE
# results$R2
# results$Beta