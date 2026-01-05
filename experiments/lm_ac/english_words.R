library("lm_ac.R")

metrics_regression <- function(y_true, y_pred) {
  ok <- !is.na(y_true) & !is.na(y_pred)
  scored <- sum(ok)
  total  <- length(y_true)
  coverage <- if (total > 0) scored / total else NA_real_
  
  if (scored == 0) {
    return(c(
      MSE      = NA_real_,
      RMSE     = NA_real_,
      MAE      = NA_real_,
      R2       = NA_real_,
      Scored   = scored,
      Total    = total,
      Coverage = coverage
    ))
  }
  
  mse <- mean((y_true[ok] - y_pred[ok])^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(y_true[ok] - y_pred[ok]))
  r2   <- {
    ybar <- mean(y_true[ok])
    1 - sum((y_true[ok] - y_pred[ok])^2) / sum((y_true[ok] - ybar)^2)
  }
  
  c(
    MSE      = mse,
    RMSE     = rmse,
    MAE      = mae,
    R2       = r2,
    Scored   = scored,
    Total    = total,
    Coverage = coverage
  )
}

format_missingness <- function(df) {
  miss <- colMeans(is.na(df)) * 100
  miss <- miss[miss > 0]
  if (!length(miss)) return("")
  paste(sprintf("%s(%.6f)", names(miss), miss), collapse = "; ")
}

# Prevent "new factor levels" in CV folds
collapse_rare_levels <- function(df, cols, min_count = 10, other_label = "Other") {
  for (col in cols) {
    if (!(col %in% names(df))) next
    
    x <- df[[col]]
    
    # Work in character to safely relabel levels
    if (is.factor(x)) x <- as.character(x)
    
    # If it's not character now, skip
    if (!is.character(x)) next
    
    freq <- table(x, useNA = "no")
    rare <- names(freq[freq < min_count])
    
    x[x %in% rare] <- other_label
    df[[col]] <- factor(x)
  }
  df
}

k_folds <- 5
seed    <- 42
set.seed(seed)

cat("\n================ English (vocab) ================\n")

english_path <- "English.csv"
if (!file.exists(english_path)) english_path <- "/cloud/project/English.csv"

english <- read.csv(
  english_path,
  na.strings = c("", "NA", "?")
)

english <- english[, !names(english) %in% c("demo", "demo_label", "n", "measure", "language", "form")]

stopifnot("vocab" %in% names(english))
english$vocab <- suppressWarnings(as.numeric(english$vocab))
if (all(is.na(english$vocab))) stop("'vocab' could not be coerced to numeric.")

# Categorical columns
categorical_cols <- c("birth_order", "ethnicity", "sex", "mom_ed")

# Ensure they are factors
for (col in intersect(categorical_cols, names(english))) {
  english[[col]] <- as.factor(english[[col]])
}

# Convert any remaining character predictors to factors
for (nm in names(english)) {
  if (is.character(english[[nm]])) {
    english[[nm]] <- factor(english[[nm]])
  }
}

# Collapse rare levels globally to avoid fold-specific "new levels" errors
english <- collapse_rare_levels(
  df          = english,
  cols        = categorical_cols,
  min_count   = 10,
  other_label = "Other"
)

# CC: Complete Cases 
cat("\n[CC] Complete Cases (k-fold CV)\n")
cc_res <- kfold_regression(
  data   = english,
  yName  = "vocab",
  k      = k_folds,
  method = "CC",
  seed   = seed
)

cat(sprintf("  Folds : %d\n", k_folds))
cat(sprintf("  MSE   : %.7f\n", cc_res$MSE_mean))
cat(sprintf("  RMSE  : %.7f\n", cc_res$RMSE_mean))
cat(sprintf("  MAE   : %.7f\n", cc_res$MAE_mean))
cat(sprintf("  R²    : %.7f\n", cc_res$R2_mean))

# AC: Available Cases
cat("\n[AC] Available Cases (k-fold CV, lm_ac)\n")
ac_res <- kfold_regression(
  data   = english,
  yName  = "vocab",
  k      = k_folds,
  method = "AC",
  seed   = seed
)

cat(sprintf("  Folds : %d\n", k_folds))
cat(sprintf("  MSE   : %.7f\n", ac_res$MSE_mean))
cat(sprintf("  RMSE  : %.7f\n", ac_res$RMSE_mean))
cat(sprintf("  MAE   : %.7f\n", ac_res$MAE_mean))
cat(sprintf("  R²    : %.7f\n", ac_res$R2_mean))

# Missingness row
miss <- format_missingness(english)
cat("\n% Average Missingness (All Data)\n")
cat("  ", miss, "\n", sep = "")

