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
    1 - sum((y_true[ok] - ybar)^2) / sum((y_true[ok] - ybar)^2)
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

k_folds <- 5
seed    <- 42
set.seed(seed)

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"
cols <- c(
  "mpg", "cylinders", "displacement", "horsepower",
  "weight", "acceleration", "model_year", "origin", "car_name"
)

mpg_raw <- read.table(
  url,
  header = FALSE,
  sep = "",
  col.names = cols,
  na.strings = "?",
  stringsAsFactors = FALSE,
  strip.white = TRUE
)

mpg_raw$horsepower <- as.numeric(mpg_raw$horsepower)

df <- subset(mpg_raw, select = -car_name)

for (nm in names(df)) {
  if (is.character(df[[nm]])) df[[nm]] <- factor(df[[nm]])
}

stopifnot(is.numeric(df$mpg))

cat("\n================ Auto MPG ================\n")

# ---- CC ----------------------------------------------------
cat("\n[CC] Complete Cases (k-fold CV)\n")
cc_res <- kfold_regression(
  data   = df,
  yName  = "mpg",
  k      = k_folds,
  method = "CC",
  seed   = seed
)
cat(sprintf("  Folds : %d\n", k_folds))
cat(sprintf("  MSE   : %.7f\n", cc_res$MSE_mean))
cat(sprintf("  RMSE  : %.7f\n", cc_res$RMSE_mean))
cat(sprintf("  MAE   : %.7f\n", cc_res$MAE_mean))
cat(sprintf("  R²    : %.7f\n", cc_res$R2_mean))

# ---- AC ----------------------------------------------------
cat("\n[AC] Available Cases (k-fold CV, lm_ac)\n")
ac_res <- kfold_regression(
  data   = df,
  yName  = "mpg",
  k      = k_folds,
  method = "AC",
  seed   = seed
)
cat(sprintf("  Folds : %d\n", k_folds))
cat(sprintf("  MSE   : %.7f\n", ac_res$MSE_mean))
cat(sprintf("  RMSE  : %.7f\n", ac_res$RMSE_mean))
cat(sprintf("  MAE   : %.7f\n", ac_res$MAE_mean))
cat(sprintf("  R²    : %.7f\n", ac_res$R2_mean))

# ---- Missingness row --------------------------------------
miss <- format_missingness(df)
cat("\n% Average Missingness (All Data)\n")
cat("  ", miss, "\n", sep = "")

