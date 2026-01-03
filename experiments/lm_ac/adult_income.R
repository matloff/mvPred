library("lm_ac.R")

compute_auc <- function(y_true, scores) {
  y_true <- as.integer(y_true)
  pos <- which(y_true == 1L)
  neg <- which(y_true == 0L)
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  r <- rank(scores, ties.method = "average")
  auc <- (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
  as.numeric(auc)
}

classification_report <- function(
    y_true,
    y_score,
    threshold = 0.5,
    header = NULL,
    quiet = FALSE
) {
  y_true <- as.integer(y_true)
  y_pred <- as.integer(y_score >= threshold)
  
  tp <- sum(y_pred == 1 & y_true == 1)
  tn <- sum(y_pred == 0 & y_true == 0)
  fp <- sum(y_pred == 1 & y_true == 0)
  fn <- sum(y_pred == 0 & y_true == 1)
  
  acc <- mean(y_pred == y_true)
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
    2 * prec * rec / (prec + rec)
  } else NA_real_
  auc  <- compute_auc(y_true, y_score)
  
  if (!quiet) {
    if (!is.null(header)) cat(header, "\n", sep = "")
    cat(sprintf("  Accuracy : %.6f\n", acc))
    cat(sprintf("  Precision: %.6f\n", prec))
    cat(sprintf("  Recall   : %.6f\n", rec))
    cat(sprintf("  F1 Score : %.6f\n", f1))
    cat(sprintf("  AUC      : %.6f\n", auc))
  }
  
  invisible(list(
    accuracy  = acc,
    precision = prec,
    recall    = rec,
    f1        = f1,
    auc       = auc
  ))
}

format_missingness <- function(df) {
  miss <- colMeans(is.na(df)) * 100
  miss <- miss[miss > 0]
  if (!length(miss)) return("")
  paste(sprintf("%s(%.6f)", names(miss), miss), collapse = "; ")
}

k_folds   <- 5
threshold <- 0.5
seed      <- 42
set.seed(seed)

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
cols <- c(
  "age","workclass","fnlwgt","education","education_num","marital_status",
  "occupation","relationship","race","sex","capital_gain","capital_loss",
  "hours_per_week","native_country","income"
)

adult <- read.table(
  url,
  sep = ",",
  header = FALSE,
  col.names = cols,
  na.strings = "?",
  strip.white = TRUE
)

adult$income_num <- ifelse(adult$income == ">50K", 1L, 0L)

adult$native_country <- trimws(adult$native_country)
adult$native_country[adult$native_country == "Holand-Netherlands"] <- NA

df <- subset(adult, select = -c(income, education, fnlwgt))

for (nm in names(df)) {
  if (is.character(df[[nm]])) df[[nm]] <- factor(df[[nm]])
}

cat("\n================ Adult Income ================\n")

# ---- CC ----------------------------------------------------
cat("\n[CC] Complete Cases (k-fold CV)\n")
cc_res <- bootstrap(
  data      = df,
  yName     = "income_num",
  k         = k_folds,
  threshold = threshold,
  method    = "CC",
  seed      = seed
)
cat(sprintf("  Folds     : %d\n", k_folds))
cat(sprintf("  Accuracy  : %.7f\n", cc_res$accuracy_mean))
cat(sprintf("  Precision : %.7f\n", cc_res$precision_mean))
cat(sprintf("  Recall    : %.7f\n", cc_res$recall_mean))
cat(sprintf("  F1        : %.7f\n", cc_res$f1_mean))
cat(sprintf("  AUC       : %.7f\n", cc_res$auc_mean))

# ---- AC ----------------------------------------------------
cat("\n[AC] Available Cases (k-fold CV, lm_ac)\n")
ac_res <- bootstrap(
  data      = df,
  yName     = "income_num",
  k         = k_folds,
  threshold = threshold,
  method    = "AC",
  seed      = seed
)
cat(sprintf("  Folds     : %d\n", k_folds))
cat(sprintf("  Accuracy  : %.7f\n", ac_res$accuracy_mean))
cat(sprintf("  Precision : %.7f\n", ac_res$precision_mean))
cat(sprintf("  Recall    : %.7f\n", ac_res$recall_mean))
cat(sprintf("  F1        : %.7f\n", ac_res$f1_mean))
cat(sprintf("  AUC       : %.7f\n", ac_res$auc_mean))

# ---- Missingness row --------------------------------------
miss <- format_missingness(df)
cat("\n% Average Missingness (All Data)\n")
cat("  ", miss, "\n", sep = "")
