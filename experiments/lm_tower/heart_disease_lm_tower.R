# Load necessary libraries
library(mvPred)
library(dplyr)
# Function to calculate missingness percentage per column
missing_percent <- function(df) {
  sapply(df, function(col) {
    sum(is.na(col)) / length(col) * 100
  })
}

# ---------------------------
# Function to compute metrics
# ---------------------------
compute_metrics <- function(true, pred_class, pred_prob) {
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

# ---------------------------
# Load Heart Disease dataset
# ---------------------------
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"
colnames <- c(
  "age", "sex", "cp", "trestbps", "chol",
  "fbs", "restecg", "thalach", "exang",
  "oldpeak", "slope", "ca", "thal", "target"
)

heart <- read.table(url, sep = ",", header = FALSE, col.names = colnames, na.strings = "?")

# Convert target to binary: 0 = no disease, 1 = disease
heart$target <- ifelse(heart$target == 0, 0, 1)

# Convert categorical columns to factors
categorical_cols <- c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal")
heart[categorical_cols] <- lapply(heart[categorical_cols], factor)

# ---------------------------
# K-fold cross-validation
# ---------------------------
set.seed(123)
k <- 5
n <- nrow(heart)
fold_ids <- sample(rep(1:k, length.out = n))

metric_results <- vector("list", k)
train_missing_list <- vector("list", k)
validation_missing_list <- vector("list", k)
discarded_validation<- vector("list", k)

for (i in 1:k) {
  cat("Processing fold", i, "of", k, "\n")
  
  train <- heart[fold_ids != i, ]
  validation <- heart[fold_ids == i, ]
  
  # Handle rare categorical levels
  for (col in categorical_cols) {
    freq <- table(train[[col]])
    rare_levels <- names(freq[freq < 5])  # Adjust threshold if needed
    train[[col]] <- as.character(train[[col]])
    train[[col]][train[[col]] %in% rare_levels] <- "Other"
    train[[col]] <- factor(train[[col]])
    
    validation[[col]] <- as.character(validation[[col]])
    validation[[col]][!validation[[col]] %in% levels(train[[col]])] <- "Other"
    validation[[col]] <- factor(validation[[col]], levels = levels(train[[col]]))
  }
  
  
  #Keep count of records that cannot be used for prediction
  removed <- sum(!complete.cases(validation))
  discarded_validation[[i]] <- removed
  
  
  # Keep complete cases for validation
  validation_complete <- validation[complete.cases(validation), ]
  
  # Store missingness for this fold
  train_missing_list[[i]] <- missing_percent(train)
  validation_missing_list[[i]] <- missing_percent(validation)
  
  cat("Missing (%) in training set for fold", i, ":\n")
  print(train_missing_list[[i]])
  
  cat("Missing (%) in validation set for fold", i, ":\n")
  print(validation_missing_list[[i]])
  
  # Fit model
  lm_tower_obj <- lm_tower(train, yName = "target")
  
  # Fit model using lm_prefill
  #lm_prefill_model <- lm_prefill(train, "target", impute_method = "mice", method = "pmm", m=5, use_dummies = FALSE)
  
  #Predict Prepare
  validation_pred <- validation_complete
  validation_pred$target <- NULL
  
  # Predict
  predictions <- predict(lm_tower_obj, newdata = validation_pred)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  true_labels <- validation_complete$target
  
  # Store metrics
  metric_results[[i]] <- compute_metrics(true_labels, predicted_class, predictions)
}

# ---------------------------
# Summarize final metrics
# ---------------------------
metrics_df <- do.call(rbind, lapply(metric_results, as.data.frame))
cat("\nAverage Metrics:\n")
print(colMeans(metrics_df, na.rm = TRUE))
# ---------------------------
# Average missingness across folds
# ---------------------------
avg_train_missing <- colMeans(do.call(rbind, train_missing_list), na.rm = TRUE)
avg_validation_missing <- colMeans(do.call(rbind, validation_missing_list), na.rm = TRUE)

cat("\nAverage missing (%) in training sets across folds:\n")
print(avg_train_missing)

cat("\nAverage missing (%) in validation sets across folds:\n")
print(avg_validation_missing)

# Convert list to numeric vector
discarded_vec <- unlist(discarded_validation)

# Number of validation rows per fold (they may differ slightly if n not divisible by k)
val_sizes <- table(fold_ids)

# Percentage removed per fold
percent_removed_per_fold <- discarded_vec / as.numeric(val_sizes) * 100

cat("\nPercentage of validation rows removed per fold:\n")
print(percent_removed_per_fold)

# Average percentage removed across all k folds
avg_percent_removed <- mean(percent_removed_per_fold)

cat("\nAverage % of validation records removed across all folds:\n")
print(avg_percent_removed)
