# Load necessary libraries
library(mvPred)
# Function to calculate missingness percentage per column
missing_percent <- function(df) {
  sapply(df, function(col) {
    sum(is.na(col)) / length(col) * 100
  })
}
# Function to compute metrics
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
  
  return(list(accuracy = accuracy, precision = precision, recall = recall, f1 = f1))
}

# Load data
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
colnames <- c(
  "age", "workclass", "fnlwgt", "education", "education_num", "marital_status",
  "occupation", "relationship", "race", "sex", "capital_gain", "capital_loss",
  "hours_per_week", "native_country", "income"
)
adult <- read.csv(url, header = FALSE, col.names = colnames,
                  na.strings = "?", strip.white = TRUE)

# Target encoding
adult$income <- factor(adult$income, levels = c("<=50K", ">50K"))
adult$income_num <- ifelse(adult$income == ">50K", 1, 0)
adult$income <- NULL

categorical_cols <- c(
  "workclass", "education", "marital_status", "occupation",
  "relationship", "race", "sex", "native_country"
)

# ---------------------------
# K-fold cross-validation
# ---------------------------
set.seed(123)
k <- 5
n <- nrow(adult)
fold_ids <- sample(rep(1:k, length.out = n))

metric_results <- vector("list", k)
train_missing_list <- vector("list", k)
validation_missing_list <- vector("list", k)
discarded_validation<- vector("list", k)

# ---------------------------
# Cross-validation loop
# ---------------------------
for (i in 1:k) {
  cat("Processing fold", i, "of", k, "\n")
  
  train <- adult[fold_ids != i, ]
  validation <- adult[fold_ids == i, ]
  
  # Handle rare categorical levels consistently
  for (col in categorical_cols) {
    train[[col]] <- as.character(train[[col]])
    validation[[col]] <- as.character(validation[[col]])
    
    freq <- table(train[[col]])
    rare_levels <- names(freq[freq < 10])
    
    train[[col]][train[[col]] %in% rare_levels] <- "Other"
    validation[[col]][validation[[col]] %in% rare_levels] <- "Other"
    
    # Force consistent levels across train and validation
    lvls <- union(unique(train[[col]]), "Other")
    train[[col]] <- factor(train[[col]], levels = lvls)
    validation[[col]] <- factor(validation[[col]], levels = lvls)
  }
  
  #Keep count of records that cannot be used for prediction
  removed <- sum(!complete.cases(validation))
  discarded_validation[[i]] <- removed
  
  # Keep complete cases
  validation_complete <- validation[complete.cases(validation), ]
  
  
  # Store missingness for this fold
  train_missing_list[[i]] <- missing_percent(train)
  validation_missing_list[[i]] <- missing_percent(validation)
  
  cat("Missing (%) in training set for fold", i, ":\n")
  print(train_missing_list[[i]])
  
  cat("Missing (%) in validation set for fold", i, ":\n")
  print(validation_missing_list[[i]])
  
  # Fit model
  lm_prefill_model <- lm_prefill(train, "income_num", 
                                 impute_method = "complete", 
                                 use_dummies = FALSE)
  
  # Predict
  predictions <- predict(lm_prefill_model, newdata = validation_complete)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  true_labels <- validation_complete$income_num
  
  # Store metrics
  metric_results[[i]] <- compute_metrics(true_labels, predicted_class, predictions)
}

# Summarize final metrics
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
