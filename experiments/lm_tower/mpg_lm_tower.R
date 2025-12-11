# Load necessary libraries
library(mvPred)

# Function to calculate missingness percentage per column
missing_percent <- function(df) {
  sapply(df, function(col) {
    sum(is.na(col)) / length(col) * 100
  })
}

# Function to compute metrics
compute_metrics <- function(true, pred) {
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


# Load data

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"

colnames <- c(
  "mpg",
  "cylinders",
  "displacement",
  "horsepower",
  "weight",
  "acceleration",
  "model_year",
  "origin",
  "car_name"
)
mpg <- read.csv(url, header = FALSE, col.names = colnames,
                na.strings = "?", strip.white = TRUE, sep = "")

# Drop categorical feature
mpg$car_name <- NULL


# Convert to numeric
mpg$cylinders <- as.numeric(as.character(mpg$cylinders))
mpg$origin <- as.numeric(as.character(mpg$origin))
mpg$model_year <- as.numeric(as.character(mpg$model_year))
mpg$horsepower <- as.numeric(as.character(mpg$horsepower))


# K-fold cross-validation
set.seed(123)
k <- 5
n <- nrow(mpg)
fold_ids <- sample(rep(1:k, length.out = n))

# Store all metrics per fold
metric_results <- vector("list", k)
train_missing_list <- vector("list", k)
validation_missing_list <- vector("list", k)
discarded_validation<- vector("list", k)

# Cross-validation loop
for (i in 1:k) {
  cat("Processing fold", i, "of", k, "\n")
  
  # training and validation sets
  train <- mpg[fold_ids != i, ]
  validation <- mpg[fold_ids == i, ]
  
  # Store missingness for this fold
  train_missing_list[[i]] <- missing_percent(train)
  validation_missing_list[[i]] <- missing_percent(validation)
  
  # Fit model
  lm_tower_obj <- lm_tower(train, yName = "mpg")
  
  #lm_prefill_model <- lm_prefill(train, "mpg",impute_method = "complete",use_dummies = FALSE)
  
  #Keep count of records that cannot be used for prediction
  removed <- sum(!complete.cases(validation))
  discarded_validation[[i]] <- removed
  validation_complete <- validation[complete.cases(validation), ]
  
  #Prep predict for tower
  validation_pred <- validation_complete
  validation_pred$mpg <- NULL
  
  # Predict
  predictions <- predict(lm_tower_obj, newdata = validation_pred)
  
  # Filter true_values to match: only keep complete cases from validation

  true_values <- validation_complete$mpg
  
  # Store all metrics
  metric_results[[i]] <- compute_metrics(true_values, predictions)
}


# Summarize final metrics
metrics_df <- do.call(rbind, lapply(metric_results, as.data.frame))
cat("\nAverage Metrics:\n")
print(colMeans(metrics_df, na.rm = TRUE))

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