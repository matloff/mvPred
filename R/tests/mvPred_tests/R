# devtools::install_github("matloff/mvPred", force=TRUE) # Uncomment this to install the package
library(mvPred)

# ============================================
# Load mvPred package
# ============================================
library(mvPred)
library(Amelia)
library(mice)
library(missForest)
library(toweranNA)
library(qeML)
library(regtools)
#install.packages("tufte")


# ============================================
# Load built-in dataset (airquality)
# ============================================
data(airquality)

# Create artificial missingness to test methods
set.seed(42)
airquality$Ozone[sample(1:nrow(airquality), 10)] <- NA
airquality$Solar.R[sample(1:nrow(airquality), 10)] <- NA

# Define holdout for methods that support holdout
n <- nrow(airquality)
holdout_vec <- rep(FALSE, n)
holdout_vec[sample(1:n, size = floor(0.20 * n))] <- TRUE

# ============================================
# Complete Cases (CC)
# ============================================
cc_result <- compCases(airquality)
lm_cc <- lm(Ozone ~ ., data = cc_result$intactData)
summary(lm_cc)

# ============================================
# Available Cases (AC)
# ============================================
model_ac <- lm_ac(airquality, "Ozone", holdout = holdout_vec)
summary(model_ac)

# ============================================
# lm_prefill with mice
# ============================================
model_prefill_mice <- lm_prefill(airquality, "Ozone",
                                 impute_method = "mice",
                                 m = 5, holdout = holdout_vec)
summary(model_prefill_mice)

# ============================================
# lm_prefill with Amelia
# ============================================
model_prefill_amelia <- lm_prefill(airquality, "Ozone",
                                   impute_method = "amelia",
                                   m = 5, holdout = holdout_vec)
summary(model_prefill_amelia)

# ============================================
# lm_prefill with missForest
# ============================================
model_prefill_mf <- lm_prefill(airquality, "Ozone",
                               impute_method = "missforest",
                               holdout = holdout_vec)
summary(model_prefill_mf)

# ============================================
# lm_tower (ToweranNA)
# ============================================
model_tower <- lm_tower(airquality, "Ozone")
summary(model_tower)

# ============================================
# qeMLna method (on a subset for simplicity)
# ============================================
# Subset for qeML example
qe_data <- airquality[, c("Ozone", "Solar.R", "Wind", "Temp")]
qe_data <- qe_data[complete.cases(qe_data), ]  # for clean example

# Use compCases with qeMLna
qe_out <- qeMLna(qe_data, "Ozone", qeMLftn = "qeLin",
                 mvFtn = "compCases",
                 retainMVFtnOut = FALSE)



# Print qeML output
print(qe_out)

# ============================================
# Predict example for AC model
# ============================================
y_pred <- predict(model_ac, newdata = model_ac$testing_data)
y_actual <- model_ac$testing_data$Ozone
rmse <- sqrt(mean((y_pred - y_actual)^2, na.rm = TRUE))

cat("\nPrediction RMSE for AC model:", rmse, "\n")


# ============================================
# lm_compare test
# ============================================
models_to_compare <- list(
  ac = model_ac,
  mice = model_prefill_mice,
  amelia = model_prefill_amelia,
  missforest = model_prefill_mf
)

compare_results <- lm_compare(models = models_to_compare, yName = "Ozone")

# Print out RMSE
cat("\n--- lm_compare RMSE ---\n")
print(compare_results$RMSE)

# Print out R²
cat("\n--- lm_compare R² ---\n")
print(compare_results$R2)

# Print out coefficients
cat("\n--- lm_compare Coefficients ---\n")
print(compare_results$Beta)

