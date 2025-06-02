# =============================
# Source lm_AC Implementation
# =============================
source("R/lm_AC.R")  # Adjust path if needed

# =============================
# Required Libraries
# =============================
library(regtools)
library(tidyr)
library(mice)
library(qeML)

# =============================
# MSE and Reporting Function
# =============================
report_test <- function(dataset_name, df, yName) {
  cat("============================================================\n")
  cat(sprintf("Testing Dataset: %s\n", dataset_name))
  cat("============================================================\n")

  # Fit standard lm using complete cases
  formula <- as.formula(paste(yName, "~ ."))
  lm_cc <- lm(formula, data = df, na.action = na.omit)

  # Fit lm_AC
  mod_ac <- lm_AC(df, yName)

  # Get effective data size
  mf_ac <- model.frame(mod_ac$formula, data = mod_ac$data, na.action = NULL)
  y_ac <- model.response(mf_ac)
  X_ac <- model.matrix(mod_ac$formula, mf_ac)
  intact_pairs <- sapply(1:ncol(X_ac), function(i) sum(!is.na(X_ac[, i]) & !is.na(y_ac)))
  cc_rows_used <- sum(complete.cases(df))
  ac_rows_effective <- mean(intact_pairs)

  cat(sprintf("[INFO] %s: CC used %d rows | AC effective rows ~%.1f\n",
              dataset_name, cc_rows_used, ac_rows_effective))

  # Evaluate on complete rows
  df_eval <- df[complete.cases(df), ]
  y_true <- df_eval[[yName]]
  mse_cc <- mean((predict(lm_cc, newdata = df_eval) - y_true)^2)
  mse_ac <- mean((predict(mod_ac, newdata = df_eval) - y_true)^2)

  cat(sprintf("[RESULT] MSE (lm): %.3f | MSE (lm_AC): %.3f\n\n", mse_cc, mse_ac))
}

# =============================
# Dataset Evaluations
# =============================

# airquality
data(airquality)
report_test("airquality", airquality[!is.na(airquality$Ozone), ], "Ozone")

# mtcars
data(mtcars)
report_test("mtcars", mtcars, "mpg")

# sleep (convert to wide format)
data(sleep)
df_sleep <- pivot_wider(sleep, names_from = group, values_from = extra, names_prefix = "group")
df_sleep <- df_sleep[!is.na(df_sleep$group2), c("group1", "group2", "ID")]
df_sleep$ID <- as.numeric(df_sleep$ID)
report_test("sleep", df_sleep, "group2")

# women
data(women)
report_test("women", women, "weight")

# longley
data(longley)
report_test("longley", longley, "Employed")

# english (qeML)
load(system.file("data", "english.rda", package = "qeML"))  # loads `english`
report_test("english (qeML)", english, "grade3")

# NHISlarge (regtools)
load(system.file("data", "NHISlarge.RData", package = "regtools"))  # loads `NHISlarge`
report_test("NHISlarge (regtools)", NHISlarge, "age_p")

# boys (mice)
data(boys)
report_test("boys (mice)", boys, "hgt")

# nhanes (mice)
data(nhanes)
report_test("nhanes (mice)", nhanes, "bmi")

# mammalsleep (mice)
data(mammalsleep)
report_test("mammalsleep (mice)", mammalsleep, "dream")