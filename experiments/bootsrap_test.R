suppressPackageStartupMessages({
  library(mvPred)
  library(Amelia)
  library(mice)
  library(missForest)
  library(toweranNA)
})

source("../R/utils.R")

set.seed(42)

# ------------------------------------------------------------
# 1) Adult Income (Classification)
# ------------------------------------------------------------
cat("\n==============================\n")
cat("Loading Adult dataset...\n")
cat("==============================\n")

url  <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
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

# Binary target
adult$income_num <- ifelse(trimws(adult$income) == ">50K", 1L, 0L)

# -----------------------------
# PREPROCESS like adult_income.R
# -----------------------------
# 1) Handle rare unseen level (becomes NA; may be dropped at scoring time)
adult$native_country[adult$native_country == "Holand-Netherlands"] <- NA

# 2) Drop columns like you did in adult_income.R
#    (avoid Amelia collinearity issues + unnecessary columns)
df_adult <- subset(adult, select = -c(income, education, fnlwgt))

# If you want to also drop native_country, uncomment:
# df_adult <- subset(adult, select = -c(income, education, fnlwgt, native_country))

cat("\nAdult rows:", nrow(df_adult), " | cols:", ncol(df_adult), "\n")

# -----------------------------
# Adult: CC (Classification)
# -----------------------------
cat("\n==============================\n")
cat("Adult Income | CC | Classification\n")
cat("==============================\n")
res_adult_cc <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "CC",
  threshold = 0.5
)
print(res_adult_cc)

# -----------------------------
# Adult: AC (Classification) 
# -----------------------------
cat("\n==============================\n")
cat("Adult Income | AC | Classification\n")
cat("==============================\n")
res_adult_ac <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "AC",
  threshold = 0.5
)
print(res_adult_ac)

# -----------------------------
# Adult: PREFILL (Classification)
# -----------------------------
cat("\n==============================\n")
cat("Adult Income | PREFILL | missforest | Classification (lm scores for now)\n")
cat("==============================\n")
res_adult_pf_mf <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "PREFILL",
  impute_method = "missforest",    
  m = 5,
  use_dummies = FALSE,
  threshold = 0.5
)
print(res_adult_pf_mf)

cat("\n==============================\n")
cat("Adult Income | PREFILL | amelia | Classification (lm scores for now)\n")
cat("==============================\n")
res_adult_pf_amelia <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "PREFILL",
  impute_method = "amelia",
  m = 5,
  use_dummies = FALSE,
  threshold = 0.5
)
print(res_adult_pf_amelia)

cat("\n==============================\n")
cat("Adult Income | PREFILL | mice | Classification (lm scores for now)\n")
cat("==============================\n")
res_adult_pf_mice <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "PREFILL",
  impute_method = "mice",
  mice_method = "pmm",
  m = 5,
  use_dummies = FALSE,
  threshold = 0.5
)
print(res_adult_pf_mice)

cat("\n==============================\n")
cat("Adult Income | PREFILL | complete | Classification (lm scores for now)\n")
cat("==============================\n")
res_adult_pf_complete <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "PREFILL",
  impute_method = "complete",
  m = 1,                 # complete doesn't need m, but keep scalar int
  use_dummies = FALSE,
  threshold = 0.5
)
print(res_adult_pf_complete)

# -----------------------------
# Adult: TOWER (Classification) 
# -----------------------------
cat("\n==============================\n")
cat("Adult Income | TOWER | Classification (lm scores for now)\n")
cat("==============================\n")
res_adult_tower <- bootstrap(
  df_adult,
  yName = "income_num",
  k = 5,
  task = "classification",
  method = "TOWER",
  threshold = 0.5,
  
  # tower args
  tower_regFtnName = "lm",
  tower_opts       = list(),
  tower_scaling    = NULL,
  tower_yesYVal    = NULL
)
print(res_adult_tower)

# ------------------------------------------------------------
# 2) mtcars mpg (Regression)
# ------------------------------------------------------------
df_mtcars <- mtcars

cat("\n==============================\n")
cat("mtcars | CC | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_cc <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "CC"
)
print(res_mtcars_cc)

cat("\n==============================\n")
cat("mtcars | AC | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_ac <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "AC"
)
print(res_mtcars_ac)

cat("\n==============================\n")
cat("mtcars | PREFILL | complete | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_pf_complete <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "PREFILL",
  impute_method = "complete",
  m = 1,                 # scalar int
  use_dummies = FALSE
)
print(res_mtcars_pf_complete)

cat("\n==============================\n")
cat("mtcars | PREFILL | missforest | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_pf_mf <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "PREFILL",
  impute_method = "missforest",
  m = 5,
  use_dummies = FALSE
)
print(res_mtcars_pf_mf)

cat("\n==============================\n")
cat("mtcars | PREFILL | mice | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_pf_mice <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "PREFILL",
  impute_method = "mice",
  mice_method = "pmm",  
  m = 5,
  use_dummies = FALSE
)
print(res_mtcars_pf_mice)

cat("\n==============================\n")
cat("mtcars | PREFILL | amelia | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_pf_amelia <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "PREFILL",
  impute_method = "amelia",
  m = 5,
  use_dummies = FALSE
)
print(res_mtcars_pf_amelia)

# -----------------------------
# mtcars: TOWER (Regression)
# -----------------------------
cat("\n==============================\n")
cat("mtcars | TOWER | Regression (mpg)\n")
cat("==============================\n")
res_mtcars_tower <- bootstrap(
  df_mtcars,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "TOWER",
  
  # tower args
  tower_regFtnName = "lm",
  tower_opts       = list(),
  tower_scaling    = NULL,
  tower_yesYVal    = NULL
)
print(res_mtcars_tower)
