suppressPackageStartupMessages({
  library(mvPred)
  library(toweranNA)
  library(VIM)
  library(qeML)
})

source("../R/lm_AC.R")
source("../R/utils.R")
set.seed(42)

# -----------------------------
# Helper: run CC vs AC vs TOWER vs PREFILL(mice)
# -----------------------------
run_cc_ac_tower_prefill <- function(
    df,
    yName,
    dataset_name,
    k = 5,
    tower_regFtnName = "lm",
    tower_opts = list(),
    tower_scaling = NULL,
    tower_yesYVal = NULL
) {
  # if (!is.data.frame(df)) df <- as.data.frame(df)
  df <- as.data.frame(df)

  if (!(yName %in% names(df))) stop(sprintf("[%s] yName '%s' not found.", dataset_name, yName))
  
  # Drop rows with missing y for ALL methods
  df <- df[!is.na(df[[yName]]), , drop = FALSE]
  
  if (nrow(df) < k) {
    return(data.frame(
      dataset = dataset_name,
      yName = yName,
      method = c("CC", "AC", "TOWER", "PREFILL_mice"),
      MSE_mean = NA_real_,
      RMSE_mean = NA_real_,
      MAE_mean = NA_real_,
      R2_mean = NA_real_,
      n = nrow(df),
      stringsAsFactors = FALSE
    ))
  }
  
  cat("\n=====================================\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Target :", yName, "\n")
  cat("Rows   :", nrow(df), " Cols:", ncol(df), "\n")
  cat("=====================================\n")
  
  # ---- CC ----
  res_cc <- bootstrap(
    data = df,
    yName = yName,
    k = k,
    task = "regression",
    method = "CC"
  )
  
  # ---- AC ----
  res_ac <- bootstrap(
    data = df,
    yName = yName,
    k = k,
    task = "regression",
    method = "AC"
  )
  
  # ---- TOWER ----
  res_tw <- bootstrap(
    data = df,
    yName = yName,
    k = k,
    task = "regression",
    method = "TOWER",
    tower_regFtnName = tower_regFtnName,
    tower_opts       = tower_opts,
    tower_scaling    = tower_scaling,
    tower_yesYVal    = tower_yesYVal
  )
  
  # ---- PREFILL (mice only) ----
  res_pf_mice <- bootstrap(
    data = df,
    yName = yName,
    k = k,
    task = "regression",
    method = "PREFILL",
    impute_method = "mice",
    m = 5,
    use_dummies = FALSE
  )
  
  out <- rbind(
    data.frame(
      dataset = dataset_name, yName = yName, method = "CC",
      MSE_mean = res_cc$MSE_mean, RMSE_mean = res_cc$RMSE_mean,
      MAE_mean = res_cc$MAE_mean, R2_mean = res_cc$R2_mean,
      n = nrow(df), stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = dataset_name, yName = yName, method = "AC",
      MSE_mean = res_ac$MSE_mean, RMSE_mean = res_ac$RMSE_mean,
      MAE_mean = res_ac$MAE_mean, R2_mean = res_ac$R2_mean,
      n = nrow(df), stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = dataset_name, yName = yName, method = "TOWER",
      MSE_mean = res_tw$MSE_mean, RMSE_mean = res_tw$RMSE_mean,
      MAE_mean = res_tw$MAE_mean, R2_mean = res_tw$R2_mean,
      n = nrow(df), stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = dataset_name, yName = yName, method = "PREFILL_mice",
      MSE_mean = res_pf_mice$MSE_mean, RMSE_mean = res_pf_mice$RMSE_mean,
      MAE_mean = res_pf_mice$MAE_mean, R2_mean = res_pf_mice$R2_mean,
      n = nrow(df), stringsAsFactors = FALSE
    )
  )
  
  print(out, row.names = FALSE)
  invisible(out)
}

# -----------------------------
# Load datasets + targets
# -----------------------------

# 1) mpg (mtcars)
auto_path <- "../data/auto-mpg.data"

auto_cols <- c(
  "mpg", "cylinders", "displacement", "horsepower",
  "weight", "acceleration", "model_year", "origin", "car_name"
)

df_auto <- read.table(
  auto_path,
  col.names = auto_cols,
  na.strings = "?",
  stringsAsFactors = FALSE
)

# Drop non-numeric ID column
df_auto$car_name <- NULL

# Coerce to numeric (horsepower has NAs)
for (nm in names(df_auto)) {
  suppressWarnings(df_auto[[nm]] <- as.numeric(df_auto[[nm]]))
}

y_mpg  <- "mpg"

# 2) tao
df_tao <- VIM::tao
y_tao  <- "Sea.Surface.Temp"

# 3) sleep
df_sleep <- VIM::sleep
y_sleep  <- "BrainWgt"

# 4) wine: ONLY points predictor, price target
df_wine <- tryCatch(VIM::wine, error = function(e) NULL)
if (is.null(df_wine)) df_wine <- tryCatch(qeML::wine, error = function(e) NULL)
if (is.null(df_wine)) stop("Could not find 'wine' in VIM or qeML.")

need_cols <- c("price", "points", "taster_twitter_handle")
miss_cols <- setdiff(need_cols, names(df_wine))
if (length(miss_cols)) stop("wine: missing columns: ", paste(miss_cols, collapse = ", "))

df_wine_pp <- data.frame(
  price = df_wine[["price"]],
  points = df_wine[["points"]],
  taster_twitter_handle = df_wine[["taster_twitter_handle"]]
)

# Ensure handle is treated as a factor (important for model.matrix / dummies)
df_wine_pp$taster_twitter_handle <- as.factor(df_wine_pp$taster_twitter_handle)

y_wine <- "price"

# ------------------------------------------------------------
# english
# ------------------------------------------------------------
load("../data/english.RData")   # loads object: english

df_english2 <- data.frame(
  age   = english[["age"]],
  birth_order = english[["birth_order"]],
  ethnicity = english[["ethnicity"]],
  sex = english[["sex"]],
  mom_ed = english[["mom_ed"]],
  vocab = english[["vocab"]]
)
y_english <- "vocab"

# ------------------------------------------------------------
# NHkids
# ------------------------------------------------------------
load("../data/NHkids.RData")   # loads object: NHkids

df_nhkids <- data.frame(
  AgeMonths = NHkids[["AgeMonths"]],
  weight    = NHkids[["Weight"]],
  height    = NHkids[["Height"]]
)
y_nhkids <- "weight"

# -----------------------------
# Run study
# -----------------------------
all_results <- rbind(
  run_cc_ac_tower_prefill(df_auto,     y_mpg,     "mpg"),
  run_cc_ac_tower_prefill(df_tao,     y_tao,     "tao"),
  run_cc_ac_tower_prefill(df_sleep,   y_sleep,   "sleep"),
  run_cc_ac_tower_prefill(df_wine_pp, y_wine,    "wine"),
  run_cc_ac_tower_prefill(df_english2,    y_english,    "english"),
  run_cc_ac_tower_prefill(df_nhkids,      y_nhkids,     "NHkids")
)


cat("\n\n====================\n")
cat("FINAL SUMMARY\n")
cat("====================\n")
print(all_results, row.names = FALSE)

# Deltas vs CC
cc_base <- subset(all_results, method == "CC", select = c("dataset","yName","n","RMSE_mean","MAE_mean","R2_mean"))
names(cc_base)[names(cc_base) %in% c("RMSE_mean","MAE_mean","R2_mean")] <- c("RMSE_cc","MAE_cc","R2_cc")

merged <- merge(all_results, cc_base, by = c("dataset","yName","n"), all.x = TRUE)
merged$RMSE_delta_vs_CC <- merged$RMSE_mean - merged$RMSE_cc
merged$MAE_delta_vs_CC  <- merged$MAE_mean  - merged$MAE_cc
merged$R2_delta_vs_CC   <- merged$R2_mean   - merged$R2_cc

cat("\n\n====================\n")
cat("DELTAS vs CC (negative RMSE/MAE is better; positive R2 is better)\n")
cat("====================\n")
print(merged[, c("dataset","yName","method","n","RMSE_delta_vs_CC","MAE_delta_vs_CC","R2_delta_vs_CC")],
      row.names = FALSE)