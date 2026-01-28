suppressPackageStartupMessages({
  library(mvPred)
  library(toweranNA)
  library(VIM)
  library(qeML)
})

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

df_mpg <- mtcars
y_mpg  <- "mpg"

df_tao <- VIM::tao
y_tao  <- "Sea.Surface.Temp"

df_sleep <- VIM::sleep
y_sleep  <- "BrainWgt"

df_wine <- tryCatch(VIM::wine, error = function(e) NULL)
if (is.null(df_wine)) df_wine <- tryCatch(qeML::wine, error = function(e) NULL)
if (is.null(df_wine)) stop("Could not find 'wine' in VIM or qeML.")

if (!("price" %in% names(df_wine)))  stop("wine: column 'price' not found.")
if (!("points" %in% names(df_wine))) stop("wine: column 'points' not found.")

df_wine_pp <- data.frame(
  price  = df_wine[["price"]],
  points = df_wine[["points"]],
  dummy  = rnorm(nrow(df_wine)) # To be removed later. Added it so I can run tower
)
y_wine <- "price"

# ------------------------------------------------------------
# Communities & Crime
# ------------------------------------------------------------
cc_data  <- "../data/communities.data"
cc_names <- "../data/communities.names"

# Parse column names
nm_lines <- readLines(cc_names, warn = FALSE)
attr_lines <- grep("^@attribute\\s+", nm_lines, value = TRUE)
col_names <- sub("^@attribute\\s+([^ ]+)\\s+.*$", "\\1", attr_lines)

# Read data
df_cc <- read.csv(
  cc_data,
  header = FALSE,
  na.strings = "?",
  stringsAsFactors = FALSE
)

# Assign names
if (ncol(df_cc) == length(col_names)) {
  names(df_cc) <- col_names
}

# Coerce to numeric
for (nm in names(df_cc)) {
  suppressWarnings(df_cc[[nm]] <- as.numeric(df_cc[[nm]]))
}

# Drop ID columns only (KEEP LEMAS)
id_cols <- intersect(
  c("state", "county", "community", "communityname", "fold"),
  names(df_cc)
)
df_commcrime <- df_cc[, setdiff(names(df_cc), id_cols), drop = FALSE]

y_commcrime <- "ViolentCrimesPerPop"

# ------------------------------------------------------------
# english
# ------------------------------------------------------------
load("../data/english.RData")   

df_english2 <- data.frame(
  age   = english[["age"]],
  vocab = english[["vocab"]],
  dummy = rnorm(nrow(english)) # To be removed later. Added it so I can run tower
)
y_english <- "vocab"

# ------------------------------------------------------------
# NHkids
# ------------------------------------------------------------
load("../data/NHkids.RData")  

df_nhkids2 <- data.frame(
  weight = NHkids[["Weight"]],
  height = NHkids[["Height"]],
  dummy  = rnorm(nrow(NHkids)) # To be removed later. Added it so I can run tower
)
y_nhkids <- "Weight"

# -----------------------------
# Run study
# -----------------------------
all_results <- rbind(
  run_cc_ac_tower_prefill(df_mpg,     y_mpg,     "mpg (mtcars)"),
  run_cc_ac_tower_prefill(df_tao,     y_tao,     "tao"),
  run_cc_ac_tower_prefill(df_sleep,   y_sleep,   "sleep"),
  run_cc_ac_tower_prefill(df_wine_pp, y_wine,    "wine (price ~ points only)"),
  run_cc_ac_tower_prefill(df_english2,    y_english,    "english (vocab ~ age)"),
  run_cc_ac_tower_prefill(df_nhkids,      y_nhkids,     "NHkids (weight)"),
  run_cc_ac_tower_prefill(df_commcrime,   y_commcrime,  "Communities & Crime")
)