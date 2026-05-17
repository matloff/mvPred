\name{lm_compare}
\alias{lm_compare}
\title{Compare Missing-Data Methods with Bootstrap Summaries}
\description{
Run \code{bootstrap()} across one or more missing-data handling methods and collect the resulting cross-validation summary metrics in a single data frame.
}
\usage{
lm_compare(
  data,
  predictors = NULL,
  yName,
  k = 5,
  task = c("regression", "classification"),
  methods = c("CC", "AC", "TOWER", "PREFILL"),
  seed = 42,
  impute_method = c("mice", "amelia", "missforest", "complete"),
  mice_method = NULL,
  m = 5,
  use_dummies = FALSE,
  tower_regFtnName = "lm",
  tower_opts = list(),
  tower_scaling = NULL,
  tower_yesYVal = NULL,
  ...
)
}
\arguments{
  \item{data}{A data frame containing predictors and the response variable.}
  \item{predictors}{Optional character vector of predictor names to keep. If \code{NULL}, all non-response columns are used.}
  \item{yName}{The name of the response column in \code{data}.}
  \item{k}{Number of folds used for cross-validation.}
  \item{task}{Whether to summarize regression or classification performance.}
  \item{methods}{Character vector of methods to evaluate. Supported values are \code{"CC"}, \code{"AC"}, \code{"TOWER"}, and \code{"PREFILL"}.}
  \item{seed}{Random seed used for reproducibility.}
  \item{impute_method}{Imputation method used when \code{"PREFILL"} is included in \code{methods}.}
  \item{mice_method}{Optional \pkg{mice} method specification passed to PREFILL runs.}
  \item{m}{Number of multiple imputations used by applicable PREFILL methods.}
  \item{use_dummies}{Logical indicating whether categorical variables should be converted to dummies for Amelia-based PREFILL runs.}
  \item{tower_regFtnName}{Regression routine passed to TOWER.}
  \item{tower_opts}{Optional list of TOWER options.}
  \item{tower_scaling}{Optional scaling object passed to TOWER.}
  \item{tower_yesYVal}{Optional flag passed to TOWER.}
  \item{...}{Additional arguments forwarded to \code{bootstrap()}.}
}
\value{
A data frame containing one row per requested method and the corresponding average performance metrics.
}
\details{
\code{lm_compare()} is intended for method comparison under cross-validation. Internally it calls \code{bootstrap()} with \code{test_split = 0}, so the reported summaries correspond to cross-validation results rather than final holdout-test performance.
}
\examples{
data("auto-mpg", package = "mvPred")

df <- auto_mpg
df$car_name <- NULL

for (nm in names(df)) {
  suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
}

lm_compare(df, yName = "mpg", k = 3, task = "regression", methods = c("CC", "AC"))
}
