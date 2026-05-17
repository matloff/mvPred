\name{bootstrap}
\alias{bootstrap}
\title{Cross-Validation for Predictive Modeling with Missing Data}
\description{
Run k-fold cross-validation for regression while applying one of the package's missing-data handling strategies.
Depending on the selected method, the function can fit complete-case, available-case, PREFILL, or TOWER-based models and return fold-level summaries.
}
\usage{
bootstrap(
  data,
  yName,
  k = 5,
  task = c("regression", "classification"),
  method = c("AC", "CC", "PREFILL", "TOWER"),
  threshold = 0.5,
  seed = 42,
  test_split = 0.2,
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
  \item{yName}{The name of the response column in \code{data}.}
  \item{k}{Number of folds used for cross-validation.}
  \item{task}{Whether to run \code{"regression"} or \code{"classification"} metrics.}
  \item{method}{Missing-data handling strategy: \code{"AC"}, \code{"CC"}, \code{"PREFILL"}, or \code{"TOWER"}.}
  \item{threshold}{Classification threshold used when \code{task = "classification"}.}
  \item{seed}{Random seed used for split reproducibility.}
  \item{test_split}{Proportion of the full dataset reserved for the final holdout test set.}
  \item{impute_method}{Imputation method used when \code{method = "PREFILL"}.}
  \item{mice_method}{Optional \pkg{mice} method specification for PREFILL runs.}
  \item{m}{Number of multiple imputations used by applicable PREFILL methods.}
  \item{use_dummies}{Logical indicating whether categorical variables should be converted to dummies for Amelia-based PREFILL runs.}
  \item{tower_regFtnName}{Regression routine passed to the TOWER method.}
  \item{tower_opts}{Optional list of control settings passed to TOWER.}
  \item{tower_scaling}{Optional scaling object passed to TOWER.}
  \item{tower_yesYVal}{Optional flag passed to TOWER.}
  \item{...}{Additional arguments forwarded to the selected fitting routine.}
}
\value{
A list containing fold-level performance summaries and fitted model objects. Depending on the selected task, the returned metrics correspond to regression or classification performance.

Common returned components include:
\itemize{
  \item \code{fold_models}: list of fitted models, one per cross-validation fold.
  \item \code{final_model}: fitted model from the final training split when \code{test_split > 0}.
  \item \code{test_results}: performance on the final holdout test set when \code{test_split > 0}.
  \item \code{training_missingness_by_fold}: missingness tables for each training fold.
  \item \code{training_missingness_average}: average missingness across training folds.
}
}
\details{
The function first optionally creates a holdout split controlled by \code{test_split}. Cross-validation is then run on the training portion using the requested missing-data handling method. For PREFILL runs, the supported imputation methods are \code{mice}, \code{Amelia}, \code{missForest}, and complete-case preprocessing.
}
\examples{
data("auto-mpg", package = "mvPred")

df <- auto_mpg
df$car_name <- NULL

for (nm in names(df)) {
  suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
}

res <- bootstrap(
  data = df,
  yName = "mpg",
  k = 3,
  task = "regression",
  method = "CC"
)

res$RMSE_mean
}
