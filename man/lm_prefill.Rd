\name{lm_prefill}
\alias{lm_prefill}
\alias{summary.lm_prefill}
\alias{predict.lm_prefill}
\title{Linear Model Fitting with Pre-Imputation of Missing Data}
\description{
  Fit linear models to data with missing values by first imputing missing data using multiple imputation methods,
  then fitting linear models to the imputed datasets. Supports multiple imputation methods including \code{mice},
  \code{Amelia}, and \code{missForest}, as well as complete-case analysis.
}
\usage{
lm_prefill(data, yName, impute_method = "mice", m = 5, use_dummies = FALSE, ...)

summary(object, ...)

predict(object, newdata, type = "response", ...)
}
\arguments{
  \item{data}{A data frame containing the variables for modeling, possibly with missing values.}
  \item{yName}{A character string specifying the name of the response variable in \code{data}.}
  \item{impute_method}{Character string specifying the imputation method to use. One of \code{"mice"} (default), \code{"amelia"}, \code{"missforest"}, or \code{"complete"}.}
  \item{m}{Integer specifying the number of multiple imputations to generate (only used for \code{"mice"} and \code{"amelia"} methods). Default is 5.}
  \item{use_dummies}{Logical indicating whether to convert categorical predictors to dummy variables before imputation (only relevant for \code{"amelia"} method). Default is \code{FALSE}.}
  \item{...}{Additional arguments passed to the imputation function (\code{mice::mice}, \code{Amelia::amelia}, or \code{missForest::missForest}) or to \code{stats::lm} for model fitting.}
  \item{object}{An object of class \code{"lm_prefill"} returned by \code{lm_prefill()}.}
  \item{newdata}{A data frame containing new data for prediction. Rows with missing values in \code{newdata} are omitted before prediction.}
  \item{type}{Type of prediction; passed to \code{stats::predict.lm}. Default is \code{"response"}.}
}
\details{
\code{lm_prefill} is a wrapper function that first imputes missing values in the dataset using one of several supported methods, then fits linear models to the imputed datasets.

Supported imputation methods:

\itemize{
  \item \code{"mice"}: Multiple Imputation by Chained Equations via the \pkg{mice} package.
  \item \code{"amelia"}: Multiple imputation using the \pkg{Amelia} package.
  \item \code{"missforest"}: Nonparametric imputation using random forests via the \pkg{missForest} package.
  \item \code{"complete"}: Complete-case analysis by omitting rows with missing values.
}

The \code{use_dummies} argument controls whether categorical variables are converted to dummy variables before imputation (only applicable for the \code{"amelia"} method). If \code{use_dummies = TRUE}, the package \pkg{regtools} is required for dummy variable creation.

Additional arguments passed via \code{...} are forwarded to the underlying imputation or modeling functions as appropriate.

The returned object is of class \code{"lm_prefill"} and contains the imputed data, fitted models, and metadata.

\code{summary.lm_prefill} returns model summaries for each imputed dataset or a single summary for complete-case or missForest methods.

\code{predict.lm_prefill} returns predictions averaged across imputations for multiple imputation methods, or direct predictions for single imputation or complete-case models.

If \code{newdata} contains missing values, those rows are omitted with a warning before prediction.
}
\value{
An object of class \code{"lm_prefill"} containing at least the following components:

\item{data}{Original input data.}
\item{yName}{Response variable name.}
\item{formula}{Model formula used for fitting.}
\item{imputed_data}{Imputed dataset(s) or imputation object.}
\item{impute_method}{Imputation method used.}
\item{fit_obj}{Fitted linear model(s).}
\item{missforest_OOBerror}{(If \code{impute_method = "missforest"}) Out-of-bag error estimate from \pkg{missForest}.}
}
\examples{
library(mice)
library(Amelia)
library(missForest)

# Example dataset with missing values
data(airquality)
airquality$Ozone[1:10] <- NA

# Fit linear model with mice imputation
lm_obj <- lm_prefill(airquality, yName = "Ozone", impute_method = "mice", m = 5)

# Summarize fitted models
summary(lm_obj)

# Predict on new data (complete cases only)
newdata <- airquality[11:20, ]
predict(lm_obj, newdata)

# Fit with Amelia and use dummy variables for categorical predictors
lm_obj2 <- lm_prefill(airquality, yName = "Ozone", impute_method = "amelia", m = 5, use_dummies = TRUE)
summary(lm_obj2)
}
\seealso{
\code{\link[mice]{mice}}, \code{\link[Amelia]{amelia}}, \code{\link[missForest]{missForest}}, \code{\link[stats]{lm}}, \code{\link[stats]{predict.lm}}, \code{\link[stats]{na.omit}}, \code{\link[regtools]{factorsToDummies}}
}
\author{
Cynthia Mascarenhas \email{cynmascarenhas@ucdavis.edu}
}
\keyword{models regression imputation missing}
