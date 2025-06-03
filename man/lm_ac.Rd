\name{lm_ac}
\alias{lm_ac}
\alias{ac_mul}
\alias{ac_vec}
\alias{summary.lm_ac}
\alias{predict.lm_ac}
\title{Linear Regression Using Available Cases Methods}
\description{
This suite of functions implements linear regression using the Available Cases (pairwise deletion) approach. 
It makes use of the avaliable cases in the dataset to fit the model and make predictions rather than deleting
any rows with missing values. The main function \code{lm_ac} fits the model even when the predictors or response contain missing values.
The helper functions compute matrix products using only complete pairs of observations. S3 methods are included for summary and prediction.
}
\usage{
lm_ac(data, yName, ...)

ac_mul(A)

ac_vec(X, y)

\method{summary}{lm_ac}(object, ...)

\method{predict}{lm_ac}(object, newdata, ...)
}
\arguments{
  \item{data}{A data.frame containing predictors and the response variable.}
  \item{yName}{The name of the predictor column in \code{data}.}
  \item{...}{Additional arguments passed to \code{model.frame} in \code{lm_ac}.}
  \item{A}{A numeric matrix or data frame to be used in \code{ac_mul} to compute ATA.}
  \item{X}{A numeric matrix of predictors for \code{ac_vec}.}
  \item{y}{A numeric vector representing the response variable for \code{ac_vec}.}
  \item{object}{An object of class \code{lm_ac} for S3 methods which is returned from calling lm_ac.}
  \item{newdata}{A data frame with predictor variables for prediction.}
}
\details{
The function \code{lm_ac} fits a linear regression model using pairwise deletion for missing values. 
Instead of removing entire rows with missing values (as in complete-case analysis), this approach computes each component 
of the normal equations (X'X and X'y) using the average of all available pairs.

\code{ac_mul} returns a symmetric matrix of average products between all pairs of columns.

\code{ac_vec} returns a vector of average products between each column of X and y.

\code{summary.lm_ac} prints the estimated coefficients.

\code{predict.lm_ac} generates predictions using the fitted model on new data.
}
\value{
\code{lm_ac} returns an object of class \code{lm_ac}, a list containing:
  \item{data}{Cleaned input data used in modeling.}
  \item{yName}{Response variable name.}
  \item{formula}{The model formula.}
  \item{fit_obj}{A list with components \code{coef} (estimated coefficients) and \code{colnames} (predictor names).}

\code{ac_mul} returns a numeric matrix.

\code{ac_vec} returns a numeric vector.

\code{summary.lm_ac} returns the coefficient vector (invisible).

\code{predict.lm_ac} returns a numeric vector of predictions.
}

\examples{
data(airquality)

# Fit model
mod <- lm_ac(airquality, "Ozone")

# Print coefficients
summary(mod)

# Predict on new data
predict(mod, newdata = airquality)
}
