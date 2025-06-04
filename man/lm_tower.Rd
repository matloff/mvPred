\name{lm_tower}
\alias{lm_tower}
\alias{summary.lm_tower}
\alias{predict.lm_tower}
\title{Fit Linear Models with Missing Data Using the Tower Method}
\description{
  \code{lm_tower} fits a linear regression model using the Tower Method implemented in the \pkg{toweranNA} package.
  This method handles missing data in both training and new datasets without explicit imputation by leveraging regression averaging.
}
\usage{
lm_tower(data, yName, regFtnName = "lm", opts = list(), scaling = NULL, yesYVal = NULL)

summary(object, ...)

predict(object, newdata, ...)
}
\arguments{
  \item{data}{A data frame containing the training data with possible missing values.}
  \item{yName}{A character string specifying the name of the response variable in \code{data}.}
  \item{regFtnName}{A character string specifying the regression function to use (default is \code{"lm"}).}
  \item{opts}{A list of options passed to \code{toweranNA::makeTower()}.}
  \item{scaling}{Optional scaling parameters passed to \code{toweranNA::makeTower()}.}
  \item{yesYVal}{Optional parameter passed to \code{toweranNA::makeTower()}.}
  \item{object}{An object of class \code{"lm_tower"} returned by \code{lm_tower()}.}
  \item{newdata}{A data frame containing new observations for prediction. May contain missing values.}
  \item{\dots}{Additional arguments (currently unused).}
}
\value{
  \code{lm_tower} returns an object of class \code{"lm_tower"} containing the fitted Tower model and related information.

  \code{summary.lm_tower} prints a brief summary of the fitted model, including the response variable, regression function, and number of complete cases used.

  \code{predict.lm_tower} returns a numeric vector of predicted values for \code{newdata}. The Tower method handles missing values internally, so predictions are returned for all rows of \code{newdata}.
}
\details{
The Tower Method implemented in the \pkg{toweranNA} package fits regression models that can handle missing data without explicit imputation. It uses regression averaging over subsets of complete cases to provide predictions even when new data contain missing values.

This function provides a convenient wrapper around \code{toweranNA::makeTower()} and associated prediction methods.

\code{lm_tower} is particularly useful when the primary goal is prediction rather than inference.
}
\examples{
library(toweranNA)

data(airquality)
airquality$Ozone[1:10] <- NA
airquality$Solar.R[5:15] <- NA

# Fit the tower model
lm_tower_obj <- lm_tower(airquality, yName = "Ozone")

# Summarize the model
summary(lm_tower_obj)

# Prepare new data with missing values
newdata <- airquality[20:30, -which(names(airquality) == "Ozone")]
newdata[1, "Solar.R"] <- NA

# Predict on new data
preds <- predict(lm_tower_obj, newdata = newdata)
print(preds)
}
\seealso{
  \code{\link[toweranNA]{makeTower}}, \code{\link[toweranNA]{predict.tower}}
}
\author{
Cynthia Mascarenhas \email{cynmascarenhas@ucdavis.edu}
}
\keyword{regression missing-data prediction}
