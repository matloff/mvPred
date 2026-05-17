\name{auto_mpg}
\alias{auto_mpg}
\alias{english}
\alias{NHkids}
\docType{data}
\title{Datasets Included in mvPred}
\description{
Datasets distributed with \pkg{mvPred} for examples, experiments, and package demonstrations.
}
\usage{
data("auto-mpg")
data(english)
data(NHkids)
}
\format{
\code{auto_mpg}: a data frame based on the Auto MPG data, including fuel-efficiency and vehicle-characteristic variables.

\code{english}: a data frame used for predictive modeling examples involving vocabulary outcomes and demographic predictors.

\code{NHkids}: a data frame used for predictive modeling examples involving child weight, height, and age in months.
}
\details{
These datasets are included to support examples and methodological comparisons within the package.
}
\examples{
data("auto-mpg", package = "mvPred")
data("english", package = "mvPred")
data("NHkids", package = "mvPred")
}
