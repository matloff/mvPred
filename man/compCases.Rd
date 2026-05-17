\name{compCases}
\alias{compCases}
\title{Extract Complete Cases from a Data Frame}
\description{
Remove rows containing missing values and report which rows were retained or discarded.
}
\usage{
compCases(data)
}
\arguments{
  \item{data}{A data frame or matrix that may contain missing values.}
}
\value{
A list with components:
\itemize{
  \item \code{intactData}: rows of \code{data} with no missing values.
  \item \code{intactRows}: logical vector indicating complete rows.
  \item \code{nonintactRows}: integer indices of rows containing at least one missing value.
}
}
\examples{
df <- data.frame(x = c(1, NA, 3), y = c(2, 4, 6))
compCases(df)
}
