\name{qeMLna}
\alias{qeMLna}
\alias{addListElement}
\title{Interface qeML Models with Missing-Value Preprocessing}
\description{
Apply a missing-value handling strategy before fitting a predictive model from the \pkg{qeML} package.
The helper \code{addListElement()} appends or replaces an element in a list.
}
\usage{
qeMLna(
  data,
  yName,
  qeMLftn,
  mvFtn,
  qeMLopts = NULL,
  mvPredOpts = NULL,
  retainMVFtnOut = TRUE,
  seed = 9999,
  holdout = 1000
)

addListElement(l, newEltName, newElt)
}
\arguments{
  \item{data}{A data frame containing predictors and the response variable.}
  \item{yName}{The name of the response column in \code{data}.}
  \item{qeMLftn}{Character string naming the predictive function from \pkg{qeML} to call.}
  \item{mvFtn}{Character string naming the missing-value handling routine to apply first.}
  \item{qeMLopts}{Optional list of arguments passed to the selected \pkg{qeML} function.}
  \item{mvPredOpts}{Optional list of arguments passed to the missing-value handling routine.}
  \item{retainMVFtnOut}{Logical indicating whether to store the missing-value routine output in the returned object.}
  \item{seed}{Random seed passed through to applicable routines.}
  \item{holdout}{Holdout size used by the predictive routine if applicable.}
  \item{l}{A list to modify. If \code{NULL}, a new list is created.}
  \item{newEltName}{Name of the element to add or replace.}
  \item{newElt}{Value of the element to add or replace.}
}
\value{
\code{qeMLna()} returns an object of class \code{"MVqeMLout"} containing the \pkg{qeML} fit and, optionally, the missing-value handling output.

\code{addListElement()} returns the modified list.
}
\examples{
\dontrun{
data("english", package = "mvPred")

eng <- english[, c("age", "birth_order", "ethnicity", "sex", "mom_ed", "vocab")]
qeMLna(eng, "vocab", "qeLin", "compCases", retainMVFtnOut = FALSE)
}

addListElement(list(a = 1), "b", 2)
}
