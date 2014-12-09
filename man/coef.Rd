\name{coef-methods}
\docType{methods}
\alias{coef,MRSP-method}
\title{Methods for Function \code{coef} in Package \pkg{MRSP}}
\description{
Extract the regression coefficients from an object of class \code{MRSP}.
}

\usage{
\S4method{coef}{MRSP}(object, type=c("original", "standardized", "prethreshold"), simplify=TRUE, ...)
}
\arguments{
  \item{object}{Object of class \code{MRSP}.}
  \item{type}{Character string specifying the type of coefficients to extract from \code{object}.}
  \item{simplify}{Logical: should the output be simplified?}
  \item{\dots}{Not used here.}
}
%\section{Methods}{
%\describe{
%\item{\code{signature(object = "MRSP"):}}{Extracts the regression coefficients of an \code{\link{MRSP-class}} object.
%}
%}}
\details{Depending on argument \code{type}, the regression coefficients
belonging to the original covariates or those belonging to standardized covariates
are returned. For \code{type="prethreshold"}, the regression coefficients prior to
(potential) thresholding (and a potential backtransformation from standardized to
original scale) are returned. \code{type} is found by a call to \code{match.arg}.

Depending on arguments \code{standardize} and \code{threshold} in the
call to \code{\link{MRSP}} that created \code{object}, the output of this function might be identical
for different choices of \code{type}.

If \code{simplify=TRUE}, the coefficient object is returned via \code{invisible}
and its actual content is printed in structured and concise fashion.
}
\value{An object of class \code{MRSP.coef}. Depending on \code{simplify}, information
about the coefficients is printed.
}
\author{Wolfgang Poessnecker}
\examples{
## see the examples in ?MRSP
}
\keyword{methods}

