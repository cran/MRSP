\name{TravelMode}
\alias{TravelMode}
\docType{data}
\title{Travel Mode Choice}
\description{Data on travel mode choice for travel between Sydney and Melbourne, Australia.}
\usage{data("TravelMode")}


\format{
A data frame containing 840 observations on 4 modes for 210 individuals.
This corresponds to so-called \bold{long} format since we have 4 rows per
individual in a situation with 4 classes/categories of the response variable \code{mode}.
 \describe{
    \item{individual}{Factor indicating individual with levels \code{1} to \code{200}.}
    \item{mode}{Factor indicating travel mode with levels
      \code{"car"}, \code{"air"}, \code{"train"}, or \code{"bus"}.}
    \item{choice}{Factor indicating choice with levels \code{"no"} and \code{"yes"}.}
    \item{wait}{Terminal waiting time, 0 for car.}
    \item{vcost}{Vehicle cost component.}
    \item{travel}{Travel time in the vehicle.}
    \item{gcost}{Generalized cost measure.}
    \item{income}{Household income.}
    \item{size}{Party size.}
  }
}
\details{
Data and description are taken from R package \bold{AER} by Christian Kleiber and Achim Zeileis.

In \code{MRSP}-terminology, variables \code{wait}, \code{vcost}, \code{travel} and \code{gcost} are
category-specific since they take different values for different categories of \code{mode}. The variables
\code{income} and \code{size} are global. In the literature on discrete choice modelling,
these two variable types are called \dQuote{alternative-specific} and \dQuote{individual-specific}, respectively.
}
\source{R package \bold{AER}, which states its source as follows:\cr

Online complements to Greene (2003).

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}
\keyword{datasets}
