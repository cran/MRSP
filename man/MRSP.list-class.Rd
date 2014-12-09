\name{MRSP.list-class}
\Rdversion{1.1}
\docType{class}
\alias{MRSP.list-class}
\alias{cv,MRSP.list-method}
\alias{extract,MRSP.list-method}
\alias{getCall,MRSP.list-method}
\alias{plot,MRSP.list-method}
\alias{refit,MRSP.list,missing-method}
\alias{select,MRSP.list-method}
\alias{show,MRSP.list-method}
\alias{update,MRSP.list-method}

\title{Class \code{"MRSP.list"}}
\description{
 An object containing the output of \code{\link{MRSP}}-calls with more than one
 tuning parameter value, i.e. model output related to a whole \dQuote{coefficient path}.
 With \code{nrlambda} denoting the number of different lambda values for which to fit
 a penalized multinomial response model, \code{MRSP.list}-objects are lists of length
 \code{nrlambda}. Each entry of this list is an object of class \code{MRSP}
 and contains the model output for one concrete lambda value. If an \code{MRSP.list}-object
 was created by function \code{\link{MRSP}}, it additionally has
 attributes \dQuote{topcall}, \dQuote{call} and \dQuote{dat}. Those attributes store, respectively, the call to
 \code{MRSP}, the call to \code{MRSP.fit} which was prepared by \code{MRSP} and
 the data object to be supplied to \code{MRSP.fit}. Note that
 \code{dat} contains the \emph{standardized} covariates if \code{standardize = TRUE}.
}
\section{Objects from the Class}{
Objects of this class are created by function \code{\link{MRSP.fit}}, which will
most likely be called internally by the end-user function \code{\link{MRSP}}.}
\section{Methods}{
  \describe{
    \item{cv}{\code{signature(object = "MRSP.list")}: \code{cv(object, k, type, parallel, cores, ...)} performs \code{k}-fold crossvalidation of the models in \code{object}.
          The output is a list whose entry \code{mean} gives the mean values of the criterion specified via argument \code{type} over the crossvalidation folds; for each model found in \code{object}.
          \code{type} can be \code{"deviance"} (the default), \code{"loglik"} or \code{"brier"}, using the deviance, loglikelihood or the Brier score. Depending on the logical argument \code{parallel} and  integer argument \code{cores},
          crossvalidation is performed either serially or with \code{cores} parallel processes.  }
    \item{extract}{\code{signature(object = "MRSP.list")}: \code{extract(object, slotname)} returns the value of the slot whose name was given as a character string
          from all elements of \code{object}. The result is coerced to a vector or matrix if appropriate. }
    \item{getCall}{\code{signature(x = "MRSP.list")}: Get the call that created the \code{MRSP.list}-object. }
    \item{plot}{\code{signature(x = "MRSP.list")}: Plot coefficient paths.  }                         %% hier noch machen
    %\item{refit}{\code{signature(object = "MRSP.list", arglist = "missing")}: Internal }
    \item{select}{\code{signature(object = "MRSP.list")}: \code{select(object, criterion=c("AIC","BIC","cv"), ...)} selects one concrete \code{MRSP} model from \code{object} based on \code{criterion}.
          For \code{criterion="cv"}, an additional, integer argument \code{k}, a logical argument \code{parallel} and an integer argument \code{cores} can be supplied, see \code{cv} above.  }                  %% noch dokumentieren
    \item{show}{\code{signature(object = "MRSP.list")}: Print some basic infos about the \code{MRSP.list} object. }
    %\item{update}{\code{signature(object = "MRSP.list")}: Not yet working }
	 }
}
\references{Tutz, G., Poessnecker, W., Uhlmann, L. (2015)
 \emph{Variable Selection in General Multinomial Logit Models} \cr
 \emph{Computational Statistics and Data Analysis, Vol. 82, 207-222.}\cr
 \url{http://www.sciencedirect.com/science/article/pii/S0167947314002709}\cr
}
\author{
Wolfgang Poessnecker
}

\seealso{\code{\linkS4class{MRSP}}}
\examples{
showClass("MRSP.list")
}
\keyword{classes}
