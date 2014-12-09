\name{MRSP-class}
\Rdversion{1.1}
\docType{class}
\alias{MRSP-class}
\alias{AIC,MRSP-method}
\alias{BIC,MRSP-method}
\alias{bootstrap,MRSP-method}
%\alias{coef,MRSP-method}
\alias{cv,MRSP-method}
\alias{extract,MRSP-method}
\alias{fitted,MRSP-method}
\alias{getCall,MRSP-method}
\alias{logLik,MRSP-method}
\alias{nobs,MRSP-method}
\alias{predict,MRSP-method}
\alias{pval,MRSP-method}
\alias{refit,MRSP,list-method}
\alias{refit,MRSP,missing-method}
\alias{residuals,MRSP-method}
\alias{se,MRSP-method}
\alias{show,MRSP-method}
\alias{summary,MRSP-method}

\title{Class \code{"MRSP"}}
\description{
An object containing the output of fitting a multinomial response model with
structured penalties for one concrete value of the tuning parameter (or one
concrete combination of all used tuning parameters).
}
\section{Objects from the Class}{
Objects of this class are created by function \code{\link{MRSP.fit}}, which will
most likely be called internally by the end-user function \code{\link{MRSP}}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"MRSP.coef"} containing the estimated
                        regression coefficients. Structurally, it is simply a list
                        of one or two matrices.}
    \item{\code{coef.stand}:}{Same as slot \code{coef}, but contains the coefficients
                              belonging to standardized covariates.}
    \item{\code{coef.pretres}:}{Same as slot \code{coef}, but contains the coefficients
                                belonging to standardized covariates before any
                                (potential) thresholding took place.}
    \item{\code{dat}:}{Data object in the form required by \code{\link{MRSP.fit}}.
                       To save memory, this slot is usually set to \code{NULL}.}
    \item{\code{x.original}:}{Original matrix containing the \code{"x"}-variables.
                       To save memory, this slot is often times set to \code{NULL}.}
    \item{\code{x.stand}:}{The \code{"x"}-matrix after standardization. }
    \item{\code{V.original}:}{The original data object containing the category-specific
                              predictors. To save memory, this slot is often times \code{NULL}.}
    \item{\code{V.stand}:}{The \code{"V"}-object after standardization.}
    \item{\code{y}:}{Response matrix in the form required by \code{\link{MRSP.fit}}.
                       To save memory, this slot is often times set to \code{NULL}.}
    \item{\code{weights}:}{Vector of observation weights.}
    \item{\code{penindex}:}{Object specifying how each covariate is penalized. See
                            the documentation of \code{\link{MRSP.fit}} for details.}
    \item{\code{grpindex}:}{Object specifying which predictors form parameter groups. See
                            the documentation of \code{\link{MRSP.fit}} for details.}
    \item{\code{penweights}:}{Object specifying the weighting of the penalty on different
                              parameters or parameter groups.  See \code{\link{MRSP.fit}}
                              for details.}
    \item{\code{guessed.active}:}{Object indexing the covariates that are found to
                                  have an effect on the response. If necessary, this
                                  refers to groups of coefficients.}
    \item{\code{guessed.active.coef}:}{An object with the same structure as slot \code{coef},
                                       indexing the atomic coefficients that are nonzero.}
    \item{\code{guessed.active.groupdiff}:}{Index of \code{"columnwise"} parameter groups
                                            that contain at least two different parameter
                                            values.}
    \item{\code{guessed.active.diff}:}{An index of nonzero pairwise differences between parameters
                                       belonging to the same covariate.}
    \item{\code{df}:}{Estimated effective degrees of freedom.}
    \item{\code{tuning}:}{List of tuning parameter values as supplied to internal
                          function \code{fista}.}
    \item{\code{lambda}:}{The lambda value, which controls the degree of penalization of most
                          traditional penalties.}
    \item{\code{lambdaR}:}{The lambda value used for ridge penalties.}
    \item{\code{lambdaF}:}{The lambda value used for fusion penalties.}
    \item{\code{fusion}:}{Either \code{FALSE} or character string specifying the type of fusion
                          that was used in fusion penalties. Note that those fusion penalties
                          are not yet supported for end-users of \code{MRSP}.}
    \item{\code{gamma}:}{A numeric  that weighs lasso penalties vs CATS lasso penalties.
                         See \code{\link{MRSP}} for details.}
    \item{\code{psi}:}{A numeric that weighs penalties on coefficients of global vs
                       category-specific predictors. See \code{\link{MRSP}} for details.}
    \item{\code{eta}:}{A \code{nobs x K} matrix of linear predictor values.}
    \item{\code{mu}:}{A \code{nobs x K} matrix of estimated probabilities for the response categories.}
    \item{\code{offset}:}{Vector of offset values that were added to the linear predictors.}
    \item{\code{residuals}:}{Object with residuals. Currently always \code{NULL}. Use method
                             \code{residuals} instead to compute residuals.}
    \item{\code{mlfit}:}{Not to be used by or of interest for end-users.}
    \item{\code{AIC}:}{The AIC of the fitted model.}
    \item{\code{BIC}:}{The BIC of the fitted model.}
    \item{\code{Brier}:}{The Brier score of the fitted model.}
    \item{\code{threshold}:}{The numeric threshold used.}
    \item{\code{refit}:}{Logical indicating whether the object results from a refitting procedure.}
    \item{\code{indg}:}{Not to be used by or of interest for end-users.}
    \item{\code{indcs}:}{Not to be used by or of interest for end-users.}
    \item{\code{model}:}{The \code{model}-object used. To save memory, this is often stored as an expression.}
    \item{\code{constr}:}{The identifiability constraint that was used. See \code{\link{MRSP}} for details.}
    \item{\code{control}:}{Object of class \code{"MRSP.control"} that contains control information.}
    \item{\code{fn.val}:}{The final value of the objective function that was minimized, i.e.
                          the negative penalized loglikelihood: \code{-loglik + lambda*penweights*penalty}.}
    \item{\code{loglik}:}{Loglikelihood value of the fitted model.}
    \item{\code{penalty}:}{Value of the penalty term for the fitted model, weighted with
                           the corresponding tuning parameters and \code{penweights}.}
    \item{\code{iter.count}:}{Number of iterations until convergence of the proximal gradient algorithm.}
    \item{\code{best.iter}:}{Iteration number with the best value of the objective function.}
    \item{\code{ridgestabil}:}{Logical indicating whether a small, untuned ridge penalty was
                               applied to all coefficients in order to stabilize otherwise
                               diverging estimates.}
    \item{\code{name}:}{A character string specifying the name and thus type of the fitted model.}
    \item{\code{fisher}:}{Fisher matrix. Currently not computed and thus always \code{NULL}.}
    \item{\code{arglist}:}{Not to be used by or of interest for end-users.}
    \item{\code{call}:}{The call to \code{\link{MRSP.fit}} that created this object.}
  }
}
\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "MRSP")}: Compute the AIC of an \code{MRSP} object.}
    \item{BIC}{\code{signature(object = "MRSP")}: Compute the BIC of an \code{MRSP} object. }
    %\item{bootstrap}{\code{signature(object = "MRSP")}: Internal }
    \item{coef}{\code{signature(object = "MRSP")}: See \code{\link{coef-methods}}. }
    %\item{cv}{\code{signature(object = "MRSP")}: See ... }                                   %% hier noch
    %\item{extract}{\code{signature(object = "MRSP")}: Included for completeness, but only really useful for MRSP-list class objects. }
    \item{fitted}{\code{signature(object = "MRSP")}: ompute fitted values, which for multinomial response correspond to class probabilities P(Y = r).
          If option \code{convert2hazard=TRUE} and a sequential model is used, \dQuote{discrete} hazard rates P(Y = r | Y >= r) are returned instead.}
    %\item{getCall}{\code{signature(x = "MRSP")}: Internal }
    \item{logLik}{\code{signature(object = "MRSP")}: Returns the loglikelihood and the (estimated, effective) degrees of freedom.  }
    \item{nobs}{\code{signature(object = "MRSP")}: Returns the number of invidual observations. }
    \item{predict}{\code{signature(object = "MRSP")}: \code{predict(object, newdata, type=c("response", "link"), ...)} predicts the response values (\code{type="response"})
          or the linear predictors (\code{type="link"}) for the observations given in \code{newdata}. Additional arguments \code{offset} and \code{weights}
          can specify offsets and weights to be used.  An argument \code{convert2hazard} can be supplied for sequential models, see \code{fitted} above. }
    %\item{pval}{\code{signature(object = "MRSP")}: Internal }
    %\item{refit}{\code{signature(object = "MRSP", arglist = "list")}: Internal }
    %\item{refit}{\code{signature(object = "MRSP", arglist = "missing")}: Internal }
    \item{residuals}{\code{signature(object = "MRSP")}: Depending on argument \code{type = c("deviance", "pearson")},
          which is matched via \code{match.arg}, deviance or pearson residuals are returned. }
    %\item{se}{\code{signature(object = "MRSP")}: Internal }
    \item{show}{\code{signature(object = "MRSP")}: Print some basic infos about the \code{MRSP} object.}
    \item{summary}{\code{signature(object = "MRSP")}: Show some slots of an \code{MRSP} object which are typically of interest.}
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

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~                                 anmerkung fuer mich selbst: hier die methoden reinschreiben, die ein eigenes Rd file bekommen. bei details darauf verweisen.
\code{\linkS4class{MRSP.list}}
}
\examples{
showClass("MRSP")
## for examples, see ?MRSP
}
\keyword{classes}
