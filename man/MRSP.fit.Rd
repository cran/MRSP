\name{MRSP.fit}
\alias{MRSP.fit}
\alias{MRSP.fit,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,list-method}
\alias{MRSP.fit,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,missing-method}
\alias{MRSP.fit,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,numeric-method}

\title{
Fitting function for multinomial response models with structured penalties
}
\description{
This function performs the actual fitting of multinomial response models with structured penalties.
Function \code{MRSP} is actually just a user-friendly wrapper that prepares calls to this function.
\code{MRSP.fit} is designed mostly for internal use and therefore not user-friendly, so that it is
highly recommended to use \code{MRSP} instead of calling \code{MRSP.fit} directly.
}
\usage{
MRSP.fit(dat, coef.init = NULL, coef.stand.init = NULL, coef.pretres.init = NULL,
         offset = NULL, weights = NULL, grpindex = NULL, penindex = NULL, lambda,
         lambdaR = lambda, lambdaF = lambda, gamma = 1, psi = 1, indg = NULL,
         indcs = NULL, model = NULL, constr = NULL, control = NULL,
         fista.control = NULL, Proximal.control = NULL, Proximal.args = NULL,
         penweights = NULL, mlfit = NULL, adaptive = FALSE, threshold = FALSE,
         refit = FALSE, fusion = FALSE, nonneg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A list that contains the data in the format required by \code{MRSP.fit}. \code{dat}
is a list with elements 'y', 'x' and, possibly, 'V'. If \code{nobs} individual
observations are available, \code{dat$y} must be a matrix of dimension \code{nobs x K}.
For \code{model=multinomlogit()}, \code{K} is equal to the number of categories of the
response variable, for\cr \code{model=sequentiallogit()}, it equals the number of response
categories minus 1. In other words: \code{K} here always refers to the column rank of the response matrix.
(Note that this is for notational convenience and in contrast to the documentation of \code{\link{MRSP}},
where \code{K} always refers to the number of categories of the multicategorical response!)
The entries of \code{dat$y} are either 0 or 1, with \code{dat$y[i,r]==1}
indicating that class r is observed for the i-th observation. For \code{model=multinomlogit()},
the rowSums of \code{dat$y} must be all 1; for \code{model=sequentiallogit()}, they must be
either 1 or 0; with a row of zeros indicating that the last category was observed.\cr

\code{dat$x} is a matrix of dimension \code{nobs x p} which contains covariates whose
value is constant across classes. They are called 'global predictors/covariates' in the
following. In the context of discrete choice modeling, they are often referred to as
'individual-specific' predictors.\cr

If available, covariates whose value varies from class to class can be included in an entry
\code{dat$V}. Such variables are called 'category-specific' in the following since their value
depends on the categories of the response variable. In the literature on discrete choice modeling,
they are often referred to as 'alternative-specific'.
These variables can either be equipped with global or with category-specific coefficients.
If a total of \code{L} category-specific variables shall be used, \code{dat$V} must be a list (!)
of length \code{K} whose elements each are matrices of dimension \code{nobs x L}.
}
  \item{coef.init}{
An optional coefficient object supplying initial coefficient values to be used. A list whose first entry is
a matrix of dimension \code{K x p}, with row \code{r} containing the coefficients for class \code{r}
and column \code{j} containing the coefficients for global predictor x_j. If category-specific
predictors are included, the second entry of \code{coef.init} is a matrix of dimension \code{K x L}
that contains the coefficients for those category-specific predictors.
}
  \item{coef.stand.init}{
Optional initial coefficient values for the standardized predictors. Same structure as \code{coef.init}.
}
  \item{coef.pretres.init}{
Optional initial coefficient values, prior to potential thresholding, for the standardized predictors. Same structure as \code{coef.init}.
}
  \item{offset}{
An optional vector or matrix of offset values to be used. Either length \code{nobs} or dimension
\code{nobs x K}.
}
  \item{weights}{
An optional vector of observation weights of length \code{nobs}.
}
  \item{grpindex}{
A list of one or two integer vectors that indicate which columns of the design matrix form a group that has
to be penalized jointly, e.g. the different dummies of a categorical predictor. The first element is the
grouping vector for x, the optional second one for V. Those columns with the same number belong to one group.
The numbers must begin with 1 and increase with every group. An example:
grpindex = list(c(1,2,3,3,4,4,4)) means that variables 1 and 2 form their own, 'scalar' group; variables 3 and 4  as well
as variables 5, 6 and 7 form multi-parameter-groups.
}
  \item{penindex}{
A list of one or two vectors which specifies the exact penalty type to use for each covariate. The first entry
specifies the penalty type for the variables in \code{dat$x}, the (optional) second entry those for the variables
in \code{dat$V}. The following penalty types are available:\cr\cr
      1:  global predictor ('x') whose coefficients shall be penalized with
          a group lasso penalty with grouping 'across' categories, i.e. CATS Lasso (see Tutz, Poessnecker and Uhlmann, 2015).\cr
      10: global predictor, unpenalized.\cr
      11: global predictor, sparse group lasso.\cr
      12: global predictor, ordinary lasso.\cr
      13: global predictor, ridge penalty. does not support penweights.\cr
      2:  category-specific predictor with global coefficient which is penalized with
          the ordinary (group-)lasso. (depending on grpindex.)\cr
      20: category-specific, unpenalized.\cr
      21: category-specific, ridge penalty.\cr
      3:  category-speficic predictor with category-specific coefficients
          that are penalized by a group lasso like in '1'.\cr
      30: category-specific with category-specific coefs, unpenalized.\cr
      31: category-specific with category-specific coefficients and sparse
          group lasso penalty.\cr
      32: cat-cat-specific, with ordinary lasso.\cr
      33: cat-cat, with ridge. does not support penweights.\cr
      4:  global predictor with global effect, penalized. (cf. the '2'-series).\cr
      40: global predictor, global effect, unpenalized.\cr
      41: global predictor, global effect, ridge penalty.\cr
      
      The '4-series' only makes sense for ordinal models!
}
  \item{lambda}{
Optional object specifying the lambda values to be used as tuning parameter(s) for
the main variable selection penalty. Either a vector or a single numeric. If missing,
a suitable grid of lambda values is computed. See arguments \code{nrlambda, lambdamin} and
\code{lambdamax} in function \code{\link{MRSP}}.
}
  \item{lambdaR}{
Lambda(s) to be used for ridge penalties. Typically, if only a ridge penalty and no other
penalty is used, one can specify the Ridge lambda via argument \code{lambda} instead.
}
  \item{lambdaF}{
Lambda(s) to be used for fusion penalties. Not available yet for end-users, but included for
compatibility with future releases of \code{MRSP}.
}
  \item{gamma}{
See argument \code{gamma} in \code{\link{MRSP}}.
}
  \item{psi}{
See argument \code{psi} in \code{\link{MRSP}}.
}
  \item{indg}{
A vector of the column indices of the category-specific variables that are equipped with
global coefficients.
}
  \item{indcs}{
A vector of the column indices of the category-specific variables that are equipped with
category-specific coefficients
}
  \item{model}{
An object of class \code{MRSP.model} that specifies the model to
be used. Currently, \code{model = multinomlogit()} and
\code{model = sequentiallogit()} are available, yielding a multinomial or
sequential logit model, respectively. Cumulative logit models will be included
in future versions of \code{MRSP}.
}
  \item{constr}{
The identifiability constraint to be used. The coefficients of
predictors which do not vary over categories (i.e. global/individual-specific
predictors) are not identifiable in (unpenalized) multinomial logit models.
If \code{constr} is an integer in [1, K], the corresponding class is used
as reference, which means that the coefficients of global predictors
for this class are set to 0. If \code{constr = "symmetric"}, a symmetric side
constraint is used, which means that all coefficients belonging to the same
global predictor sum to zero. If \code{constr = "none"}, no constraint is used
for penalized parameter groups and identifiability is ensured by the penalty
term (see Friedman, Hastie and Tibshirani, 2010.)
For ordinal regression, \code{constr} must take value \code{"none"}.
If left unspecified, a symmetric side constraint is used for multinomial and
no constraint for ordinal models.
}
  \item{control}{
An object of class \code{MRSP.control} that stores control
information. It's slots \code{max.iter} and \code{rel.tol} specify the max
number of iterations and the relative change in penalized log-likelihood
values that indicates convergence. The other slots should not be changed
unless by experienced users.
}
  \item{fista.control}{
An object of class \code{fista.control} that contains control information for
the FISTA algorithm that is internally used to compute numerical estimates. Not
intended for end-users!
}
  \item{Proximal.control, Proximal.args}{
Arguments to be passed to the proximal gradient algorithm. Not intended for
end-users!
}
  \item{penweights}{
An optional list containing weights for the various penalty terms of different
coefficients or coefficient groups. Assuming that category-specific covariates
are present, the first element of penweights is a list of length two, with the
first element of \code{penweights[[1]]} being a numeric of length \code{p} that
contains the weights for the CATS penalty on the group of coefficients of the
global covariates for the response classes. The second element of \code{penweights[[1]]}
is a numeric of length \code{L} with the group penalty weights for the
category-specific variables.
The second element of penweights is again a list of length two. The first element
of \code{penweights[[2]]} is a \code{K x p} matrix with penalty weights for
unstructured lasso penalties on atomic coefficients belonging to  global
predictors. The second element of \code{penweights[[2]]} is a \code{K x L}
matrix with penalty weights for unstructured lasso penalties on category-specific
covariates.
}
  \item{mlfit}{
A list that contains information about the ML or 'pseudo-ML' coefficients of the
specified model. It must contain at least one entry called 'coef.stand' that has
the same structure as the coefficient object (see \code{coef.init}). The
value of those parameters must be known to compute the effective degrees of
freedom of penalized parameter estimates with grouped penalties.
}
\item{adaptive}{
Should adaptive weights be used? Use \code{adaptive="ML"} to
   obtain the traditional adaptive weights proposed in the literature. Using
   \code{adaptive = TRUE} computes the penalized estimator with whatever penalty is
   specified and no adaptive weights and computes adaptive weights from the output
   of this penalized model. The final output is the computed with those
   adaptive weights. It is strongly recommended to prefer \code{adaptive="ML"}
   over \code{adaptive=TRUE}.}
  \item{threshold}{If \code{TRUE}, the coefficients will be thresholded with an
   appropriate threshold value. You can also specify an explicit nonnegative
   value to be used as the threshold.}
  \item{refit}{Should refitting be performed? If \code{TRUE}, the model is first fit
   traditionally, and then refitted on the active set found by this first fit.
   This can improve variable selection, but tends to be rather slow and time-consuming.}
  \item{fusion}{If fusion penalties are used, this specifies the type of fusion.
   Not yet supported for end-users of \code{MRSP}, but included
   for compatibility with future releases of \code{MRSP}.}
  \item{nonneg}{If \code{TRUE}, all coefficients are restricted to be nonnegative.}
  \item{\dots}{
Further arguments or objects to be passed to \code{MRSP.fit}.
}
}
\details{
This function does the actual work of fitting multinomial response models with
structured penalties. It is intended mainly for internal use. The main purpose
of function \code{\link{MRSP}} is to provide a user-friendly wrapper that prepares
and evaluates a call to \code{MRSP.fit}.

}
\value{
Depending on \code{nrlambda}, either an object of class \code{MRSP} or of
 class \code{MRSP.list}, which are lists of length
 \code{nrlambda} whose elements are \code{MRSP} objects.
}
\references{Tutz, G., Poessnecker, W., Uhlmann, L. (2015)
 \emph{Variable Selection in General Multinomial Logit Models} \cr
 \emph{Computational Statistics and Data Analysis, Vol. 82, 207-222.}\cr
 \url{http://www.sciencedirect.com/science/article/pii/S0167947314002709}\cr
}
\author{
Wolfgang Poessnecker
}
\keyword{CATS Lasso}

