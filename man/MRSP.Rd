\name{MRSP}
\alias{MRSP}
\alias{MRSP,ANY-method}

\title{Multinomial Response Models with Structured Penalties}
\description{Fit models with multinomial response (including ordinal regression)
 with structured penalties.}
\usage{
MRSP(formula, data, class.names = NULL, model = multinomlogit(), constr = NULL,
     offset = NULL, weights = NULL, penweights = NULL, standardize = TRUE,
     nrlambda = 50, lambdamin = 0.01, lambdamax = NULL, control = NULL,
     penalty = TRUE, group.classes = TRUE, group.dummies = TRUE,
     sparse.groups = FALSE, adaptive = FALSE, threshold = FALSE, refit = FALSE,
     lambda, lambdaR = 0, lambdaF = 0, gamma = 1, psi = 1, fusion = FALSE,
     nonneg = FALSE, y = NULL, X = NULL, Z = NULL, penindex = NULL,
     grpindex = NULL, mlfit = NULL, perform.fit = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A symbolic description of the model to be fitted. The left-hand
    side specifies the response, which must be a categorical variable with \code{K} 
   categories. The right-hand-side specifies the covariate structure. See details below.}
  \item{data}{A data frame containing the variables in \code{formula}. It must 
   be in \bold{long} format. This means that \code{K} rows are required for each 
   individual observation. For details and a possibility to convert to this 
   format, see \code{\link[mlogit]{mlogit.data}} from the \pkg{mlogit} package. 
   The response variable must be either a 0-1-vector or a logical 
   vector, with 1 or TRUE indicating the observed class/category/alternative.}
  \item{class.names}{An optional character vector of length \code{K} specifying 
   the names of the response classes/categories.}
  \item{model}{An object of class \code{MRSP.model} that specifies the model to 
   be used. Currently, \code{model = multinomlogit()} and 
   \code{model = sequentiallogit()} are available, yielding a multinomial or
   sequential logit model, respectively. Cumulative logit models will be included 
   in future versions of \code{MRSP}.}
  \item{constr}{The identifiability constraint to be used. The coefficients of
   predictors which do not vary over categories (i.e. global/individual-specific 
   predictors) are not identifiable in (unpenalized) multinomial logit models. 
   Either an integer in \code{[1,K]} or \code{"symmetric"} or \code{"none"}. See
   details below.}
  \item{offset}{An optional vector of offsets. Must be of appropriate length.}
  \item{weights}{An optional vector of observations weights. Must be of 
   appropriate length.}
  \item{penweights}{An object containing weights that modify the penalty terms 
   on different parameters. See \code{MRSP.fit} for details.}
  \item{standardize}{If \code{TRUE}, predictors are mean-centered and standardized to 
   unit variance. The reported coefficients, by default, correspond to the 
   original covariates.}
  \item{nrlambda, lambdamin, lambdamax}{If a sequence of \code{lambda} values is
   not specified explicitly via argument \code{lambda}, \code{MRSP} computes a 
   suitable grid of length \code{nrlambda}, ranging from \code{lambdamin} to 
   \code{lambdamax}. If \code{lambdamax} is missing, \code{MRSP} tries to find
   a suitable value for it.}
  \item{control}{An object of class \code{MRSP.control} that stores control 
   information. It's slots \code{max.iter} and \code{rel.tol} specify the max
   number of iterations and the relative change in penalized log-likelihood 
   values that indicates convergence. The other slots should not be changed 
   unless by experienced users.}
  \item{penalty}{Specifies the general type of penalty to be applied. FALSE 
   means no penalty is used. \code{TRUE} means that a lasso-type penalty is applied.
   \code{"Ridge"} applies a ridge penalty. Arguments \code{group.classes}, 
   \code{group.dummies} and \code{sparse.groups} have no effect unless 
   \code{penalty = TRUE}.}
  \item{group.classes}{If \code{TRUE}, lasso-type penalties will be grouped across all 
   coefficients belonging to the same covariate. This corresponds to the 
   \sQuote{CATS Lasso} penalty proposed in Tutz, Poessnecker and Uhlmann (2015).}   
  \item{group.dummies}{If entries on the rhs of \code{formula} refer to factor
   variables with more than 2 levels, several dummy variables will enter the 
   model that are related to the same actual covariate. Setting this argument to
   \code{TRUE} will yield penalty terms that treat all corresponding coefficients as 
   one parameter group. This corresponds to the original \sQuote{Group Lasso}
   of Yuan and Lin (2006).}
  \item{sparse.groups, gamma}{If \code{TRUE}, parameter groups will also be penalized 
   by an L1 penalty on top of the unsquared L2 penalty. If the L2 penalty uses a 
   tuning parameter of value lambda, the L1 penalty will use tuning parameter 
   lambda*gamma.}
  \item{adaptive}{Should adaptive weights be used? Use \code{adaptive="ML"} to 
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
   This can improve variable selection, but tends to be rather slow and memory-consuming.}
  \item{lambda}{Tuning parameter(s) to be used. Can be a nonnegative scalar or 
   vector. If missing, \code{MRSP} computes a suitable grid of lambda values, 
   see also \code{nrlambda}.}
  \item{lambdaR}{Lambda value(s) for ridge penalties. Same structure as 
   \code{lambda}.}
  \item{lambdaF}{Lambda values(s) for fusion penalties. Same structure as 
   \code{lambda}. Not yet supported for end-users of \code{MRSP}, but included 
   for compatibility with future releases of \code{MRSP}.}
  \item{psi}{A numeric that balances the weighting of penalties on global and 
   class-/category-/alternative-specific predictors (if present) in multinomial 
   logit models. Penalties on global predictors are weighted with \code{psi},
   penalties on class-specific predictors are weighted with \code{(2-psi)}.}
  \item{fusion}{If fusion penalties are used, this specifies the type of fusion.
   Not yet supported for end-users of \code{MRSP}, but included 
   for compatibility with future releases of \code{MRSP}.}
  \item{nonneg}{If \code{TRUE}, all coefficients are restricted to be nonnegative.}
  \item{y, X, Z, penindex, grpindex, mlfit}{These optional arguments allow to 
   supply the corresponding objects directly to \code{MRSP.fit}. It is highly 
   recommended not to use those arguments and to use \code{formula} instead for
   the model specification.}
  \item{perform.fit}{If TRUE, the model is fitted. If FALSE, function \code{MRSP}
   prepares the call to \code{MRSP.fit} and returns a list from which this call 
   can be accessed.}
  \item{\dots}{Further arguments to be passed.}
}  
\details{For \code{model = multinomlogit()}, a formula of the form
 \dQuote{\code{Y ~ x | z1 | z2}} yields linear predictors
 \deqn{\eta_r = \beta_{0r} + x^T\beta_r + z_{1_r}^T\gamma + z_{2_r}^T\delta_r}{%, \qquad r=1,\ldots,K,}{%
       eta_r = beta_0r  +  x ' beta_r  +  z1_r ' gamma  +  z2_r ' delta_r}%,         r=1,...,K,}
 for \eqn{r=1,\ldots,K,}
 which are connected to probabilities by the multinomial logit link,
 also known as softmax function:
 \deqn{P(Y = r | x, z_1, z_2) = \frac{\exp(\eta_{r})}{\sum_{s=1}^K \exp(\eta_{s})}}{%
       P(Y = r | x, z1, z2) = exp(eta_r) / sum[ exp(eta_s)]_{s=1,...,K}.}
       
 This means that the \code{x}-variables have global values that are class-/category-/alternative-unspecific.
 The coefficients belonging to those global variables are not identifiable in the
 generic form of the multinomial logit model as given above. Therefore, an identifiability
 constraint can be specified with argument \code{constr}: By setting \code{constr = r} (with \eqn{r \in [1,K]}{r in [1,K]}),
 category \code{r} is chosen as reference category. Technically, this means setting \eqn{\beta_r = 0}{beta_r = 0}.
 Alternatively, \code{constr="symmetric"} specifies a so-called symmetric side constraint, which technically
 means imposing that
 \deqn{\sum_{s=1}^K \beta_{sj} = 0}{% \qquad \forall j=1,\ldots,p.}{%
       sum[ beta_sj ]_{s=1,...,K} = 0,}%          j=1,...,p.}
 for all \eqn{j=1,\ldots,p.}
 If \code{constr = "none"}, no constraint is used for penalized parameter groups and
 identifiability is ensured by the penalty term (see Friedman, Hastie and Tibshirani, 2010).
 Coefficients of an unpenalized \code{x}-variable are subject to a symmetric side constraint in this case. \cr
 The \code{z1}- and \code{z2}-variables are class-/category-/alternative-specific. The \code{z1}-variables
 have a global effect while the category-specific \code{z2}-variables are equipped with coefficients
 that are also category-specific. Note that no identifiability constraints are required for category-specific variables.\cr\cr
 
 
 
 For \code{model = sequentiallogit()}, a formula of the form \dQuote{\code{Y ~ x1 | x2}}
 yields linear predictors of the form
 \deqn{\eta_r = \beta_{0r} + x_1^T\alpha + x_2^T\beta_r}{%, \qquad r=1,\ldots,K-1,}{%
       eta_r = beta_0r  +  x1 ' alpha  +  x2 ' beta_r}%,       r=1,...,K-1,}
 for \eqn{r=1,\ldots,K-1,}
 which are connected to conditional probabilities via the logit link in the following way:
 \deqn{P(Y = r | Y \ge r, x_1, x_2) = \frac{\exp(\eta_r)}{1 + \exp(\eta_r)}}{%,\qquad r=1,\ldots,K-1}{%
       P(Y = r | Y >= r, x1, x2) = exp(eta_r) / (1 + exp(eta_r))}%,     r=1,...,K-1.}
 for \eqn{r=1,\ldots,K-1.}\cr
 Note that in the sequential case, slot \code{"mu"} of an \code{MRSP-class} object
 contains the unconditional class probabilities P(Y = r). If you want to get the
 \dQuote{discrete} hazard rates P(Y = r | Y >= r), use methods \code{fitted} or
 \code{predict} with argument \code{convert2hazard=TRUE}.}
\value{Depending on \code{nrlambda}, either an object of class \code{MRSP} or of
 class \code{MRSP.list}, which are lists of length
 \code{nrlambda} whose elements are \code{MRSP} objects. Additionally, the 
 attributes \dQuote{topcall}, \dQuote{call} and \dQuote{dat} store, respectively, the call to
 \code{MRSP}, the call to \code{MRSP.fit} which was prepared by \code{MRSP} and 
 the data object to be supplied to \code{MRSP.fit}. Note that
 \code{dat} contains the \emph{standardized} covariates if \code{standardize = TRUE}.}
\author{Wolfgang Poessnecker}
\keyword{Regularization}
\keyword{Penalization}
\keyword{Group Lasso}
\keyword{Multinomial Logit Model}
\keyword{CATS Lasso} 
\references{Tutz, G., Poessnecker, W. and Uhlmann, L. (2015):
 \emph{Variable Selection in General Multinomial Logit Models} \cr
 \emph{Computational Statistics and Data Analysis, Vol. 82, 207-222.}\cr
 \url{http://www.sciencedirect.com/science/article/pii/S0167947314002709}\cr\cr
 Friedman, J., Hastie, T. and Tibshirani, R. (2010):
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent},   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software, Vol. 33(1), 1-22.}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr\cr
 Yuan, M. and Lin, Y. (2006):
 \emph{Model selection and estimation in regression with grouped variables}\cr
 \emph{Journal of the Royal Statistical Society Series B, Vol. 68(1), 49-67.}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

## load data
data(TravelMode, package="MRSP")
## bring the response variable to the form required by MRSP
TravelMode$choice <- ifelse(TravelMode$choice=="yes",1,0)

## construct a list of fitted models for different lambda values. 
## income is a global predictor, wait is a class-specific predictor with global
## coefficients, vcost and travel are specified as class-specific predictors 
## with class-specific coefficients. The fourth category ("car") is chosen as 
## reference. 
fit <- MRSP(choice~income|wait|vcost+travel, data=TravelMode, constr=4,
            class.names=levels(TravelMode$mode), lambdamax=150, nrlambda=10,
            group.classes=TRUE, sparse.groups=FALSE, adaptive="ML")
            
fit

## slots can be extracted from all elements via function 'extract':
BICs <- extract(fit, "BIC")
lambdagrid <- extract(fit, "lambda")


## to select a concrete lambda/model, one can use function 'select'. Here, we
## chose the best model according to its AIC value
bestfit <- select(fit, "AIC")
bestfit

## some methods:
summary(bestfit)
BIC(bestfit)
fitted(bestfit)[1:6,]
bestfit@coef
## get the coefficients belonging to standardized predictors:
coef(bestfit, type="stand")
residuals(bestfit)
predict(bestfit, newdata = TravelMode[1:40,c(4,5,6,8)])

## plot some coefficient paths:
par(mfrow=c(1,2))
## you can either specify the number of the variable...
plot(fit, 2, legendpars=list(x="bottomright"))
## ... or its name as a character string. lcex is the cex parameter for legends.
## set it to zero to disable legend plotting. 
plot(fit,"travel", lambda = bestfit@lambda, lcex=0)

}
