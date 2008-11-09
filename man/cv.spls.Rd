\name{cv.spls}
\alias{cv.spls}
\title{Compute and plot the cross-validated error for SPLS regression}
\description{
 Draw the heatmap-type plot of
 the v-fold cross-validated mean square prediction error
 and return the optimal eta (thresholding parameter)
 and K (number of hidden components).
}
\usage{
cv.spls( x, y, fold=10, K, eta, kappa=0.5,
        select="pls2", fit="simpls",
        scale.x=TRUE, scale.y=FALSE )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector or matrix of responses. }
  \item{fold}{ Number of cross-validation folds. Default is 10-folds. }
  \item{K}{ Number of hidden components. }
  \item{eta}{ Thresholding parameter. \code{eta} should be between 0 and 1. }
  \item{kappa}{ Parameter to control the effect of
     the concavity of the objective function
     and the closeness of the original and surrogate direction vectors.
     \code{kappa} is relevant only when the responses are multivariate.
     \code{kappa} should be between 0 and 0.5. Default is 0.5. }
  \item{select}{ PLS algorithm for variable selection.
                    Alternatives are \code{"pls2"} or \code{"simpls"}.
                    Default is \code{"pls2"}.}
  \item{fit}{ PLS algorithm for model fitting. Alternatives are
                    \code{"kernelpls"}, \code{"widekernelpls"},
                    \code{"simpls"}, or \code{"oscorespls"}.
                    Default is \code{"simpls"}.}
  \item{scale.x}{ Scale the predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{scale.y}{ Scale the responses by dividing each response variable
                by its sample standard deviation? }
}
\value{
Invisibly returns a list with components:
  \item{mspemat}{ Matrix of the cross-validated mean squared prediction error.
                Rows correspond to \code{eta} and
                columns correspond to the number of components (\code{K}). }
  \item{etaopt}{ Optimal \code{eta}. }
  \item{Kopt}{ Optimal \code{K}. }
}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ print, plot, predict, and coef methods for spls. }
\examples{
data(yeast)
set.seed(1)
# MSPE plot. eta is searched between 0.1 and 0.9 and
# the number of hidden components is searched between 1 and 10
cv <- cv.spls( yeast$x, yeast$y, K = c(1:10), eta = seq(0.1,0.9,0.1) )
# Optimal eta and K
cv$etaopt
cv$Kopt
}
\keyword{regression}
\keyword{multivariate}
