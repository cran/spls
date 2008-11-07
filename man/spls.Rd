\name{spls}
\alias{spls}
\title{Fit SPLS regression models}
\description{
 Fit a SPLS regression.
}
\usage{
spls( x, y, K, eta, kappa=0.5, select="pls2", fit="simpls",
    scale=TRUE, center=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100)
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector or matrix of responses. }
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
  \item{scale}{ Scale the predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{center}{ Center the predictors? }
  \item{scale.y}{ Scale the responses by dividing each response variable
                by its sample standard deviation? }
  \item{eps}{ An effective zero. Default is 1e-4. }
  \item{maxstep}{ Maximum number of iterations when fitting direction vectors.
                    Default is 100. }
}
\details{
The SPLS method is described in detail in Chun and Keles (2007).
SPLS directly imposes sparsity on the dimension reduction step of PLS
in order to achieve accurate prediction and variable selection simultaneously.
The option \code{select} refers to the PLS algorithm for variable selection.
The option \code{fit} refers to the PLS algorithm for model fitting
and \code{spls} utilizes the algorithms offered by the \pkg{pls} package for this purpose.
See the help files of the function \code{plsr} in the \pkg{pls} package for more detail.
The user should install the \pkg{pls} package before using \pkg{spls} functions.
The choices for \code{select} and \code{fit} are independent.
}
\value{A \pkg{spls} object is returned.
print, plot, predict, coef, ci.spls, coefplot.spls methods use this object.}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ print, plot, predict, coef, ci.spls, coefplot.spls methods for spls. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
print(f)
# Print out coefficients
coef.f <- coef(f)
coef.f[,1]
# Coefficient path plot
plot( f, yvar=1 )
x11()
# Coefficient plot of the selected variables
coefplot.spls( f, xvar=c(1:4) )
}
\keyword{regression}
\keyword{multivariate}
