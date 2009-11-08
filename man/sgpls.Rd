\name{sgpls}
\alias{sgpls}
\title{Fit SGPLS classification models}
\description{
 Fit a SGPLS classification model.
}
\usage{
sgpls( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector of class indices. }
  \item{K}{ Number of hidden components. }
  \item{eta}{ Thresholding parameter. \code{eta} should be between 0 and 1. }
  \item{scale.x}{ Scale predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{eps}{ An effective zero for change in estimates. Default is 1e-5. }
  \item{denom.eps}{ An effective zero for denominators. Default is 1e-20. }
  \item{zero.eps}{ An effective zero for success probabilities. Default is 1e-5. }
  \item{maxstep}{ Maximum number of Newton-Raphson iterations.
                    Default is 100. }
  \item{br}{ Apply Firth's bias reduction procedure? }
  \item{ftype}{ Type of Firth's bias reduction procedure.
     Alternatives are \code{"iden"} (the approximated version)
     or \code{"hat"} (the original version).
     Default is \code{"iden"}.}
}
\details{
The SGPLS method is described in detail in Chung and Keles (2009).
SGPLS provides PLS-based classification with variable selection,
by incorporating sparse partial least squares (SPLS) proposed in Chun and Keles (2009)
into a generalized linear model (GLM) framework.
\code{y} is assumed to have numerical values, 0, 1, ..., G,
where G is the number of classes subtracted by one.
}
\value{A \code{sgpls} object is returned.
print, predict, coef methods use this object.}
\references{
Chung, D. and Keles, S. (2009). "Sparse partial least squares classification
  for high dimensional data"
(\url{http://www.stat.wisc.edu/~keles/Papers/C_SPLS.pdf}).
  
Chun, H. and Keles, S. (2009). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
  To appear in \emph{Journal of the Royal Statistical Society - Series B}
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{print.sgpls}}, \code{\link{predict.sgpls}}, and \code{\link{coef.sgpls}}. }
\examples{
data(prostate)
# SGPLS with eta=0.6 & 3 hidden components
f <- sgpls( prostate$x, prostate$y, K=3, eta=0.6, scale.x=FALSE )
print(f)
# Print out coefficients
coef.f <- coef(f)
coef.f[ coef.f!=0, ]
}
\keyword{models}
\keyword{multivariate}
