\name{correct.spls}
\alias{correct.spls}
\title{Correct the initial SPLS coefficient estimates
        based on the bootstrapped confidence intervals}
\description{
 Correct the initial SPLS coefficient estimates of the selected predictors
 based on the bootstrapped confidence intervals and draw
 heatmap-type plots of the original and the corrected coefficient estimates.
}
\usage{
correct.spls( object, plot.it=TRUE )
}
\arguments{
  \item{object}{ An object obtained from the function \code{ci.spls}. }
  \item{plot.it}{ Draw the heatmap-type plots of
                    the original coefficient estimates
                    and the corrected coefficient estimates? }
}
\details{
The set of the selected variables is updated by setting
the coefficients with zero-containing confidence intervals to zero.
}
\value{ Invisibly returns a matrix of corrected coefficient estimates.
}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ ci.spls. }
\examples{
data(mice)
# SPLS with eta=0.6 & 1 latent components
f <- spls( mice$x, mice$y, K=1, eta=0.6 )
# Calculate the confidence intervals of coefficients
ci.f <- ci.spls(f)
# Corrected coefficient estimates
cf <- correct.spls( ci.f )
cf[20,1:5]
}
\keyword{regression}
\keyword{multivariate}
