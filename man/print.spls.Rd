\name{print.spls}
\alias{print.spls}
\title{Print fits of a SPLS object}
\description{
 Print out SPLS fits, the number and the list of selected predictors.
}
\usage{
\method{print}{spls}( x, ... )
}
\arguments{
  \item{x}{ A fitted SPLS object.}
  \item{...}{ Additonal arguments for generic \code{print}.}
}
\value{NULL.}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection"
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ plot, predict, and coef methods for spls. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
print(f)
}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
