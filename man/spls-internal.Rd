\name{spls-internal}
\alias{heatmap.spls}
\alias{spls.dv}
\alias{ust}
\alias{correctp}
\title{Internal SPLS functions}
\description{
 Internal SPLS functions.
}
\usage{
heatmap.spls( mat, coln=16, as='n', ... )
spls.dv( x, y, eta, kappa, eps, maxstep )
ust( b, eta )
correctp( x, y, eta, K, kappa, select, fit )
}
\details{
These are not to be called by the user.
}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\keyword{ internal }
