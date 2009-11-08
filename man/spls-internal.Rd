\name{spls-internal}
\alias{heatmap.spls}
\alias{spls.dv}
\alias{ust}
\alias{correctp}
\alias{cv.split}
\alias{wpls}
\alias{sgpls.binary}
\alias{sgpls.multi}
\alias{cv.sgpls.binary}
\alias{cv.sgpls.multi}
\title{Internal SPLS functions}
\description{
 Internal SPLS functions.
}
\usage{
heatmap.spls( mat, coln=16, as='n', ... )
spls.dv( Z, eta, kappa, eps, maxstep )
ust( b, eta )
correctp( x, y, eta, K, kappa, select, fit )
cv.split( y, fold )
wpls( x, y, V, K=ncol(x), type="pls1",
        center.x=TRUE, scale.x=FALSE )
sgpls.binary( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
sgpls.multi( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
cv.sgpls.binary( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
    br=TRUE, ftype='iden' )
cv.sgpls.multi( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
    br=TRUE, ftype='iden' )
}
\details{
These are not to be called by the user.
}
\references{
Chung, D. and Keles, S. (2009). "Sparse partial least squares classification
  for high dimensional data"
(\url{http://www.stat.wisc.edu/~keles/Papers/C_SPLS.pdf}).
  
Chun, H. and Keles, S. (2009). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
  To appear in \emph{Journal of the Royal Statistical Society - Series B}
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\keyword{ internal }
