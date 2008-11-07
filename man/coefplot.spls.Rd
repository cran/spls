\name{coefplot.spls}
\alias{coefplot.spls}
\title{Plot the estimated coefficients of the SPLS object}
\description{
 Plot the estimated coefficients of the selected predictors in the SPLS object.
}
\usage{
coefplot.spls( object, nwin=c(2,2),
            xvar=c(1:length(object$A)), ylimit=NA )
}
\arguments{
  \item{object}{ A fitted SPLS object. }
  \item{nwin}{ Vector of the numbers of rows and columns in a plotting area.
                Default is two rows and two columns, i.e. four plots.}
  \item{xvar}{ Index of the variables to be plotted among the set of the selected predictors.
                Default is to plot the coefficients of all the selected predictors.}
  \item{ylimit}{ Range of the y axis (the coefficients) in the plot.
                If \code{ylimit} is not specified, the y axis of the plot has the range
                between the minimum and the maximum of all coefficient estimates. }
}
\details{
This plot is useful for visualizing the coefficient estimates of a variable
for different responses. Hence, the function is applicable
only with multivariate response SPLS.
}
\value{NULL.}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}). }
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ plot, ci.spls, and correct.spls methods for spls. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
# Draw the estimated coefficient plot of the first four variables
# among the selected predictors
coefplot.spls( f, xvar=c(1:4), nwin=c(2,2) )
}
\keyword{ hplot }
