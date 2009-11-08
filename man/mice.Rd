\name{mice}
\docType{data}
\alias{mice}
\title{Mice Dataset}
\description{
  This is the Mice dataset used in Chun and Keles (2009).
}
\usage{ data(mice) }
\format{
  A list with two components:
  \describe{
  \item{x}{ Marker map data. A matrix with 60 rows and 145 columns.}
  \item{y}{ Gene expression data. A matrix with 60 rows and 83 columns.}
  }
}
\details{
  The Mice dataset was published by Lan et al. (2006). Matrix \code{x} is
  the marker map consisting of 145 microsatellite markers from 19 non-sex mouse chromosomes. 
  Matrix \code{y} is gene expression measurements of the 83 transcripts
  from liver tissues of 60 mice. This group of the 83 transcripts is one of the clusters
  analyzed by Chun and Keles (2009). See Chun and Keles (2009) for more details.                                  
}
\source{
Lan, H., M. Chen, J. B. Flowers, B. S. Yandell, D. S. Stapleton, C. M. Mata,
E. T-K Mui, M. T. Flowers, K. L. Schueler, K. F. Manly, R. W. Williams,
C. Kendziorski, and A. D. Attie (2006). "Combined expression
trait correlations and expression quantitative trait locus mapping",
\emph{PLoS Genetics}, 2, e6.}
\references{
Chun, H. and Keles, S. (2009). "Expression quantitative trait loci mapping 
with multivariate sparse partial least squares regression", \emph{Genetics}, 182, pp. 79--90.
}
\examples{
data(mice)
mice$x[1:5,1:5]
mice$y[1:5,1:5]
}
\keyword{datasets}
