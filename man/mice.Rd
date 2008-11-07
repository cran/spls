\name{mice}
\docType{data}
\alias{mice}
\title{Mice Dataset}
\description{
  This is the Mice dataset used in Chun and Keles (2008).
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
  The Mice dataset was published by Lan et al. (2006). The \code{x} matrix is
  the marker map consisting of 145 microsatellite markers from 19 non-sex mouse chromosomes. 
  The \code{y} matrix is the gene expression measurements of the 83 transcripts
  from liver tissues of 60 mice. This group of the 83 transcripts is one of the clusters
  analyzed by Chun and Keles (2008). See Chun and Keles (2008) for more details.                                  
}
\source{
Lan, H., M. Chen, J. B. Flowers, B. S. Yandell, D. S. Stapleton, C. M. Mata,
E. T-K Mui, M. T. Flowers, K. L. Schueler, K. F. Manly, R. W. Williams,
C. Kendziorski, and A. D. Attie (2006). "Combined expression
trait correlations and expression quantitative trait locus mapping",
\emph{PLoS Genetics}, 2, e6.}
\references{
Chun, H. and Keles, S. (2008). "Expression quantitative trait loci mapping 
with multivariate sparse partial least squares regression",
(\url{http://www.stat.wisc.edu/~keles/Papers/chun_keles_eQTL_SPLS_submit.pdf}).
}
\examples{
data(mice)
mice$x[1:5,1:5]
mice$y[1:5,1:5]
}
\keyword{datasets}
