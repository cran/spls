\name{yeast}
\docType{data}
\alias{yeast}
\title{Yeast Cell Cycle Dataset}
\description{
  This is the Yeast Cell Cycle dataset used in Chun and Keles (2007).
}
\usage{ data(yeast) }
\format{
  A list with two components:
  \describe{
  \item{x}{ ChIP-chip data. A matrix with 542 rows and 106 columns.}
  \item{y}{ Cell cycle gene expression data.
                A matrix with 542 rows and 18 columns.}
  }
}
\details{
  The \code{y} matrix is the cell cycle gene expression data (Spellman et al., 1998)
  of 542 genes from an \eqn{\alpha} factor based experiment.
  Each column corresponds to mRNA levels
  measured at every 7 minutes during 119 minutes (a total of 18 measurements).
  The \code{x} matrix is the chromatin immunoprecipitation on chip (ChIP-chip) data of
  Lee et al. (2002) and it contains the binding information for 106
  transcription factors. See Chun and Keles (2007) for more details. 
}
\source{
Lee, T. I., N. J. Rinaldi, F. Robert, D. T. Odom, Z. Bar-Joseph, G. K. Gerber, N. M. Hannett,
C. T. Harbison, C. M. Thomson, I. Simon, J. Zeitlinger, E. G. Jennings, H. L. Murray, D.
B. Gordon, B. Ren, J. J. Wyrick, J.-B. Tagne, T. L. Volkert, E. Fraenkel, D. K. Gifford,
and R. A. Young (2002). "Transcriptional regulatory networks in \emph{Saccharomyces cerevisiae}",
\emph{Science}, 298, pp. 799--804.

Spellman, P. T., G. Sherlock, M. Q. Zhang, V. R. Iyer, K. Anders, M. B. Eisen, P. O. Brown, D.
Botstein, and B. Futcher (1998). "Comprehensive identification of cell cycle-regulated genes of
the yeast \emph{Saccharomyces cerevisiae} by microarray hydrization",
\emph{Molecular Biology of the Cell}, 9, pp. 3273--3279.
}
\references{
Chun, H. and Keles, S. (2007). "Sparse partial least squares
  for simultaneous dimension reduction and variable selection",
(\url{http://www.stat.wisc.edu/~keles/Papers/SPLS_Nov07.pdf}).
}
\examples{
data(yeast)
yeast$x[1:5,1:5]
yeast$y[1:5,1:5]
}
\keyword{datasets}
