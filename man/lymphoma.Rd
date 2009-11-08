\name{lymphoma}
\docType{data}
\alias{lymphoma}
\title{Lymphoma Gene Expression Dataset}
\description{
  This is the Lymphoma Gene Expression dataset used in Chung and Keles (2009).
}
\usage{ data(lymphoma) }
\format{
  A list with two components:
  \describe{
  \item{x}{ Gene expression data. A matrix with 62 rows and 4026 columns.}
  \item{y}{ Class index. A vector with 62 elements.}
  }
}
\details{
  The lymphoma dataset consists of 42 samples of diffuse large B-cell lymphoma (DLBCL),
  9 samples of follicular lymphoma (FL),
  and 11 samples of chronic lymphocytic leukemia (CLL).  
  DBLCL, FL, and CLL classes are coded in 0, 1, and 2, respectively, in \code{y} vector.
  Matrix \code{x} is gene expression data and
  arrays were normalized, imputed, log transformed, and standardized 
  to zero mean and unit variance across genes as described
  in Dettling (2004) and Dettling and Beuhlmann (2002).
  See Chung and Keles (2009) for more details. 
}
\source{
Alizadeh, A., Eisen, M. B., Davis, R. E., Ma, C., Lossos, I. S., Rosenwald, A., Boldrick, J. C.,
Sabet, H., Tran, T., Yu, X., Powell, J. I., Yang, L., Marti, G. E., Moore, T., Jr., J. H., Lu,
L., Lewis, D. B., Tibshirani, R., Sherlock, G., Chan, W. C., Greiner, T. C., Weisenburger,
D. D., Armitage, J. O., Warnke, R., Levy, R., Wilson, W., Grever, M. R., Byrd, J. C.,
Botstein, D., Brown, P. O., , and Staudt, L. M. (2000). "Distinct types of diffuse large
B-cell lymphoma identified by gene expression profiling", \emph{Nature}, 403, pp. 503--511.
}
\references{
Chung, D. and Keles, S. (2009). "Sparse partial least squares classification
  for high dimensional data"
(\url{http://www.stat.wisc.edu/~keles/Papers/C_SPLS.pdf}). 
  
Dettling, M. (2004). "BagBoosting for tumor classification with gene expression data",
\emph{Bioinformatics}, 20, pp. 3583--3593.

Dettling, M. and Beuhlmann, P. (2002). "Supervised clustering of genes", \emph{Genome Biology},
3, pp. research0069.1--0069.15.
}
\examples{
data(lymphoma)
lymphoma$x[1:5,1:5]
lymphoma$y
}
\keyword{datasets}
