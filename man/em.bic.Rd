\name{em.bic}
\alias{em.bic}
\title{
  BIC for Functional Mixture Gaussian Models
}
\description{
  Compute the BIC (Bayesian Information Criterion) for functional
  mixture Gaussian models  given the negative loglikelihood, the dimension of the data,
  and the trace of the smoothing matrix. 
}
\usage{
em.bic(likelihood, trc, n, m)
}
\arguments{  
  \item{likelihood}{
    The negative loglikelihood for a data set with respect to the functional mixture model.
  }
    \item{trc}{
    The trace of the smoothing matrix.
  } 
 
  \item{n}{
    The number of the functional data use to compute
    \code{loglik}.
  }
  \item{m}{
    The number of the repeated measurements in the data used to compute \code{loglik}.
  }
}
\value{
  The BIC or Bayesian Information Criterion for the given input arguments.
}
\section{References}{ 
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{em.clust}}.
}
\examples{


}
\keyword{cluster}
