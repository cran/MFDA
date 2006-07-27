\name{compute.reject}
\alias{compute.reject}
\title{
  The rejection controlled EM.
}
\description{
  Algorithm to estimate the weight in M-step.  
}
\usage{
compute.reject(tk,k,thres)
}
\arguments{ 
  \item{tk}{
    A matrix obtained from E-step. 
  }
  \item{k}{
    The kth Gaussian component.
  }
  \item{thres}{ 
   The threshold for  rejection-controlled EM. 		
  }
 
}
\value{
  A list including the following components: 
  \item{em.select}{
    A vector.   
  }
  \item{sem.select}{
    A vector.
 }
}

\section{References}{
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
}
\examples{

}
\keyword{cluster}

