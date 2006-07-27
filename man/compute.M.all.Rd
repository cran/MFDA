\name{compute.M.all}
\alias{compute.M.all}
\title{
  M-step in the EM algorithm for functional mixture models.
}
\description{
  Maximization step in the EM algorithm for functional mixture models.
}
\usage{
compute.M.all(x, weight, my.label)
}
\arguments{

  \item{x}{
    A matrix, or data frame of observations, rows correspond to observations and
    columns correspond to variables.
  }
  \item{weight}{
    A matrix whose \code{[i,k]}th entry is the
    conditional probability of the ith observation belonging to
    the \emph{k}th component of the mixture.  
  }
  \item{my.label}{
    A list whose \emph{k}th entry is the
    indicator value indicating whether the ith functional observation is 
    participant of M-step for cluster k. 
  }
  
  
}
\value{
  A list including the following components: 
  
  \item{mu}{
    A matrix whose kth column is the
    mean of the \emph{k}th component of the mixture model.
  }
  \item{zeta}{
    A numerical vector specifying the estimate of cluster precision parameters. 
  }
  \item{varht}{
    A numerical vector specifying the estimate of cluster error variance. 
  }
  \item{mu}{
    A matrix specifying the estimate of cluster mean profile. 
  }
  \item{zeta}{
    A numerical vector specifying the estimate of cluster precision parameters. 
  }
  \item{trc}{
    A numerical vector whose kth entry specifying the trace of smoothing matrix of cluster k. 
  }
}
\section{References}{
    Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{Estep.tk}}
}
\examples{
}
\keyword{cluster}
% docclass is function
