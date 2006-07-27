\name{compute.center}
\alias{compute.center}
\title{
  Update parameters for cluster k in M-step for functional mixture models.
}
\description{
  Update  and trace of the smoothing matrix for cluster k in the M step of 
  EM algorithm for functional mixture models. 
}
\usage{
compute.center(x,weight)
}
\arguments{
 
  \item{x}{
    A matrix of obervations with each row representing the observation and each column representing the experiment at each time point. 
  }
  \item{weight}{
    A matrix with each element is the condtional posterior probability that functional observation belongs to cluster k.  
  }
 
}
\value{
  A list including the following components: 
  \item{mu}{
    A vector whose \emph{k}th row represent the center mean curve for the kth cluster 
  }
  \item{zeta}{
     A vector with each element representing the random effect $b_i$ for the kth cluster.
  }
  \item{varht}{
    A vector with each element representing the estimated error variance.
  }
  \item{trc}{
    A vector with each element representing the trace of the smoothing matrix for the kth cluster.		
}

}

\section{References}{
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{compute.center}}
}
\examples{

}
\keyword{cluster}

