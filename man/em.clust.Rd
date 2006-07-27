\name{em.clust}
\alias{em.clust}
\title{
  Model-Based Clustering for a single Markov chain
}
\description{
  EM algorithm
  for functional Gaussian mixture models.
}
\usage{
em.clust(x,clust,mu,zeta,varht,thres,iter.max,alpha)
}
\arguments{
  \item{x}{
    A numeric matrix, or data frame of observations, rows
    correspond to functional observations and columns correspond 
    to the number of repeated measurements. 
  }
  \item{clust}{
    An integer vector specifying the initial clustering membership. 
  }
  \item{mu}{
    A matrix specifying the initial estimate of cluster mean profile. 
  }
  \item{zeta}{
    A numerical vector specifying the initial estimate of cluster precision parameters. 
  }
  \item{varht}{
    A numerical vector specifying the initial estimate of cluster error variance. 
  }
  \item{thres}{
    A threshold value for rejection-controlled step.  
  }

  \item{iter.max}{
    An integer limit on the number of EM iterations. 
  }
  \item{alpha}{
    The prior for the cluster proportions p.
  }

}
\value{
  A list contains estimated cluster membership, estimated cluster mean profile,
  Bayesian Information Criterion, negative loglikelihood.
  
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
% docclass is function
