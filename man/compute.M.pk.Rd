\name{compute.M.pk}
\alias{compute.M.pk}
\title{
  Compute the cluster proportion.
}
\description{
  Compute cluster proportion in Gaussian mixture model.
}
\usage{
  compute.M.pk(tk, nclust, alpha)
}
\arguments{
 
  \item{tk}{
    An n by k matrix obtained from Estep.tk.  
  }
  \item{nclust}{
   The number of cluster.
  }
  \item{alpha}{
   The parameters of a Dirichilet distribution.
 }
}
\value{
  A list including the following components: 
  \item{p_k}{
    A vector of updated cluster proportion.
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


