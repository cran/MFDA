\name{Estep.tk}
\alias{Estep.tk}
\title{
  Compute the cluster proportion  in E-step for functional mixture models.
}
\description{
  Compute the cluster proportion in the expectation step of EM algorithm for functional
  mixture models. 
}
\usage{
Estep.tk(x, old, nclust)
}
\arguments{
 
  \item{x}{
    A numeric matrix, or data frame of observations, rows
    correspond to functional observations and columns correspond 
    to the number of repeated measurements. 
    
  }
  \item{old}{
    The list of output from maximization step of EM algorithm.
  }
  \item{nclust}{
    The number of clusters (componenets in mixture model).
  }
}
\value{
  A list including the following components: 
\item{clust}{
     A vector whose \emph{i}th entry is the indicator value 
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  }
  \item{tk}{
    A matrix whose \code{[i,k]}th entry is the conditional probability
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  } 
  \item{loglikelihood}{
    The negative logliklihood for the data in the mixture model. 
  }

}

\section{References}{
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{Estep.tik}}
}
\examples{

}
\keyword{cluster}
% docclass is function

