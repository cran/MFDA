\name{Estep.tik}
\alias{Estep.tik}
\title{
  Compute the conditional probability of subject i in cluster k in E-step for functional mixture models.
}
\description{
  Compute the conditional probability of subject i in cluster k in the expectation step of EM algorithm for functional
  mixture models. 
}
\usage{
Estep.tik(xx.i, old)
}
\arguments{
 
  \item{xx.i}{
    A numeric vector corresponds 
    to a functional observation. 
  }
  \item{old}{
    The list of output from maximization step of EM algorithm.
  }
 
}
\value{
  A list including the following components: 
  \item{clust}{
    A vector whose \emph{i}th entry is the indicator value
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  }
  \item{t.ik}{
    A vector whose \emph{k}th entry is the conditional probability
    of the \emph{i}th observation belonging to the \emph{k}th component
    of the mixture.   
  }
  \item{loglike}{
    The negative logliklihood for the data in the mixture model. 
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

\keyword{cluster}
% docclass is function

