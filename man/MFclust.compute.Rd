\name{MFclust.compute}
\alias{MFclust.compute}
\title{
  Computational Component in Model-Based Functional Data Clustering
}
\description{
  The computation componenet for clustering via rejection-controlled 
  EM initialized by kmeans clustering for
  functional Gaussian mixture models. The number of clusters and the
  clustering model is chosen to minimize the BIC. 
}
\usage{
MFclust.compute(data, minG, maxG, nchain,thres,iter.max,my.alpha, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to time.
  }
  \item{minG}{
    An integer vector specifying the minimum number of mixture components 
    (clusters) to be considered.
    The default is \code{1} component.
  }
  \item{maxG}{
    An integer vector specifying the maximum number of mixture components 
    (clusters) to be considered.
    The default is \code{9} components.
  }
    \item{nchain}{
    An integer specifying the number of Markov chains in RCEM
    The default is \code{5} chains.
  }
   \item{thres}{
    A number between 0 and 1 specifying the threshold value of rejection step in RCEM.
    The default is \code{0.5}.
  }
    \item{iter.max}{
    An iteger specifying the maximum number of iteration in RCEM.
    The default is \code{10}.
  }
  \item{my.alpha}{
    The prior for the cluster proportion.
    The default is \code{1}.
  }
  \item{\dots}{
    The arguments to be part of the function.
  }
}
\value{
  A list representing the best model (according to BIC) for the given range
  of numbers of clusters. 

}

\section{References}{
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{MFclust}}
}
\examples{

}
\keyword{cluster}

