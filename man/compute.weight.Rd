\name{compute.weight}
\alias{compute.weight}
\title{
  Compute the weight for the penalized likelihood in M-step.
}
\description{
  Compute the weight for the penalized likelihood in M-step.
}
\usage{
  compute.weight(tk,thres)
}
\arguments{
 
  \item{tk}{
    A n by K matrix obtained from Estep.tk.  
  }

  \item{thres}{
    A scaler specifying the threshold for the rejection step in RCEM.
 }
}
\value{
  A list including the following components: 
  \item{weight}{
    A matrix of the weight.
  }	
 \item{my.label}{
    A list indicating whether the ith functional observation is a participant in penalized
    likelihood estimation in M-step. 
  }	

}
\section{References}{
   Ma, P., Castillo-Davis, C., Zhong, W., and Liu, J. S. (2006)
   A data-driven clustering method for time course gene expression data,  
   Nucleic Acids Research, 34 (4), 1261-1269. 
}
\seealso{
  \code{\link{compute.M.pk}}
}
\examples{

}
\keyword{cluster}
% docclass is function

