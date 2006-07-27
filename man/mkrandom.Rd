\name{mkrandom}
\alias{mkrandom}
\title{
    Generating Random Effects for cluster model
}
\description{
    Generate entries representing random effect matrix 
}
\usage{
mkrandom(formula, data)
}
\arguments{
    \item{formula}{Symbolic description of the random effects.}
    \item{data}{Data frame.}
}
\details{
    This function generates random effects terms from simple grouping
    variables, for use in nonparametric mixed-effect models as described
    in Gu and Ma (2003a, b).
}
\value{
    A list of three components.
    \item{z}{Z matrix.}
    \item{sigma}{Sigma matrix to be evaluated through
        \code{sigma$fun(para,sigma$env)}.}
    \item{init}{Initial parameter values.}
}
\note{
    This program intrinsic function for MFclust. User cannot specify the parameters.
}
\author{pingma@uiuc.edu}
\references{}
\examples{}
\keyword{models}
\keyword{regression}
