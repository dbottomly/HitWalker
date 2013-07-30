\name{rwrParams-class}
\Rdversion{1.1}
\docType{class}
\alias{rwrParams-class}
\alias{rwrParams}
\alias{getConvergeMaxIt,rwrParams-method}
\alias{getConvergeMaxIt}
\alias{getConvergeThresh,rwrParams-method}
\alias{getConvergeThresh}
\alias{getRestartProb,rwrParams-method}
\alias{getRestartProb}
\alias{verbLevel,rwrParams-method}
\alias{show,rwrParams-method}
\alias{makeRwrParams}

\title{Class \code{"rwrParams"}}
\description{
A class providing control over the parameters of the random walk with restarts (RWR) algorithm
}
\section{Objects from the Class}{
\describe{
Objects can be created by calls of the form \code{makeRwrParams(c.val, threshold, maxit, verbose)}. Or equivalently through \code{new{"rwrParams", ...}}
}
}
\section{Slots}{
  \describe{
    \item{\code{c.val}:}{Object of class \code{"numeric"} containing the probability of restart.}
    \item{\code{threshold}:}{Object of class \code{"numeric"} containing the threshold to judge convergence.}
    \item{\code{maxit}:}{Object of class \code{"numeric"} containing the maximum number of iterations.}
    \item{\code{verbose}:}{Object of class \code{"logical"} specifying whether to print additional information}
  }
}
\section{Methods}{
  \describe{
    \item{getConvergeMaxIt}{\code{signature(obj = "rwrParams")}}
    \item{getConvergeThresh}{\code{signature(obj = "rwrParams")}}
    \item{getRestartProb}{\code{signature(obj = "rwrParams")}}
    \item{verbLevel}{\code{signature(obj = "rwrParams")}}
	 }
}
\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{random.walk}}
}
\examples{
rwr.param <- makeRwrParams(c.val=.3, threshold=1e-10, maxit=100, verbose=TRUE)
show(rwr.param)
}
\keyword{classes}
