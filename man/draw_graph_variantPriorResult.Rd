\name{draw.graph.variantPriorResult}
\alias{draw.graph.variantPriorResult}

\title{
Plot a \code{variantPriorResult} object
}
\description{
Plot a subset of a \code{variantPriorResult} object based on the shortest paths between hits and variants.  This function wraps the workhorse \code{make.target.graph} and \code{draw.graph} functions.
}
\usage{
draw.graph.variantPriorResult(var.obj, graph.params)
}

\arguments{
  \item{var.obj}{
	A \code{variantPriorResult} object
}
  \item{graph.params}{
	A \code{graphDispParams} object
}
}
\details{
For the top ranked variants and hits, the shortest path is found between the hits and variants as well as between the hits and between the variants.
Currently only the 1st order interactions are kept for variant-variant paths and the shortest paths are unweighted.  See \code{make.target.graph} for details.
Graphs are drawn using the \code{Rgraphviz} package in conjunction with the \code{graph} package, see \code{draw.graph} for details.
}
\value{
None
}
\author{
Daniel Bottomly
}

\seealso{
\code{\link{draw.graph}},\code{\link{make.target.graph}}, \code{\linkS4class{variantPriorResult}}, \code{\linkS4class{graphDispParams}}
}
\examples{
set.seed(123)
data(prior_obj)
graph.params <- makeGraphDispParams(file.name=character())
draw.graph.variantPriorResult(prior.obj, graph.params)
}

\keyword{dplot}
