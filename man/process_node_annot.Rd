\name{process.node.annot}
\alias{process.node.annot}

\title{
Apply a function to split annotation
}
\description{
Applies a specified function (in the \code{graphDispParams} object) to each section of the annotation \code{data.frame} after being split by the summary ID.  This facilitates
mapping metadata supplied in the annotation \code{data.frame} to nodes.
}
\usage{
process.node.annot(obj, targ.nodes, graph.params, func.name)
}

\arguments{
  \item{obj}{
	A \code{variantPriorResult} object
}
  \item{targ.nodes}{
	The node names of the graph to be plotted.  Should correspond to the summary IDs.
}
  \item{graph.params}{
	A \code{graphDispParams} object
}
  \item{func.name}{
	The name of the function-calling method to apply, the resulting function should take a \code{data.frame} as the sole argument
}
}

\value{
\describe{
\item{process.node.annot}{Returns a \code{character} vector named by the appropriate summary ID}
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{variantPriorResult}}, \code{\linkS4class{graphDispParams}}, \code{\linkS4class{annotatedIGraph}}, \code{\link{make.target.graph}}
}
\examples{
set.seed(123)
data(prior_obj, package="HitWalker")

test.symbs <- sample(annotation(getGraph(prior.obj))$symbol, 10)

graph.params <- makeGraphDispParams(file.name=character())

process.node.annot(obj=prior.obj, targ.nodes=test.symbs, graph.params=graph.params, func.name=shapeFunc)

}
\keyword{utilities}
