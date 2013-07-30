\name{loadGraph}
\alias{loadGraph}
\title{
Load an annotatedIGraph object
}
\description{
Load an \code{annotatedIGraph} object and verify that it is consistent with the supplied \code{priorDbParams} object.
}
\usage{
loadGraph(file.name, param.obj)
}
\arguments{
 \item{file.name}{
	A \code{character} vector of length 1 indicating the file containing the \code{annotatedIGraph} object created by the \code{saveGraph} method
}
 \item{param.obj}{
	A \code{priorDbParams} object
}
}
\value{
\describe{
An \code{annotatedIGraph} object
}
}
\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{igraph}}, \code{\link{get.protein.protein.graph}}
}

\examples{
set.seed(123)

if (require(RSQLite) && require(HitWalkerData))
{
	data(params)
	graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
}


}
\keyword{functions}
