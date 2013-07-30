\name{make.target.graph}
\alias{make.target.graph}
\alias{get.graph.sp}
\title{
Creating shortest path subgraphs
}
\description{
The \code{make.target.graph} function uses \code{get.graph.sp} to create subgraphs corresponding to the shortest paths between hits or between variants and hits or between variants.  These subgraphs are subsequently joined together to form the graph structure that will be plotted.
}
\usage{
make.target.graph(use.graph, hit.nodes, mut.nodes)
get.graph.sp(use.graph, from.nodes, to.nodes, path.len.restrict)
}
\arguments{

  \item{use.graph}{
	An \code{igraph} object.
}
  \item{hit.nodes}{
	A \code{character} vector of nodes corresponding to hits.
}
  \item{mut.nodes}{
	A \code{character} vector of nodes corresponding to variant nodes.
}
  \item{from.nodes}{
	A \code{character} vector of nodes corresponding to the starting point of the shortest path.
}

  \item{to.nodes}{
	A \code{character} vector of nodes corresponding to the ending point of the shortest path (or NULL).
}

  \item{path.len.restrict}{
	Maximum number of nodes on the shortest path.
}

}

\value{
\describe{
\item{make.target.graph}{returns a \code{graphNEL} object.}
\item{get.graph.sp}{returns a matrix containing two columns indicating the nodes that have an edge.}
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\link{igraph}}, \code{\linkS4class{graphNEL}}
}
\examples{
set.seed(123)

if (require(RSQLite) && require(HitWalkerData))
{
  data(params)
  graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
  graph.params <- makeGraphDispParams(file.name=character())
  nodes <- annotation(graph.obj)$protein
  samp.nodes <- sample(nodes, 6)
  res.graph <- make.target.graph(use.graph=getGraph(graph.obj), hit.nodes=samp.nodes[1:3], mut.nodes=samp.nodes[4:6])
  
  if (require(Rgraphviz)) plot(res.graph)
  
  res.graph.sp <- get.graph.sp(getGraph(graph.obj), from.nodes=samp.nodes[1:3], to.nodes=samp.nodes[4:6], path.len.restrict=NULL)
  trans.graph.sp <- matrix(V(getGraph(graph.obj))$name[as.integer(res.graph.sp)], ncol=2, nrow=nrow(res.graph.sp))
  head(trans.graph.sp)
}
}


\keyword{utilities}
