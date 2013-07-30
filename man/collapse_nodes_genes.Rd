\name{collapse.nodes.genes}
\alias{collapse.nodes.genes}

\title{
Collapse nodes into a single representative nodes.
}
\description{
Rename the nodes of a \code{graphNEL} object according to a \code{character} vector named by current node names with elements indicating the new node names.  Duplicates are
summarized to be a single node containing all edges and indicated by the maximum edge weights
}
\usage{
collapse.nodes.genes(targ.graph, node.names)
}

\arguments{
  \item{targ.graph}{
	A \code{graphNEL} object.
}
  \item{node.names}{
	A named \code{character} vector.
}
}

\value{
A \code{graphNEL} object.
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{graphNEL}}
}
\examples{
set.seed(123)
rand.graph <- randomGraph(V=sample(colors()[1:10]), M=1:4, p=.4)

node.vec <- nodes(rand.graph)
names(node.vec) <- node.vec

node.vec[sample.int(length(node.vec), 3)] <- "hello"

new.graph <- collapse.nodes.genes(rand.graph, node.vec)

edgeWeights(new.graph)["hello"]
edgeWeights(rand.graph)[names(node.vec)[node.vec == "hello"]]

if (require(Rgraphviz))
{
	par(mfcol=c(2,1))
	plot(rand.graph)
	plot(new.graph)
}
}
\keyword{manip}
