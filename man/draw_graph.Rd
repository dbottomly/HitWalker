\name{draw.graph}
\alias{draw.graph}

\title{
Low-level wrapper for drawing a graph based on supplied metadata. 
}
\description{
A wrapper function for processing node-specific metadata and drawing the resulting graph.
}
\usage{
draw.graph(targ.graph, var.obj, node.names, graph.params)
}

\arguments{

  \item{targ.graph}{
	A \code{graphNEL} object
}
  \item{var.obj}{
	A \code{variantPriorResult} object
}
  \item{node.names}{
	A named \code{character} vector that maps the original graph node names to the summary IDs
}
  \item{graph.params}{
	A \code{graphDispParams} object
}
}

\value{
None
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{graphDispParams}}, \code{\linkS4class{variantPriorResult}}, \code{\linkS4class{graphNEL}}
}
\examples{
	
	set.seed(123)
	if (require(HitWalkerData) && require(Rgraphviz))
	{
	  data(prior_obj)
	  param.obj <- prior.obj@param.obj
	  graph.params <- makeGraphDispParams(file.name=character())
	  graph.obj <- loadGraph(graph.file.path(), param.obj)
	  use.graph <- getGraph(graph.obj)
  
	  ##make an induced subgraph of the shortest path between two proteins containing variants for this sample. 	
	  annot <- annotation(graph.obj)
	  node.names <- sapply(split(annot, annot$protein), "[[", "symbol")
	  prots.in.graph <- intersect(getVarDta(prior.obj)$protein, annotation(graph.obj)$protein)
	  var.graph <- induced.subgraph(use.graph, get.shortest.paths(use.graph, V(use.graph)[prots.in.graph[1]], V(use.graph)[prots.in.graph[2]])[[1]])
	  
	  ##make a version of graph.nel without duplicated edges and gene symbol names
	  graph.nel <- igraph.to.graphNEL(var.graph)
	  edge.l <- lapply(edgeWeights(graph.nel), function(x){new.x=x[!duplicated(names(x))];list(edges=as.character(node.names[names(new.x)]), weights=as.numeric(new.x))})
	  names(edge.l) <- as.character(node.names[names(edge.l)])
	  graph.nel <- new("graphNEL", nodes=names(edge.l), edgeL=edge.l, edgemode="undirected")
	  
	  draw.graph(graph.nel, prior.obj,node.names ,graph.params)
	}
}

\keyword{utilities}
