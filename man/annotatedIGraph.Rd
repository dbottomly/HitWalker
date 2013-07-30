\name{annotatedIGraph-class}
\Rdversion{1.1}
\docType{class}
\alias{annotatedIGraph-class}
\alias{annotatedIGraph}
\alias{annotation,annotatedIGraph-method}
\alias{annotation<-,annotatedIGraph-method}
\alias{annotation}
\alias{annotation<-}
\alias{getGraph,annotatedIGraph-method}
\alias{getGraph}
\alias{saveGraph,annotatedIGraph-method}
\alias{saveGraph}
\alias{show,annotatedIGraph-method}

\title{Class \code{"annotatedIGraph"}}
\description{
This class encapsulates the \code{igraph} class from the igraph package providing a \code{data.frame} of annotation for the nodes that can be accessed and modified in a manner consistent with the rest of the package.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("annotatedIGraph", graph, annotation)}.
}

\section{Slots}{
  \describe{
    \item{\code{graph}:}{Object of class \code{"igraph"}}
    \item{\code{annotation}:}{Object of class \code{"data.frame"} with one row per node indicating the annotation for the graph.  For instance in a PPI graph each node would be mapped to a protein as well as a gene ID.  Also contained in this object would be graph-specific metadata which would be used in downstream functions for edge/node drawing characteristics.}
  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(obj = "annotatedIGraph")}: Retrieve the annotation \code{data.frame} }
    \item{annotation<-}{\code{signature(obj = "annotatedIGraph")}: Replace the annotation \code{data.frame} }
    \item{getGraph}{\code{signature(obj = "annotatedIGraph")}: Retrieve the \code{igraph} object}
    \item{saveGraph}{\code{signature(obj = "annotatedIGraph")}: Save an \code{annotatedIGraph} object to disk }
	 }
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{loadGraph}}, \code{\link{igraph}}, \code{\link{get.protein.protein.graph}}
}
\examples{
stopifnot(require(igraph))
set.seed(123)
sim.graph <- erdos.renyi.game(n=10, p.or.m=.5, type="gnp", directed=FALSE, loops=FALSE)
sim.dta <- data.frame(node.name=paste("node", 1:10, sep="_"), colors=sample(colors(), 10))

annot.graph <- new("annotatedIGraph", graph=sim.graph, annotation=sim.dta)
\dontrun{plot(getGraph(annot.graph))}
head(annotation(annot.graph))

}
\keyword{classes}
