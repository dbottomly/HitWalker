\name{get.protein.protein.graph}
\alias{get.protein.protein.graph}

\title{
Retrieve a graph from a datasource.
}
\description{
Retrieves a graph from a datasource converting it to an \code{igraph} object as well as retrieving accompanying annotation.
}
\usage{
get.protein.protein.graph(db.con, param.obj, get.annotation = TRUE, on.symbol.protein.dup="remove.dup")
}

\arguments{
  \item{db.con}{
	An object inheriting from \code{DBIConnection} (e.g. MySQLConnection)
}
  \item{param.obj}{
	A \code{priorDbParams} object
}
  \item{get.annotation}{
	A \code{logical} value specifying whether to retrieve annotation
}
  \item{on.symbol.protein.dup}{
	This controls how duplicated proteins are handled in this program.  Currently, the default and only (but not optimal) option is "remove.dup".  This
	removes the duplicate entries keeping only one for each protein.
}
}
\details{
The graph has to be in form usable by \code{graph.data.frame} function in the \code{igraph} package.  The graph will be retrieved using the \code{getMatrixQuery} method of \code{priorDbParams}.
If \code{get.annotation} is TRUE, the \code{getAnnotQuery} method will be used to retrieve annotation.  If \code{nodeWeightCol} is specified, then the specified values will be multiplied
with the entire set of edges for each gene and the edge scores will be modified accordingly.

}

\value{
An \code{annotatedIGraph} object
}

\author{
Daniel Bottomly
}
\note{
Edge re-weighting through a specified \code{nodeWeightCol} should be considered highly experimental and is subject to change.
}


\seealso{
\code{\link{graph.data.frame}}, \code{\linkS4class{priorDbParams}}, \code{\link{dbConnect}}
}
\examples{

if (require(HitWalkerData) && require(RSQLite))
{
  data(params)
  db.con <- dbConnect("SQLite", annot.db.path())
  annotQuery(examp.prior.param) <- function(sample.id.list, param.obj){"SELECT * FROM annotation"}
  graph.obj <- get.protein.protein.graph(db.con, examp.prior.param, get.annotation=TRUE)
  \dontrun{saveGraph(graph.obj, file="graph.RData")}
}

}

\keyword{utilities}
