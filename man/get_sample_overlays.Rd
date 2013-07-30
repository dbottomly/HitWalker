\name{get.sample.overlays}
\alias{get.sample.overlays}
\title{
Retrieve sample-specific metadata
}
\description{
This function retrieves sample-level gene annotation (e.g. differential expression, methylation status) and returns a revised \code{annotatedIGraph} object with the
sample-specific data merged in with the existing annotation.
}
\usage{
get.sample.overlays(db.con, param.obj, graph.obj, sample.id.list)
}
\arguments{
  \item{db.con}{
	An object inheriting from \code{DBIConnection} representing the datasource
}
  \item{param.obj}{
	A \code{priorDbParams} object
}
  \item{graph.obj}{
	An \code{annotatedIGraph} object
}
  \item{sample.id.list}{
	A list with three elements: 'recon' is the human readable name, 'ind' is the database index value and 'num.pats' is the number of samples in the database.
}
}

\value{
An \code{annotatedIGraph}{containing the sample-specific data merged in with the exising annotation supplied in \code{graph.obj}.}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{annotatedIGraph}}, \code{\linkS4class{priorDbParams}}
}
\examples{
if (require(RSQLite) && require(HitWalkerData))
{
	set.seed(123)
	data(params)
	db.con <- dbConnect("SQLite", hitwalker.db.path())
	graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
	sample.id.list <- reconcile.sample.name("08-00102", db.con, examp.prior.param)
	new.graph.obj <- get.sample.overlays(db.con, examp.prior.param, graph.obj, sample.id.list)
	##Note that the graph.obj will look the same in this example.  See the Hitwalker_Add_Metadata vignette for an actual example
}
}

\keyword{methods}
