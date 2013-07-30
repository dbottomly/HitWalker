\name{run.prioritization.patient}
\alias{run.prioritization.patient}

\title{
Run RWR prioritization for a given sample
}
\description{
A wrapper function that runs the prioritization workflow for a given sample.
}
\usage{
run.prioritization.patient(db.con, sample.id, graph.obj, param.obj)
}

\arguments{
  \item{db.con}{
	An object inheriting from \code{DBIConnection} to connect to the datasource.
}
  \item{sample.id}{
	The sample ID that should be uniquely identifiable even if not an exact match to the datasource.
}
  \item{graph.obj}{
	An \code{annotatedIGraph} object	
}
  \item{param.obj}{
	A \code{priorDbParams} object
}
}

\value{
  \describe{
 A \code{variantPriorResult}{object containing the results of the prioritization}
 }
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{annotatedIGraph}}, \code{\linkS4class{priorDbParams}}, \code{\linkS4class{variantPriorResult}}
}
\examples{

if (require(RSQLite) && require(HitWalkerData))
{
  data(params)
  db.con <- dbConnect("SQLite", hitwalker.db.path())
  graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
  test.out <- run.prioritization.patient(db.con=db.con, sample.id="08-00102", graph.obj=graph.obj, param.obj=examp.prior.param)
  head(summary(test.out))
}

}

\keyword{utilities}
