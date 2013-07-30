\name{get.hit.prots}
\alias{get.hit.prots}

\title{
Retrieve functional assay hits from datasource
}
\description{
This function retrieves the hit data from the specified connection and summarizes the scores as specified in the function retrieved through the \code{getScoreSummaryFunc} method of \code{priorDbParams}.
}
\usage{
get.hit.prots(sample.id.list, db.con, param.obj)
}

\arguments{
  \item{sample.id.list}{
	A list with three elements: recon is the human readable name, ind is the database index value and num.pats is the number of samples in the database.
}
  \item{db.con}{
	An object inherited from \code{DBIConnection} which represents the datasource
}
  \item{param.obj}{
	A \code{priorDbParams} object
}
}

\note{
  The score summary function should take a \code{data.frame} as the argument and return a modified \code{data.frame} containing an additional column named as in the \code{getCombScoreName} method of \code{priorDbParams}.	
}

\value{
A list with three elements:
\describe{
\item{protein.scores}{A vector containing scores for each of the specified summary IDs values with corresponding names.}
\item{query.dta}{The \code{data.frame} retrieved from the hit query with the additional combined hit score column.}
\item{summary.id}{The summary ID column to use for hit data downstream.}
}
}
\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}
}
\examples{
set.seed(123)
if (require(RSQLite) && require(HitWalkerData))
{
	data(params)
	db.con <- dbConnect("SQLite", hitwalker.db.path())
	sample.id.list <- reconcile.sample.name("08-00102", db.con, examp.prior.param)
	hit.list <- get.hit.prots(sample.id.list, db.con, examp.prior.param)
	dbDisconnect(db.con)
}
}

\keyword{methods}
