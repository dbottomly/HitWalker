\name{retrieve.param.query}
\alias{retrieve.param.query}

\title{
Send and retrieve the result of an SQL query
}
\description{
A function that sends and retrieves queries to the specified datasources and optionally performs sanity checks in terms of returned column names.  If \code{can.fail} is TRUE, an empty \code{data.frame} will be returned if a query does not return a result.  Otherwise program execution is halted.
}
\usage{
retrieve.param.query(sample.id.list, db.con, param.obj, query.func, exp.col = NULL, can.fail=FALSE)
}

\arguments{
  \item{sample.id.list}{
	Either a \code{sample.id.list} style list or NULL.  A \code{sample.id.list} is defined as a list with three elements: \code{recon} is the human readable name, \code{ind} is the database index value and \code{num.pats} is the number of samples in the database.
}
  \item{db.con}{
	An object inherited from \code{DBIConnection} that provides access to the datasource.
}
  \item{param.obj}{
	A \code{priorDbParams} object.
}
  \item{query.func}{
	A method generating the function that generates the SQL query to be sent.  Needs to take the \code{sample.id.list} and \code{priorDbParams} object as arguments.
}
  \item{exp.col}{
	A \code{character} vector containing the expected columns or NULL if none.
}
  \item{can.fail}{
	A \code{logical} value.
}
}
\value{
\describe{
A \code{data.frame}{containing the results of the query.}
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{dbGetQuery}}
}
\examples{

if(require(RSQLite) && require(HitWalkerData))
{
  data(params)
  db.con <- dbConnect("SQLite", hitwalker.db.path())
  sample.id.list <- reconcile.sample.name("08-00102", db.con, examp.prior.param)
  query.func <- hitQuery
  query.res <- retrieve.param.query(sample.id.list, db.con, examp.prior.param, query.func, exp.col = NULL, can.fail=FALSE)
  head(query.res)
}


}
\keyword{utilities}
