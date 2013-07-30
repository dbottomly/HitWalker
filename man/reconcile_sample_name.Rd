\name{reconcile.sample.name}
\alias{reconcile.sample.name}

\title{
Partial name matching in datasource.
}
\description{
Provides a simple mechanism to decipher a user-specified sample ID based on the specified columns and table in the \code{priorDbParams} object.
Matching is performed in the sense of \code{grep} so at least a unique substring should be specified in \code{sample.id}.
}
\usage{
reconcile.sample.name(sample.id, db.con, param.obj)
}
\arguments{
  \item{sample.id}{
	A \code{character} vector of length 1 containing the value to use for the sample ID matching.
}
  \item{db.con}{
	An object that inherits from \code{DBIConnection}, for example a \code{MySQLConnection} object.
}
  \item{param.obj}{
	A \code{priorDbParams} object.
}
}

\value{
A list containing three elements:
\describe{
\item{recon}{The reconciled sample ID, usually a character string}
\item{ind}{The unique index of the sample ID, usually an integer}
\item{num.pats}{The total number of samples in the datasource}
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{dbConnect}}
}
\examples{
if(require(RSQLite) && require(HitWalkerData))
{
  data(params)
  db.con <- dbConnect("SQLite", hitwalker.db.path())
  sample.id.list <- reconcile.sample.name("08-00102", db.con, examp.prior.param)
  dbDisconnect(db.con)
}


}
\keyword{utilities}
