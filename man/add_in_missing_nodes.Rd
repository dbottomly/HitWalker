\name{add.in.missing.nodes}
\alias{add.in.missing.nodes}

\title{
Fill in missing values in a character vector.
}
\description{
Fills in any missing values from a named \code{character} vector (prot.list) relative to another \code{character} vector (targ.nodes) potentially containing additional names with a supplied 
\code{character} string serving as a default.
}
\usage{
add.in.missing.nodes(prot.list, targ.nodes, default.val)
}

\arguments{
  \item{prot.list}{Named \code{character} vector.}
  \item{targ.nodes}{\code{character} vector of names.}
  \item{default.val}{Value to assign each added vector element.}
}

\value{
\code{character} vector of length >= the supplied prot.list with added elements equal to default.val
}
\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{variantPriorResult}}, \code{\link{draw.graph}}
}
\examples{
set.seed(123)

vector.one <- rep("string", 10)
names(vector.one) <- paste("name", 1:10, sep="_")
vector.two <- paste("name", 1:15, sep="_")

add.in.missing.nodes(vector.one, vector.two, "default")
}

\keyword{utilities}
