\name{protein.list.to.gene}
\alias{protein.list.to.gene}

\title{
Summarize named character vector based on supplied mappings
}
\description{
Renames a character vector based on a supplied mapping between symbols in the form of a named \code{character} vector.  If the renamed vector contains duplicated names, 
and \code{collapse.char} is present, the elements are split by the character and the longest split element is kept as long as the shorter split elements are subsets of it.  If \code{collapse.char} is NULL, the elements have to be equivalent.  If these requirements are not met an error is thrown. Only the names in common between the two vectors are considered.
}
\usage{
protein.list.to.gene(prot.list, node.names, collapse.char=NULL)
}

\arguments{
  \item{prot.list}{
	A named \code{character} vector with names corresponding to the names of \code{node.names}.
}
  \item{node.names}{
	A named \code{character} vector.
}
  \item{collapse.char}{
	Either NULL or a \code{character}.
}
}

\value{
\describe{
A potentially shorter version of \code{prot.list} with new names.
}
}
\author{
Daniel Bottomly
}

\examples{
set.seed(123)

vector.1 <- paste("value", 1:10, sep=".")
names(vector.1) <- paste("sample", 1:10, sep=".")

vector.2 <- as.character(sapply(1:5, function(x) rep(paste("new.sample", x, sep="."),2)))
names(vector.2) <-  paste("sample", 1:10, sep=".")

\dontrun{protein.list.to.gene(prot.list=vector.1, node.names=vector.2)}##Not run; should fail as there are different values in the elements of vector.1 to be summarized

vector.1 <- as.character(sapply(1:5, function(x) rep(paste("value", x, sep="."),2)))
names(vector.1) <- paste("sample", 1:10, sep=".")

protein.list.to.gene(prot.list=vector.1, node.names=vector.2)

vector.1[1:2] <- paste(paste("value", 1:3, sep="."), collapse=",")

protein.list.to.gene(prot.list=vector.1, node.names=vector.2, collapse.char=",")

vector.1[2] <- "value.1"
protein.list.to.gene(prot.list=vector.1, node.names=vector.2, collapse.char=",")

vector.1[2] <- paste(paste("value", 4:6, sep="."), collapse=",")
\dontrun{protein.list.to.gene(prot.list=vector.1, node.names=vector.2, collapse.char=",")}##should fail as the elements are the same length but not the same

vector.1[2] <- "value.4"
\dontrun{protein.list.to.gene(prot.list=vector.1, node.names=vector.2, collapse.char=",")}##the value for vector.1 is not a subset of vector.2

}

\keyword{utilities}
