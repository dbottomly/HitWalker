\name{prior.obj}
\alias{prior.obj}
\docType{data}
\title{
A sample \code{variantPriorResult} object
}
\description{
This is a \code{variantPriorResult} object for a single leukemia patient sample demonstrated in Bottomly et al. submitted.
}
\usage{data(prior_obj)}
\source{
The PPI was obtained from STRING (http://string-db.org/newstring_download/protein.links.detailed.v9.0.txt.gz on 1-30-2012) and human proteins interactions with scores > 400 were kept.
}
\examples{
data(prior_obj)
head(summary(prior.obj))
}
\keyword{datasets}

