\name{run.protocol}
\alias{run.protocol}
\alias{random.walk}
\title{
Prioritizing variants using random walk with restarts
}
\description{
The \code{random.walk} function is the workhorse function for generating association scores in this implementation.  The \code{run.protocol} function wraps \code{random.walk} and specifically determines the association scores for the sample variants and ranks them for downstream analysis.
}
\usage{
run.protocol(graph.sp.mat, seed.prots, samp.vars, param.obj)
random.walk(graph.sp.mat, seed.prots, rwr.params)
}
\arguments{
  \item{graph.sp.mat}{
	A \code{dgCMatrix} from the Matrix package.  Probably other types of sparse matrices would work as well, though their use is untested.
}
  \item{seed.prots}{
	A numeric vector of hit scores (termed seeds) corresponding to nodes.  The vector needs to be named by node names in \code{graph.sp.mat}.
}
  \item{samp.vars}{
	A \code{data.frame} of variants.  A column corresponding to the node IDs is necessary.
}
  \item{param.obj}{
	A \code{priorDbParams} object
}
  \item{rwr.params}{
	A \code{rwrParams} object
}
}
\value{
For run.protocol a list with 4 elements:
\describe{
\item{\code{res.dta}:}{A \code{data.frame} containing the prioritization results}
\item{\code{summary.id}:}{The column name of the node IDs}
\item{\code{rank.col}:}{The column of \code{res.dta} containing the rank}
\item{\code{type.id}:}{The column of \code{res.dta} containing whether the node ID corresponds to a seed (i.e. hit) or query (i.e. variant)}
}
For random.walk a list with 3 elements:
\describe{
\item{\code{prox.vector}:}{A named ordered \code{numeric} vector containing the proximity scores for all nodes in the graph}
\item{\code{seed.prots}:}{The names of the hits that were present in the graph}
\item{\code{prot.weights}:}{The corresponding scores of the values in seed.prots}
}
}

\note{
The RWR implementation was based off of the Matlab code of Erten et al. 2011 (http://compbio.case.edu/dada/).
}

\references{
S. Erten, G. Bebek, R. Ewing and M. Koyuturk. DADA: Degree-aware algorithms for network-based disease gene prioritization. BMC BioData Mining, 4:19, 2011. 
}
\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\linkS4class{rwrParams}}, \code{\linkS4class{dgCMatrix}}
}
\examples{
data(prior_obj, package="HitWalker")

use.graph <- getGraph(prior.obj)
param.obj <- getParameters(prior.obj)
hit.dta <- getHitDta(prior.obj)

hit.prot.list <- split(hit.dta, hit.dta$protein)
seed.prots <- sapply(hit.prot.list, "[[", "sum.score")

graph.sp.mat <- transformGraph(param.obj)(use.graph)
samp.vars <- getVarDta(prior.obj)

run.res <- run.protocol(graph.sp.mat, seed.prots, samp.vars, param.obj)

rwr.res <- random.walk(graph.sp.mat, seed.prots, rwrParamsObj(param.obj))
}

\keyword{utilities}
