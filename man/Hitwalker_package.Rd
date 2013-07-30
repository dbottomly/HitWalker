\name{HitWalker-package}
\alias{HitWalker-package}
\alias{HitWalker}
\docType{package}
\title{
Network-based variant prioriziation using functional assay scores
}
\description{
HitWalker provides the framework for prioritizing sample variants detected using an exome capture (or other method) relative to functional assay results based on protein-protein interaction data. In addition to prioritization, it also provides a simple visualization facility.  It uses the R/DBI package and therefore requires a supported datasource as well as appropriate queries. 
}
\details{
\tabular{ll}{
Package: \tab HitWalker\cr
Version: \tab 0.99.0\cr
Date: \tab 2012-11-29\cr
Depends: \tab R (>= 2.13.0), igraph, DBI, Matrix, methods, graph, Rgraphviz\cr
Suggests: \tab RMySQL, RSQLite, biomaRt, cgdsr, reshape2, Streamer, HitWalkerData\cr
License: \tab GPL-3\cr
URL: \tab www.biodevlab.org/software.html\cr
BugReports: \tab \cr
Built: \tab R 2.15.1; ; 2012-11-29 22:13:03 UTC; unix\cr
}

Index:
\preformatted{
HitWalker-package                                   Network-based variant prioriziation using functional assay scores
add.in.missing.nodes                                Fill in missing values in a character vector
annotatedIGraph-class                               Class '"annotatedIGraph"'
collapse.nodes.genes                                Collapse nodes into a single representative nodes
draw.graph                                          Low-level wrapper for drawing a graph based on supplied metadata.
draw.graph.variantPriorResult                       Plot a 'variantPriorResult' object
get.hit.prots                                       Retrieve functional assay hits from datasource
get.protein.protein.graph                           Retrieve a graph from a datasource
get.sample.overlays                                 Retrieve sample-specific metadata
get.sample.variants                                 Retrieving and filtering variants
graphDispParams-class                               Class '"graphDispParams"'
loadGraph                                           Load an annotatedIGraph object
make.target.graph                                   Creating shortest path subgraphs
prior.obj                                           A sample 'variantPriorResult' object
priorDbParams-class                                 Class '"priorDbParams"'
process.node.annot                                  Apply a function to split annotation
protein.list.to.gene                                Summarize named character vector based on supplied mappings
reconcile.sample.name                               Partial name matching in datasource
retrieve.param.query                                Send and retrieve the result of an SQL query
run.prioritization.patient                          Run RWR prioritization for a given sample
run.protocol                                        Prioritizing variants using random walk with restarts
rwrParams-class                                     Class '"rwrParams"'
variantFilter-class                                 Class '"variantFilter"'
variantPriorResult-class                            Class '"variantPriorResult"'
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{HitWalker} \tab  (source, pdf)\cr
\code{Hitwalker_Add_Metadata} \tab (source, pdf)\cr
\code{Hitwalker_Create_DB} \tab (source, pdf)\cr
\code{Hitwalker_glio} \tab (source, pdf)\cr
}
}
\author{
Daniel Bottomly

Maintainer: Daniel Bottomly <bottomly@ohsu.edu>
}
\keyword{ package }
\examples{

if (require(RSQLite) && require(HitWalkerData))
{
    data(params)
    
    db.con <- dbConnect("SQLite", hitwalker.db.path())
    
    graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
        
    test.out <- run.prioritization.patient(db.con, "08-00102", graph.obj, examp.prior.param)
    
    stopifnot(dbDisconnect(db.con))
}

}
