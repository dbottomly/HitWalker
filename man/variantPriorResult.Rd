\name{variantPriorResult-class}
\Rdversion{1.1}
\docType{class}
\alias{variantPriorResult-class}
\alias{variantPriorResult}
\alias{getBestVarText,variantPriorResult-method}
\alias{getBestVarText}
\alias{getBordColors,variantPriorResult-method}
\alias{getBordColors}
\alias{getFillColors,variantPriorResult-method}
\alias{getFillColors}
\alias{getGraph,variantPriorResult-method}
\alias{getHitDta,variantPriorResult-method}
\alias{getHitDta}
\alias{getHitSummaryID,variantPriorResult-method}
\alias{getHitSummaryID}
\alias{getNodeShapes,variantPriorResult-method}
\alias{getNodeShapes}
\alias{getNodeStyles,variantPriorResult-method}
\alias{getNodeStyles}
\alias{getParameters,variantPriorResult-method}
\alias{getParameters}
\alias{getProteinScores,variantPriorResult-method}
\alias{getProteinScores}
\alias{getProxDta,variantPriorResult-method}
\alias{getProxDta}
\alias{getProxRankCol,variantPriorResult-method}
\alias{getProxRankCol}
\alias{getProxSummaryID,variantPriorResult-method}
\alias{getProxSummaryID}
\alias{getProxTypeCol,variantPriorResult-method}
\alias{getProxTypeCol}
\alias{getRankedQueries,variantPriorResult-method}
\alias{getRankedQueries}
\alias{getSampleName,variantPriorResult-method}
\alias{getSampleName}
\alias{getVarDta,variantPriorResult-method}
\alias{getVarDta}
\alias{getVarSummaryID,variantPriorResult-method}
\alias{getVarSummaryID}
\alias{show,variantPriorResult-method}
\alias{summary,variantPriorResult-method}
\alias{makeVariantPriorResult}
\alias{plot}

\title{Class \code{"variantPriorResult"}}
\description{
Contains the results from a variant prioritization run for a single sample for additional summarization and/or analysis.
}
\section{Objects from the Class}{
\describe{
Objects can be created by calls of the form \code{makeVariantPriorResult(prox.dta, graph, sample, var.dta, hit.dta, param.obj, prox.summary.id, 
									var.summary.id, hit.summary.id, prox.rank.col, prox.type.col)}.
Or equivalently through \code{new("variantPriorResult", ...)}

\item{The end user should not need to instantiate an object of this class directly, instead one will be created from a call to \code{run.prioritization.patient}.}{}
}
}
\section{Slots}{
  \describe{
    \item{\code{prox.dta}:}{Object of class \code{"data.frame"} containing the results from a RWR run produced from \code{run.protocol}.}
    \item{\code{graph}:}{Object of class \code{"annotatedIGraph"}}
    \item{\code{var.dta}:}{Object of class \code{"data.frame"} as created by a call to \code{get.sample.variants}.}
    \item{\code{hit.dta}:}{Object of class \code{"data.frame"} as created by a call to \code{get.hit.prots}.}
    \item{\code{sample}:}{Object of class \code{"character"} containing the sample name}
    \item{\code{param.obj}:}{Object of class \code{"priorDbParams"}}
    \item{\code{prox.summary.id}:}{Object of class \code{"character"} containing the summary ID used in creating the \code{prox.dta} \code{data.frame}.}
    \item{\code{var.summary.id}:}{Object of class \code{"character"} containing the summary ID used in creating the variant \code{data.frame}.}
    \item{\code{hit.summary.id}:}{Object of class \code{"character"} containing the summary ID used in creating the hit \code{data.frame}.}
    \item{\code{prox.rank.col}:}{Object of class \code{"character"} containing the column in \code{prox.dta} containing the rank of the nodes.}
    \item{\code{prox.type.col}:}{Object of class \code{"character"} containing the column in \code{prox.dta} containing the type of node (i.e. seed vs query) as named in \code{priorDbParams}.}
  }
}
\section{Methods}{
  \describe{

    Specific methods used for plotting, all of these have three additional arguments that need to be supplied:
    	
	node.names: A \code{character} vector of summary IDs named by the core IDs
    	targ.nodes: The node names of the subgraph created from \code{make.target.graph}
    	graph.param: A \code{graphDispParams} object

	All return named \code{character} vectors.

    \item{getBestVarText}{\code{signature(obj = "variantPriorResult")}: Retrieve node text for the representative variant of each node.  The variant is chosen for the 
									variant summary ID based on patient cohort frequency and quality filter as specified in the 
									variantFilter object.  The information that is displayed is controlled through the 
									\code{label.func.list} slot in \code{graphDispParams}.}
    \item{getBordColors}{\code{signature(obj = "variantPriorResult")}: 	Retrieve border colors for the nodes.  By default borders are black, see \code{getFillColors} below for 									an example of how to color borders based on node-specific data in the supplied objects.}
    \item{getFillColors}{\code{signature(obj = "variantPriorResult")}:  Retrieve the fill colors for the nodes.  By default nodes are colored based on whether they are a seed 
									or query with coloring based on multiple categories being supported.}
    \item{getNodeShapes}{\code{signature(obj = "variantPriorResult")}: 	Retrieve node shapes as specified in the \code{getShapeFunc} method of \code{graphDispParams} and the 
								       	\code{process.node.annot} function.}
    \item{getNodeStyles}{\code{signature(obj = "variantPriorResult")}:  Retrieve node border styles (as in the lty graphical parameter) as specified in the \code{getStyleFunc} 									method of \code{graphDispParams} and the \code{process.node.annot} function.}

    \item{General methods that access the specified slots}{}
    \item{getGraph}{\code{signature(obj = "variantPriorResult")}}
    \item{getHitDta}{\code{signature(obj = "variantPriorResult")}}
    \item{getHitSummaryID}{\code{signature(obj = "variantPriorResult")}}
    \item{getParameters}{\code{signature(obj = "variantPriorResult")}}
    \item{getProxDta}{\code{signature(obj = "variantPriorResult")}}
    \item{getProxRankCol}{\code{signature(obj = "variantPriorResult")}}
    \item{getProxSummaryID}{\code{signature(obj = "variantPriorResult")}}
    \item{getProxTypeCol}{\code{signature(obj = "variantPriorResult")}}
    \item{getSampleName}{\code{signature(obj = "variantPriorResult")}}
    \item{getVarDta}{\code{signature(obj = "variantPriorResult")}}
    \item{getVarSummaryID}{\code{signature(obj = "variantPriorResult")}}
    
    \item{Methods that summarize data}{}
    \item{summary}{\code{signature(object = "variantPriorResult")}: Produces a merged \code{data.frame} summary of the prioritization.  Desired columns can be specified as arguments.}
    \item{plot}{\code{signature(object="variantPriorResult")}: Produces a plot of the results. In order to modify the graphical parameters, a \code{graphDispParams} object and/or values corresponding to the slots of the object can be supplied}
    \item{getProteinScores}{\code{signature(obj = "variantPriorResult")}: Retrieves a \code{data.frame} of the hit scores used for prioritization ordered by score.  There are 										two columns, the hit summary ID and summary score column specified in \code{priorDbParams}.}	
    \item{getRankedQueries}{\code{signature(obj = "variantPriorResult")}: Retrieves an ordered \code{data.frame} of prioritization ranks.  There are two columns, 
									the prioritization summary ID and the rank column specified in \code{priorDbParams}.} 
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{graphDispParams}}, \code{\linkS4class{priorDbParams}}, \code{\link{make.target.graph}}, \code{\link{process.node.annot}}
}
\examples{

if(require(RSQLite) && require(HitWalkerData))
{
	data(params)
	db.con <- dbConnect("SQLite", hitwalker.db.path())
	graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
	test.out <- run.prioritization.patient(db.con=db.con, sample.id="08-00102", graph.obj=graph.obj, param.obj=examp.prior.param)
	head(summary(test.out))
}


}
\keyword{classes}
