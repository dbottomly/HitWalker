\name{priorDbParams-class}
\Rdversion{1.1}
\docType{class}
\alias{priorDbParams-class}
\alias{priorDbParams}
\alias{patVarID,priorDbParams-method}
\alias{patVarID<-,priorDbParams-method}
\alias{patVarID}
\alias{patVarID<-}
\alias{additionalIDs,priorDbParams-method}
\alias{additionalIDs<-,priorDbParams-method}
\alias{additionalIDs}
\alias{additionalIDs<-}
\alias{annotQuery,priorDbParams-method}
\alias{annotQuery<-,priorDbParams-method}
\alias{annotQuery}
\alias{annotQuery<-}
\alias{combScoreName,priorDbParams-method}
\alias{combScoreName<-,priorDbParams-method}
\alias{combScoreName}
\alias{combScoreName<-}
\alias{coreID,priorDbParams-method}
\alias{coreID<-,priorDbParams-method}
\alias{coreID}
\alias{coreID<-}
\alias{getFilterBestCategory,priorDbParams-method}
\alias{getFilterBestCategory}
\alias{getFilterCategory,priorDbParams-method}
\alias{getFilterCategory}
\alias{getFilterCollapseCol,priorDbParams-method}
\alias{getFilterCollapseCol}
\alias{getFilterIgnoreCol,priorDbParams-method}
\alias{getFilterIgnoreCol}
\alias{getFilterLabels,priorDbParams-method}
\alias{getFilterLabels}
\alias{getFilterSep,priorDbParams-method}
\alias{getFilterSep}
\alias{getFilterValues,priorDbParams-method}
\alias{getFilterValues}
\alias{hitQuery,priorDbParams-method}
\alias{hitQuery<-,priorDbParams-method}
\alias{hitQuery}
\alias{hitQuery<-}
\alias{idTable,priorDbParams-method}
\alias{idTable<-,priorDbParams-method}
\alias{idTable}
\alias{idTable<-}
\alias{matrixQuery,priorDbParams-method}
\alias{matrixQuery<-,priorDbParams-method}
\alias{matrixQuery}
\alias{matrixQuery<-}
\alias{patientFilterCol,priorDbParams-method}
\alias{patientFilterCol<-,priorDbParams-method}
\alias{patientFilterCol}
\alias{patientFilterCol<-}
\alias{patientFreqCol,priorDbParams-method}
\alias{patientFreqCol<-,priorDbParams-method}
\alias{patientFreqCol}
\alias{patientFreqCol<-}
\alias{patientOverlay,priorDbParams-method}
\alias{patientOverlay<-,priorDbParams-method}
\alias{patientOverlay}
\alias{patientOverlay<-}
\alias{possibleIdCols,priorDbParams-method}
\alias{possibleIdCols<-,priorDbParams-method}
\alias{possibleIdCols}
\alias{possibleIdCols<-}
\alias{queryLabel,priorDbParams-method}
\alias{queryLabel<-,priorDbParams-method}
\alias{queryLabel}
\alias{queryLabel<-}
\alias{reconciledIdCol,priorDbParams-method}
\alias{reconciledIdCol<-,priorDbParams-method}
\alias{reconciledIdCol}
\alias{reconciledIdCol<-}
\alias{rwrParamsObj,priorDbParams-method}
\alias{rwrParamsObj<-,priorDbParams-method}
\alias{rwrParamsObj}
\alias{rwrParamsObj<-}
\alias{sampleIdCol,priorDbParams-method}
\alias{sampleIdCol<-,priorDbParams-method}
\alias{sampleIdCol}
\alias{sampleIdCol<-}
\alias{seedLabel,priorDbParams-method}
\alias{seedLabel<-,priorDbParams-method}
\alias{seedLabel}
\alias{seedLabel<-}
\alias{getSummarizationLevel,priorDbParams-method}
\alias{getSummarizationLevel}
\alias{summaryID,priorDbParams-method}
\alias{summaryID<-,priorDbParams-method}
\alias{summaryID}
\alias{summaryID<-}
\alias{summaryMap,priorDbParams-method}
\alias{summaryMap<-,priorDbParams-method}
\alias{summaryMap}
\alias{summaryMap<-}
\alias{variantFilterObj,priorDbParams-method}
\alias{variantFilterObj<-,priorDbParams-method}
\alias{variantFilterObj}
\alias{variantFilterObj<-}
\alias{variantQuery,priorDbParams-method}
\alias{variantQuery<-,priorDbParams-method}
\alias{variantQuery}
\alias{variantQuery<-}
\alias{variantTypeCol,priorDbParams-method}
\alias{variantTypeCol<-,priorDbParams-method}
\alias{variantTypeCol}
\alias{variantTypeCol<-}
\alias{nodeWeightCol,priorDbParams-method}
\alias{nodeWeightCol}
\alias{show,priorDbParams-method}
\alias{transformGraph,priorDbParams-method}
\alias{transformGraph<-,priorDbParams-method}
\alias{transformGraph}
\alias{transformGraph<-}
\alias{verbLevel,priorDbParams-method}
\alias{verbLevel<-,priorDbParams-method}
\alias{verbLevel}
\alias{verbLevel<-}
\alias{queryParams,priorDbParams-method}
\alias{queryParams<-,priorDbParams-method}
\alias{queryParams}
\alias{queryParams<-}
\alias{scoreSummaryFunc,priorDbParams-method}
\alias{scoreSummaryFunc<-,priorDbParams-method}
\alias{scoreSummaryFunc}
\alias{scoreSummaryFunc<-}

\title{Class \code{"priorDbParams"}}
\description{
The main parameter class for HitWalker that controls how the prioritization via RWR is performed. This class will facilitate adapting this package for alternative uses and/or databases though once set appropriately can be mostly left alone by end users. Not all slots and retrieval methods will be applicable to all deployments and can be removed or addition ones added.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("priorDbParams", ...)}.
With the slots described in the below section as the arguments in lieu of '....'
}
\section{Slots}{
  \describe{
    \item{\code{id.table}:}{Object of class \code{"character"} Name of database table containing the sample information }
    \item{\code{possible.id.cols}:}{Object of class \code{"character"} Columns containing the possible ID values supplied by the user }
    \item{\code{reconciled.id.col}:}{Object of class \code{"character"} If a match is found, the column containing the human interpretable name }
    \item{\code{sample.id.col}:}{Object of class \code{"character"} If a match is found, the column containing the database ID }
    \item{\code{variant.query}:}{Object of class \code{"function"} Function that returns a suitable query for retrieving variant information }
    \item{\code{query.params}:}{Object of class \code{"list"} A list containing additional input needed by any functions with access to a \code{priorDbParams} object} 
    \item{\code{filter.values}:}{Object of class \code{"variantFilter"} A variant filter object}
    \item{\code{matrix.query}:}{Object of class \code{"function"} Function resulting in an SQL query used to create the association graph }
    \item{\code{hit.query}:}{Object of class \code{"function"}  Function resulting in an SQL query used to retrieve hit (e.g. functional assay) data}
    \item{\code{graph.trans.func}:}{Object of class \code{"function"} Function used to transform an association matrix into one of types used for input into an RWR algorithm }
    \item{\code{annot.query.func}:}{Object of class \code{"function"} Function resulting in an SQL query used to retrieve annotation for an \code{annotatedIGraph} object }
    \item{\code{score.summary.func}:}{Object of class \code{"function"} Function used to combine multiple functional hit score columns into a summary score}
    \item{\code{node.weight.col}:}{Object of class \code{"character"} Column from the resulting association graph (in FT form) from \code{matrix.query} used as a weight }
    \item{\code{patient.overlay.func}:}{Object of class \code{"function"} Function resulting in an SQL query retrieving the patient-specific metadata }
    \item{\code{seed.label}:}{Object of class \code{"character"} Value indicating how the functional assay hits should be referred to in the downstream functions }
    \item{\code{query.label}:}{Object of class \code{"character"} Value indicating how the variants should be referred to in the downstream functions }
    \item{\code{core.id}:}{Object of class \code{"character"} Columns in the SQL-derived \code{data.frames} that refer to the core IDs (i.e. those corresponding to nodes in the graph) }
    \item{\code{summary.id}:}{Object of class \code{"character"} Columns in the SQL-derived \code{data.frames} that refer to the summary IDs (e.g. gene symbol or Ensembl gene ID) }
    \item{\code{func.summary.map}:}{Object of class \code{"character"} A named \code{character} vector indicating whether the specified function uses a core or summary ID }
    \item{\code{hit.comb.score.name}:}{Object of class \code{"character"} Column that will be added to the hit SQL query result data.frame, indicating the summarized score if applicable }
    \item{\code{pat.freq.col}:}{Object of class \code{"character"} Column indicating the number of patients with a particular allele for a given variant in the variant SQL result \code{data.frame}}
    \item{\code{filter.col}:}{Object of class \code{"character"} Column in the variant SQL result \code{data.frame} that indicates a QA/QC filter for the given variant }
    \item{\code{variant.type.col}:}{Object of class \code{"character"} Column in the variant SQL result \code{data.frame} that indicates the type of variant (e.g. SNV, INDEL) }
    \item{\code{rwr.params}:}{Object of class \code{"rwrParams"} Parameter object controlling the RWR algorithm }
    \item{\code{additional.annot.ids}:}{Object of class \code{"character"} Additional IDs provided in the \code{annotatedIGraph} object that should be used in downstream functions }
    \item{\code{pat.var.id}:}{Object of class \code{"character"} The column in the retrieved variant \code{data.frame} that uniquely identifies the variants for each patient}  
    \item{\code{verbose}:}{Object of class \code{"logical"} Indicates whether additional output should be printed in downstream functions}
}
}
\section{Methods}{
    
    Convenience methods that retrieve data from the indicated slots of the \code{variantFilter} object
    contained in a \code{priorDbParams} object.
    \describe{
    \item{getFilterBestCategory}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterCategory}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterCollapseCol}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterIgnoreCol}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterLabels}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterSep}{\code{signature(obj = "priorDbParams")}}
    \item{getFilterValues}{\code{signature(obj = "priorDbParams")}}
    }
    Retrieval and replacemement methods for the indicated slots for a \code{priorDbParams} object.
    \describe{
    \item{patVarID}{\code{signature(obj = "priorDbParams")}}
    \item{patVarID<-}{\code{signature(obj = "priorDbParams")}}
    \item{additionalIDs}{\code{signature(obj = "priorDbParams")}}
    \item{additionalIDs<-}{\code{signature(obj = "priorDbParams")}}
    \item{annotQuery}{\code{signature(obj = "priorDbParams")}}
    \item{annotQuery<-}{\code{signature(obj = "priorDbParams")}}
    \item{combScoreName}{\code{signature(obj = "priorDbParams")}}
    \item{combScoreName<-}{\code{signature(obj = "priorDbParams")}}
    \item{coreID}{\code{signature(obj = "priorDbParams")}}
    \item{coreID<-}{\code{signature(obj = "priorDbParams")}}
    \item{hitQuery}{\code{signature(obj = "priorDbParams")}}
    \item{hitQuery<-}{\code{signature(obj = "priorDbParams")}}
    \item{idTable}{\code{signature(obj = "priorDbParams")}}
    \item{idTable<-}{\code{signature(obj = "priorDbParams")}}
    \item{matrixQuery}{\code{signature(obj = "priorDbParams")}}
    \item{matrixQuery<-}{\code{signature(obj = "priorDbParams")}}
    \item{patientFilterCol}{\code{signature(obj = "priorDbParams")}}
    \item{patientFilterCol<-}{\code{signature(obj = "priorDbParams")}}
    \item{patientFreqCol}{\code{signature(obj = "priorDbParams")}}
    \item{patientFreqCol<-}{\code{signature(obj = "priorDbParams")}}
    \item{patientOverlay}{\code{signature(obj = "priorDbParams")}}
    \item{patientOverlay<-}{\code{signature(obj = "priorDbParams")}}
    \item{possibleIdCols}{\code{signature(obj = "priorDbParams")}}
    \item{possibleIdCols<-}{\code{signature(obj = "priorDbParams")}}
    \item{queryLabel}{\code{signature(obj = "priorDbParams")}}
    \item{queryLabel<-}{\code{signature(obj = "priorDbParams")}}
    \item{reconciledIdCol}{\code{signature(obj = "priorDbParams")}}
    \item{reconciledIdCol<-}{\code{signature(obj = "priorDbParams")}}
    \item{rwrParamsObj}{\code{signature(obj = "priorDbParams")}}
    \item{rwrParamsObj<-}{\code{signature(obj = "priorDbParams")}}
    \item{sampleIdCol}{\code{signature(obj = "priorDbParams")}}
    \item{sampleIdCol<-}{\code{signature(obj = "priorDbParams")}}
    \item{seedLabel}{\code{signature(obj = "priorDbParams")}}
    \item{seedLabel<-}{\code{signature(obj = "priorDbParams")}}
    \item{getSummarizationLevel}{\code{signature(obj = "priorDbParams")}}
    \item{summaryID}{\code{signature(obj = "priorDbParams")}}
    \item{summaryID<-}{\code{signature(obj = "priorDbParams")}}
    \item{summaryMap}{\code{signature(obj = "priorDbParams")}}
    \item{summaryMap<-}{\code{signature(obj = "priorDbParams")}}
    \item{variantFilterObj}{\code{signature(obj = "priorDbParams")}}
    \item{variantFilterObj<-}{\code{signature(obj = "priorDbParams")}}
    \item{variantQuery}{\code{signature(obj = "priorDbParams")}}
    \item{variantQuery<-}{\code{signature(obj = "priorDbParams")}}
    \item{variantTypeCol}{\code{signature(obj = "priorDbParams")}}
    \item{variantTypeCol<-}{\code{signature(obj = "priorDbParams")}}
    \item{nodeWeightCol}{\code{signature(obj = "priorDbParams")}}
    \item{show}{\code{signature(object = "priorDbParams")}}
    \item{transformGraph}{\code{signature(obj = "priorDbParams")}}
     \item{transformGraph<-}{\code{signature(obj = "priorDbParams")}}
    \item{verbLevel}{\code{signature(obj = "priorDbParams")}}
    \item{verbLevel<-}{\code{signature(obj = "priorDbParams")}}
    \item{queryParams}{\code{signature(obj = "priorDbParams")}}
    \item{queryParams<-}{\code{signature(obj = "priorDbParams")}}
    \item{scoreSummaryFunc}{\code{signature(obj = "priorDbParams")}}
    \item{scoreSummaryFunc<-}{\code{signature(obj = "priorDbParams")}}
}
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{variantFilter}}, \code{\linkS4class{rwrParams}}, \code{\link{run.prioritization.patient}}, \code{\linkS4class{variantPriorResult}}
}
\examples{
param.obj <- HitWalker:::testPriorDbParams()
}
\keyword{classes}
