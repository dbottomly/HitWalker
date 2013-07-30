setClass(Class="variantFilter", representation=list(filter.dta="data.frame", filter.sep="character", filter.collapse.col="character", filter.ignore.col="character", filter.best.category="character"),
         prototype=prototype(filter.dta=data.frame(), filter.sep=";", filter.collapse.col="filter", filter.ignore.col="variant_qual_id", filter.best.category="OK"))

#note that this assumes that the filters are supplied in an ordinal fashion, best to worst
makeVariantFilter <- function(filter.dta, filter.sep=";", filter.collapse.col="filter", filter.ignore.col="variant_qual_id", filter.best.category="OK")
{    
    if (missing(filter.dta) || is.null(filter.dta))
    {
        filter.dta <- data.frame()
    }
    else if (class(filter.dta) != "data.frame" || all(c("category", "label", "type") %in% colnames(filter.dta)) == FALSE)
    {
        stop("ERROR: filter.dta needs to be a data.frame containing three columns: category, label and type")
    }
    
    return(new("variantFilter", filter.dta=filter.dta, filter.sep=filter.sep, filter.collapse.col=filter.collapse.col, filter.ignore.col=filter.ignore.col, filter.best.category=filter.best.category))
}

setMethod("show", signature("variantFilter"), function(object)
          {
            cat(paste("variantFilter containing", nrow(object@filter.dta), "filters\n"))
          })

setGeneric("filterDta", def=function(obj,...) standardGeneric("filterDta"))
setMethod("filterDta", signature("variantFilter"), function(obj)
          {
                return(obj@filter.dta)
          })

setGeneric("filterLabels", def=function(obj,...) standardGeneric("filterLabels"))
setMethod("filterLabels", signature("variantFilter"), function(obj)
          {
                filter.dta <- filterDta(obj)
                filter.vec <- filter.dta$type
                names(filter.vec) <- filter.dta$label
                return(filter.vec)
          })

setGeneric("filterSep", def=function(obj,...) standardGeneric("filterSep"))
setMethod("filterSep", signature("variantFilter"), function(obj)
          {
                return(obj@filter.sep)
          })

setGeneric("filterCollapseCol", def=function(obj,...) standardGeneric("filterCollapseCol"))
setMethod("filterCollapseCol", signature("variantFilter"), function(obj)
          {
                return(obj@filter.collapse.col)
          })

setGeneric("filterIgnoreCol", def=function(obj,...) standardGeneric("filterIgnoreCol"))
setMethod("filterIgnoreCol", signature("variantFilter"), function(obj)
          {
                return(obj@filter.ignore.col)
          })

setGeneric("filterCategory", def=function(obj,...) standardGeneric("filterCategory"))
setMethod("filterCategory", signature("variantFilter"), function(obj)
          {
                filter.dta <- filterDta(obj)
                
                ret.vec <- filter.dta$category
                names(ret.vec) <- filter.dta$label
                
                return(ret.vec)
          })

setGeneric("filterBestCategory", def=function(obj,...) standardGeneric("filterBestCategory"))
setMethod("filterBestCategory", signature("variantFilter"), function(obj)
          {
                return(obj@filter.best.category)
          })


setClass(Class="rwrParams", representation=list(c.val="numeric", threshold="numeric", maxit="numeric", verbose="logical"), prototype=prototype(c.val=.3, threshold=1e-10, maxit=100, verbose=FALSE))

setGeneric("getRestartProb", def=function(obj,...) standardGeneric("getRestartProb"))
setMethod("getRestartProb", signature("rwrParams"), function(obj)
          {
                return(obj@c.val)
          })

setGeneric("getConvergeThresh", def=function(obj,...) standardGeneric("getConvergeThresh"))
setMethod("getConvergeThresh", signature("rwrParams"), function(obj)
          {
                return(obj@threshold)
          })

setGeneric("getConvergeMaxIt", def=function(obj,...) standardGeneric("getConvergeMaxIt"))
setMethod("getConvergeMaxIt", signature("rwrParams"), function(obj)
          {
                return(obj@maxit)
          })
setGeneric("verbLevel", def=function(obj,...) standardGeneric("verbLevel"))
setMethod("verbLevel", signature("rwrParams"), function(obj)
          {
                return(obj@verbose)
          })

setMethod("show", signature("rwrParams"), function(object)
          {
                cat("rwrParams object:\n")
                for (slot.val in slotNames(object))
                {
                    cat(paste(slot.val,": ", slot(object, slot.val),  "\n", sep=""))
                }
          })


defaultRwrParams <- function()
{
    return(makeRwrParams(c.val=.3, threshold=1e-10, maxit=100, verbose=FALSE))
}

makeRwrParams <- function(c.val=.3, threshold=1e-10, maxit=100, verbose=FALSE)
{
    return(new("rwrParams", c.val=c.val, threshold=threshold, maxit=maxit, verbose=verbose))
}

base.query.func <- function(sample.id.list, param.obj){}

defaultPriorDbParams <- function()
{
    query.param.list <- list(pat.freq=.5, af.freq=.01, db.only.pat.num=1)
    
    return(new("priorDbParams", id.table="sample_experiment", possible.id.cols=c("core_id", "recon_name"), reconciled.id.col="recon_name",
                                                    sample.id.col="sample_id", variant.query=defaultVariantQuery, query.params=query.param.list,filter.values=defaultVariantFilter(),
                                                    matrix.query=defaultMatrixQuery, hit.query=defaultHitQuery, graph.trans.func=defaultGraphTransFunc, annot.query.func=defaultAnnotQuery, score.summary.func=default.score.summary.func,
                                                    node.weight.col=character(), patient.overlay.func=defaultPatientOverlay, seed.label="SEED", query.label="QUERY", core.id="protein", summary.id="symbol",
                                                    func.summary.map=defaultFuncSummaryMap(), hit.comb.score.name="sum.score", pat.freq.col="pat_freq", filter.col="filter",
                                                    variant.type.col="variant_type", rwr.params=defaultRwrParams(), additional.annot.ids="refseq_id", pat.var.id="pat_var_id", verbose=FALSE))
}

defaultLLSPriorDbParams <- function()
{
    query.param.list <- list(pat.freq=.5, af.freq=.01, db.only.pat.num=1)
    
    return(new("priorDbParams", id.table="sample_experiment", possible.id.cols=c("core_id", "recon_name"), reconciled.id.col="recon_name",
                                                    sample.id.col="sample_id", variant.query=defaultLLSVariantQuery, query.params=query.param.list,filter.values=defaultVariantFilter(),
                                                    matrix.query=defaultLLSMatrixQuery, hit.query=defaultLLSHitQuery, graph.trans.func=defaultGraphTransFunc, annot.query.func=defaultLLSAnnotQuery, score.summary.func=default.score.summary.func,
                                                    node.weight.col=character(), patient.overlay.func=defaultPatientOverlay, seed.label="SEED", query.label="QUERY", core.id="gene", summary.id="symbol",
                                                    func.summary.map=defaultFuncSummaryMap(), hit.comb.score.name="sum.score", pat.freq.col="pat_freq", filter.col="filter",
                                                    variant.type.col="variant_type", rwr.params=defaultRwrParams(), additional.annot.ids="refseq_id", pat.var.id="pat_var_id", verbose=FALSE))
}

testPriorDbParams <- function()
{
    query.param.list <- list(pat.freq=.5, af.freq=.01, db.only.pat.num=1)
    
    return(new("priorDbParams", id.table="sample", possible.id.cols=c("core_id", "recon_name"), reconciled.id.col="recon_name",
                                                    sample.id.col="sample_id", variant.query=testVariantQuery, query.params=query.param.list,filter.values=defaultVariantFilter(),
                                                    matrix.query=defaultMatrixQuery, hit.query=testHitQuery, graph.trans.func=defaultGraphTransFunc, annot.query.func=defaultAnnotQuery, score.summary.func=default.score.summary.func,
                                                    node.weight.col=character(), patient.overlay.func=testPatientOverlay, seed.label="SEED", query.label="QUERY", core.id="protein", summary.id="symbol",
                                                    func.summary.map=defaultFuncSummaryMap(), hit.comb.score.name="sum.score", pat.freq.col="pat_freq", filter.col="filter",
                                                    variant.type.col="variant_type", rwr.params=defaultRwrParams(), additional.annot.ids="refseq_id", pat.var.id="pat_var_id", verbose=FALSE))
}

basePriorDbParams <- function()
{
    return(new("priorDbParams", id.table="sample", possible.id.cols="recon_name", reconciled.id.col="recon_name",
                                                    sample.id.col="sample_id", variant.query=testVariantQuery, query.params=list(),filter.values=makeVariantFilter(),
                                                    matrix.query=defaultMatrixQuery, hit.query=testHitQuery, graph.trans.func=defaultGraphTransFunc, annot.query.func=defaultAnnotQuery,
                                                    score.summary.func=base.score.summary.func,node.weight.col=character(), patient.overlay.func=base.query.func, seed.label="SEED",
                                                    query.label="QUERY", core.id="protein", summary.id="symbol",func.summary.map=defaultFuncSummaryMap(), hit.comb.score.name="sum.score",
                                                    pat.freq.col="pat_freq", filter.col="filter",variant.type.col="variant_type", rwr.params=defaultRwrParams(), additional.annot.ids=character(),
                                                    pat.var.id="pat_var_id", verbose=FALSE))   
}

setClass(Class="priorDbParams", representation=list(id.table="character", possible.id.cols="character", reconciled.id.col="character", sample.id.col="character",
                                                    variant.query="function", query.params="list", filter.values="variantFilter", matrix.query="function",
                                                    hit.query="function", graph.trans.func="function", annot.query.func="function",score.summary.func="function",
                                                    node.weight.col="character", patient.overlay.func="function", seed.label="character", query.label="character", core.id="character",
                                                    summary.id="character", func.summary.map="character", hit.comb.score.name="character", pat.freq.col="character",
                                                    filter.col="character", variant.type.col="character", rwr.params="rwrParams", additional.annot.ids="character", pat.var.id="character", verbose="logical"),
                                prototype=prototype())


setMethod("show", signature("priorDbParams"), function(object)
          {
            cat(paste("priorDbParams object\n\n"))
            
            for (i in slotNames(object))
            {
                  if (is.function(slot(object, i))==FALSE && is.vector(slot(object, i)))
                  {
                           cat(paste(i, ":", paste(slot(object, i), collapse=", "), "\n"))
                  }
            }
          })


setMethod("verbLevel", signature("priorDbParams"), function(obj)
          {
                return(obj@verbose)
          })

setGeneric("verbLevel<-", def=function(obj, value) standardGeneric("verbLevel<-"))
setReplaceMethod("verbLevel", signature("priorDbParams"), function(obj, value)
          {
                slot(obj, "verbose") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("queryParams", def=function(obj,...) standardGeneric("queryParams"))
setMethod("queryParams", signature("priorDbParams"), function(obj)
          {
                return(obj@query.params)
          })
setGeneric("queryParams<-", def=function(obj, value) standardGeneric("queryParams<-"))
setReplaceMethod("queryParams", signature("priorDbParams"), function(obj, value)
          {
                slot(obj, "query.params") <- value
                validObject(obj)
                return(obj)
          })


setGeneric("scoreSummaryFunc", def=function(obj,...) standardGeneric("scoreSummaryFunc"))
setMethod("scoreSummaryFunc", signature("priorDbParams"), function(obj)
          {
                return(obj@score.summary.func)
          })

setGeneric("scoreSummaryFunc<-", def=function(obj, value) standardGeneric("scoreSummaryFunc<-"))
setReplaceMethod("scoreSummaryFunc", signature("priorDbParams"), function(obj, value)
          {
                slot(obj, "score.summary.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("patVarID", def=function(obj,...) standardGeneric("patVarID"))
setMethod("patVarID", signature("priorDbParams"), function(obj)
          {
                return(obj@pat.var.id)
          })

setGeneric("patVarID<-", def=function(obj, value) standardGeneric("patVarID<-"))
setReplaceMethod("patVarID", signature("priorDbParams"), function(obj, value)
          {
                slot(obj, "pat.var.id") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("additionalIDs", def=function(obj,...) standardGeneric("additionalIDs"))
setMethod("additionalIDs", signature("priorDbParams"), function(obj)
          {
                return(obj@additional.annot.ids)
          })

setGeneric("additionalIDs<-", def=function(obj, value) standardGeneric("additionalIDs<-"))
setReplaceMethod("additionalIDs", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "additional.annot.ids") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("rwrParamsObj", def=function(obj,...) standardGeneric("rwrParamsObj"))
setMethod("rwrParamsObj", signature("priorDbParams"), function(obj)
          {
                return(obj@rwr.params)
          })
setGeneric("rwrParamsObj<-", def=function(obj, value) standardGeneric("rwrParamsObj<-"))
setReplaceMethod("rwrParamsObj", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "rwr.params") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("variantTypeCol", def=function(obj,...) standardGeneric("variantTypeCol"))
setMethod("variantTypeCol", signature("priorDbParams"), function(obj)
          {
                return(obj@variant.type.col)
          })

setGeneric("variantTypeCol<-", def=function(obj, value) standardGeneric("variantTypeCol<-"))
setReplaceMethod("variantTypeCol", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "variant.type.col") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("patientFreqCol", def=function(obj,...) standardGeneric("patientFreqCol"))
setMethod("patientFreqCol", signature("priorDbParams"), function(obj)
          {
                return(obj@pat.freq.col)
          })

setGeneric("patientFreqCol<-", def=function(obj, value) standardGeneric("patientFreqCol<-"))
setReplaceMethod("patientFreqCol", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "pat.freq.col") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("patientFilterCol", def=function(obj,...) standardGeneric("patientFilterCol"))
setMethod("patientFilterCol", signature("priorDbParams"), function(obj)
          {
                return(obj@filter.col)
          })

setGeneric("patientFilterCol<-", def=function(obj, value) standardGeneric("patientFilterCol<-"))
setReplaceMethod("patientFilterCol", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "filter.col") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("combScoreName", def=function(obj,...) standardGeneric("combScoreName"))
setMethod("combScoreName", signature("priorDbParams"), function(obj)
          {
                return(obj@hit.comb.score.name)
          })

setGeneric("combScoreName<-", def=function(obj, value) standardGeneric("combScoreName<-"))
setReplaceMethod("combScoreName", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "hit.comb.score.name") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("getSummarizationLevel", def=function(obj,...) standardGeneric("getSummarizationLevel"))
setMethod("getSummarizationLevel", signature("priorDbParams"), function(obj, func.name)
          {
                sum.vec <- summaryMap(obj)
                
                if (func.name %in% names(sum.vec) == FALSE)
                {
                    stop("ERROR: Function name is not found in func.summary.map")
                }
                
                summary.name <- slot(obj, sum.vec[func.name])
                
                return(summary.name)
          })

setGeneric("summaryMap", def=function(obj,...) standardGeneric("summaryMap"))
setMethod("summaryMap", signature("priorDbParams"), function(obj)
          {
                return(obj@func.summary.map)
          })

setGeneric("summaryMap<-", def=function(obj, value) standardGeneric("summaryMap<-"))
setReplaceMethod("summaryMap", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "func.summary.map") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("coreID", def=function(obj,...) standardGeneric("coreID"))
setMethod("coreID", signature("priorDbParams"), function(obj)
          {
                return(obj@core.id)
          })

setGeneric("coreID<-", def=function(obj, value) standardGeneric("coreID<-"))
setReplaceMethod("coreID", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "core.id") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("summaryID", def=function(obj,...) standardGeneric("summaryID"))
setMethod("summaryID", signature("priorDbParams"), function(obj)
          {
                return(obj@summary.id)
          })

setGeneric("summaryID<-", def=function(obj, value) standardGeneric("summaryID<-"))
setReplaceMethod("summaryID", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "summary.id") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("seedLabel", def=function(obj,...) standardGeneric("seedLabel"))
setMethod("seedLabel", signature("priorDbParams"), function(obj)
          {
                return(obj@seed.label)
          })

setGeneric("seedLabel<-", def=function(obj, value) standardGeneric("seedLabel<-"))
setReplaceMethod("seedLabel", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "seed.label") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("queryLabel", def=function(obj,...) standardGeneric("queryLabel"))
setMethod("queryLabel", signature("priorDbParams"), function(obj)
          {
                return(obj@query.label)
          })

setGeneric("queryLabel<-", def=function(obj, value) standardGeneric("queryLabel<-"))
setReplaceMethod("queryLabel", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "query.label") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("idTable", def=function(obj,...) standardGeneric("idTable"))
setMethod("idTable", signature("priorDbParams"), function(obj)
          {
                return(obj@id.table)
          })

setGeneric("idTable<-", def=function(obj, value) standardGeneric("idTable<-"))
setReplaceMethod("idTable", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "id.table") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("possibleIdCols", def=function(obj,...) standardGeneric("possibleIdCols"))
setMethod("possibleIdCols", signature("priorDbParams"), function(obj)
          {
                return(obj@possible.id.cols)
          })

setGeneric("possibleIdCols<-", def=function(obj, value) standardGeneric("possibleIdCols<-"))
setReplaceMethod("possibleIdCols", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "possible.id.cols") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("reconciledIdCol", def=function(obj,...) standardGeneric("reconciledIdCol"))
setMethod("reconciledIdCol", signature("priorDbParams"), function(obj)
          {
                return(obj@reconciled.id.col)
          })

setGeneric("reconciledIdCol<-", def=function(obj, value) standardGeneric("reconciledIdCol<-"))
setReplaceMethod("reconciledIdCol", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "reconciled.id.col") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("sampleIdCol", def=function(obj,...) standardGeneric("sampleIdCol"))
setMethod("sampleIdCol", signature("priorDbParams"), function(obj)
          {
                return(obj@sample.id.col)
          })

setGeneric("sampleIdCol<-", def=function(obj, value) standardGeneric("sampleIdCol<-"))
setReplaceMethod("sampleIdCol", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "sample.id.col") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("variantQuery", def=function(obj, ...) standardGeneric("variantQuery"))
setMethod("variantQuery", signature("priorDbParams"), function(obj)
          {
                return(obj@variant.query)
          })

setGeneric("variantQuery<-", def=function(obj, value) standardGeneric("variantQuery<-"))
setReplaceMethod("variantQuery", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "variant.query") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("hitQuery", def=function(obj, ...) standardGeneric("hitQuery"))
setMethod("hitQuery", signature("priorDbParams"), function(obj)
          {
                return(obj@hit.query)
          })

setGeneric("hitQuery<-", def=function(obj, value) standardGeneric("hitQuery<-"))
setReplaceMethod("hitQuery", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "hit.query") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("patientOverlay", def=function(obj, ...) standardGeneric("patientOverlay"))
setMethod("patientOverlay", signature("priorDbParams"), function(obj)
          {
                return(obj@patient.overlay.func)
          })

setGeneric("patientOverlay<-", def=function(obj, value) standardGeneric("patientOverlay<-"))
setReplaceMethod("patientOverlay", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "patient.overlay.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("matrixQuery", def=function(obj, ...) standardGeneric("matrixQuery"))
setMethod("matrixQuery", signature("priorDbParams"), function(obj)
          {
                return(obj@matrix.query)
          })

setGeneric("matrixQuery<-", def=function(obj, value) standardGeneric("matrixQuery<-"))
setReplaceMethod("matrixQuery", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "matrix.query") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("transformGraph", def=function(obj, ...) standardGeneric("transformGraph"))
setMethod("transformGraph", signature("priorDbParams"), function(obj)
          {
                return(obj@graph.trans.func)
          })

setGeneric("transformGraph<-", def=function(obj, value) standardGeneric("transformGraph<-"))
setReplaceMethod("transformGraph", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "graph.trans.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("annotQuery", def=function(obj, ...) standardGeneric("annotQuery"))
setMethod("annotQuery", signature("priorDbParams"), function(obj)
          {
                return(obj@annot.query.func)
          })

setGeneric("annotQuery<-", def=function(obj, value) standardGeneric("annotQuery<-"))
setReplaceMethod("annotQuery", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "annot.query.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("variantFilterObj", def=function(obj, ...) standardGeneric("variantFilterObj"))
setMethod("variantFilterObj", signature("priorDbParams"), function(obj)
          {
                return(obj@filter.values)
          })

setGeneric("variantFilterObj<-", def=function(obj, value) standardGeneric("variantFilterObj<-"))
setReplaceMethod("variantFilterObj", signature("priorDbParams"), function(obj, value)
         {
                slot(obj, "filter.values") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("getFilterValues", def=function(obj, ...) standardGeneric("getFilterValues"))
setMethod("getFilterValues", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterLabels(var.filt))
          })

setGeneric("getFilterSep", def=function(obj, ...) standardGeneric("getFilterSep"))
setMethod("getFilterSep", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterSep(var.filt))
          })
setGeneric("getFilterCollapseCol", def=function(obj, ...) standardGeneric("getFilterCollapseCol"))
setMethod("getFilterCollapseCol", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterCollapseCol(var.filt))
          })

setGeneric("getFilterIgnoreCol", def=function(obj, ...) standardGeneric("getFilterIgnoreCol"))
setMethod("getFilterIgnoreCol", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterIgnoreCol(var.filt))
          })

setGeneric("getFilterLabels", def=function(obj, ...) standardGeneric("getFilterLabels"))
setMethod("getFilterLabels", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(names(filterLabels(var.filt)))
          })

setGeneric("getFilterCategory", def=function(obj, ...) standardGeneric("getFilterCategory"))
setMethod("getFilterCategory", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterCategory(var.filt))
          })

setGeneric("getFilterBestCategory", def=function(obj, ...) standardGeneric("getFilterBestCategory"))
setMethod("getFilterBestCategory", signature("priorDbParams"), function(obj)
          {
                var.filt <- variantFilterObj(obj)
                
                return(filterBestCategory(var.filt))
          })

setGeneric("nodeWeightCol", def=function(obj, ...) standardGeneric("nodeWeightCol"))
setMethod("nodeWeightCol", signature("priorDbParams"), function(obj)
          {
                return(obj@node.weight.col)
          })

