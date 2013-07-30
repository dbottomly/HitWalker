setClass(Class="variantPriorResult", representation=list(prox.dta="data.frame", graph="annotatedIGraph", var.dta="data.frame", hit.dta="data.frame",
                                                        sample="character", param.obj="priorDbParams", prox.summary.id="character", var.summary.id="character",
                                                        hit.summary.id="character", prox.rank.col="character", prox.type.col="character"),
                                                        prototype=prototype())

makeVariantPriorResult <- function(prox.dta, graph, sample, var.dta, hit.dta, param.obj,prox.summary.id, var.summary.id, hit.summary.id, prox.rank.col, prox.type.col)
{
    return(new("variantPriorResult", prox.dta=prox.dta, graph=graph, var.dta=var.dta, hit.dta=hit.dta, sample=sample, param.obj=param.obj, prox.summary.id=prox.summary.id,
               var.summary.id=var.summary.id, hit.summary.id=hit.summary.id, prox.rank.col=prox.rank.col, prox.type.col=prox.type.col))
}

setMethod("show", signature("variantPriorResult"), function(object)
          {
            cat(paste("Variant Prioritization Result object for sample:",object@sample,"\n"))
          })

setGeneric("getVarSummaryID", def=function(obj, ...) standardGeneric("getVarSummaryID"))
setMethod("getVarSummaryID", signature("variantPriorResult"), function(obj)
          {
                return(obj@var.summary.id)
          })

setGeneric("getProxTypeCol", def=function(obj, ...) standardGeneric("getProxTypeCol"))
setMethod("getProxTypeCol", signature("variantPriorResult"), function(obj)
          {
                return(obj@prox.type.col)
          })

setGeneric("getHitSummaryID", def=function(obj, ...) standardGeneric("getHitSummaryID"))
setMethod("getHitSummaryID", signature("variantPriorResult"), function(obj)
          {
                return(obj@hit.summary.id)
          })

setGeneric("getProxSummaryID", def=function(obj, ...) standardGeneric("getProxSummaryID"))
setMethod("getProxSummaryID", signature("variantPriorResult"), function(obj)
          {
                return(obj@prox.summary.id)
          })

setGeneric("getProxRankCol", def=function(obj, ...) standardGeneric("getProxRankCol"))
setMethod("getProxRankCol", signature("variantPriorResult"), function(obj)
          {
                return(obj@prox.rank.col)
          })

setGeneric("getSampleName", def=function(obj, ...) standardGeneric("getSampleName"))
setMethod("getSampleName", signature("variantPriorResult"), function(obj)
          {
                return(obj@sample)
          })

setGeneric("getParameters", def=function(obj, ...) standardGeneric("getParameters"))
setMethod("getParameters", signature("variantPriorResult"), function(obj)
          {
                return(obj@param.obj)
          })

setGeneric("getRankedQueries", def=function(obj, ...) standardGeneric("getRankedQueries"))
setMethod("getRankedQueries", signature("variantPriorResult"), function(obj, in.graph=FALSE)
          {
                prox.dta <- getProxDta(obj)
                prox.dta <- prox.dta[prox.dta$type == queryLabel(getParameters(obj)),]
                
                col.name <- getProxSummaryID(obj)
                rank.col <- getProxRankCol(obj)
                
                if (in.graph == TRUE)
                {
                    annot <- annotation(getGraph(obj))
                    prox.dta <- prox.dta[prox.dta[,col.name] %in% annot[,col.name],]
                }
                
                prox.dta <- prox.dta[order(as.numeric(prox.dta[,rank.col]), decreasing=FALSE),c(col.name, rank.col)]
                
                return(prox.dta)
          })

setGeneric("getProteinScores", def=function(obj, ...) standardGeneric("getProteinScores"))
setMethod("getProteinScores", signature("variantPriorResult"), function(obj, in.graph=FALSE)
          {
                hit.dta <- getHitDta(obj)
                col.name <- getHitSummaryID(obj)
                sum.score.name <- combScoreName(getParameters(obj))
                hit.dta <- hit.dta[,c(col.name, sum.score.name)]
                
                if (in.graph == TRUE)
                {
                    annot <- annotation(getGraph(obj))
                    
                    hit.dta <- hit.dta[hit.dta[,col.name] %in% annot[,col.name],]
                }
                
                hit.dta <- hit.dta[order(as.numeric(hit.dta[,sum.score.name]), decreasing=TRUE),]
                
                return(hit.dta)
          })

setGeneric("getProxDta", def=function(obj, ...) standardGeneric("getProxDta"))
setMethod("getProxDta", signature("variantPriorResult"), function(obj)
          {
                prox.dta <- obj@prox.dta
                return(prox.dta)
          })

setGeneric("getHitDta", def=function(obj, ...) standardGeneric("getHitDta"))
setMethod("getHitDta", signature("variantPriorResult"), function(obj)
          {
                hit.dta <- obj@hit.dta
                return(hit.dta)
          })

setGeneric("getVarDta", def=function(obj, ...) standardGeneric("getVarDta"))
setMethod("getVarDta", signature("variantPriorResult"), function(obj)
          {
                var.dta <- obj@var.dta
                return(var.dta)
          })

setGeneric("getFillColors", def=function(obj, ...) standardGeneric("getFillColors"))
setMethod("getFillColors", signature("variantPriorResult"), function(obj, node.names, targ.nodes, graph.param)
          {
                prox.dta <- getProxDta(obj)
                hit.dta <- getHitDta(obj)
                
                seed.label <- seedLabel(getParameters(obj))
                query.label <- queryLabel(getParameters(obj))
                
                hit.colors <- hitColors(graph.param)
                query.color <- queryColor(graph.param)
                
                collapse.char <- collapseChar(graph.param)
                
                hit.summary.id <- getHitSummaryID(obj)
                prox.summary.id <- getProxSummaryID(obj)
                
                prox.type <- getProxTypeCol(obj)
                
                if (hit.summary.id != prox.summary.id)
                {
                    stop("ERROR: hit.summary.id needs to equal prox.summary.id for getFillColors")
                }
                
                #if no name is present it assumes the same color for everything
                if (is.null(names(hit.colors)) && length(hit.colors) == 1)
                {
                    hit.colors <- c(is_hit=hit.colors)
                    stopifnot("is_hit" %in% names(hit.dta) == FALSE)
                    hit.dta$is_hit <- 1
                }
                else if (is.null(names(hit.colors)) && length(hit.colors) > 1)
                {
                    default.fail("ERROR: Need to either specify names for hit.colors or supply a vector of length 1")
                }
                else if (all(names(hit.colors) %in% colnames(hit.dta)) == FALSE)
                {
                    default.fail("ERROR: Columns specified in hit.colors should be present in the hit results")
                }
                
                split.hit <- split(hit.dta[,c(hit.summary.id, names(hit.colors))], hit.dta[,hit.summary.id])
                split.prox <- split(prox.dta, prox.dta[,prox.summary.id])
                
                fill.vec <- sapply(names(split.prox), function(x)
                       {
                            ret.col <- vector(mode="character", length=0)
                        
                            if (any(split.prox[[x]][,prox.type] == seed.label))
                            {
                                #check against split.hit for type assignment
                                temp.hit <- split.hit[[x]]
                                if (nrow(temp.hit) != 1)
                                {
                                    default.fail("ERROR: Unexpected number of hits")
                                }
                                
                                ret.col <- hit.colors[unlist(temp.hit[,names(hit.colors)]) == 1]
                                
                            }
                            if (any(split.prox[[x]][,prox.type] == query.label))
                            {
                                ret.col <- append(ret.col, query.color)
                            }
                            
                            return(paste(ret.col, collapse=collapse.char))
                       })
                
                #rename the protein to gene symbols
                cor.fill.vec <- protein.list.to.gene(fill.vec, node.names, collapse.char)

                #add in the missing nodes from targ.nodes
                cor.fill.vec <- add.in.missing.nodes(cor.fill.vec, targ.nodes, default.val=defaultColor(graph.param))

                return(cor.fill.vec)
          })

setGeneric("getBordColors", def=function(obj, ...) standardGeneric("getBordColors"))
setMethod("getBordColors", signature("variantPriorResult"), function(obj, node.names, targ.nodes, graph.param)
          {
                seed.label <- seedLabel(getParameters(obj))
                query.label <- queryLabel(getParameters(obj))
                
                query.color <- queryColor(graph.param)
                query.border <- defaultBorder(graph.param)
            
                prox.dta <- getProxDta(obj)
                prox.summary.id <- getProxSummaryID(obj)
                collapse.char <- collapseChar(graph.param)
                
                split.prox <- split(prox.dta, prox.dta[,prox.summary.id])
                
                bord.vec <- sapply(split.prox, function(x)
                       {

                            return(defaultBorder(graph.param))
                       })
                
                cor.bord.vec <- protein.list.to.gene(bord.vec, node.names, collapse.char)
                cor.bord.vec <- add.in.missing.nodes(cor.bord.vec, targ.nodes, default.val=defaultBorder(graph.param))
                
                return(cor.bord.vec)  
        })

#revise for new filtering code....
setGeneric("getBestVarText", def=function(obj, ...) standardGeneric("getBestVarText"))
setMethod("getBestVarText", signature("variantPriorResult"), function(obj, node.names, targ.nodes, graph.params)
          {
                var.dta <- getVarDta(obj)
                prox.dta <- getProxDta(obj)
                
                params <- getParameters(obj)
                filter.vals <- getFilterLabels(params)
                filter.cat <- getFilterCategory(params)
                best.filt <- getFilterBestCategory(params)
                pat.freq.col <- patientFreqCol(params)#pat_freq
                filter.col <- patientFilterCol(params)#filter
                label.sep <- labelSep(graph.params)
                label.func.list <- labelFuncList(graph.params)
                
                prox.summary.id <- getProxSummaryID(obj)
                var.summary.id <- getVarSummaryID(obj)
                rank.col <- getProxRankCol(obj)
                
                sub.prox <- prox.dta[prox.dta$type == queryLabel(getParameters(obj)),]
                split.prox <- split(sub.prox, sub.prox[,prox.summary.id])
                
                split.var <- split(var.dta, var.dta[,var.summary.id])
                
                text.vec <- sapply(names(split.prox), function(x)
                       {
                            #need to figure out why there are some duplications...
                            nd.x <- split.var[[x]][!duplicated(split.var[[x]]),]
                            
                            if (nrow(nd.x) > 1)
                            {
                                #choose representative based on minimum number of patients observed
                                nd.x <- nd.x[nd.x[,pat.freq.col] == min(nd.x[,pat.freq.col]),]
                                
                                if (nrow(nd.x) > 1)
                                {
                                    #then choose only the highest quality variant(s) assuming the filter values are supplied in an ordinal form
                                    nd.x$ord.filts <- factor(nd.x[,filter.col], levels=filter.vals, labels=filter.vals)
                                    nd.x <- nd.x[as.integer(nd.x$ord.filts) == min(as.integer(nd.x$ord.filts)),]
                                    
                                    if (nrow(nd.x) > 1)
                                    {
                                        #ok, choose one at random for now...
                                        nd.x <- nd.x[sample.int(n=nrow(nd.x), size=1),]
                                    }
                                }
                            }
                            
                            label.vec <- sapply(label.func.list, function(func)
                                                {
                                                    func(nd.x, obj, split.prox[[x]])
                                                })

                            return(paste(label.vec, collapse=label.sep))
                       })
                
                cor.text.vec <- protein.list.to.gene(text.vec, node.names)
                cor.text.vec <- add.in.missing.nodes(cor.text.vec, targ.nodes, default.val="")
                
                node.names <- names(cor.text.vec)
                
                cor.text.vec <- paste(names(cor.text.vec), cor.text.vec, sep=nameTextSep(graph.params))
                names(cor.text.vec) <- node.names
                
                return(cor.text.vec)
          })

setGeneric("getNodeShapes", def=function(obj, ...) standardGeneric("getNodeShapes"))
setMethod("getNodeShapes", signature("variantPriorResult"), function(obj, node.names, targ.nodes, graph.params)
          {
                return(process.node.annot(obj, targ.nodes, graph.params, shapeFunc))
          })

setGeneric("getNodeStyles", def=function(obj, ...) standardGeneric("getNodeStyles"))
setMethod("getNodeStyles", signature("variantPriorResult"), function(obj, node.names, targ.nodes, graph.params)
          {
                return(process.node.annot(obj,targ.nodes, graph.params, styleFunc))
          })

setGeneric("summary", def=function(object, ...) standardGeneric("summary"))
setMethod("summary", signature("variantPriorResult"), function(object, keep.cols=NULL,filter.category.col="quality.category")
          {
                #create a writeable data.frame representation of the result object for import into Excel
                
                params <- getParameters(object)
                
                prox.dta <- getProxDta(object)
                prox.dta <- prox.dta[prox.dta$type == queryLabel(params),]
                
                core.id <- coreID(params)
                summary.id <- summaryID(params)
                
                rank.col <- getProxRankCol(object)
                
                var.dta <- getVarDta(object)
                #may be able to remove this once the view is cleaned up...
                var.dta <- var.dta[! duplicated(var.dta),]
                
                hit.dta <- getHitDta(object)
                
                #add in the summary ids through merging
                annot.dta <- annotation(getGraph(object))
                
                merged.var <- merge(var.dta, prox.dta, by=core.id, all=FALSE, incomparables=NA, sort=FALSE)
                
                #natural join in this case
                merged.var <- merge(merged.var, annot.dta, all.x=TRUE, all.y=FALSE, incomparables=NA, sort=FALSE)
                
                merged.var <- merge(merged.var, hit.dta, by=core.id, all.x=TRUE, all.y=FALSE, incomparables=NA, sort=FALSE)
                
                merged.var <- merged.var[order(as.numeric(merged.var[,rank.col]), decreasing=FALSE),]
                
                var.filt <- variantFilterObj(params)
                
                if (nrow(var.filt@filter.dta) > 0)
                {
                    filt.col <- patientFilterCol(params)
                
                    var.categ <- filterCategory(var.filt)
                    
                    if ((is.null(filter.category.col) || is.na(filter.category.col)) == FALSE)
                    {
                        merged.var <- cbind(merged.var, as.character(var.categ[merged.var[,filt.col]]))
                        names(merged.var)[ncol(merged.var)] <- filter.category.col
                    }
                }
                
                if (missing(keep.cols) || is.null(keep.cols) || is.na(keep.cols))
                {
                    keep.cols <- names(merged.var)
                }
                
                merged.var <- merged.var[,keep.cols]
                
                return(merged.var)
          })

#setGeneric("getGraph", def=function(obj, ...) standardGeneric("getGraph"))
setMethod("getGraph", signature("variantPriorResult"), function(obj)
          {
                return(obj@graph)
          })

setMethod("plot", signature("variantPriorResult"), function(x, y, ...)
          {
            if (missing(y) || is.null(y))
            {
                y <- makeGraphDispParams(file.name=character())
            }
            
            add.args <- list(...)
            valid.args <- formals(makeGraphDispParams)
            
            if (length(add.args) > 0 && (is.null(names(add.args))==TRUE || all(names(add.args) %in% names(valid.args)) == FALSE))
            {
                default.fail("ERROR: In addition to x and optional y, please supply arguments that would be valid to a call to 'makeGraphDispParams'")
            }
            else if (length(add.args) > 0)
            {
                for (i in names(add.args))
                {
                    slot(y, i, check=TRUE) <- add.args[[i]]
                }
            }
            
            draw.graph.variantPriorResult(var.obj=x, graph.params=y)
          })

