setOldClass("igraph")

setClass(Class="annotatedIGraph", representation=list(graph="igraph", annotation="data.frame"), prototype=prototype())

setMethod("show", signature("annotatedIGraph"), function(object)
          {
                cat("An annotatedIGraph Object:\n")
                show(getGraph(object))
          })

setGeneric("getGraph", def=function(obj,...) standardGeneric("getGraph"))
setMethod("getGraph", signature("annotatedIGraph"), function(obj)
          {
                return(obj@graph)
          })

setGeneric("annotation", def=function(obj,...) standardGeneric("annotation"))
setMethod("annotation", signature("annotatedIGraph"), function(obj)
          {
                return(obj@annotation)
          })

setGeneric("annotation<-", def=function(obj, value) standardGeneric("annotation<-"))
setReplaceMethod("annotation", signature("annotatedIGraph"), function(obj, value)
          {
                slot(obj, "annotation") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("saveGraph", def=function(obj,...) standardGeneric("saveGraph"))
setMethod("saveGraph", signature("annotatedIGraph"), function(obj, file.name="graph.RData")
          {
                save(obj, file=file.name)
          })

loadGraph <- function(file.name, param.obj)
{
    if (!file.exists(file.name))
    {
        stop("ERROR file does not exist")
    }
    
    graph.name <- load(file.name)
    
    if (length(graph.name) == 1)
    {
        graph <- get(graph.name)
        
        if (class(graph) != "annotatedIGraph")
        {
            stop("ERROR: Unexpected class for loaded object")
        }
    }
    else
    {
        stop("ERROR: Expecting only a single object")
    }
    
    #simple sanity check of object

    if (coreID(param.obj) %in% names(annotation(graph)) && summaryID(param.obj) %in% names(annotation(graph)))
    {
        return(graph)
    }
    else
    {
        stop("ERROR: param.obj is not consistent with loaded annotatedIGraph object")
    }
}

