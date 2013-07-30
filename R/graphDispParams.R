basic.plot.image <- function(open.plot, draw.node, graph.params, pat.name)
{
    plotTargetGraph(open.plot, drawNode=draw.node, xpd=TRUE, graph.params=graph.params, main=pat.name)
}

setClass(Class="graphDispParams", representation=list(shape.mapping.func="function", style.mapping.func="function",
                                                      file.name="character", width="numeric", height="numeric", hit.colors="character", query.color="character",
                                                      default.color="character", collapse.char="character", default.border="character", node.comb.func="function",
                                                      label.sep="character", max.plot.vars="numeric", max.plot.hits="numeric", plot.func="function", label.func.list="list",
                                                      name.text.sep="character"),
         prototype=prototype(shape.mapping.func=defaultShapeMappingFunc, style.mapping.func=defaultStyleMappingFunc,
                             file.name=character(), width=10, height=10, hit.colors=c(is_sirna="red", is_score="green"), query.color="lightblue", collapse.char=",",
                             default.color="white", default.border="black", node.comb.func=defaultNodeCombFunc, label.sep="$", max.plot.vars=3, max.plot.hits=3,
                             plot.func=basic.plot.image, label.func.list=defaultLabelFuncList(), name.text.sep=""))

makeGraphDispParams <- function(shape.mapping.func=defaultShapeMappingFunc, style.mapping.func=defaultStyleMappingFunc, file.name=character(), width=10, height=10,
                                label.sep="$", max.plot.vars=3, max.plot.hits=3, plot.func=basic.plot.image, label.func.list=basicLabelFuncList(), name.text.sep="", node.comb.func=defaultNodeCombFunc,
                                hit.colors=c(is_sirna="red", is_score="green"), query.color="lightblue", collapse.char=",", default.color="white", default.border="black")
{
    return(new("graphDispParams", shape.mapping.func=shape.mapping.func, style.mapping.func=style.mapping.func, file.name=file.name, width=width, height=height, label.sep=label.sep, max.plot.vars=max.plot.vars,
               max.plot.hits=max.plot.hits, plot.func=plot.func, label.func.list=label.func.list, name.text.sep=name.text.sep, node.comb.func=node.comb.func, hit.colors=hit.colors,
               query.color=query.color, collapse.char=collapse.char, default.color=default.color, default.border=default.border))
}


setMethod("show", signature("graphDispParams"), function(object)
          {
                cat("graphDispParams object\n")
          })

setGeneric("nameTextSep", def=function(obj, ...) standardGeneric("nameTextSep"))
setMethod("nameTextSep", signature("graphDispParams"), function(obj)
          {
                return(obj@name.text.sep)
          })
setGeneric("nameTextSep<-", def=function(obj, value) standardGeneric("nameTextSep<-"))
setReplaceMethod("nameTextSep", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "name.text.sep") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("labelFuncList", def=function(obj, ...) standardGeneric("labelFuncList"))
setMethod("labelFuncList", signature("graphDispParams"), function(obj)
          {
                return(obj@label.func.list)
          })
setGeneric("labelFuncList<-", def=function(obj, value) standardGeneric("labelFuncList<-"))
setReplaceMethod("labelFuncList", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "label.func.list") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("plotFunc", def=function(obj, ...) standardGeneric("plotFunc"))
setMethod("plotFunc", signature("graphDispParams"), function(obj)
          {
                return(obj@plot.func)
          })
setGeneric("plotFunc<-", def=function(obj, value) standardGeneric("plotFunc<-"))
setReplaceMethod("plotFunc", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "plot.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("maxPlotVars", def=function(obj, ...) standardGeneric("maxPlotVars"))
setMethod("maxPlotVars", signature("graphDispParams"), function(obj)
          {
                return(obj@max.plot.vars)
          })

setGeneric("maxPlotVars<-", def=function(obj, value) standardGeneric("maxPlotVars<-"))
setReplaceMethod("maxPlotVars", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "max.plot.vars") <- value
                validObject(obj)
                return(obj)
          })


setGeneric("maxPlotHits", def=function(obj, ...) standardGeneric("maxPlotHits"))
setMethod("maxPlotHits", signature("graphDispParams"), function(obj)
          {
                return(obj@max.plot.hits)
          })

setGeneric("maxPlotHits<-", def=function(obj, value) standardGeneric("maxPlotHits<-"))
setReplaceMethod("maxPlotHits", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "max.plot.hits") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("labelSep", def=function(obj, ...) standardGeneric("labelSep"))
setMethod("labelSep", signature("graphDispParams"), function(obj)
          {
                return(obj@label.sep)
          })

setGeneric("labelSep<-", def=function(obj, value) standardGeneric("labelSep<-"))
setReplaceMethod("labelSep", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "label.sep") <- value
                validObject(obj)
                return(obj)
          })


setGeneric("nodeCombinationFunc", def=function(obj, ...) standardGeneric("nodeCombinationFunc"))
setMethod("nodeCombinationFunc", signature("graphDispParams"), function(obj)
          {
                return(obj@node.comb.func)
          })

setGeneric("nodeCombinationFunc<-", def=function(obj, value) standardGeneric("nodeCombinationFunc<-"))
setReplaceMethod("nodeCombinationFunc", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "node.comb.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("defaultBorder", def=function(obj, ...) standardGeneric("defaultBorder"))
setMethod("defaultBorder", signature("graphDispParams"), function(obj)
          {
                return(obj@default.border)
          })
setGeneric("defaultBorder<-", def=function(obj, value) standardGeneric("defaultBorder<-"))
setReplaceMethod("defaultBorder", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "default.border") <- value
                validObject(obj)
                return(obj)
          })


setGeneric("defaultColor", def=function(obj, ...) standardGeneric("defaultColor"))
setMethod("defaultColor", signature("graphDispParams"), function(obj)
          {
                return(obj@default.color)
          })
setGeneric("defaultColor<-", def=function(obj, value) standardGeneric("defaultColor<-"))
setReplaceMethod("defaultColor", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "default.color") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("collapseChar", def=function(obj, ...) standardGeneric("collapseChar"))
setMethod("collapseChar", signature("graphDispParams"), function(obj)
          {
                return(obj@collapse.char)
          })

setGeneric("collapseChar<-", def=function(obj, value) standardGeneric("collapseChar<-"))
setReplaceMethod("collapseChar", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "collapse.char") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("hitColors", def=function(obj, ...) standardGeneric("hitColors"))
setMethod("hitColors", signature("graphDispParams"), function(obj)
          {
                return(obj@hit.colors)
          })

setGeneric("hitColors<-", def=function(obj, value) standardGeneric("hitColors<-"))
setReplaceMethod("hitColors", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "hit.colors") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("queryColor", def=function(obj, ...) standardGeneric("queryColor"))
setMethod("queryColor", signature("graphDispParams"), function(obj)
          {
                return(obj@query.color)
          })

setGeneric("queryColor<-", def=function(obj, value) standardGeneric("queryColor<-"))
setReplaceMethod("queryColor", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "query.color") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("shapeFunc", def=function(obj, ...) standardGeneric("shapeFunc"))
setMethod("shapeFunc", signature("graphDispParams"), function(obj)
          {
                return(obj@shape.mapping.func)
          })

setGeneric("shapeFunc<-", def=function(obj, value) standardGeneric("shapeFunc<-"))
setReplaceMethod("shapeFunc", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "shape.mapping.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("styleFunc", def=function(obj, ...) standardGeneric("styleFunc"))
setMethod("styleFunc", signature("graphDispParams"), function(obj)
          {
                return(obj@style.mapping.func)
          })

setGeneric("styleFunc<-", def=function(obj, value) standardGeneric("styleFunc<-"))
setReplaceMethod("styleFunc", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "style.mapping.func") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("fileName", def=function(obj, ...) standardGeneric("fileName"))
setMethod("fileName", signature("graphDispParams"), function(obj)
          {
                return(obj@file.name)
          })

setGeneric("fileName<-", def=function(obj, value) standardGeneric("fileName<-"))
setReplaceMethod("fileName", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "file.name") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("height", def=function(obj, ...) standardGeneric("height"))
setMethod("height", signature("graphDispParams"), function(obj)
          {
                return(obj@height)
          })

setGeneric("height<-", def=function(obj, value) standardGeneric("height<-"))
setReplaceMethod("height", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "height") <- value
                validObject(obj)
                return(obj)
          })

setGeneric("width", def=function(obj, ...) standardGeneric("width"))
setMethod("width", signature("graphDispParams"), function(obj)
          {
                return(obj@width)
          })

setGeneric("width<-", def=function(obj, value) standardGeneric("width<-"))
setReplaceMethod("width", signature("graphDispParams"), function(obj, value)
          {
                slot(obj, "width") <- value
                validObject(obj)
                return(obj)
          })



