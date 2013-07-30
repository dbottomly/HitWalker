\name{graphDispParams-class}
\Rdversion{1.1}
\docType{class}
\alias{graphDispParams-class}
\alias{graphDispParams}
\alias{fileName,graphDispParams-method}
\alias{fileName<-,graphDispParams-method}
\alias{fileName}
\alias{fileName<-}
\alias{collapseChar,graphDispParams-method}
\alias{collapseChar<-,graphDispParams-method}
\alias{collapseChar}
\alias{collapseChar<-}
\alias{defaultBorder,graphDispParams-method}
\alias{defaultBorder<-,graphDispParams-method}
\alias{defaultBorder}
\alias{defaultBorder<-}
\alias{defaultColor,graphDispParams-method}
\alias{defaultColor<-,graphDispParams-method}
\alias{defaultColor}
\alias{defaultColor<-}
\alias{hitColors,graphDispParams-method}
\alias{hitColors<-,graphDispParams-method}
\alias{hitColors}
\alias{hitColors<-}
\alias{labelFuncList,graphDispParams-method}
\alias{labelFuncList<-,graphDispParams-method}
\alias{labelFuncList}
\alias{labelFuncList<-}
\alias{labelSep,graphDispParams-method}
\alias{labelSep<-,graphDispParams-method}
\alias{labelSep}
\alias{labelSep<-}
\alias{maxPlotHits,graphDispParams-method}
\alias{maxPlotHits<-,graphDispParams-method}
\alias{maxPlotHits}
\alias{maxPlotHits<-}
\alias{maxPlotVars,graphDispParams-method}
\alias{maxPlotVars<-,graphDispParams-method}
\alias{maxPlotVars}
\alias{maxPlotVars<-}
\alias{nameTextSep,graphDispParams-method}
\alias{nameTextSep<-,graphDispParams-method}
\alias{nameTextSep}
\alias{nameTextSep<-}
\alias{nodeCombinationFunc,graphDispParams-method}
\alias{nodeCombinationFunc<-,graphDispParams-method}
\alias{nodeCombinationFunc}
\alias{nodeCombinationFunc<-}
\alias{plotFunc,graphDispParams-method}
\alias{plotFunc<-,graphDispParams-method}
\alias{plotFunc}
\alias{plotFunc<-}
\alias{queryColor,graphDispParams-method}
\alias{queryColor<-,graphDispParams-method}
\alias{queryColor}
\alias{queryColor<-}
\alias{shapeFunc,graphDispParams-method}
\alias{shapeFunc<-,graphDispParams-method}
\alias{shapeFunc}
\alias{shapeFunc<-}
\alias{styleFunc,graphDispParams-method}
\alias{styleFunc<-,graphDispParams-method}
\alias{styleFunc}
\alias{styleFunc<-}
\alias{height,graphDispParams-method}
\alias{height<-,graphDispParams-method}
\alias{height}
\alias{height<-}
\alias{show,graphDispParams-method}
\alias{width,graphDispParams-method}
\alias{width<-,graphDispParams-method}
\alias{width}
\alias{width<-}
\alias{makeGraphDispParams}

\title{Class \code{"graphDispParams"}}
\description{
A parameter class controlling aspects of plot creation and display.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{makeGraphDispParams(shape.mapping.func, style.mapping.func, file.name, width, height,label.sep, max.plot.vars, max.plot.hits, plot.func, label.func.list, name.text.sep, node.comb.func, hit.colors,query.color, collapse.char, default.color, default.border)}.  Or equivalently through new("graphDispParams", ...)

}
\section{Slots}{
  \describe{
    \item{\code{shape.mapping.func}:}{Object of class \code{"function"} Maps a \code{data.frame} split by the summary ID to a node shape.}
    \item{\code{style.mapping.func}:}{Object of class \code{"function"}  Maps a \code{data.frame} split by the summary ID to a style (i.e. line type).}
    \item{\code{file.name}:}{Object of class \code{"character"} A  vector specifying the file name that the plot will be written to.  Specify \code{character(0)} for on screen display.}
    \item{\code{width}:}{Object of class \code{"numeric"} Width of the plot in inches (if applicable)}
    \item{\code{height}:}{Object of class \code{"numeric"} Height of the plot in inches (if applicable)}
    \item{\code{hit.colors}:}{Object of class \code{"character"} A vector containing colors to be used for the nodes of the hit(s).  The names should indicate the columns of the hit data.frame.}
    \item{\code{query.color}:}{Object of class \code{"character"} A vector of length 1 containing the color of the variant nodes.}
    \item{\code{default.color}:}{Object of class \code{"character"} A vector of length 1 containing the default node color.}
    \item{\code{collapse.char}:}{Object of class \code{"character"} A vector of length 1 containing the character delimiting multiple colors if a node falls into multiple color categories (e.g. has a variant and is a hit).}
    \item{\code{default.border}:}{Object of class \code{"character"} A vector of length 1 containing the default border color.}
    \item{\code{node.comb.func}:}{Object of class \code{"function"} A function that allows additional logic to be used in terms of combinations of node attributes when determining the final plot node attributes.}
    \item{\code{label.sep}:}{Object of class \code{"character"} A vector specifying the character used to break up labels onto different lines for each node.}
    \item{\code{max.plot.vars}:}{Object of class \code{"numeric"} The maximum number of variant nodes to consider for the visualization.}
    \item{\code{max.plot.hits}:}{Object of class \code{"numeric"} The maximum number of hit nodes to consider for the visualization.}
    \item{\code{plot.func}:}{Object of class \code{"function"} Function used for plotting (see HitWalker:::basic.plot.image).}
    \item{\code{label.func.list}:}{Object of class \code{"list"} A list that specifies the functions and ordering of the node labels.}
    \item{\code{name.text.sep}:}{Object of class \code{"character"} A vector specifying the character used to seperate the node name from the node labels.}
  }
}
\section{Methods}{
  \describe{
    Simple slot accession methods:
    \item{fileName}{\code{signature(obj = "graphDispParams")}}
    \item{fileName<-}{\code{signature(obj = "graphDispParams")}}
    \item{collapseChar}{\code{signature(obj = "graphDispParams")}}
    \item{collapseChar<-}{\code{signature(obj = "graphDispParams")}}
    \item{defaultBorder}{\code{signature(obj = "graphDispParams")}}
    \item{defaultBorder<-}{\code{signature(obj = "graphDispParams")}}
    \item{defaultColor}{\code{signature(obj = "graphDispParams")}}
    \item{defaultColor<-}{\code{signature(obj = "graphDispParams")}}
    \item{hitColors}{\code{signature(obj = "graphDispParams")}}
    \item{hitColors<-}{\code{signature(obj = "graphDispParams")}}
    \item{labelFuncList}{\code{signature(obj = "graphDispParams")}}
    \item{labelFuncList<-}{\code{signature(obj = "graphDispParams")}}
    \item{labelSep}{\code{signature(obj = "graphDispParams")}}
    \item{labelSep<-}{\code{signature(obj = "graphDispParams")}}
    \item{maxPlotHits}{\code{signature(obj = "graphDispParams")}}
    \item{maxPlotHits<-}{\code{signature(obj = "graphDispParams")}}
    \item{maxPlotVars}{\code{signature(obj = "graphDispParams")}}
    \item{maxPlotVars<-}{\code{signature(obj = "graphDispParams")}}
    \item{nameTextSep}{\code{signature(obj = "graphDispParams")}}
    \item{nameTextSep<-}{\code{signature(obj = "graphDispParams")}}
    \item{nodeCombinationFunc}{\code{signature(obj = "graphDispParams")}}
    \item{nodeCombinationFunc<-}{\code{signature(obj = "graphDispParams")}}
    \item{plotFunc}{\code{signature(obj = "graphDispParams")}}
    \item{plotFunc<-}{\code{signature(obj = "graphDispParams")}}
    \item{queryColor}{\code{signature(obj = "graphDispParams")}}
    \item{queryColor<-}{\code{signature(obj = "graphDispParams")}}
    \item{shapeFunc}{\code{signature(obj = "graphDispParams")}}
    \item{shapeFunc<-}{\code{signature(obj = "graphDispParams")}}
    \item{styleFunc}{\code{signature(obj = "graphDispParams")}}
    \item{styleFunc<-}{\code{signature(obj = "graphDispParams")}}
    \item{height}{\code{signature(obj = "graphDispParams")}}
    \item{height<-}{\code{signature(obj = "graphDispParams")}}
    \item{width}{\code{signature(obj = "graphDispParams")}}
    \item{width<-}{\code{signature(obj = "graphDispParams")}}
	 }
}
\author{
Daniel Bottomly
}

\seealso{

\code{\link{draw.graph.variantPriorResult}}, \code{\linkS4class{variantPriorResult}}
}
\examples{
data(prior_obj, package="HitWalker")

graph.params <- makeGraphDispParams(file.name=character())
    
draw.graph.variantPriorResult(var.obj=prior.obj, graph.params)
}
\keyword{classes}
