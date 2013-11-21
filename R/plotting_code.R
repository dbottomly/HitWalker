#won't export these functions

check.label.params <- function(nd.x, obj, sub.prox)
{
    if ((class(obj) == "variantPriorResult" && is.data.frame(nd.x) && is.data.frame(sub.prox)) == FALSE)
    {
        default.fail("ERROR: Invalid arguments supplied to label function")
    }
}

check.col.data.frame <- function(value, dta)
{
    if (value %in% colnames(dta) == FALSE)
    {
        default.fail(paste("ERROR:",value,"is not in supplied data.frame"))
    }
}


lls.plot.image <- function(open.plot, draw.node, graph.params, pat.name)
{
    layout(rbind(c(2,3), c(1,1)), widths=c(.25,.75), heights=c(.15, .85))
    plot(open.plot, drawNode=draw.node, xpd=TRUE, graph.params=graph.params)
    par(mar=rep(0,4))
    plot.new()
    legend("left", legend=c("siRNA", "Gene Target", "Mutation", "Not assay hit", "Not on SeqCap"), pch=c(19,19,19, NA_integer_, NA_integer_),
           col=c(hitColors(graph.params)[c("is_sirna", "is_score")], queryColor(graph.params), "black", "black"),
           lty=c("blank", "blank", "blank", "dashed", "dotted"))
    plot.new()
    text(x=grconvertX(.35, from="npc", to="user"), y=grconvertY(.7, from="npc", to="user"), paste("Sample:", pat.name), cex=2)
    text(x=grconvertX(.35, from="npc", to="user"), y=grconvertY(.5, from="npc", to="user"), "Circles indicate genes and contain the gene name", cex=1)
    text(x=grconvertX(.35, from="npc", to="user"), y=grconvertY(.4, from="npc", to="user"), "Genes with variants also indicate rank, frequency and (W) if variant is of lower quality", cex=1)
    text(x=grconvertX(.35, from="npc", to="user"), y=grconvertY(.3, from="npc", to="user"), "Triangles indicate gene over/under expression", cex=1)
    text(x=grconvertX(.35, from="npc", to="user"), y=grconvertY(.2, from="npc", to="user"), "Half triangles indicate exon dropout events", cex=1)
}

#adapted from igraph package
igraph.to.graphNEL.slim <- function (graph)
{
    if (!is.igraph(graph)) {
        stop("Not an igraph graph")
    }
    require(graph)
    if ("name" %in% list.vertex.attributes(graph) && is.character(V(graph)$name)) {
        name <- V(graph)$name
    }
    else {
        name <- as.character(seq(vcount(graph)))
    }
    edgemode <- if (is.directed(graph))
        "directed"
    else "undirected"
    if ("weight" %in% list.edge.attributes(graph) && is.numeric(E(graph)$weight)) {
        al <- get.adjedgelist(graph, "out")
        for (i in seq(along = al)) {
            edges <- get.edges(graph, al[[i]])
            edges <- ifelse(edges[, 2] == i, edges[, 1], edges[, 2])
            weights <- E(graph)$weight[al[[i]]]

            dup.edges <- duplicated(edges)
            edges <- edges[!dup.edges]
            weights <- weights[!dup.edges]

            #DWB was here...
            #al[[i]] <- list(edges = edges, weights = weights)
            al[[i]] <- list(edges = V(graph)[edges]$name, weights = weights)
        }
    }
    else {
        al <- get.adjlist(graph, "out")
        al <- lapply(al, function(x) list(edges = x))
    }
    names(al) <- name
    res <- new("graphNEL", nodes = name, edgeL = al, edgemode = edgemode)
    res
}

#adapted from Rgraphviz packages

setGeneric("plotTargetGraph", def=function(x,y,...) standardGeneric("plotTargetGraph"))
setMethod("plotTargetGraph", signature("Ragraph"),
          function (x, y, ...) 
{
    .local <- function (x, y, edgeAttrs = list(), ..., main = NULL, 
        cex.main = NULL, col.main = "black", sub = NULL, cex.sub = NULL, 
        col.sub = "black", drawNode = drawAgNode, xlab, ylab, graph.params=NULL) 
    {
        if (missing(y)) 
            y <- x@layoutType
        x <- graphLayout(x, y)
        plot.new()
        old.mai = par(mai = 0.01 + c(0.83 * (!is.null(sub)), 
            0, 0.83 * (!is.null(main)), 0))
        on.exit(par(mai = old.mai), add = TRUE)
        ur <- upRight(boundBox(x))
        bl <- botLeft(boundBox(x))
        if (x@bg != "") 
            par(bg = x@bg)
        if (x@fg != "") 
            par(fg = x@fg)
        plot.window(xlim = c(getX(bl), getX(ur)), ylim = c(getY(bl), 
            getY(ur)), log = "", asp = NA, ...)
        xy <- xy.coords(NA, NA)
        
        plot.xy(xy, type = "n", ...)
        if (!missing(xlab) && !missing(ylab)) 
            stop("Arguments 'xlab' and 'ylab' are not handled.")
        if (!is.null(sub) || !is.null(main)) 
            title(main, sub, cex.main = cex.main, col.main = col.main, 
                cex.sub = cex.sub, col.sub = col.sub)
        agn <- AgNode(x)
        nodeDims <- sapply(agn, function(n) {
            c(getNodeRW(n) + getNodeLW(n), getNodeHeight(n))
        })
        strDims <- sapply(agn, function(n) {
            s <- labelText(txtLabel(n))
            #changes from DWB
            
            if (length(s) == 0) {
                rv <- c(strwidth(" "), strheight(" "))
            }
            else {
                split.s <- strsplit(s, labelSep(graph.params), fixed=TRUE)[[1]]
                #rv <- c(strwidth(s) * 1.1, strheight(s) * 1.4)
                rv <- c(strwidth(split.s[1]) * 1.1, strheight(split.s[1]) * 1.4)
            }
            return(rv)
        })
        cex <- min(nodeDims/strDims)
        if (is.finite(cex) && cex > 0) {
            old.cex <- par(cex = cex)
            on.exit(par(cex = old.cex), add = TRUE)
        }
        if (length(drawNode) == 1) {
            lapply(agn, drawNode)
        }else {
            if (length(drawNode) == length(AgNode(x))) {
                for (i in seq(along = drawNode)) {
                  drawNode[[i]](agn[[i]], graph.params=graph.params)
                }
            }
            else {
                stop(paste("Length of the drawNode parameter is ", 
                  length(drawNode), ", it must be either length 1 or the number of nodes.", 
                  sep = ""))
            }
        }
        arrowLen <- par("pin")[1]/diff(par("usr")[1:2]) * min(nodeDims)/pi
        lapply(AgEdge(x), lines, len = arrowLen, edgemode = edgemode, 
            ...)
        invisible(x)
    }
    .local(x, y, ...)
}
)
                                                            
drawAgNodeTG <- function (node, ur, attrs, radConv, graph.params) 
{
    
    nodeCenter <- getNodeCenter(node)
    nodeX <- getX(nodeCenter)
    nodeY <- getY(nodeCenter)
    lw <- getNodeLW(node)
    rw <- getNodeRW(node)
    rad <- (lw + rw)/2
    height <- getNodeHeight(node)
    fg <- color(node)
    style <- style(node)
    shape <- shape(node)
    if (shape == "") 
        shape <- "ellipse"
    if (fg == "") 
        fg <- "black"
    bg <- fillcolor(node)
    if (bg == "") {
        bg <- "transparent"
    }
    
    color.sep <- collapseChar(graph.params)
    
    bg <- unlist(strsplit(bg, color.sep))
    
    lwd <- 1
    lty <- style
    
    label.sep <- labelSep(graph.params)
    
    #browser()
    switch(shape,
           ellipse = {ellipse(x = nodeX, 
                y = nodeY, height = height, width = rad * 2, fg = fg, 
                bg = bg, lwd=lwd, lty=lty); draw.multiline.text(txtLabel(node), nodeX, nodeY, height, label.sep=label.sep)},
           up_std={draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="up_std");
                draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type=NULL, label.sep=label.sep)},
           down_std={draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="down_std");
                    draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type="down_std", label.sep=label.sep)},
           up_half={draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="up_half");
                    draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type="up_half", label.sep=label.sep)},
           down_half={draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="down_half");
                    draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type="down_half", label.sep=label.sep)},
           both_down={draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="down_std");
                    draw.triangle(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="down_half");
                    draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type="down_half", label.sep=label.sep)},
            rectangle={draw.rectangle(nodeX, nodeY, rw, lw, bg, fg, lwd, lty);
            draw.multiline.text(txtLabel(node), nodeX, nodeY, height, type=NULL, label.sep=label.sep)}
    , stop("Unimplemented node shape: ", shape(node)))
}

#taken from Rgraphviz ellipse, but with bigger borders...
ellipse <- function (x, y, width, height = width, theta = 2 * pi, npoints = 100, 
    fg = par("fg"), bg = par("bg"), lwd=1, lty=1) 
{
    a <- width/2
    b <- height/2
    xcoord <- seq(-a, a, length = npoints)
    ycoord.neg <- sqrt(b^2 * (1 - (xcoord)^2/a^2))
    ycoord.pos <- -sqrt(b^2 * (1 - (xcoord)^2/a^2))
    xx <- c(xcoord, xcoord[npoints:1])
    yy <- c(ycoord.neg, ycoord.pos)
    x.theta <- xx * cos(2 * pi - theta) + yy * sin(2 * pi - theta) + x
    y.theta <- yy * cos(2 * pi - theta) - xx * sin(2 * pi - theta) + y
    
    if (length(bg) > 1)
    {
        multi.shading(cur.points.x=x.theta, cur.points.y=y.theta, bg=bg, width=width)
        
        polygon(x.theta, y.theta, density = NA, border = fg, col = "transparent", lwd=lwd, lty=lty)
    }
    else if (length(bg) == 1)
    {
        polygon(x.theta, y.theta, density = NA, border = fg, col = bg, lwd=lwd, lty=lty)
        #width/72
    }
    else
    {
        stop("ERROR: Unexpected length of bg")
    }
}

multi.shading <- function(cur.points.x, cur.points.y, bg, width)
{
    for (i in 1:(length(bg)-1))
    {
        mid.diffs <- cur.points.x - (min(cur.points.x) + (width/length(bg)))
        diff.pos <- which(abs(mid.diffs) == min(abs(mid.diffs)))
        if (length(diff.pos) != 2)
        {
            #if there are ties as is often the case in bg of length 2, take the first
            if (length(diff.pos) == 4 && (diff.pos[1] + 1) == diff.pos[2] && (diff.pos[3] + 1) == diff.pos[4])
            {
                diff.pos <- diff.pos[c(1, 3)]
            }
            else
            {
                default.fail("ERROR: diff.pos not equal to 2")
            }
        }
        #stopifnot(length(diff.pos) == 2)
        rm.sect <- diff.pos[1]:diff.pos[2]
            
        polygon(cur.points.x[-rm.sect], cur.points.y[-rm.sect], density = NA, border = bg[i], col = bg[i])
            
        cur.points.x <- cur.points.x[rm.sect]
        cur.points.y <- cur.points.y[rm.sect]
    }
   
    polygon(cur.points.x, cur.points.y, density = NA, border = bg[length(bg)], col = bg[length(bg)])
}

draw.rectangle <- function(nodeX, nodeY, rw, lw, bg, fg, lwd, lty)
{
   
    if (length(bg) == 1)
    {
        rect(nodeX - lw, nodeY - (height/2), nodeX + rw, nodeY + (height/2), col = bg, border = fg,  lty = lty, lwd = lwd)
    }
    else
    {
        stop("ERROR: multi bg coloring not yet supported for the rectangle")
    }
    
}

draw.triangle <- function(nodeX, nodeY, height, rw, lw, bg, fg, lwd, lty, type="up_std")
{
    
    if (type == "up_std")
    {
        x.coords <- c(nodeX - lw,nodeX,nodeX + rw)
        y.coords <- c(nodeY-(height/2), nodeY+(height/2) , nodeY-(height/2))
    }
    else if (type == "down_std")
    {
        x.coords=c(nodeX - lw, nodeX, nodeX + rw)
        y.coords=c(nodeY+(height/2), nodeY-(height/2) , nodeY+(height/2))
    }
    else if (type == "up_half")
    {
        x.coords <- c(nodeX - lw, nodeX - (lw/2), nodeX ,nodeX + rw, nodeX + rw)
        y.coords <- c(nodeY-(height/2), nodeY+(height/2), nodeY ,nodeY,nodeY-(height/2))
    }
    else if (type == "down_half")
    {
        x.coords <- c(nodeX - lw, nodeX - (lw/2), nodeX ,nodeX + rw, nodeX + rw)
        y.coords <- c(nodeY+(height/2), nodeY-(height/2), nodeY ,nodeY,nodeY+(height/2))
    }
    else
    {
        stop("ERROR: Unknown type of triangle")
    }
    
    if (length(bg) > 1)
    {
        #need to interpolate some points to make this work..
        if (type %in% c("down_std", "up_std"))
        {
            main.approx <- approx(x.coords, y.coords, n=1000)
            last.approx <- approx(x.coords[c(length(x.coords), 1)], y.coords[c(length(y.coords), 1)], n=1000)
            last.x.ord <- order(last.approx$x, decreasing=TRUE)
            
            use.coords.x <- c(main.approx$x, last.approx$x[last.x.ord])
            use.coords.y <- c(main.approx$y, last.approx$y[last.x.ord])
        }
        else
        {
            main.approx <- approx(x.coords[1:4], y.coords[1:4], n=1000)
            last.approx.1 <- approx(y.coords[4:5], x.coords[4:5], n=1000, method="constant")
            last.approx.2 <- approx(x.coords[c(5, 1)], y.coords[c(5,1)], n=1000)
            
            #have to flip the x and y to perform the interpolation correctly
            last.approx.1.temp <- list(x=last.approx.1$y, y=last.approx.1$x)
            last.x.1.ord <- order(last.approx.1.temp$x, decreasing=TRUE)
            last.x.2.ord <- order(last.approx.2$x, decreasing=TRUE)
            
            use.coords.x <- c(main.approx$x, last.approx.1.temp$x[last.x.1.ord], last.approx.2$x[last.x.2.ord])
            use.coords.y <- c(main.approx$y, last.approx.1.temp$y[last.x.1.ord], last.approx.2$y[last.x.2.ord])
            
            #stop("ERROR: Unimplemented shading procedure")
        }
       
        multi.shading(cur.points.x=use.coords.x, cur.points.y=use.coords.y, bg=bg, width=lw+rw)
        polygon(x.coords, y.coords, density = NA, border = fg, col = "transparent", lwd=lwd, lty=lty)
    }
    else if (length(bg) == 1)
    {
        polygon(x=x.coords, y=y.coords, col=bg, border=fg, lwd=lwd, lty=lty)
    }
    else
    {
        stop("ERROR: Unexpected length of bg")
    }
    
}

draw.multiline.text <- function(use.text, nodeX, nodeY, height, type=NULL, label.sep="$")
{
    if (class(use.text) == "character")
    {
        text.lines <- strsplit(use.text, label.sep, fixed=TRUE)[[1]]
    }
    else
    {
        text.lines <- strsplit(labelText(use.text), label.sep, fixed=TRUE)[[1]]
    }
    
    text.lines <- text.lines[text.lines != ""]
    
    if (length(text.lines) > 1)
    {
        if (is.null(type))
        {
            height.len <- (height/length(text.lines))/1.2
            for (i in 1:length(text.lines))
            {   
                text(labels=text.lines[i], x= nodeX, y=(nodeY+((height/2) - (i*height.len))), cex=.8)
            }
        }
        else if (type == "up_std")
        {
            for (i in 1:length(text.lines))
            {
                use.y <- (nodeY - (i*(strheight(text.lines[i], cex=.8)*1.5)))
                text(labels=text.lines[i], x= nodeX, y=use.y, cex=.8)
            }
        }
        else if (type == "down_std")
        {
            for (i in 1:length(text.lines))
            {
                use.y <- (nodeY - (i*(strheight(text.lines[i], cex=.8)*1.5))) + (height/2)
                text(labels=text.lines[i], x= nodeX, y=use.y, cex=.8)
            }
        }
        else if (type == "down_half")
        {
            text.lines <- rev(text.lines)
            for (i in 1:length(text.lines))
            {
                use.y <- ((nodeY - (height/8)) + (i*(strheight(text.lines[i], cex=.8)*1.5)))
                text(labels=text.lines[i], x= nodeX, y=use.y, cex=.8)
            }
        }
        else if (type == "up_half")
        {
            for (i in 1:length(text.lines))
            {
                use.y <- ((nodeY + (height/8)) - (i*(strheight(text.lines[i], cex=.8)*1.5)))
                text(labels=text.lines[i], x= nodeX, y=use.y, cex=.8)
            }
        }
        else
        {
            stop("ERROR: Unimplemented shape")
        }
        
    }
    else
    {
        if (is.null(type))
        {
            text(labels=text.lines, x= nodeX, y=nodeY, cex=1)
        }
        else if (type %in% c("up_std", "up_half"))
        {
            text(labels=text.lines, x= nodeX, y=nodeY-(height/4), cex=1)
        }
        else if (type %in% c("down_std", "down_half"))
        {
            text(labels=text.lines, x= nodeX, y=nodeY+(height/8), cex=1)
        }
    }
}

plot.target.graph <- function(target.graph, node.attrs, draw.node, var.obj, graph.params)
{   
    file.name <- fileName(graph.params)
    width <- width(graph.params)
    height <- height(graph.params)
    pat.name <- getSampleName(var.obj)
    
    eAttrs <- list()
    
    all.edge.weights <- unlist(edgeWeights(target.graph))
    
    #make prettier to look at
    
    all.edge.weights <- trunc(all.edge.weights*100)
    
    names(all.edge.weights) <- sub("\\.", "~", names(all.edge.weights))
    
    #a test agopen to make sure the right edge names are found
    #test.ag <- agopen(target.graph, name="test_graph", recipEdges="combined")
    #issue found when testing on R-2.15.1 on Windows as use.edge.names was seen to be a list.
    #edgeNames does not seem to work for the Ragraph class
    #use.edge.names <- edgeNames(test.ag)
    
    #so the above was commented out and below statement was used as a workaround and probably better way to deal with
    #any duplicated edges per the Rgraphviz manual per the 'combined' value for recipEdges
    use.edge.names <- setdiff(seq(along=all.edge.weights), removedEdges(target.graph))
    
    keep.edge.weights <- all.edge.weights[use.edge.names]
    
    eAttrs$label <- as.character(keep.edge.weights)
    names(eAttrs$label) <- names(keep.edge.weights)
    
    stopifnot((2*length(eAttrs$label)) == length(unlist(graph::edges(target.graph))))
    
    #adjust height and width
    
    #by default:

    #$node$height
    #[1] "0.5"
    #
    #$node$width
    #[1] "0.75"
    #
    #
    node.attrs$height <- rep(.5*1.8, length(nodes(target.graph)))
    names(node.attrs$height) <- nodes(target.graph)
    
    node.attrs$width <- rep(.75*1.5, length(nodes(target.graph)))
    names(node.attrs$width) <- nodes(target.graph)
    open.plot <- agopen(target.graph, name="", layoutType="dot", nodeAttrs=node.attrs, edgeAttrs=eAttrs, recipEdges="combined")
   
    #add in the styles manually for some reason...
    
    nodes <- AgNode(open.plot)
    for (i in 1:length(nodes))
    {
        nodes[[i]]@style <- as.character(node.attrs$style[nodes[[i]]@name])
    }
    
    AgNode(open.plot) <- nodes
    
    if (length(file.name) == 0)
    {
        plotFunc(graph.params)(open.plot, draw.node, graph.params, pat.name)
    }
    else
    {
        pdf (file=file.name, width=width, height=height)
        plotFunc(graph.params)(open.plot, draw.node, graph.params, pat.name)
        dev.off()
    }
    
}

#will export these functions

process.node.annot <- function(obj, targ.nodes, graph.params, func.name)
{
    if ((class(obj) == "variantPriorResult" && is.character(targ.nodes) && class(graph.params) == "graphDispParams" && is.function(func.name)) == FALSE)
    {
        stop("ERROR: Invalid arguments supplied to process.node.annot")
    }
    
    graph.annot <- annotation(getGraph(obj))
    
    summary.col <- summaryID(getParameters(obj))
    
    sub.graph.annot <- graph.annot[graph.annot[,summary.col] %in% targ.nodes,]
    split.annot <- split(sub.graph.annot, sub.graph.annot[,summary.col])
    
    node2annot <- sapply(split.annot, function(x)
                       {    
                            return(func.name(graph.params)(x))
                       })
                
    return(node2annot)
}

add.in.missing.nodes <- function(prot.list, targ.nodes, default.val)
{
    if ((is.character(prot.list) && is.character(targ.nodes) && is.character(default.val) && is.null(names(prot.list)) == FALSE && length(default.val) == 1) == FALSE)
    {
        stop("ERROR: Invalid arguments supplied to add.in.missing.nodes")
    }
    
    missing.nodes <- setdiff(targ.nodes, names(prot.list))
    
    if (length(missing.nodes) > 0)
    {
        subs.val <- rep(default.val, length(missing.nodes))
        names(subs.val) <- missing.nodes
        
        return(c(prot.list, subs.val))
    }
    else
    {
        return(prot.list)
    }
}


protein.list.to.gene <- function(prot.list, node.names, collapse.char=NULL)
{
    if (((is.character(prot.list) && is.character(node.names)) && (is.null(names(prot.list)) == FALSE && is.null(names(node.names)) == FALSE)) == FALSE)
    {
        default.fail("ERROR: Invalid arguments supplied to protein.list.to.gene")
    }
    
    common.prots <- intersect(names(prot.list), names(node.names))
    
    if (length(common.prots) != length(prot.list))
    {
        default.fail("ERROR: Protein name cannot be found in node.names")
    }
    
    if ((is.null(collapse.char) || (is.character(collapse.char) && length(collapse.char) == 1)) == FALSE)
    {
        default.fail("ERROR: collapse.char should either be NULL or a single character string")
    }
    
    sub.node.names <- node.names[common.prots]
    
    prot.list <- prot.list[names(sub.node.names)]
    
    names(prot.list) <- sub.node.names
    
    if (any(duplicated(names(prot.list))))
    {
        #fixed this 7-25-2012 to allow the selective inclusion of multiple protein-gene relationships
        
        dup.names <- unique(names(prot.list)[duplicated(names(prot.list))])
        
        #fixed this 7-29-2012 to allow the duplicates to be subsets of each other instead of requiring equivalence
        ok.names <- sapply(dup.names, function(x)
                           {
                                dup.elems <- prot.list[names(prot.list) == x]
                                if (is.null(collapse.char))
                                {
                                    if (length(unique(dup.elems)) > 1)
                                    {
                                        default.fail("ERROR: Invalid duplicated symbol for protein name found")
                                    }
                                    
                                    return(as.character(dup.elems[1]))
                                }
                                else
                                {
                                    split.dup.elems <- strsplit(dup.elems, collapse.char)
                                    elem.len <- sapply(split.dup.elems, length)
                                    
                                    is.max <- elem.len == max(elem.len)
                                    
                                    if (sum(is.max) > 1)
                                    {
                                        #they need to be equaivalent or an error is thrown
                                        if (all(sapply(split.dup.elems[is.max], function(y) y == split.dup.elems[is.max][[1]])) == FALSE)
                                        {
                                            default.fail("ERROR: Multiple non-equivalent duplicated items detected")
                                        }
                                        #if this goes through, choose one arbitrarily to be the max to compare the others to
                                    }
                                   
                                    which.is.max <- which(is.max)[1]
                                    
                                    if (all(sapply(split.dup.elems, function(y) all(y %in% split.dup.elems[[which.is.max]]))) == FALSE)
                                    {
                                        default.fail("ERROR: Invalid duplicated symbol for protein name found")
                                    }
                                    
                                    return(paste(split.dup.elems[[which.is.max]], collapse=collapse.char))
                                }
                           })
        
        prot.list <- prot.list[!duplicated(names(prot.list))]
        prot.list[names(ok.names)] <- as.character(ok.names)
    }
    
    return(prot.list)
}

draw.graph.variantPriorResult <- function(var.obj, graph.params)
{
    if (class(var.obj) != "variantPriorResult")
    {
        default.fail("ERROR: var.obj needs to be of class variantPriorResult")
    }
   
    if (class(graph.params) != "graphDispParams")
    {
        default.fail("ERROR: graph.params needs to be of class graphDispParams")
    }
    
    #extract the parameter object

    params <- getParameters(var.obj)
    
    ranked.query <- getRankedQueries(var.obj, in.graph=TRUE)
    
    queries.to.graph <- ranked.query[1:min(maxPlotVars(graph.params), nrow(ranked.query)), getProxSummaryID(var.obj)]
    
    hit.scores <- getProteinScores(var.obj, in.graph=TRUE)
    
    hits.to.graph <- hit.scores[1:min(maxPlotHits(graph.params), nrow(hit.scores)), getHitSummaryID(var.obj)]
    
    use.graph <- getGraph(var.obj)
    
    targ.graph <- make.target.graph(use.graph=getGraph(use.graph), hit.nodes=hits.to.graph, mut.nodes=queries.to.graph)
    
    #both label graph and collapse nodes by gene name and reassign edges if necessary

    graph.annot <- annotation(getGraph(var.obj))

    split.annots <- split(graph.annot, graph.annot[,coreID(params)])
    
    node.names <- sapply(split.annots, function(x) x[,summaryID(params)])
    
    #check for annotation consistency
    stopifnot(class(node.names) == "character")

    targ.graph <- collapse.nodes.genes(targ.graph, node.names)
    
    draw.graph(targ.graph, var.obj, node.names, graph.params)
}

##new function for computing the k-shortest paths based on implementation in the kBestShortestPaths package
#here, and probably in make.target.graph use.graph should actually be the var.obj  object
make.target.graph.kbest <- function(var.obj, hit.nodes, mut.nodes, kbest=5)
{
    if (class(var.obj) != "variantPriorResult")
    {
        default.fail("ERROR: use.graph must be an object of class igraph")
    }
    
    if (is.character(hit.nodes)==FALSE || is.character(mut.nodes) == FALSE)
    {
        default.fail("ERROR: hit.nodes and mut.nodes need to be character vectors")
    }
    
    use.graph <- getGraph(var.obj)
    
    graph.verts <- V(getGraph(use.graph))
    
    hits.in.graph <- intersect(hit.nodes, graph.verts$name)
    muts.in.graph <- intersect(mut.nodes, graph.verts$name)
    
    cat("Found", length(hits.in.graph), "target nodes and", length(muts.in.graph), "mut nodes \n")
    
    cat("Finding target--mutation distances\n")
    
    #need to adapt me to run through multiple sources and destinations maybe in the form of a data.frame
    #kb.res <- kBestShortestPathsIgraph(g=getGraph(use.graph), source.node, dest.node, best.count=kbest)
    
    
}

#started as a response to a reviewer.  Seems to be working at the moment
##code adapted from kBestShortestPaths (which appears to be broken at the moment)
##source.node and dest.node are character vectorsŒ
##need to also add kBestShortestPaths to the DESCRIPTION and NAMESPACE files and allow for incorporation of C code somehow
kBestShortestPathsIgraph <- function(g, source.node, dest.node, best.count)
{
    #edge.list is really a matrix...
    edge.list <- get.edgelist(g)
    edge.list <- edge.list[!duplicated(edge.list),]
    
    #also the reverse
    e.names <- c(sprintf("%s~%s", edge.list[,1], edge.list[,2]), sprintf("%s~%s", edge.list[,2], edge.list[,1]))
    
    e.weights <- rep(1, length(e.names))
    
    node.names <- V(g)$name
    
    output <- .Call("bozghale",
              node.names,
              length(node.names),
              e.names,
              e.weights,
              length(e.names),
              as.integer(1),
              source.node, 
              dest.node,
              as.integer(best.count), PACKAGE="kBestShortestPaths");
  return (new('Paths', best.paths=output$best.paths, path.scores=output$path.weights));
}


draw.graph <- function(targ.graph, var.obj, node.names, graph.params)
{
    if(class(targ.graph) != "graphNEL")
    {
        default.fail("ERROR: targ.graph needs to be of class graphNEL")
    }
    
    if (class(var.obj) != "variantPriorResult")
    {
        default.fail("ERROR: var.obj needs to be of class variantPriorResult")
    }
    
    if (class(graph.params) != "graphDispParams")
    {
        default.fail("ERROR: graph.params needs to be of class graphDispParams")
    }
    
    if ((is.character(node.names) && is.null(names(node.names)) == FALSE) == FALSE)
    {
        default.fail("ERROR: node.names needs to be a named character vector")
    }
    
    fill.colors <- getFillColors(var.obj, node.names, targ.nodes=nodes(targ.graph), graph.param=graph.params)
    bord.colors <- getBordColors(var.obj, node.names, targ.nodes=nodes(targ.graph), graph.param=graph.params)    

    paste.names <- getBestVarText(var.obj, node.names, targ.nodes=nodes(targ.graph), graph.params=graph.params)

    node2shape <- getNodeShapes(var.obj, node.names, targ.nodes=nodes(targ.graph), graph.params=graph.params)

    node2style <- getNodeStyles(var.obj, node.names, targ.nodes=nodes(targ.graph), graph.params=graph.params)

    #this function allows the various vectors to be modified based on the values of other vectors as specified by the user
    node.list <- nodeCombinationFunc(graph.params)(list(label=paste.names[nodes(targ.graph)], fillcolor=fill.colors[nodes(targ.graph)],
                                                      color=bord.colors[nodes(targ.graph)], shape=node2shape[nodes(targ.graph)],
                                                      style=node2style[nodes(targ.graph)]), var.obj, graph.params)
    
    node.attrs <- makeNodeAttrs(targ.graph, label=node.list$label, fillcolor=node.list$fillcolor, color=node.list$color, shape=node.list$shape, style=node.list$style)
    
    node.draw.list <- lapply(1:length(nodes(targ.graph)), function(x) drawAgNodeTG)
    names(node.draw.list) <- nodes(targ.graph)
    
    plot.target.graph(target.graph=targ.graph, node.attrs=node.attrs, node.draw.list, var.obj, graph.params)
}

make.target.graph <- function(use.graph, hit.nodes, mut.nodes)
{
    if (class(use.graph) != "igraph")
    {
        default.fail("ERROR: use.graph must be an object of class igraph")
    }
    
    if (is.character(hit.nodes)==FALSE || is.character(mut.nodes) == FALSE)
    {
        default.fail("ERROR: hit.nodes and mut.nodes need to be character vectors")
    }
    
    graph.verts <- V(use.graph)
    
    hits.in.graph <- intersect(hit.nodes, graph.verts$name)
    muts.in.graph <- intersect(mut.nodes, graph.verts$name)
    
    cat("Found", length(hits.in.graph), "target nodes and", length(muts.in.graph), "mut nodes \n")
    
    cat("Finding target--mutation distances\n")
    
    all.paths <- get.graph.sp(use.graph, from.nodes=muts.in.graph, to.nodes=hits.in.graph, path.len.restrict=NULL)
    
    #find hit-hit distances

    cat("Finding target--target distances\n")
    
    if (length(hits.in.graph) > 1)
    {
        target.paths <- get.graph.sp(use.graph, from.nodes=hits.in.graph, to.nodes=NULL,path.len.restrict=NULL)
        all.paths <- rbind(all.paths, target.paths)
    }

    #if there are not too many edges, connect the mutants

    #only connect the mutants if they are first order interactions...
    
    cat("Finding mut--mut distances\n")
    mut.paths <- get.graph.sp(use.graph, from.nodes=muts.in.graph, to.nodes=NULL, path.len.restrict=2)
    
    all.paths <- rbind(all.paths, mut.paths)

    #make a graph from the sp list

    cat("Making subgraph\n")
    
    all.paths <- all.paths[!duplicated(all.paths),]
    
    eids <- sapply(1:nrow(all.paths), function(y)
           {
                as.numeric(E(use.graph, P=as.numeric(all.paths[y,])))
           })
    
    ret.graph <- subgraph.edges(graph=use.graph, eids=eids)
    
    cat("Converting to graphNEL\n")
    
    ret.graph.nel <- igraph.to.graphNEL.slim(ret.graph)
    
    return(ret.graph.nel)
}

get.graph.sp <- function(use.graph, from.nodes, to.nodes=NULL, path.len.restrict=NULL)
{
    if (class(use.graph) != "igraph")
    {
        default.fail("ERROR: use.graph needs to be an igraph object")
    }
    
    if (is.character(from.nodes) == FALSE)
    {
        default.fail("ERROR: from.nodes needs to be a character vector")
    }
    
    if ((is.null(to.nodes) || is.character(to.nodes)) == FALSE)
    {
        default.fail("ERROR: to.nodes needs to be NULL or a character vector")
    }
    
    if ((is.null(path.len.restrict) || is.numeric(path.len.restrict)) == FALSE)
    {
        default.fail("ERROR: path.len.restrict needs to be NULL or a numeric vector")
    }
    
    graph.verts <- V(use.graph)
    graph.edges <- E(use.graph)$weight
    
    if(missing(to.nodes) || is.null(to.nodes) || is.na(to.nodes))
    {
        node.combn <- combn(from.nodes, 2)
        node.grid <- data.frame(t(node.combn), stringsAsFactors=FALSE)
        names(node.grid) <- c("from.nodes", "to.nodes")
    }
    else
    {
        node.grid <- expand.grid(from.nodes=from.nodes, to.nodes=to.nodes)
        node.grid <- node.grid[from.nodes != node.grid$to.nodes,]
    }
    
    hit.paths <- mapply(function(x, y)
                        {
                            #print(x)
                            #browser()
                            paths <- get.all.shortest.paths(use.graph, from=graph.verts[x], to = graph.verts[y], mode = "all", weights=NA)$res
                            
                            paths <- paths[!duplicated(paths)]
                            
                            #exclude situations where a protein is both a hit and mut
                            keep.paths <- sapply(paths, function(y) length(y) > 1)
                            
                            if (is.null(path.len.restrict) || is.na(path.len.restrict))
                            {
                                keep.paths <- keep.paths
                            }
                            else
                            {
                                keep.paths <- keep.paths & sapply(paths, function(y) length(y) <= path.len.restrict)
                            }
                            
                            if (sum(keep.paths) > 0)
                            {
                                paths <- paths[keep.paths]
                                
                                #assuming there are more than one valid shortest path, select the one crudely with the highest confidence
                                path.dist <- sapply(paths, function(z)
                                                    {
                                                        sum(E(use.graph, path=as.numeric(z))$weight)
                                                    })
                                
                                ret.paths <- lapply(paths[which.max(path.dist)], function(y)
                                                   {
                                                        
                                                        edge.mat <- cbind(y[seq(1, length(y)-1)], y[seq(2, length(y))])
                                                        
                                                        return(edge.mat)
                                                   })
                                
                                return(ret.paths)
                            }
                            else
                            {
                                return(NA)
                            }
                            
                        }, as.character(node.grid$from.nodes), as.character(node.grid$to.nodes), USE.NAMES=FALSE)
    
    #encode the list of lists into a better structure
    
    if (class(hit.paths) != "list")
    {
        return(matrix(nrow=0, ncol=2))
    }
    else
    {
        hit.path.dta <- do.call("rbind", hit.paths)

        return(hit.path.dta[complete.cases(hit.path.dta),])
    }
    
}


