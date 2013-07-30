#non-exported helper functions

default.fail <- function(err.str)
{
    stop(err.str)
}

check.sample.id.list <- function(sample.id.list)
{
    if ((is.null(names(sample.id.list)) == FALSE && length(sample.id.list) == 3) == FALSE)
    {
        default.fail("ERROR: sample.id.list needs to have three elements and be named")
    }
    
    if (all(c("recon", "ind", "num.pats") %in% names(sample.id.list)) == FALSE)
    {
        default.fail("ERROR: recon, ind and num.pats need to be names in sample.id.list")
    }
    
    if (with(sample.id.list, (length(recon) == 1 && length(ind) == 1 && length(num.pats) == 1 && is.numeric(num.pats)) == FALSE))
    {
        default.fail("ERROR: mispecified elements in sample.id.list")
    }
}

check.db.params <- function(db.con, param.obj)
{
    if (inherits(db.con, "DBIConnection") == FALSE)
    {
        default.fail("ERROR: Please supply an object inherited from DBIConnection (e.g. MySQLConnection) for db.con")
    }
    
    if (class(param.obj) != "priorDbParams")
    {
        default.fail("ERROR: Please supply a priorDbParams object for param.obj")
    }
}

#exported main functions

reconcile.sample.name <- function(sample.id, db.con, param.obj)
{
    if ((is.character(sample.id) && length(sample.id) == 1) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector of length 1 for sample.id")
    }
    
    if (idTable(param.obj) %in% dbListTables(db.con) == FALSE)
    {
        default.fail(paste("ERROR: idTable", idTable(param.obj), "is not a valid table name"))
    }
    
    check.db.params(db.con, param.obj) 
    
    db.samps <- dbReadTable(db.con, idTable(param.obj))
    
    lo.cols <- setdiff(c(reconciledIdCol(param.obj), sampleIdCol(param.obj)), colnames(db.samps))
    
    if (length(lo.cols) > 0)
    {
        default.fail(paste("ERROR: Cannot find columns", paste(lo.cols, collapse=","), "in sample table"))
    }
    
    matches <- sapply(possibleIdCols(param.obj), function(x)
           {
                if (x %in% colnames(db.samps) == FALSE) default.fail(paste("ERROR:", x, "is not found in the sample table"))
                grep.match <- grep(sample.id, db.samps[,x])
                return(ifelse(length(grep.match)==1, grep.match, NA))
           })
    
    good.match <- matches[is.na(matches) == FALSE]
    
    if ((length(good.match) == 1) || (length(good.match) > 1 & length(unique(good.match)) == 1))
    {
        good.match <- good.match[1]
        
        recon <- db.samps[good.match,reconciledIdCol(param.obj)]
        ind <- db.samps[good.match,sampleIdCol(param.obj)]
        
        return(list(recon=recon, ind=ind, num.pats=length(unique(db.samps[,sampleIdCol(param.obj)]))))
    }
    else
    {
        stop(paste("ERROR: No valid unique matches were found for sample.id:", sample.id))
    }
}

get.protein.protein.graph <- function(db.con, param.obj, get.annotation=TRUE, on.symbol.protein.dup="remove.dup")
{
    check.db.params(db.con, param.obj)
    
    on.symbol.protein.dup <- match.arg(on.symbol.protein.dup)
    
    if (is.logical(get.annotation) == FALSE && length(get.annotation) == 1)
    {
        default.fail("ERROR: Need to supply a single boolean value for get.annotation")
    }
    
    db.node1 <- paste(coreID(param.obj), "1", sep="")
    db.node2 <- paste(coreID(param.obj), "2", sep="")
    
    cat("Retrieving graph structure\n")
    param.query <- retrieve.param.query(sample.id.list=NULL, db.con, param.obj, matrixQuery)
    
    if (get.annotation == TRUE)
    {
        cat("Retrieving gene annotation\n")
        annot.query <- retrieve.param.query(sample.id.list=NULL, db.con, param.obj, annotQuery)
        
        graph.prots <- unlist(param.query[,c(db.node1, db.node2)])
        
        cat(paste("There are", length(unique(graph.prots)), "nodes in the graph\n"))
        cat(paste("There are", length(unique(annot.query[,coreID(param.obj)])), "nodes in the annotation\n"))
        
        common.prots <- intersect(graph.prots, annot.query[,coreID(param.obj)])
        
        cat(paste("There are", length(common.prots), "nodes consistent between the graph and annotation\n"))
        
        param.query <- param.query[ param.query[,db.node1] %in% common.prots & param.query[,db.node2] %in% common.prots,]
        
        #retrieve the nodes where connections exist between valid nodes and make sure the annotation is consistent

        new.graph.prots <- unlist(param.query[,c(db.node1, db.node2)])
        
        annot <- annot.query[annot.query[,coreID(param.obj)] %in% new.graph.prots,]
        
        cat(paste("After removing the resulting invalid edges, creating graph with", length(unique(new.graph.prots)), "nodes\n"))
        
        #if available weight all the edges associated with specified nodes and specifed weights
        #this annotation is provided in the annotation part of the annotatedIGraph 
    
        node.weight.col <- nodeWeightCol(param.obj)
        
        if (length(node.weight.col) == 1 && node.weight.col %in% names(annot))
        {
            cat("Reweighting specified nodes\n")
           
           #assuming an undirected network
           
           param.query <- merge(param.query, annot[,c("protein", node.weight.col)], by.x=db.node1, by.y="protein", all=TRUE, incomparables=NA, sort=FALSE)
           colnames(param.query)[colnames(param.query) == node.weight.col] <- "from.weight"
           
           param.query <- merge(param.query, annot[,c("protein", node.weight.col)], by.x=db.node2, by.y="protein", all=TRUE, incomparables=NA, sort=FALSE)
           colnames(param.query)[colnames(param.query) == node.weight.col] <- "to.weight"
           
           param.query$score <- with(param.query, weight*from.weight*to.weight)
           
           param.query <- param.query[,-which(colnames(param.query) %in% c("from.weight", "to.weight"))]
        }
    }
    else
    {
        annot <- data.frame() 
    }
    
    #remove duplications in graph
    
    graph.dups <- duplicated(param.query[c(db.node1, db.node2)])
    
    if (sum(graph.dups) > 0)
    {
        cat("Found duplicated edges, choosing representative edge with highest score\n")
        
        all.edge <- paste(param.query[,db.node1], param.query[,db.node2])
        
        which.dup <- all.edge[duplicated(all.edge)]
        
        nd.graph <- param.query[all.edge %in% which.dup == FALSE,]
        
        dup.graph <- param.query[all.edge %in% which.dup == TRUE,]
        
        split.dup <- split(dup.graph, paste(dup.graph[,db.node1], dup.graph[,db.node2]))
        
        recon.dup <- do.call("rbind", lapply(split.dup, function(x)
               {
                    sum.x <- x[x$weight == max(x$weight),]
                    
                    if (nrow(sum.x) > 1)
                    {
                        return(sum.x[1,])
                    }
                    else if (nrow(sum.x) == 1)
                    {
                        return(sum.x)
                    }
                    else
                    {
                        stop("ERROR: Unexpected value as a result of the summarization of duplicated elements")
                    }
               }))
        
        param.query <- rbind(nd.graph, recon.dup)
        
        stopifnot(sum(duplicated(param.query[c(db.node1, db.node2)])) == 0)
    }
    
    use.graph <- graph.data.frame(param.query, directed=FALSE)
    
    if (nrow(annot) > length(V(use.graph)))
    {
        if (sum(duplicated(annot[,coreID(param.obj)])) > 0)
        {
            cat(paste("Found protein duplication, attempting to correct using the", on.symbol.protein.dup, "function\n"))
            
            annot <- switch(on.symbol.protein.dup, remove.dup=remove.annot.dup(annot, coreID(param.obj)))
        }
        
        if (nrow(annot) != length(V(use.graph)))
        {
            default.fail("ERROR: Corrected annotation is still does not have an equivalent number of proteins to the graph")
        }
        
        cat("Correction appears to have been successful\n")
    }
    else if (nrow(annot) < length(V(use.graph)))
    {
        default.fail("ERROR: There are fewer proteins in the annotation than in the graph")
    }
    
    return(new("annotatedIGraph", graph=use.graph, annotation=annot))
}


remove.annot.dup <- function(annot, core.id)
{
    return(annot[!duplicated(annot[,core.id]),])
}

get.sample.variants <- function(sample.id.list, db.con, param.obj)
{
    check.db.params(db.con, param.obj)
    check.sample.id.list(sample.id.list)
    
    variant.type.col <- variantTypeCol(param.obj)#variant_type
    verbose <- verbLevel(param.obj)
    
    pat.var.id <- patVarID(param.obj)
    
    summarization.col <- getSummarizationLevel(param.obj, as.character(sys.call()[1]))
    
    #note here that the point of defining a summarization.col here is to enforce that a column with a specified name is present
    
    query.res <- retrieve.param.query(sample.id.list, db.con, param.obj, variantQuery, exp.col=summarization.col)
    
    filter.values <- getFilterValues(param.obj)
    filter.sep <- getFilterSep(param.obj)
    collapse.col <- getFilterCollapseCol(param.obj)
    ignore.cols <- getFilterIgnoreCol(param.obj)
    
    if (is.null(filter.values))
    {
        return(list(var.dta=query.res, summary.id=summarization.col))
    }
    
    if (variant.type.col %in% colnames(query.res) == FALSE)
    {
        default.fail("ERROR: the value specified for variant.type.col in the supplied priorDbParams object cannot be found in the input")
    }
    
    if (collapse.col %in% colnames(query.res) ==FALSE)
    {
        default.fail("ERROR: the value specified for collapse.col in the supplied variantFilter object cannot be found in the input")
    }
    
    if (is.null(ignore.cols) == FALSE && ignore.cols %in% colnames(query.res) ==FALSE)
    {
        default.fail("ERROR: the value specified for ignore.cols in the supplied variantFilter object cannot be found in the input")
    }
    
    if (verbose == TRUE) cat("Applying additional filters\n")
    
    res.dta <- data.frame()
    
    for (i in unique(filter.values))
    {
        use.filts <- filter.values[filter.values == i]
        filt.list <- strsplit(names(use.filts), filter.sep)
        use.types <- unlist(strsplit(i, filter.sep))
        
        if (all(use.types %in% query.res[,variant.type.col]) == FALSE)
        {
            diff.types <- setdiff(use.types, query.res[,variant.type.col])
            #changed this from error to warning as it shouldn't fail if only one type or variant or the other is found
            warning(paste("Unexpected variant type(s) supplied to variantFilter:", paste(diff.types, collapse=",")))
        }
        
        res.dta <- rbind(res.dta, filter.variant.annots(query.res[query.res[,variant.type.col] %in% use.types,], filt.list,
                                                        pat.var.id, filter.sep, collapse.col, ignore.cols))
    }

    return(list(var.dta=res.dta, summary.id=summarization.col))
}

filter.variant.annots <- function(var.dta, filter.list, var.id.col, filter.sep, collapse.col="filter", ignore.cols="variant_qual_id")
{
    if (class(var.dta) != "data.frame")
    {
        default.fail("ERROR: Please supply a valid data.frame to var.dta")
    }
    
    if (class(filter.list) != "list")
    {
        default.fail("ERROR: Please supply a valid list to filter.list")
    }
    
    if ((is.character(collapse.col) && length(collapse.col) == 1 && collapse.col %in% names(var.dta)) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector of length one that is one of the columns in var.dta")
    }
    
    #ignore.cols can be missing
    if ((is.null(ignore.cols) || (is.character(ignore.cols) && all(ignore.cols %in% names(var.dta)))) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector for ignore.cols or set it to NULL")
    }
    
    if ((class(var.id.col) == "character" && length(var.id.col) == 1 && var.id.col %in% names(var.dta)) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector of length 1 with that indicates a column name of var.dta")
    }
    
    if ((class(filter.sep) == "character" && length(filter.sep) == 1) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector of length 1 for filter.sep")
    }
    
    filt.len <- sapply(filter.list, length)
    
    split.var <- split(var.dta, var.dta[,var.id.col])
    
    #attempt to filter in a more controlled manner, also attempt to preserve the names provided in the param.obj as they will be used downstream
    
    should.keep <- sapply(split.var, function(x)
                          {
                            which.hit <- which(sapply(filter.list, function(y) length(intersect(x[,collapse.col], y)) == length(union(x[,collapse.col], y))))
                            
                            if (length(which.hit) == 1)
                            {
                                return(paste(filter.list[[which.hit]], collapse=filter.sep))
                            }
                            else if (length(which.hit) == 0)
                            {
                                return(NA)
                            }
                            else
                            {
                                default.fail("ERROR: Unexpected number of matches")
                            }
                          })
    
    which.should.keep <- which(is.na(should.keep) == FALSE)
    
    keep.var <- split.var[which.should.keep]
    
    if (length(keep.var) > 0)
    {
        res.var.list <- lapply(1:length(keep.var), function(x)
           {
                temp.x <- keep.var[[x]]
                sub.x <- temp.x[,-which(colnames(temp.x) %in% c(collapse.col, ignore.cols))]
                sub.x <- sub.x[!duplicated(sub.x),]
                sub.x <- cbind(sub.x, should.keep[which.should.keep][x])
                colnames(sub.x)[ncol(sub.x)] <- collapse.col
                return(sub.x)
           })
    
        res.var.dta <- data.frame(do.call("rbind", res.var.list), stringsAsFactors=FALSE)
    }
    else
    {
        res.var.dta <- data.frame()
    }
    
    return(res.var.dta)
    
}

get.hit.prots <- function(sample.id.list, db.con, param.obj)
{
    check.db.params(db.con, param.obj)
    check.sample.id.list(sample.id.list)
    
    summary.col <- getSummarizationLevel(param.obj, as.character(sys.call()[1]))
    comb.score.name <- combScoreName(param.obj)
    
    query.res <- retrieve.param.query(sample.id.list, db.con, param.obj, hitQuery, exp.col=summary.col)
    
    sum.query.res <- scoreSummaryFunc(param.obj)(query.res, comb.score.name)
    
    if (comb.score.name %in% colnames(sum.query.res) == FALSE)
    {
        default.fail("ERROR: Need to supply a score.summary.func that returns a valid comb.score.name")
    }
    
    #assuming that the inputs will  uniformly be a named vector for the hit proteins
    hit.prot.list <- split(sum.query.res, sum.query.res[,summary.col])
    hit.prot.vec <- sapply(hit.prot.list, "[[", comb.score.name)
    
    return(list(protein.scores=hit.prot.vec, query.dta = sum.query.res, summary.id=summary.col))
}

get.sample.overlays <- function(db.con, param.obj, graph.obj, sample.id.list)
{
    if (class(param.obj) != "priorDbParams")
    {
        default.fail("ERROR: Please supply a priorDbParams object for param.obj")
    }
    
    if (class(graph.obj) != "annotatedIGraph")
    {
        default.fail("ERROR: Please supply an annotatedIGraph object for graph.obj")
    }
    
    check.sample.id.list(sample.id.list)
    
    summarization.col <- getSummarizationLevel(param.obj, as.character(sys.call()[1]))
    
    query.res <- retrieve.param.query(sample.id.list, db.con, param.obj, patientOverlay, exp.col=summarization.col, can.fail=TRUE)
    
    if (nrow(query.res) == 0)
    {
        return(graph.obj)
    }
    else
    {
        cur.annot <- annotation(graph.obj)
        
        stopifnot(summarization.col %in% colnames(cur.annot))
        
        merged.annot <- merge(cur.annot, query.res, by=summarization.col, all.x=TRUE, all.y=FALSE, incomparables=NA, sort=FALSE)
       
        new.cols <- setdiff(colnames(merged.annot), colnames(cur.annot))
        
        for (i in new.cols)
        {
            merged.annot[is.na(merged.annot[,i]),i] <- 0
        }
        
        annotation(graph.obj) <- merged.annot
        
        return(graph.obj)
    }
}

retrieve.param.query <- function(sample.id.list, db.con, param.obj, query.func, exp.col=NULL, can.fail=FALSE)
{
    check.db.params(db.con, param.obj)
    
    if (missing(sample.id.list) || is.null(sample.id.list) == FALSE)
    {
        check.sample.id.list(sample.id.list)
    }
    
    if (is.function(query.func) == FALSE)
    {
        default.fail("ERROR: Please supply a function for query.func")
    }
    
    if ((is.null(exp.col) || (is.character(exp.col) && length(exp.col) > 0)) == FALSE)
    {
        default.fail("ERROR: Please supply a character vector of length > 0 or a NULL for exp.col")
    }
    
    if ((is.logical(can.fail) && length(can.fail) == 1) == FALSE)
    {
        default.fail("ERROR: Please supply a TRUE or FALSE argument for can.fail")
    }
    
   #need to return a data.frame if no valid function (or empty function) is specfied

    verbose <- verbLevel(param.obj)
    
    query.str <- query.func(param.obj)(sample.id.list, param.obj)
    
    if (is.null(query.str))
    {
        return(data.frame())
    }
    
    if (verbose == TRUE) cat(paste(query.str, "\n"))
    
    if (verbose == TRUE) cat("Sending query...\n")
    
    query.res <- dbGetQuery(db.con, query.str)
    
    if (nrow(query.res) == 0)
    {
        if (can.fail == TRUE)
        {
            return(data.frame())
        }
        else
        {
            default.fail("ERROR: query appeared to fail")
        }
    }
    
    if (is.null(exp.col) == FALSE && all(exp.col %in% colnames(query.res)) == FALSE)
    {
        diff.cols <- setdiff(exp.col, colnames(query.res))
        stop("ERROR: Expected column(s) in result not found:", paste(diff.cols, collapse=","))
    }
    
    if (verbose == TRUE) cat("Retrieved query\n")
    
    return(query.res)
}

run.prioritization.patient <- function(db.con, sample.id, graph.obj, param.obj)
{
    #check db.con

    check.db.params(db.con, param.obj)
    
    if(is.character(sample.id) == FALSE || length(sample.id) == 0)
    {
        default.fail("ERROR: Please provide a character string for sample.id")
    }
    
    if (class(graph.obj) != "annotatedIGraph")
    {
        default.fail("ERROR: graph.obj should be an annotatedIGraph object")
    }

    #first task is attempt to reconcile the sample.id with database

    sample.id.list <- reconcile.sample.name(sample.id, db.con, param.obj)
    
    #retrieve sample specific overlay data and add to graph.obj

    graph.obj <- get.sample.overlays(db.con, param.obj, graph.obj, sample.id.list)
    
    #now that we have a decent sample id, find all the variants/transcripts associated with that sample

    samp.vars <- get.sample.variants(sample.id.list, db.con, param.obj)
    
    hit.prots <- get.hit.prots(sample.id.list, db.con, param.obj)
    
    graph.mat <- transformGraph(param.obj)(graph.obj, param.obj)
    
    prox.dta <- run.protocol(graph.mat, hit.prots$protein.scores, samp.vars$var.dta, param.obj)
    
    var.obj <- makeVariantPriorResult(prox.dta=prox.dta$res.dta, graph=graph.obj, sample=sample.id.list$recon, var.dta=samp.vars$var.dta, hit.dta=hit.prots$query.dta, param.obj=param.obj,
                                      prox.summary.id=prox.dta$summary.id, var.summary.id=samp.vars$summary.id, hit.summary.id=hit.prots$summary.id, prox.rank.col=prox.dta$rank.col,
                                      prox.type.col=prox.dta$type.id)
    
    return(var.obj)
    
}

run.protocol <- function(graph.sp.mat, seed.prots, samp.vars, param.obj)
{
    if (class(graph.sp.mat) != "dgCMatrix")
    {
        default.fail("ERROR: Please provide an object of class dgCMatrix")
    }
    
    if((is.numeric(seed.prots) && length(seed.prots) > 0 && is.null(names(seed.prots)) == FALSE) == FALSE)
    {
        default.fail("ERROR: seed.prots needs to be a named numeric vector")
    }
    
    if (class(param.obj) != "priorDbParams")
    {
        default.fail("ERROR: param.obj needs to be of class priorDbParams")
    }
    
    summary.col <- getSummarizationLevel(param.obj, as.character(sys.call()[1]))

    if ((is.data.frame(samp.vars) && summary.col %in% names(samp.vars)) == FALSE)
    {
        default.fail(paste("ERROR: samp.vars needs to be a data.frame with a column name:", summary.col))   
    }
    
    pat.prots <- unique(samp.vars[,summary.col])
    
    ranked.list.np <- random.walk(graph.sp.mat=graph.sp.mat, seed.prots=seed.prots, rwr.params=rwrParamsObj(param.obj))
    
    #sanity check to make sure the seed genes in the graph are the top ranked genes if all are equally weighted...
    
    if (all(ranked.list.np$prot.weights == ranked.list.np$prot.weights[1]))
    {
        stopifnot(all(ranked.list.np$seed.prots %in% names(ranked.list.np$prox.vector)[1:length(ranked.list.np$seed.prots)]))
    }
    
    #order pat.prots by rank

    use.pats.prots <- intersect(pat.prots, names(ranked.list.np$prox.vector))
    
    ord.pat.prots <- sort(ranked.list.np$prox.vector[use.pats.prots], decreasing=TRUE)
    
    all.id <- c(ranked.list.np$seed.prots, names(ord.pat.prots))
    
    #make a pretty table of the results
    res.dta <- data.frame(use.id=all.id, type=c(rep(seedLabel(param.obj), length(ranked.list.np$seed.prots)), rep(queryLabel(param.obj), length(ord.pat.prots))),
                          assoc.score=ranked.list.np$prox.vector[all.id], type.rank=c(rep(NA, length(ranked.list.np$seed.prots)),
                                                                                      1:length(ord.pat.prots)), stringsAsFactors=FALSE)
    
    colnames(res.dta)[colnames(res.dta) == "use.id"] <- summary.col
    
    return(list(res.dta=res.dta, summary.id=summary.col, rank.col="type.rank", type.id="type"))
}

#here seed.prots is a named vector containing the weights
random.walk <- function(graph.sp.mat, seed.prots, rwr.params)
{
    if (class(graph.sp.mat) != "dgCMatrix")
    {
        default.fail("ERROR: Please provide an object of class dgCMatrix")
    }
    
    if((is.numeric(seed.prots) && length(seed.prots) > 0 && is.null(names(seed.prots)) == FALSE) == FALSE)
    {
        default.fail("ERROR: seed.prots needs to be a named numeric vector")
    }
    
    if (class(rwr.params) != "rwrParams")
    {
        default.fail("ERROR: rwr.params needs to be of class rwrParams")
    }
    
    residue <- 1
    iter <- 1
    
    #c=.3, threshold=1e-10, maxit=100
    c.val <- getRestartProb(rwr.params)
    threshold <- getConvergeThresh(rwr.params)
    maxit <- getConvergeMaxIt(rwr.params)
    verbose <- verbLevel(rwr.params)
    
    #probability of restart...
    prox.vector <- rep(0, nrow(graph.sp.mat))
    names(prox.vector) <- rownames(graph.sp.mat)
    
    seed.genes.in.graph <- intersect(names(seed.prots), names(prox.vector))
    
    if (length(seed.genes.in.graph) == 0)
    {
        default.fail("ERROR: None of the seeds were found in the supplied graph")
    }
    
    if (verbose == TRUE) cat (paste("Using", length(seed.genes.in.graph), "seed genes\n"))
    
    #make sure the supplied values sum to 1

    if(sum(seed.prots[seed.genes.in.graph]) != 1)
    {
        if (verbose == TRUE) cat("Sum of selected seed.prots not equal to one, renormalizing...\n")
        seed.prots[seed.genes.in.graph] <- seed.prots[seed.genes.in.graph]/sum(seed.prots[seed.genes.in.graph])
    }
    
    prox.vector[seed.genes.in.graph] <- seed.prots[seed.genes.in.graph] #if the weights are being uniformly applied... 1/length(seed.genes.in.graph)
    
    restart.vector <- prox.vector
    
    while (residue > threshold && iter < maxit)
    {
        if (verbose == TRUE) cat (paste("iter:", iter, "residue:", residue, "\n"))
        old.prox.vector <-  prox.vector
        prox.vector <- as.numeric((1-c.val)*graph.sp.mat %*% prox.vector + (c.val*restart.vector))
        residue <- as.numeric(dist(rbind(prox.vector, old.prox.vector), method="euclidean"))
        iter <- iter + 1; 
    }
    
    if (iter == maxit)
    {
        warning(paste("Seed prots:", paste(seed.prots, collapse=","), "did not converge in 100 iterations..."))
    }
    
    names(prox.vector) <- rownames(graph.sp.mat)
    
    ord.prox.vector <- sort(prox.vector, decreasing=TRUE)
    
    return(list(prox.vector=ord.prox.vector, seed.prots=seed.genes.in.graph, prot.weights=seed.prots[seed.genes.in.graph]))
}

collapse.nodes.genes <- function(targ.graph, node.names)
{
    if (class(targ.graph) != "graphNEL")
    {
        default.fail("ERROR: targ.graph needs to be a graphNEL object")
    }
    
    if ((is.character(node.names) && is.null(names(node.names)) == FALSE) == FALSE)
    {
        default.fail("ERROR: node.names needs to be a named character vector")
    }
    
    edge.weights.pr <- edgeWeights(targ.graph)

    new.edge.list <- lapply(edge.weights.pr, function(x)
                            {
                                names(x) <- as.character(node.names[names(x)])
                                #deal with duplicate edge names here
                                new.vec <- sapply(split(x, names(x)), max)
                                
                                return(list(edges=names(new.vec), weights=as.numeric(new.vec)))
                            })
    names(new.edge.list) <- node.names[names(new.edge.list)]
    
    if (any(is.na(names(new.edge.list))))
    {
        default.fail("ERROR: Node named cannot be mapped to annotation, check graph creation")
    }
    
    which.multi.gene <- which(duplicated(names(new.edge.list)))
    
    if (length(which.multi.gene) > 0)
    {
        multi.name <- unique(names(new.edge.list)[which.multi.gene])
        
        fixed.dup.list <- lapply(multi.name, function(x)
               {
                    ap.list <- new.edge.list[names(new.edge.list) == x]
                    
                    ret.list <- list(edges = as.character(unlist(lapply(ap.list, function(x) x$edges))),
                                     weights = as.numeric(unlist(lapply(ap.list, function(x) x$weights))))
                    return(ret.list)
               })
            
         #only keep the maximum score for each edge if it is a duplicate

        fixed.dup.list.cor <- lapply(fixed.dup.list, function(x)
                                     {
                                        bound.edges <- data.frame(edges=x$edges, weights=x$weights, stringsAsFactors=FALSE)
                                        split.edges <- split(bound.edges, bound.edges$edges)
                                        max.weights <- sapply(split.edges, function(x) max(x$weights))
                                        
                                        return(list(edges=names(max.weights), weights=as.numeric(max.weights)))
                                     })
        
        names(fixed.dup.list.cor) <- multi.name
         #remove the names, then re-append when clean
        new.edge.list <- new.edge.list[-which(names(new.edge.list) %in% multi.name)]
        
        new.edge.list <- append(new.edge.list, fixed.dup.list.cor)
    }
    
    new.graph <-  new("graphNEL", nodes=names(new.edge.list), edgeL=new.edge.list, edgemode="undirected")
    
    return(new.graph)
}

