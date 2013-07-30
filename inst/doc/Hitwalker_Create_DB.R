### R code from vignette source 'Hitwalker_Create_DB.Rnw'

###################################################
### code chunk number 1: Hitwalker_Create_DB.Rnw:18-81 (eval = FALSE)
###################################################
## 
## stopifnot(require(Streamer))
## stopifnot(require(RSQLite))
## options(stringsAsFactors=FALSE)
## 
## #we will hold onto the header and use it at each iteration
## tab.head <- strsplit(readLines("protein.links.detailed.v9.0.txt.gz", n=1), "\\s")[[1]]
## rt.prod <- ReadTableProducer("protein.links.detailed.v9.0.txt.gz", nrows=10000,
##     header=FALSE, sep="")  
## 
## hu.tab <- data.frame()
## 
## iter <- 0
## 
## while(length(cur.tab <- yield(rt.prod)))
## {
##     if (iter %% 10000 == 0)
##     {
##         print (iter)
##     }
##     
##     #remove the header if its present in the data
##     if (all(as.character(unlist(cur.tab[1,])) == tab.head))
##     {
##         cur.tab <- cur.tab[-1,]
##     }
##     
##     names(cur.tab) <- tab.head
##     
##     cur.tab.hu <- cur.tab[grepl("^9606\\.", cur.tab$protein1) | grepl("^9606\\.",
##         cur.tab$protein2),]
##     
##     if (nrow(cur.tab.hu) > 0)
##     {
##         if (nrow(hu.tab) == 0)
##         {
##             hu.tab <- cur.tab.hu
##         }
##         else
##         {
##             hu.tab <- rbind(hu.tab, cur.tab.hu)
##         }
##     }
##     
##     iter <- iter + 1
## }
## 
## close(rt.prod)
## 
## #check to ensure that all are ensembl human proteins
## 
## stopifnot(sum(grepl("^9606\\.ENSP\\d+", hu.tab$protein1) & grepl("^9606\\.ENSP\\d+",
##     hu.tab$protein2)) == nrow(hu.tab))
## 
## #ok, then strip off the organism IDs
## 
## hu.tab$protein1 <- sub("9606\\.", "", hu.tab$protein1)
## hu.tab$protein2 <- sub("9606\\.", "", hu.tab$protein2)
## 
## db.con <- dbConnect("SQLite", "ens_b57_string_v9_annot.db")
## 
## dbWriteTable(db.con, "string_v9_db", hu.tab, row.names=FALSE)
## 


###################################################
### code chunk number 2: Hitwalker_Create_DB.Rnw:90-133 (eval = FALSE)
###################################################
## 
## #we are reading the entire file in, though a Streamer based approach similar to above
## #could be used.
## ens.gtf <- read.delim("Homo_sapiens.GRCh37.57.gtf.gz", sep="\t", header=FALSE)
## 
## #reduce the number of rows to just those potentially of interest
## sub.gtf <- subset(ens.gtf, V3 == "CDS")
## split.annots <- strsplit(sub.gtf$V9, ";")
## 
## #retrieve the annotations that we need
## annot.mat <- sapply(split.annots, function(x)
##        {
##             ens.gene <- grep("gene_id", x)
##             ens.trans <- grep("transcript_id", x)
##             ens.name <- grep("gene_name", x)
##             ens.prot <- grep("protein_id", x)
##             
##             stopifnot(length(ens.gene) == 1 && length(ens.trans) == 1 &&
##                 length(ens.name) == 1 && length(ens.prot) == 1)
##             
##             return(x[c(ens.name, ens.gene, ens.trans, ens.prot)])
##        })
## 
## annot.mat <- t(annot.mat)
## 
## #strip out all the identifier strings and spaces
## anno.subs <- c("\\s+gene_name\\s+", "\\s+gene_id\\s+", "\\s+transcript_id\\s+",
##     "\\s+protein_id\\s+")
## 
## for (i in 1:length(anno.subs))
## {
##     annot.mat[,i] <- sub(anno.subs[i], "", annot.mat[,i])
## }
## 
## #rename and write to the database
## colnames(annot.mat) <- c("symbol", "gene", "transcript", "protein")
## 
## annot.mat <- annot.mat[!duplicated(annot.mat),]
## 
## dbWriteTable(db.con, "annotation",as.data.frame(annot.mat), row.names=FALSE)
## 
## dbDisconnect(db.con)
## 


###################################################
### code chunk number 3: Hitwalker_Create_DB.Rnw:140-185
###################################################

stopifnot(require(HitWalker))
stopifnot(require(RSQLite))
stopifnot(require(HitWalkerData))

#first we adjust some of the queries to be used in forming the annotatedIGraph object
test.parm <- HitWalker:::basePriorDbParams()

#note that in order for get.protein.protein.graph to be used, both the protein1, protein2
#and weight columns need to be specified.
matrixQuery(test.parm) <- function(sample.id.list = NULL, param.obj) 
{
    query <- "SELECT protein1, protein2, combined_score/1000.0 AS
        weight FROM string_v9_db WHERE combined_score > 400;"
    return(query)
}

annotQuery(test.parm) <- function(sample.id.list = NULL, param.obj)
{
    query <- "SELECT * FROM annotation;"
}

scoreSummaryFunc(test.parm) <- HitWalker:::default.score.summary.func

variantFilterObj(test.parm) <- HitWalker:::defaultVariantFilter()

#This is where the ens_b57_string_v9_annot.db lives in
#HitWalkerData
db.con <- dbConnect("SQLite", annot.db.path())

anno.igraph <- get.protein.protein.graph(db.con, test.parm, get.annotation=TRUE,
    on.symbol.protein.dup="remove.dup")

dbDisconnect(db.con)

var.hit.db <- dbConnect("SQLite", hitwalker.db.path())

test.out <- run.prioritization.patient(var.hit.db, "08-00102", anno.igraph, test.parm)

dbDisconnect(var.hit.db)

graph.params <- makeGraphDispParams(file.name = character())

plot(test.out, graph.params)
    


