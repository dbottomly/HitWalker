### R code from vignette source 'Hitwalker_Add_Metadata.Rnw'

###################################################
### code chunk number 1: Hitwalker_Add_Metadata.Rnw:19-82
###################################################

stopifnot(require(HitWalker))
stopifnot(require(RSQLite))
stopifnot(require(HitWalkerData))

#Copy the database to the current directory so it is not overwritten each time.

local.annot.db <- file.path(getwd(), basename(annot.db.path()))

stopifnot(file.copy(from=annot.db.path(), to=local.annot.db))

graph.con <- dbConnect("SQLite", local.annot.db)

#Here we add an additional table which contains two types of metadata.  The column
#'cap_probes' indicates the presence of sequence capture probes
#The column 'assay_targets' indicates the genes that are targets of functional assays.  

symbol <- c('EPHA4','ZAK','PIK3CB','AC084035.1','CBL','FRK','FYN','MAP2K6','KHDRBS1',
    'PTEN','PRKCE','MAPK14','IRS4','STAT5A','RAC1','DDR1','SRC','JAK3','FLT3')
cap_probes <- c(1,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,1,1)
assay_targets <- c(1,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1)

fixed.annot <- data.frame(symbol=symbol, cap_probes=cap_probes,
    assay_targets=assay_targets)

dbWriteTable(graph.con, "fixed_annots", fixed.annot, row.names=FALSE, overwrite=TRUE)

#Put together a priorDbParams object from the basic one that is supplied

test.parm <- HitWalker:::basePriorDbParams()
    
matrixQuery(test.parm) <- function(sample.id.list = NULL, param.obj) 
{
    query <- "SELECT protein1, protein2, combined_score/1000.0 AS weight FROM
        string_v9_db WHERE combined_score > 400;"
    return(query)
}

#now we can modify the annotation retrieval query to use this new table

annotQuery(test.parm) <- function(sample.id.list = NULL, param.obj)
{
    query <- "SELECT symbol, gene, transcript, protein, IFNULL(cap_probes,0) AS
        cap_probes, IFNULL(assay_targets,0) AS assay_targets FROM annotation
        LEFT OUTER JOIN fixed_annots USING (symbol);"
}

scoreSummaryFunc(test.parm) <- HitWalker:::default.score.summary.func
variantFilterObj(test.parm) <- HitWalker:::defaultVariantFilter()

anno.igraph <- get.protein.protein.graph(graph.con, test.parm, get.annotation=TRUE,
    on.symbol.protein.dup="remove.dup")

head(annotation(anno.igraph))

#Alternatively, this data could be added to the anno.igraph through modification of the
#data.frame retrieved through the annotation method. It is safer to go through the
#database at this time though.

dbDisconnect(graph.con)

stopifnot(file.remove(local.annot.db))



###################################################
### code chunk number 2: Hitwalker_Add_Metadata.Rnw:90-102
###################################################

#First prioritize the variants as usual
var.hit.db <- dbConnect("SQLite", hitwalker.db.path())
test.out <- run.prioritization.patient(var.hit.db, "08-00102", anno.igraph, test.parm)
dbDisconnect(var.hit.db)

#The metadata is perpetuated to the data.frame resulting from the summary method

summary(test.out)[1:5,c("gene", "symbol", "type.rank", "cap_probes", "assay_targets")]

#If we plot it using the default parameters nothing special happens



###################################################
### code chunk number 3: Hitwalker_Add_Metadata.Rnw:107-112
###################################################

graph.params <- makeGraphDispParams(file.name = character())
    
plot(test.out, graph.params)



###################################################
### code chunk number 4: Hitwalker_Add_Metadata.Rnw:115-145
###################################################

#However, if we modify the mapping functions, we can change the shape and/or borders

styleFunc(graph.params) <- function(dta)
{
    if (all(c("cap_probes", "assay_targets") %in% colnames(dta)) == 
        FALSE) {
        return("solid")
    }
    else if (all(dta$cap_probes == 1 & dta$assay_targets == 0))
    {
        return("solid")
    }
    else if (all((dta$cap_probes %in% 0:1) & (dta$assay_targets == 1)))
    {
        return("longdash")
    }
    else if (all((dta$cap_probes == 0) & (dta$assay_targets == 0)))
    {
        return("dotted")
    }
    else
    {
        warning("Unexpected style categories")
        return("twodash")
    }
}

plot(test.out, graph.params)



###################################################
### code chunk number 5: Hitwalker_Add_Metadata.Rnw:161-200
###################################################

#Note that we add the sample-specific variables into the same database as the rest of the
#sample data as it will be retrieved through the same connection.
#Again we will make a local copy of the database

local.hitwalker.db <- file.path(getwd(), basename(hitwalker.db.path()))

stopifnot(file.copy(from=hitwalker.db.path(), to=local.hitwalker.db))

var.hit.db <- dbConnect("SQLite", local.hitwalker.db)

symbol <- c('EPHA4','ZAK','PIK3CB','AC084035.1','CBL','FRK','FYN','MAP2K6','KHDRBS1',
    'PTEN','PRKCE','MAPK14','IRS4','STAT5A','RAC1','DDR1','SRC','JAK3','FLT3')
triangle <- c(1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1)
half_tri <- c(1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1)

samp.annot <- data.frame(recon_name="08-00102", symbol=symbol, triangle=triangle,
    half_tri=half_tri)

dbWriteTable(var.hit.db, "samp_overlay", samp.annot, row.names=FALSE, overwrite=TRUE)

#Although not necessary in this case, this would typically be a situation where it would
#be necessary to use the parameters of this function such as sample.id.list.
patientOverlay(test.parm) <- function(sample.id.list, param.obj)
{
    query <- paste("SELECT symbol, triangle, half_tri FROM samp_overlay WHERE
        recon_name = '", sample.id.list$recon, "'", sep="")
}

test.out <- run.prioritization.patient(var.hit.db, "08-00102", anno.igraph, test.parm)

#The new annotation is present after prioritization
annotation(getGraph(test.out))[1:5,c("protein", "cap_probes", "assay_targets", "triangle",
    "half_tri")]

dbDisconnect(var.hit.db)

stopifnot(file.remove(local.hitwalker.db))



###################################################
### code chunk number 6: Hitwalker_Add_Metadata.Rnw:209-240
###################################################

shapeFunc(graph.params) <- function(dta)
{
    if (all(c("triangle", "half_tri") %in% colnames(dta)) == FALSE) {
        return("ellipse")
    }
    else if (all(dta$triangle == 1 & dta$half_tri == 0))
    {
        return("down_std")
    }
    else if (all(dta$triangle == 0 & dta$half_tri == 1))
    {
        return("down_half")
    }
    else if (all(dta$triangle == 1 & dta$half_tri == 1))
    {
        return("both_down")
    }
    else if (all(dta$triangle == 0 & dta$half_tri == 0))
    {
        return("ellipse")
    }
    else
    {
        warning("Unexpected gene categories")
        return("rectangle")
    }
}

plot(test.out, graph.params)



