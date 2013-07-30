### R code from vignette source 'Hitwalker_glio.Rnw'

###################################################
### code chunk number 1: Hitwalker_glio.Rnw:24-55 (eval = FALSE)
###################################################
## 
## stopifnot(require(cgdsr))
## stopifnot(require(biomaRt))
## 
## canc.genes <- read.csv("Table_1_full_2012-03-15.csv", header=TRUE, stringsAsFactors=FALSE)
## 
## #connect to the server
## cgds.con <- CGDS("http://www.cbioportal.org/public-portal/")
## 
## #first retrieve the mutation data. Can retrieve based on Entrez IDs per
## #http://www.cbioportal.org/public-portal/web_api.jsp (which cgdsr queries)
## 
## gbm.muts <- getMutationData(cgds.con, caseList="gbm_tcga_pub_all",
##     geneticProfile="gbm_tcga_pub_mutations", genes=canc.genes$GeneID)
## 
## #next we will retrieve data from the gene expression analysis (as an example)
##     #these appear to be Zscores
## 
## gbm.expr <- getProfileData(cgds.con, caseList="gbm_tcga_pub_all",
##     geneticProfile="gbm_tcga_pub_mrna", genes=canc.genes$GeneID)
##     
## #the other piece is to be able to reconcile the Entrez Gene IDs with protein IDs,
## #which in this case are Ensembl protein IDs.  One way of doing this is through
## #the biomaRt package which can provide the appropriate conversions
## 
## ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
## 
## entrez.to.ens <- getBM(mart=ensembl, attributes=c("hgnc_symbol", "entrezgene",
##     "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
##     filters="entrezgene", values=canc.genes$GeneID)
## 


###################################################
### code chunk number 2: Hitwalker_glio.Rnw:65-96
###################################################

stopifnot(require(reshape2))
stopifnot(require(RSQLite))
stopifnot(require(HitWalker))
stopifnot(require(HitWalkerData))

data(glio)

common.samps <- intersect(gbm.muts$case_id, gsub("\\.", "-", rownames(gbm.expr)))

samp.db <- data.frame(sample_id=1:length(common.samps), recon_name=common.samps,
    alt_name=make.names(common.samps))

#Note that alt_name is the conversion to R safe names

head(samp.db)

db.con <- dbConnect("SQLite", "TCGA_glio_mut_exp.db")

dbWriteTable(db.con, "sample", samp.db, row.names=FALSE)

#If we retrieve a base priorDbParams object we can modify it appropriately as
#we go using the numerous replacement methods.

base.parm <- HitWalker:::basePriorDbParams()

idTable(base.parm) <- "sample"
possibleIdCols(base.parm) <- c("recon_name", "alt_name")
reconciledIdCol(base.parm) <- "recon_name"
sampleIdCol(base.parm) <- "sample_id"



###################################################
### code chunk number 3: Hitwalker_glio.Rnw:106-147
###################################################

#Specify how the annotatedIGraph object is to be constructed 
matrixQuery(base.parm) <- function(sample.id.list = NULL, param.obj) 
{
    query <- "SELECT protein1, protein2, combined_score/1000.0 AS
        weight FROM string_v9_db WHERE combined_score > 400;"
    return(query)
}

annotQuery(base.parm) <- function(sample.id.list = NULL, param.obj)
{
    query <- "SELECT * FROM annotation;"
}

graph.con <- dbConnect("SQLite", annot.db.path())

graph.obj <- get.protein.protein.graph(graph.con, base.parm, get.annotation=TRUE,
    on.symbol.protein.dup="remove.dup")

dbDisconnect(graph.con)

kept.ens.prots <- subset(entrez.to.ens, ensembl_peptide_id %in%
    annotation(graph.obj)$protein)

gbm.muts.merge <- merge(gbm.muts, kept.ens.prots, by.x="entrez_gene_id",
    by.y="entrezgene", all=FALSE, incomparables=NA, sort=FALSE)

dbWriteTable(db.con, "variant", gbm.muts.merge, row.names=FALSE)

queryParams(base.parm) <- list()

variantQuery(base.parm) <- function(sample.id.list, param.obj)
    {
        return(paste("SELECT gene_symbol AS symbol, ensembl_peptide_id AS protein,
            case_id, mutation_type, validation_status, amino_acid_change,
            functional_impact_score FROM variant WHERE case_id = '",
            sample.id.list$recon,"'", sep=""))
    }

variantFilterObj(base.parm) <- makeVariantFilter()



###################################################
### code chunk number 4: Hitwalker_glio.Rnw:159-190
###################################################

gbm.expr.df <- melt(as.matrix(gbm.expr))
colnames(gbm.expr.df) <- c("sample", "symbol", "Zscore")

#Remove NAs and convert them to 0's
gbm.expr.df$Zscore[is.na(gbm.expr.df$Zscore)] <- 0

gbm.expr.df$symbol <- as.character(gbm.expr.df$symbol)

merged.gbm.expr <- merge(gbm.expr.df, kept.ens.prots, by.x="symbol", by.y="hgnc_symbol",
    all=FALSE, incomparables=NA, sort=FALSE)

dbWriteTable(db.con, "gene_hit", merged.gbm.expr, row.names=FALSE)

#Here we are defining genes with absolute Zscores greater than 2 to be hits
hitQuery(base.parm) <- function(sample.id.list, param.obj)
{
    return(paste("SELECT sample, Zscore,ensembl_peptide_id AS protein FROM
        gene_hit WHERE Sample = '", make.names(sample.id.list$recon) ,
        "' AND ABS(Zscore) > 2", sep=""))
}

#We also modify this function to ensure that the scores are computed from the right column
#Will take the absolute values of the Zscores as directionality is not really taken into
#account using the RWR framework
scoreSummaryFunc(base.parm) <- function (dta, summary.col) 
{
    dta[, summary.col] <- abs(dta$Zscore)
    return(dta)
}



###################################################
### code chunk number 5: Hitwalker_glio.Rnw:196-221
###################################################

int.pats <- dbGetQuery(db.con, "SELECT sample, Zscore FROM gene_hit JOIN sample ON
    (sample=alt_name) GROUP BY sample HAVING MAX(ABS(Zscore)) > 6")

int.pats

samp.res.1 <- run.prioritization.patient(db.con, int.pats$sample[1], graph.obj, base.parm)
samp.res.2 <- run.prioritization.patient(db.con, int.pats$sample[2], graph.obj, base.parm)

dbDisconnect(db.con)

summary(samp.res.1)[,c("protein", "symbol", "mutation_type", "amino_acid_change",
    "type.rank")]
summary(samp.res.2)[,c("protein", "symbol", "mutation_type", "amino_acid_change",
    "type.rank")]

#Make a basic plot but need to change the colors attributed to the gene hits
#Also, as we just want a node label indicating ranks, we have to adjust the labelFuncList
#method

graph.params <- makeGraphDispParams(file.name = character())
hitColors(graph.params) <- "red"
labelFuncList(graph.params) <- list(HitWalker:::no.label, HitWalker:::no.label,
    HitWalker:::rank.label)



###################################################
### code chunk number 6: Hitwalker_glio.Rnw:226-229
###################################################

plot(samp.res.1, graph.params)



###################################################
### code chunk number 7: Hitwalker_glio.Rnw:232-235
###################################################

plot(samp.res.2, graph.params)



