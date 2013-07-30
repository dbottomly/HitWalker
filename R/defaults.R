

#Please note that the test* functions are minimalistic toy examples, the default versions below are more realistic though not perfect
baseTestQuery <- function(table.name)
{
    query <- paste("SELECT * FROM", table.name)
    return(query)
}

testHitQuery <- function(sample.id.list, param.obj)
{
    return(baseTestQuery("hit"))
}

testVariantQuery <- function(sample.id.list, param.obj)
{
    return(baseTestQuery("variant"))
}

testPatientOverlay <- function(sample.id.list, param.obj)
{
    return(baseTestQuery("samp_overlay"))
}

#Need to modify the query in some respect maybe something like:
#SELECT * FROM xref join gene on xref_id = display_xref_id join transcript using (gene_id) join translation using (transcript_id) join translation_stable_id using (translation_id) where display_label = 'JAK3';
defaultHitQuery <- function(sample.id.list, param.obj)
{
    sample.id <- sample.id.list$ind
    query <- paste("SELECT stable_id AS protein, score, hit, (CASE WHEN score > 0 THEN 1 ELSE 0 END) AS is_score,(CASE WHEN hit = 1 THEN 1 ELSE 0 END) AS is_sirna FROM hit_score JOIN sample_experiment USING (recon_name) JOIN xref ON (xref_gene=display_label) JOIN gene ON (xref_id = display_xref_id)
                        JOIN transcript USING (gene_id) JOIN translation USING (transcript_id) JOIN translation_stable_id USING (translation_id) WHERE (sample_id=",sample.id, ") AND
                        (score > 0 OR hit = 1);", sep="")
    
    return(gsub("\\n\\s+", " ", query))
}

defaultLLSHitQuery <- function (sample.id.list, param.obj) 
    {
        sample.id <- sample.id.list$ind
        query <- paste("SELECT stable_id AS gene, score, hit, (CASE WHEN score > 0 THEN 1 ELSE 0 END) AS is_score,(CASE WHEN hit = 1 THEN 1 ELSE 0 END) AS is_sirna
                       FROM hit_score JOIN sample_experiment USING (recon_name) JOIN xref ON (xref_gene=display_label) JOIN gene ON (xref_id = display_xref_id)
                       JOIN gene_stable_id USING (gene_id) WHERE (sample_id=", sample.id, ") AND (score > 0 OR hit = 1);", sep = "")
        return(gsub("\\n\\s+", " ", query))
    }

#note that an additional table was made:

#create table protein_refseq as (select translation_stable_id.stable_id AS protein,group_concat(dbprimary_acc separator ',')
#AS refseq_id from transcript join object_xref on transcript_id = ensembl_id join xref using (xref_id) join translation using (transcript_id)
#join translation_stable_id using (translation_id) where ensembl_object_type = 'Transcript' and external_db_id = 1800 and info_type = 'DIRECT' group by stable_id);
#ALTER TABLE `lls_variants`.`protein_refseq` ADD INDEX `protein_index`(`protein`);

#also the view main_chrs was created to keep this analysis to the main chromosomes removing alternative haplotypes

defaultAnnotQuery <- function(sample.id.list, param.obj)
{
    query <- "SELECT translation_stable_id.stable_id AS protein, gene_stable_id.stable_id AS gene, display_label AS symbol, refseq_id,
    IFNULL(on_seq_cap, 0) AS on_seq_cap, IFNULL(on_sirna, 0) AS on_sirna, IFNULL(gene_target, 0) AS gene_target FROM translation_stable_id LEFT OUTER JOIN protein_refseq ON stable_id = protein
    JOIN translation USING (translation_id) JOIN transcript USING (transcript_id) JOIN gene USING (gene_id)
    JOIN gene_stable_id USING (gene_id) JOIN xref ON gene.display_xref_id = xref_id LEFT OUTER JOIN overlay_data ON gene_stable_id.stable_id = ens_id
    JOIN main_chrs ON gene.seq_region_id = main_chrs.seq_region_id WHERE gene.biotype = 'protein_coding';"
    
    return(gsub("\\n\\s+", " ", query))
}

defaultLLSAnnotQuery <- function (sample.id.list, param.obj) 
{
    query <- "SELECT gene_stable_id.stable_id AS gene, display_label AS symbol, IFNULL(on_seq_cap, 0) AS on_seq_cap, IFNULL(on_sirna, 0) AS on_sirna,
    IFNULL(gene_target, 0) AS gene_target FROM gene JOIN gene_stable_id USING (gene_id) JOIN xref ON gene.display_xref_id = xref_id LEFT OUTER JOIN
    overlay_data ON gene_stable_id.stable_id = ens_id JOIN main_chrs ON gene.seq_region_id = main_chrs.seq_region_id WHERE gene.biotype = 'protein_coding';"
    return(gsub("\\n\\s+", " ", query))
}

#default view creation code for reference
#CREATE VIEW effect_to_protein AS SELECT variant_id, variant_allele_id, trans_change, aa_change, cds_change, coding_disrupt, transcript_stable_id.stable_id
#AS transcript, translation_stable_id.stable_id AS protein, gene_stable_id.stable_id AS gene FROM allele_effect join trans_var using (trans_var_id) join variant using (variant_id)
#join variant_alleles using (variant_allele_id, variant_id) join transcript using (transcript_id) join transcript_stable_id using (transcript_id)
#join translation using (transcript_id) join translation_stable_id using (translation_id) join gene_stable_id using (gene_id) where coding_disrupt in
#('non-synonymous', 'frameshift', 'non.frameshift', 'stop.impacting');
    
defaultVariantQuery <- function(sample.id.list, param.obj)
{   
    if (missing(sample.id.list) || is.null(sample.id.list) || is.na(sample.id.list))
    {
        stop("ERROR: Need to specify sample.id.list")
    }
    
    if (missing(param.obj) || is.null(param.obj))
    {
        stop("ERROR: Need to specify param.obj")
    }
    
    query.param.list <- queryParams(param.obj)
    
    pat.num <- round(query.param.list$pat.freq * sample.id.list$num.pats)
    af.freq <- query.param.list$af.freq
    db.only.pat.num <- query.param.list$db.only.pat.num
    sample.id <- sample.id.list$ind

    main.query <- paste("SELECT * FROM sample JOIN variant_alleles USING (variant_id) JOIN variant_quality USING (pat_var_id) JOIN
                        allele_ext USING (variant_allele_id, variant_id) JOIN effect_to_protein USING (variant_id, variant_allele_id)
                        WHERE sample_id = ",sample.id," AND genotype_2 = alt_allele", sep="")
    
    filter.query <- paste("(pat_freq < ",pat.num,") AND (((dbsnp_id IS NULL) AND (in_1kg = 0)) OR
                          ((in_1kg = 1) AND (af < ",af.freq,")) OR (dbsnp_id IS NOT NULL AND in_1kg = 0 AND pat_freq = ",db.only.pat.num,"));", sep="")
    
    return(gsub("\\n\\s+", " ", paste(main.query, filter.query, sep=" AND ")))
}

defaultLLSVariantQuery <- function(sample.id.list, param.obj)
{   
        if (missing(sample.id.list) || is.null(sample.id.list) || is.na(sample.id.list))
        {
            stop("ERROR: Need to specify sample.id.list")
        }
        
        if (missing(param.obj) || is.null(param.obj))
        {
            stop("ERROR: Need to specify param.obj")
        }
        
        query.param.list <- queryParams(param.obj)
        
        pat.num <- round(query.param.list$pat.freq * sample.id.list$num.pats)
        af.freq <- query.param.list$af.freq
        db.only.pat.num <- query.param.list$db.only.pat.num
        sample.id <- sample.id.list$ind
    
        main.query <- paste("SELECT * FROM sample JOIN variant_alleles USING (variant_id) JOIN variant_quality USING (pat_var_id) JOIN
                            allele_ext USING (variant_allele_id, variant_id) JOIN effect_to_protein USING (variant_id, variant_allele_id)
                            LEFT OUTER JOIN protein_refseq USING (protein) WHERE sample_id = ",sample.id," AND genotype_2 = alt_allele", sep="")
        
        filter.query <- paste("(pat_freq < ",pat.num,") AND (((dbsnp_id IS NULL) AND (in_1kg = 0)) OR
                              ((in_1kg = 1) AND (af < ",af.freq,")) OR (dbsnp_id IS NOT NULL AND in_1kg = 0 AND pat_freq = ",db.only.pat.num,"));", sep="")
        
        return(gsub("\\n\\s+", " ", paste(main.query, filter.query, sep=" AND ")))
}

defaultPatientOverlay <- function(sample.id.list, param.obj)
{
    if (missing(sample.id.list) || is.null(sample.id.list) || is.na(sample.id.list))
    {
        stop("ERROR: Need to specify sample.id.list")
    }
    
    if (missing(param.obj) || is.null(param.obj))
    {
        stop("ERROR: Need to specify param.obj")
    }
    
    query <- paste("SELECT DISTINCT display_label AS symbol, gene_outlier, gene_drop, (ABS(gene_outlier) + ABS(gene_drop)) AS gene_sum FROM
                   overlay_data_patient JOIN gene_stable_id ON ens_id = stable_id JOIN gene USING (gene_id) JOIN xref ON display_xref_id = xref_id
                   WHERE recon_name ='", sample.id.list$recon,"' GROUP BY symbol HAVING gene_sum = MAX(gene_sum);", sep="")
    
    return(gsub("\\n\\s+", " ", query))
}

defaultGraphTransFunc <- function(use.graph, param.obj)
{
    if (class(use.graph) != "annotatedIGraph")
    {
        stop("ERROR: An annotatedIGraph object needs to be provided to defaultGraphTransFunc")
    }
    
    act.graph <- getGraph(use.graph)
    
    temp.graph <- -graph.laplacian(act.graph, normalized=TRUE, sparse=TRUE)
    diag(temp.graph) <- 0
    
    return(temp.graph)
}

defaultFuncSummaryMap <- function()
{
    funcs <- c(get.sample.overlays="summary.id", get.sample.variants="core.id", get.hit.prots="core.id", run.protocol="core.id")
    
    return(funcs)
}

defaultRwrParams <- function()
{
    return(makeRwrParams(c.val=.3, threshold=1e-10, maxit=100))
}

defaultNodeCombFunc <- function(node.list, var.obj, graph.params)
{
    #here just adjust the line styles based on whether an actual siRNA hit was observed

    style.nodes <- node.list$style
    color.nodes <- node.list$fillcolor
    
    collapse.char <- collapseChar(graph.params)
    hit.colors <- hitColors(graph.params)
    
    stopifnot(all(names(style.nodes) == names(color.nodes)))
    
    style.nodes <- mapply(function(style, color)
                          {
                            if (style != "solid" & any(color %in% hit.colors))
                            {
                                return("solid")
                            }
                            else
                            {
                                return(style)
                            }
                            
                            },style.nodes, strsplit(color.nodes, collapse.char))
    
    node.list$style <- style.nodes
    
    return(node.list)
}

#These two functions could probably be prettied up a little and maybe generalized
defaultShapeMappingFunc <- function(dta)
{
    if (all(c("gene_outlier", "gene_drop") %in% colnames(dta)) == FALSE)
    {
        return("ellipse")
    }
    
    if (all(dta$gene_outlier > 0 & dta$gene_drop == 0))
    {
        return("up_std")
    }
    else if (all(dta$gene_outlier < 0 & dta$gene_drop == 0))
    {
        return("down_std")
    }
    else if (all(dta$gene_outlier == 0 & dta$gene_drop != 0))
    {
        return("down_half")
    }
    else if (all(dta$gene_outlier < 0 & dta$gene_drop != 0))
    {
        return("both_down")
    }
    else if (all(dta$gene_outlier == 0 & dta$gene_drop == 0))
    {
        return("ellipse")
    }
    else
    {
        warning("Unexpected gene categories")
        return("rectangle")
    }
}

defaultStyleMappingFunc <- function(dta)
{
    if (all(c("on_seq_cap", "on_sirna", "gene_target") %in% colnames(dta)) == FALSE)
    {
        return("solid")
    }
    
    if (all(dta$on_seq_cap == 1 & dta$on_sirna == 0 & dta$gene_target == 0))
    {
        return("solid")
    }
    else if (all((dta$on_seq_cap == 0) & (dta$on_sirna == 1 | dta$gene_target == 1)))
    {
        return("longdash")
    }
    else if (all((dta$on_seq_cap == 1) & (dta$on_sirna == 1 | dta$gene_target == 1)))
    {
        return("longdash")
    }
    else if (all(dta$on_seq_cap == 0 & dta$on_sirna == 0 & dta$gene_target == 0))
    {
        return("dotted")
    }
    else
    {
        warning("Unexpected style categories")
        return("twodash")
    }
}

defaultLabelFuncList <- function()
{
    return(list(qa.qc.label, freq.label, rank.label))
}

basicLabelFuncList <- function()
{
    return(list(no.label, freq.label, rank.label))
}

no.label <- function(nd.x, obj, sub.prox)
{
    return("")
}

qa.qc.label <- function(nd.x, obj, sub.prox)
{
    check.label.params(nd.x, obj, sub.prox)
    
    obj.parameters <- getParameters(obj)
    
    filter.cat <- getFilterCategory(obj.parameters)
    best.filt <- getFilterBestCategory(obj.parameters)
    filter.col <- getFilterCollapseCol(obj.parameters)
    
    check.col.data.frame(filter.col, nd.x)
    
    use.cats <- filter.cat[as.character(nd.x[,filter.col])]
    
    qa.qc.label <- ifelse(use.cats == best.filt, "", paste("(",substr(use.cats, 1, 1),")", sep=""))
    
    return(qa.qc.label)
}

freq.label <- function(nd.x, obj, sub.prox)
{
    check.label.params(nd.x, obj, sub.prox)
    pat.freq.col <- patientFreqCol(getParameters(obj))
    
    check.col.data.frame(pat.freq.col, nd.x)
    
    freq.label <- paste("F:", nd.x[,pat.freq.col])
    return(freq.label)
}

rank.label <- function(nd.x, obj, sub.prox)
{
    check.label.params(nd.x, obj, sub.prox)
    rank.col <- getProxRankCol(obj)
    
    check.col.data.frame(rank.col, sub.prox)
    
    rank.label <- paste("R:", na.omit(sub.prox[,rank.col]))
}

base.score.summary.func <- function(dta, summary.col)
{
    dta[,summary.col] <- dta$score
    return(dta)
}

default.score.summary.func <- function(dta, summary.col)
{
    dta[,summary.col] <- dta$score
    
    if (all(dta$sum.score == 0))
    {
        dta[,summary.col] <- 1
    }
    else
    {
        dta[,summary.col] <- ifelse(dta$hit == 1, max(dta$score), dta[,summary.col])
    }
    
    return(dta)
}

default.summary.cols <- c("symbol", "transcript", "protein", "refseq_id", "chr", "pos", "ref", "alt", "pat_freq",
                                                                                  "dbsnp_id", "in_1kg", "af", "trans_change", "cds_change", "aa_change",
                                                                                  "coding_disrupt", "assoc.score", "type.rank", "filter", "quality.category")

defaultVariantFilter <- function()
{
    category <- c("OK", "OK", "WARN", "WARN", "WARN", "WARN")
    label <- c("PASS", "StrandBiasFilter", "GATKStandardFilter", "HARD_TO_VALIDATE", "HARD_TO_VALIDATE;GATKStandardFilter", "GATKStandardFilter;StrandBiasFilter")
    type <- c("SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "INDEL")
    
    filter.dta <- data.frame(category, label, type, stringsAsFactors=FALSE)
    
    return(makeVariantFilter(filter.dta, filter.sep=";", filter.collapse.col="filter", filter.ignore.col="variant_qual_id", filter.best.category="OK"))
}

defaultMatrixQuery <- function(sample.id.list=NULL, param.obj)
{
    query <- "SELECT protein1, protein2, combined_score/1000 AS weight from protein_links WHERE combined_score > 400;"
    
    return(query)
}

defaultLLSMatrixQuery <- function (sample.id.list = NULL, param.obj) 
{
    query <- "SELECT gene1, gene2, combined_score/1000 AS weight from gene_links WHERE combined_score > 400;"
    return(query)
}
