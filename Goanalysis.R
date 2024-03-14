BiocManager::install("GOstats")
BiocManager::install("org.At.tair.db")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library("AnnotationDbi")
library("GOstats")
library("clusterProfiler")
library("GenomicFeatures")
library("RColorBrewer")
library("Rgraphviz")
library("org.At.tair.db")
library("dplyr")
library("pathview")
help('clusterProfiler')
help('goplot')


keytypes(txdb)
keytypes(org.At.tair.db)
setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker/bammerge/')
##############################################################################################################################################
HDLmerge_entrez <- HDL_annot$geneId

# Return the gene symbol for the set of Entrez IDs
HDLmerge_annot_db <- AnnotationDbi::select(org.At.tair.db,
                                           keys = HDLmerge_entrez,
                                           columns = c("GENENAME"),
                                           keytype = "TAIR")

# Change IDs to character type to merge
HDLmerge_annot_db$TAIR <- as.character(HDLmerge_annot_db$TAIR)

# Write to file
HDL_annot %>% left_join(HDLmerge_annot_db, by=c("geneId"="TAIR")) %>% write.table(file="HDLmerge_GO.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis: BP
HDLmerge_ego <- enrichGO(gene = HDLmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
HDLmerge_cluster_summary <- data.frame(HDLmerge_ego)
write.csv(HDLmerge_cluster_summary, "HDLmerge_BP.csv")

# Dotplot visualization
dotplot(HDLmerge_ego,  x = "Count",
        color = "qvalue",
        showCategory = 20,
        font.size = 10)
help(enrichGO)
# Run GO enrichment analysis: ALL
HDLmerge_all <- enrichGO(gene = HDLmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
HDLmerge_cluster_summary_ALL <- data.frame(HDLmerge_all)
write.csv(HDLmerge_cluster_summary_ALL, "HDLmerge_ALL.csv")

# Dotplot visualization
dotplot(HDLmerge_all,  x = "Count",
             color = "qvalue",
             showCategory = 10,
             font.size = 8, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

#################################################################################################################################
# Get the entrez IDs
mhd_wt_entrez <- mhd_wt_bamcommon_annot$geneId

# Return the gene symbol for the set of Entrez IDs
mhd_wt_annot_db <- AnnotationDbi::select(org.At.tair.db,
                                         keys = mhd_wt_entrez,
                                         columns = c("GENENAME"),
                                         keytype = "TAIR")
View(mhd_wt_annot_db)

# Change IDs to character type to merge
mhd_wt_annot_db$TAIR <- as.character(mhd_wt_annot_db$TAIR)

# Write to file
mhd_wt_bamcommon_annot %>% left_join(mhd_wt_annot_db, by=c("geneId"="TAIR")) %>% write.table(file="mhd_wt_bamcommon_GO.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis: BH refers to the Benjamini-Hochberg method
mhd_wt_bamcommon_ego <- enrichGO(gene = mhd_wt_entrez, 
                keyType = "TAIR", 
                OrgDb = org.At.tair.db, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Output results from GO analysis to a table
mhd_wt_bamcommon_cluster_summary <- data.frame(mhd_wt_bamcommon_ego)
write.csv(mhd_wt_bamcommon_cluster_summary, "mhd_wt_bamcommon_ALL.csv")

# Dotplot visualization
dotplot(mhd_wt_bamcommon_ego,  x = "Count",
        color = "qvalue",
        showCategory = 50,
        font.size = 5)
##############################################################################################################################################
# Get the entrez IDs
wtmerge_entrez <- wtmerge_in2_peaks_annot$geneId

# Return the gene symbol for the set of Entrez IDs
wtmerge_annot_db <- AnnotationDbi::select(org.At.tair.db,
                                          keys = wtmerge_entrez,
                                          columns = c("GENENAME"),
                                          keytype = "TAIR")

# Change IDs to character type to merge
wtmerge_annot_db$TAIR <- as.character(wtmerge_annot_db$TAIR)

# Write to file
wtmerge_in2_peaks_annot %>% left_join(wtmerge_annot_db, by=c("geneId"="TAIR")) %>% write.table(file="wtmerge_GO.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis: BP
wtmerge_ego <- enrichGO(gene = wtmerge_entrez, 
                        keyType = "TAIR", 
                        OrgDb = org.At.tair.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
# Output results from GO analysis to a table
wtmerge_cluster_summary <- data.frame(wtmerge_ego)
write.csv(wtmerge_cluster_summary, "wtmerge_BP.csv")

goplot(wtmerge_ego, font.size = 1, showCategory = 10,
       color = "qvalue",
       layout = "sugiyama",
       geom = "label")

# Dotplot visualization
dotplot(wtmerge_ego,  x = "Count",
        color = "qvalue",
        showCategory = 20,
        font.size = 10)

# Run GO enrichment analysis: ALL
wtmerge_all <- enrichGO(gene = wtmerge_entrez, 
                        keyType = "TAIR", 
                        OrgDb = org.At.tair.db, 
                        ont = "ALL", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
# Output results from GO analysis to a table
wtmerge_cluster_summary_all <- data.frame(wtmerge_all)
write.csv(wtmerge_cluster_summary_all, "wtmerge_ALL.csv")

# Dotplot visualization
wt<-dotplot(wtmerge_all,  x = "Count",
        color = "qvalue",
        showCategory = 10,
        font.size = 10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
##############################################################################################################################################
# Get the entrez IDs
mhdmerge_entrez <- mhdmerge_in2_peaks_annot$geneId

# Return the gene symbol for the set of Entrez IDs
mhdmerge_annot_db <- AnnotationDbi::select(org.At.tair.db,
                                           keys = mhdmerge_entrez,
                                           columns = c("GENENAME"),
                                           keytype = "TAIR")

# Change IDs to character type to merge
mhdmerge_annot_db$TAIR <- as.character(mhdmerge_annot_db$TAIR)

# Write to file
mhdmerge_in2_peaks_annot %>% left_join(mhdmerge_annot_db, by=c("geneId"="TAIR")) %>% write.table(file="mhdmerge_GO.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis: BP
mhdmerge_ego <- enrichGO(gene = mhdmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
mhdmerge_cluster_summary <- data.frame(mhdmerge_ego)
write.csv(mhdmerge_cluster_summary, "mhdmerge_BP.csv")

# Dotplot visualization
dotplot(mhdmerge_ego,  x = "Count",
        color = "qvalue",
        showCategory = 20,
        font.size = 10)
help(enrichGO)
# Run GO enrichment analysis: ALL
mhdmerge_all <- enrichGO(gene = mhdmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
mhdmerge_cluster_summary_ALL <- data.frame(mhdmerge_all)
write.csv(mhdmerge_cluster_summary_ALL, "mhdmerge_ALL.csv")

# Dotplot visualization
mhd<-dotplot(mhdmerge_all,  x = "Count",
        color = "qvalue",
        showCategory = 20,
        font.size = 10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

###################################################################FDR#########################################################################
setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker/bammerge/FDR')
# Get the entrez IDs
mhdmerge_entrez <- mhdmerge_in2_peaks_annot$geneId

# Return the gene symbol for the set of Entrez IDs
mhdmerge_annot_db <- AnnotationDbi::select(org.At.tair.db,
                                           keys = mhdmerge_entrez,
                                           columns = c("GENENAME"),
                                           keytype = "TAIR")

# Change IDs to character type to merge
mhdmerge_annot_db$TAIR <- as.character(mhdmerge_annot_db$TAIR)

# Write to file
mhdmerge_in2_peaks_annot %>% left_join(mhdmerge_annot_db, by=c("geneId"="TAIR")) %>% write.table(file="mhdmerge_GO.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis: BP
mhdmerge_ego <- enrichGO(gene = mhdmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "BP", 
                         pAdjustMethod = "fdr", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
mhdmerge_cluster_summary <- data.frame(mhdmerge_ego)
write.csv(mhdmerge_cluster_summary, "mhdmerge_BP.csv")

# Dotplot visualization
dotplot(mhdmerge_ego,  x = "Count",
        color = "p.adjust",
        showCategory = 20,
        font.size = 10)
help(enrichGO)
# Run GO enrichment analysis: ALL
mhdmerge_all <- enrichGO(gene = mhdmerge_entrez, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "ALL", 
                         pAdjustMethod = "fdr", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
# Output results from GO analysis to a table
mhdmerge_cluster_summary_ALL <- data.frame(mhdmerge_all)
write.csv(mhdmerge_cluster_summary_ALL, "mhdmerge_ALL.csv")

# Dotplot visualization
dotplot(mhdmerge_all,  x = "Count",
        color = "p.adjust",
        showCategory = 20,
        font.size = 10)

##############################################################################################################################################
# Dotplot visualization: combined


dotplot(mhdmerge_all,  x = "Count",
        color = "qvalue",
        showCategory = 20,
        font.size = 10, split="ONTOLOGY")
##############################################################################################################################################
getwd()
setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker/bammerge')
list.files()
keytypes("org.At.tair.db")
###########################################################################################################################
# Get GOstats results for GO categories and KEGG pathways
calcGOstats <- function(universeIDs, targetIDs, outGO, outKEGG)
{
  # universeIDs : universe ids
  # targetIDs : target ids
  # outFile : file with significant GO categories
  # outKEGG : significant KEGG pathways
  # Specify parameters for hypergeometric test against GO categories biological process
  paramsBP <- new("mhd_wt_bamcommon_annot",
                  geneIds=targetIDs,
                  universeGeneIds=universeIDs,
                  annotation="org.At.tair.db",
                  ontology="ALL",
                  pvalueCutoff=0.05,
                  conditional=FALSE,
                  testDirection="over")
  
  # Run enrichment against Go categories - Biological process
  hgOver <- hyperGTest(paramsBP)
  GOtab <- summary(hgOver)
  write.table(GOtab, file=outGO, sep="\t", col.names=T, row.names=F, quote=F)
  
  paramsKEGG <- new("KEGGHyperGParams",
                    geneIds=targetIDs,
                    universeGeneIds=universeIDs,
                    annotation="org.At.tair.db",
                    pvalueCutoff=0.1,
                    categoryName="KEGG",
                    testDirection="over")
  
  # Run enrichment against Go categories - Biological process
  hgOverKEGG <- hyperGTest(paramsKEGG)
  KEGGtab <- summary(hgOverKEGG)
  
  # Retrieve gene ids found in enriched pathways
  keggMaps <- lapply(KEGGtab$KEGGID, function(x) unlist(mget(x, org.At.tairPATH2TAIR)))
  
  # subset the selected genes based on those in the mappings
  keggSelected <- lapply(keggMaps, function(x) targetIDs[targetIDs %in% x])
  
  # join together with a semicolon to make up the last column
  KEGGtab$inGenes <- unlist(lapply(keggSelected, function(x) paste(x, collapse=";")))
  write.table(KEGGtab, file=outKEGG, sep="\t", col.names=T, row.names=F, quote=F)
}

############################################################################################################################
# Calculate for Crysp response
# Read the results file
res <- read.table("Crysp_1h_vs_PBS_1h_annot.txt", sep="\t", header=T, fill = T, quote = "")
universeIDs <- as.character(res$gene_id)
# Up-regulated in Crysp
up <- res[res$prob > 0.95 & res$M > 0,]
upIDs <- as.character(unique(up$gene_id))
upIDs <- upIDs[!is.na(upIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = upIDs, outGO = "UP_in_Crysp_GO.txt", outKEGG = "UP_in_Crysp_KEGG.txt" ) 

# Down-regulated in Crysp
down <- res[res$prob > 0.95 & res$M < 0,]
downIDs <- as.character(unique(down$gene_id))
downIDs <- downIDs[!is.na(downIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = downIDs, outGO = "DOWN_in_Crysp_GO.txt", outKEGG = "DOWN_in_Crysp_KEGG.txt" ) 

# all significant in Crysp
all <- res[res$prob > 0.95,]
allIDs <- as.character(unique(all$gene_id))
allIDs <- allIDs[!is.na(allIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = allIDs, outGO = "AllSig_in_Crysp_GO.txt", outKEGG = "AllSig_in_Crysp_KEGG.txt" ) 

###############################################################################################################################
# Calculate for Gyr response
# Read the results file
res <- read.table("Gyr_1h_vs_PBS_1h_annot.txt", sep="\t", header=T, fill = T, quote = "")
universeIDs <- as.character(res$gene_id)
# Up-regulated in Gyr
up <- res[res$prob > 0.95 & res$M > 0,]
upIDs <- as.character(unique(up$gene_id))
upIDs <- upIDs[!is.na(upIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = upIDs, outGO = "UP_in_Gyr_GO.txt", outKEGG = "UP_in_Gyr_KEGG.txt" ) 

# Down-regulated in Gyr
down <- res[res$prob > 0.95 & res$M < 0,]
downIDs <- as.character(unique(down$gene_id))
downIDs <- downIDs[!is.na(downIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = downIDs, outGO = "DOWN_in_Gyr_GO.txt", outKEGG = "DOWN_in_Gyr_KEGG.txt" ) 

# all significant in Gyr
all <- res[res$prob > 0.95,]
allIDs <- as.character(unique(all$gene_id))
allIDs <- allIDs[!is.na(allIDs)]
calcGOstats(universeIDs = universeIDs, targetIDs = allIDs, outGO = "AllSig_in_Gyr_GO.txt", outKEGG = "AllSig_in_Gyr_KEGG.txt" ) 

###########################################################################################################################
# draw pathway
pathway_ids <- c("04626", "00592", "04075", "04144", "04141", "00040", "00270", "04070")

# read expression data from Crysp response
crysp_exp <- read.table("expression_data/expressions_vsd.txt", header = T, sep = "\t")
crysp_exp_diff <- crysp_exp[,c(3,4)] - crysp_exp[,c(1,2)]
crysp_exp_diff <- as.matrix(crysp_exp_diff)

# read expression data from Gyr response
gyr_exp <- read.table("../../expression_data/expressions_vsd.txt", header = T, sep = "\t")
gyr_exp_diff <- gyr_exp[,c(5,6)] - gyr_exp[,c(1,2)]
gyr_exp_diff <- as.matrix(gyr_exp_diff)

setwd("../Gyr_pathway_images/")
sapply(pathway_ids, function(pid) pathview(gene.data = gyr_exp_diff, 
                                           pathway.id = pid, 
                                           species = "ath",
                                           kegg.native=T, 
                                           sign.pos="bottomleft",
                                           gene.annotpkg="org.At.tair.db",
                                           out.suffix = pid,
                                           gene.idtype = "KEGG"))
