BiocManager::install("AnnotationDbi")
BiocManager::install("org.At.tair.db")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library("AnnotationDbi")
library("clusterProfiler")
library("org.At.tair.db")
library("dplyr")
library("pathview")


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

