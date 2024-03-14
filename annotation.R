setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker/bammerge/bammerge_comp_HDLchip/')
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF(file="TAIR10_GFF3_genes_transposons.gff")
View(txdb)

#
wt_in2 <- readPeakFile("wt_in2.HDL.bed", header=F)
wt_in2peakAnnoList <- annotatePeak(wt_in2, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
wt_in2peakAnnoList

write.table(wt_in2peakAnnoList,"wt_in2.HDL_annot.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(wt_in2peakAnnoList)

#
mhd_in2 <- readPeakFile("mhd_in2.HDL.bed", header=F)
mhd_in2peakAnnoList <- annotatePeak(mhd_in2, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhd_in2peakAnnoList

write.table(mhd_in2peakAnnoList,"mhd_in2.HDL_annot.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(mhd_in2peakAnnoList)

#
wtmerge_in2 <- readPeakFile("wtmerge_in2_peaks.bed", header=F)
wtmerge_in2peakAnnoList <- annotatePeak(wtmerge_in2, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
wtmerge_in2peakAnnoList

write.table(wtmerge_in2peakAnnoList,"wtmerge_in2_annot.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(wtmerge_in2peakAnnoList)

# tss-500, 500
wtmerge_in2 <- readPeakFile("wtmerge_in2_peaks.bed", header=F)
wtmerge_in2peakAnnoList <- annotatePeak(wtmerge_in2, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
wtmerge_in2peakAnnoList

write.table(wtmerge_in2peakAnnoList,"wtmerge_in2_annot.tss500.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(wtmerge_in2peakAnnoList)

#
mhdmerge_in2 <- readPeakFile("mhdmerge_in2_peaks.bed", header=F)
mhdmerge_in2peakAnnoList <- annotatePeak(mhdmerge_in2, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhdmerge_in2peakAnnoList

write.table(mhdmerge_in2peakAnnoList,"mhdmerge_in2_annot.tss500.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(mhdmerge_in2peakAnnoList)

#
HDLmerge <- readPeakFile("HDLmerge_peaks.bed", header=F)
HDLmergepeakAnnoList <- annotatePeak(HDLmerge, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
HDLmergepeakAnnoList

write.table(HDLmergepeakAnnoList,"HDLmerge_annot.tss500.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(HDLmergepeakAnnoList)
HDL_annot <- data.frame(HDLmergepeakAnnoList@anno)
