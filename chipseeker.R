getwd()

#install arabidopsis snnotation library
BiocManager::install("ChIPseeker")
BiocManager::install("clusterProfiler")
BiocManager::install("GenomicFeatures", force = TRUE)
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart51")

keytypes(txdb)
View(txdb)

setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker')
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
library(TxDb.Athaliana.BioMart.plantsmart51)
library("org.At.tair.db")
library(ggplot2)

txdb <- makeTxDbFromGFF(file="TAIR10_GFF3_genes_transposons.gff")
View(txdb)

hdl_diffbound <- readPeakFile("hdl_diffbound.bed", header=F)
hdl_diffboundAnnoList <- annotatePeak(hdl_diffbound, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
hdl_diffboundAnnoList
hdl_diffbound_annot <- data.frame(hdl_diffboundAnnoList@anno)

write.table(hdl_diffboundAnnoList,"hdl_diffbound_annot500.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(hdl_diffboundAnnoList)

##############################################################################################################################
setwd('/data/2023_DAPseq/HDL_DAP-seq/11chipseeker/bammerge')
mhd_wt_bamcommon <- readPeakFile("mhd_wt_bamcommon.bed", header=F)
mhd_wt_bamcommon
mhd_wt_bamcommonAnnoList <- annotatePeak(mhd_wt_bamcommon, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhd_wt_bamcommonAnnoList
mhd_wt_bamcommon_annot <- data.frame(mhd_wt_bamcommonAnnoList@anno)
mhd_wt_bamcommon_annot
write.table(mhd_wt_bamcommonAnnoList,"mhd_wt_bamcommon_annot500.csv", sep="\t", col.names=T, row.names = F)
plotAnnoPie(mhd_wt_bamcommonAnnoList)

#ChIP peaks coverage plot
covplot(mhd_wt_bamcommon, weightCol="V5")
#Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrix <- getTagMatrix(mhd_wt_bamcommon, windows=promoter)
tagHeatmap(tagMatrix) +
  scale_fill_distiller(palette = "RdYlGn")
plotPeakProf2(peak = mhd_wt_bamcommon, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb,ignore_strand = F)

#plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2),
#conf = 0.95, by = "gene", type = "body", nbin = 800,
#TxDb = txdb, weightCol = "V5",ignore_strand = F)

upsetplot(mhd_wt_bamcommonAnnoList)

plotDistToTSS(mhd_wt_bamcommonAnnoList,
              title="Distribution of transcription factor-binding loci relative to TSS")


wtmerge_in2_peaks <- readPeakFile("wtmerge_in2_peaks.bed", header=F)
wtmerge_in2_peaksAnnoList <- annotatePeak(wtmerge_in2_peaks, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
wtmerge_in2_peaksAnnoList
wtmerge_in2_peaks_annot <- data.frame(wtmerge_in2_peaksAnnoList@anno)
write.table(wtmerge_in2_peaksAnnoList,"wtmerge_in2_peaks_annot500.csv", sep="\t", col.names=T, row.names = F)
plotAnnoPie(wtmerge_in2_peaksAnnoList)

mhdmerge_in2_peaks <- readPeakFile("mhdmerge_in2_peaks.bed", header=F)
mhdmerge_in2_peaksAnnoList <- annotatePeak(mhdmerge_in2_peaks, TxDb=txdb, tssRegion=c(-500, 500), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhdmerge_in2_peaksAnnoList
mhdmerge_in2_peaks_annot <- data.frame(mhdmerge_in2_peaksAnnoList@anno)
write.table(mhdmerge_in2_peaksAnnoList,"mhdmerge_in2_peaks_annot500.csv", sep="\t", col.names=T, row.names = F)
plotAnnoPie(mhdmerge_in2_peaksAnnoList)

################################################################################################################################

setwd('/data/2023_dapseq/7chipseeker/broadpeak_merged_3_replicates')
# annotation of peaks in all mhd replicates: merged
mhdpeak <- readPeakFile("mhd_mergedpeak.bed", header=F)
mhdmergepeakAnnoList <- annotatePeak(mhdpeak, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhdmergepeakAnnoList
mhd_annot <- data.frame(mhdmergepeakAnnoList@anno)

write.table(mhdmergepeakAnnoList,"mhd_mergedpeak_annot.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(mhdmergepeakAnnoList)
vennpie(mhdmergepeakAnnoList)

# annotation of peaks in all hdl replicates: merged
hdlmergepeak <- readPeakFile("hdl_mergedpeak.bed", header=F)
hdlmergepeakpeakAnnoList <- annotatePeak(hdlmergepeak, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
hdlmergepeakpeakAnnoList
hdl_annot <- data.frame(hdlmergepeakpeakAnnoList@anno)

write.table(hdlmergepeakpeakAnnoList,"hdl_mergedpeak_annot.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(hdlmergepeakpeakAnnoList)

#  annotation of peaks common in hdl and mhd: 3 replicates peaks merged
merged123peakmhd_hdl <- readPeakFile("merged123peakmhd_hdl.bed", header=F)
merged123peakmhd_hdlpeakAnnoList <- annotatePeak(merged123peakmhd_hdl, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
merged123peakmhd_hdlpeakAnnoList
merged123peakmhd_hdl_annot <- data.frame(merged123peakmhd_hdlpeakAnnoList@anno)

write.table(merged123peakmhd_hdlpeakAnnoList,"merged123peakmhd_hdl_annotatpeaks.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(merged123peakmhd_hdlpeakAnnoList)

# annotation of unique peaks in hdl: 3 replicates peaks merged
merged123peakuniq_hdl <- readPeakFile("merged123peakuniq_hdl.bed", header=F)
merged123peakuniq_hdlpeakAnnoList <- annotatePeak(merged123peakuniq_hdl, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
merged123peakuniq_hdlpeakAnnoList
merged123peakuniq_hdl_annot <- data.frame(merged123peakuniq_hdlpeakAnnoList@anno)

write.table(merged123peakuniq_hdlpeakAnnoList,"merged123peakuniq_hdl_annotatpeaks.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(merged123peakuniq_hdlpeakAnnoList)

# annotation of unique peaks in mhd: 3 replicates peaks merged
merged123peakuniq_mhd <- readPeakFile("merged123peakuniq_mhd.bed", header=F)
merged123peakuniq_mhdpeakAnnoList <- annotatePeak(merged123peakuniq_mhd, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
merged123peakuniq_mhdpeakAnnoList
merged123peakuniq_mhd_annot <- data.frame(merged123peakuniq_mhdpeakAnnoList@anno)

write.table(merged123peakuniq_mhdpeakAnnoList,"merged123peakuniq_mhd_annotatpeaks.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(merged123peakuniq_mhdpeakAnnoList)

setwd('/data/2023_dapseq/7chipseeker/peak_called_after_bammerge_3_rep')
# annotation of unique peaks in hdl after peak calling with merged bams
uniq_hdl_bammerge <- readPeakFile("uniq_hdl_bammerge.bed", header=F)
uniq_hdl_bammergepeakAnnoList <- annotatePeak(uniq_hdl_bammerge, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
uniq_hdl_bammergepeakAnnoList
uniq_hdl_bam <- data.frame(uniq_hdl_bammergepeakAnnoList@anno)

write.table(uniq_hdl_bammergepeakAnnoList,"uniq_hdl_bammerge_annotatpeaks.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(uniq_hdl_bammergepeakAnnoList)

# annotation of common peaks in mhd after peak calling with merged bams
mhd_merge_bam_peaks <- readPeakFile("mhd_merge_bam_peaks.bed", header=F)
mhd_merge_bam_peakspeakAnnoList <- annotatePeak(mhd_merge_bam_peaks, TxDb=txdb, tssRegion=c(-3000, 3000), assignGenomicAnnotation = TRUE, verbose=FALSE)
mhd_merge_bam_peakspeakAnnoList
mhd_merge_bam <- data.frame(mhd_merge_bam_peak_peaksAnnoList@anno)

write.table(mhd_merge_bam_peakspeakAnnoList,"mhd_merge_bam_peaks_annotatpeaks.csv", sep="\t", col.names=T, row.names = F)

plotAnnoPie(mhd_merge_bam_peakspeakAnnoList)
