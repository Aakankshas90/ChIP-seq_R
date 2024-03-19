BiocManager::install("ChIPQC")

## Load libraries
library(ChIPQC)

## Load sample data
samples_hdl <- read.csv("meta/hdl_dap.csv")
names(samples_hdl)
samples_hdl

## Create ChIPQC object
chipObj <- ChIPQC(samples_hdl)

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: hdlvsmhd", reportFolder="ChIPQCreport")

# Generating a QC report for a ChIP-seq single sample
hdl1bam = ChIPQCsample("data/bowtie/hdl1_sort.bam")
hdl1bam
ChIPQCreport(hdl1bam, reportName="hdl1bam", reportFolder="ChIPQCreport/hdl1bam")

hdl2bam = ChIPQCsample("data/bowtie/hdl2_sort.bam")
ChIPQCreport(hdl2bam, reportName="hdl2bam", reportFolder="ChIPQCreport/hdl2bam")

hdl3bam = ChIPQCsample("data/bowtie/hdl3_sort.bam")
hdl3bam
ChIPQCreport(hdl3bam, reportName="hdl3bam", reportFolder="ChIPQCreport/hdl3bam")

hdl4bam = ChIPQCsample("data/bowtie/hdl4_sort.bam")
ChIPQCreport(hdl4bam, reportName="hdl4bam", reportFolder="ChIPQCreport/hdl4bam")

mhd1bam = ChIPQCsample("data/bowtie/mhd1_sort.bam")
ChIPQCreport(mhd1bam, reportName="mhd1bam", reportFolder="ChIPQCreport/mhd1bam")

mhd2bam = ChIPQCsample("data/bowtie/mhd2_sort.bam")
ChIPQCreport(mhd2bam, reportName="mhd2bam", reportFolder="ChIPQCreport/mhd2bam")

mhd3bam = ChIPQCsample("data/bowtie/mhd3_sort.bam")
ChIPQCreport(mhd3bam, reportName="mhd3bam", reportFolder="ChIPQCreport/mhd3bam")

mhd4bam = ChIPQCsample("data/bowtie/mhd4_sort.bam")
ChIPQCreport(mhd4bam, reportName="mhd4bam", reportFolder="ChIPQCreport/mhd4bam")

gst2bam = ChIPQCsample("data/bowtie/gst2_sort.bam")
ChIPQCreport(gst2bam, reportName="gst2bam", reportFolder="ChIPQCreport/gst2bam")

# Constructing a ChIPQCexperiment object (exclude Chr C & M)
hdlvmhd = ChIPQC(samples_hdl, consensus=TRUE, bCount=TRUE, summits=250, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))
hdlvmhd

plotCoverageHist(hdlvmhd)
plotPrincomp(hdlvmhd,attributes=c("Condition"))

