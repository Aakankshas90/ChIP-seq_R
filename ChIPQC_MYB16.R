setwd("/data/2023_DAPseq/MYB16_DAPseq/chipQC")

library("ChIPQC")
samples_MYB16 <- read.csv("MYB16.csv")
names(samples_MYB16)
samples_MYB16
MYB16Obj <- ChIPQC(samples_MYB16)
MYB16Obj

ChIPQCreport(MYB16Obj, reportName="ChIP QC report: MYB16", reportFolder="ChIPQCreport")

