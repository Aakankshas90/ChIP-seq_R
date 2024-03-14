#GO dotplot
setwd("D:/Hamsa/")
library(ggplot2)
theme_classic()
datad=read.csv("D:/Hamsa/downreg.csv")


head(datad)
hdown <-ggplot(datad, aes(x = Sample, y = Description, size =Number)) +
  geom_point(aes(color= log10(FDR))) + scale_x_discrete(limits = c('C2-C2W', 'E2-E2W', 'Z2-Z2W', 'C4-C4W', 'E4-E4W', 'Z4-Z4W'))

downreg<- hdown + theme(axis.text = element_text(colour = "Black", size = rel(0.5)))

downreg

ggsave("downreg.png", units="in", width=6, height=6, dpi=300)