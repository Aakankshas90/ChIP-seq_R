getwd()
library(ggplot2)
library(VennDiagram)

overrideTriple=T
draw.triple.venn(area1=1920, area2=8140, area3=2634, n12=1092, n23=1139, n13=295, n123=236, 
                 category = c("hdl1", "hdl2", "hdl3"), lty = rep("solid", 3), 
                 col = c("cyan", "pink", "violet"), fill = c("cyan", "pink", "violet"), 
                 alpha = rep(0.7, 3), label.col = rep("black", 7), cex = rep(1, 7),
                 fontface = rep("plain", 7), fontfamily = rep("Calibri", 7),
                 cat.pos = c(-40, 40, 180),cat.dist = c(0.05, 0.05, 0.025),
                 cat.col = rep("black", 3), cat.cex = rep(1, 3), cat.fontface = rep("plain", 3),
                 cat.fontfamily = rep("serif", 3), cat.just = list(c(0.5, 1), c(0.5, 1),c(0.5, 0)),
                 cat.default.pos = "outer")


draw.triple.venn(area1=419, area2=3729, area3=75, n12=124, n23=49, n13=33, n123=30, 
                 category = c("mhd1", "mhd2", "mhd3"), lty = rep("solid", 3), 
                 col = c("cyan", "pink", "violet"), fill = c("cyan", "pink", "violet"), 
                 alpha = rep(0.8, 3), label.col = rep("black", 7), cex = rep(0.5, 7), 
                 fontface = rep("plain", 7), fontfamily = rep("Calibri", 7),
                 cat.pos = c(-40, 40, 180), cat.dist = c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
                 cat.cex = rep(1, 3), cat.fontface = rep("plain", 3),
                 cat.fontfamily = rep("serif", 3), cat.just =
                   list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
                 = "outer")

draw.pairwise.venn(area1=6478, area2=73, cross.area=21, category = c("HDL ampDAP", "HDL ChIP"), 
                   col = c("cyan", "pink"), fill = c("cyan", "pink"))
                   
