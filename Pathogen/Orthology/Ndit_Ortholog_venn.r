#!/usr/bin/Rscript

# Plot a 5-way Venn diagram from a tab delimited file containing a matrix showing
 # presence /absence of orthogroups within 5 isolates.

 # This is intended to be used on the output of the orthoMCL pipeline following
 # building of the matrix using:
 # ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2tab.py

# This script requires the optparse R package. This can be downloaded by opening
# R and running the following command:
# install.packages("optparse",repos="http://cran.uk.r-project.org")
# When given the option, install this package to a local library.

# The script also requires the colorspace package. This can be downloaded by
# opening R and running the following command:
# install.packages(colorspace)

#get config options
library(optparse)
library(colorspace)
library(VennDiagram, lib.loc="/home/armita/R-packages/")
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of presence of orthogroups"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
#    make_option("--maxrf", type="double", default=0.2, help="max rf to consider as linked"),
#    make_option("--minlod", type="double", default=20.0, help="min LOD to consider as linked")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

orthotabs <-data.frame()
orthotabs <- read.table(f)
df1 <- t(orthotabs)
summary(df1)

area1=sum(df1[, "Ag02", "Ag04", "Ag05", "Ag06", "H199", "R09", "R602", "R603", "A11A", "A11B", "A11C", "Ag08", "A09A", "R37", "R39", "R41", "R42", "R45", "199R"])
area2=sum(df1[, "BGV3", "P112","OPC3"])
area3=sum(df1[, "R68])
area4=sum(df1[, "ND8", "ND9"])
area5=sum(df1[, "R305","R324"])

#print(area1, area2, area3, area4, area5)

colname1 <- paste("NorthWest")
colname2 <- paste("Spain")
colname3 <- paste("Italy")
colname4 <- paste("Brazil")
colname5 <- paste("NewZealand")

# Full labels
# label1 <- paste(colname1, ' (20869 genes in ', area1, ' ortholog groups)', sep="" )
# label2 <- paste(colname2, ' (17787 genes in ', area2, ' ortholog groups)', sep="" )
# label3 <- paste(colname3, ' (20822 genes in ', area3, ' ortholog groups)', sep="" )
# label4 <- paste(colname4, ' (19805 genes in ', area4, ' ortholog groups)', sep="" )
# label5 <- paste(colname5, ' (26584 genes in ', area5, ' ortholog groups)', sep="" )

# Species abreviation labels
#label1 <- paste(colname1, sep="" )
#label2 <- paste(colname2, sep="" )
#label3 <- paste(colname3, sep="" )
#label4 <- paste(colname4, sep="" )
#label5 <- paste(colname5, sep="" )

# Full species name labels
label1 <- paste("NorthWest2")
label2 <- paste("Spain2")
label3 <- paste("Italy2")
label4 <- paste("Brazil2")
label5 <- paste("NewZealand2")

# No label
#label1 <- paste("", sep="" )
#label2 <- paste("", sep="" )
#label3 <- paste("", sep="" )
#label4 <- paste("", sep="" )
#label5 <- paste("", sep="" )

n12=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1))
n13=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Italy"] == 1))
n14=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Brazil"] == 1))
n15=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"NewZealand"] == 1))
n23=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Italy"] == 1))
n24=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Brazil"] == 1))
n25=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"NewZealand"] == 1))
n34=nrow(subset(df1, df1[,"Italy"] == 1 & df1[,"Brazil"] == 1))
n35=nrow(subset(df1, df1[,"Italy"] == 1 & df1[,"NewZealand"] == 1))
n45=nrow(subset(df1, df1[,"Italy"] == 1 & df1[,"NewZealand"] == 1))
n123=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Italy"] == 1))
n124=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Brazil"] == 1))
n125=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"NewZealand"] == 1))
n134=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1))
n135=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Italy"] == 1 & df1[,"NewZealand"] == 1))
n145=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n234=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1))
n235=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"NewZealand"] == 1))
n245=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n345=nrow(subset(df1, df1[,"Italy"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n1234=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1))
n1235=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"NewZealand"] == 1))
n1245=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n1345=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n2345=nrow(subset(df1, df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))
n12345=nrow(subset(df1, df1[,"NorthWest"] == 1 & df1[,"Spain"] == 1 & df1[,"Italy"] == 1 & df1[,"Brazil"] == 1 & df1[,"NewZealand"] == 1))

summary(n12)
summary(n123)
summary(n1234)
summary(n12345)

pdf(o)
draw.quintuple.venn(
  area1, area2, area3, area4, area5,
  n12, n13, n14, n15, n23, n24, n25, n34, n35, n45,
  n123, n124, n125, n134, n135, n145, n234, n235, n245, n345,
  n1234, n1235, n1245, n1345, n2345,
  n12345,
  category = c(label1, label2, label3, label4, label5),
  lwd = rep(2, 5),
	lty = rep("solid", 5),
  col = rep("black", 5),
  fill = c(rainbow_hcl(5)),
  alpha = rep(0.5, 5),
  label.col = rep("black", 31),
  cex = rep(1, 31),
  fontface = rep("plain", 31),
  fontfamily = rep("serif", 31),
  cat.pos = c(0, 310, 215, 145, 60),
  cat.dist = rep(0.25, 5),
  cat.col = rep("black", 5),
  cat.cex = rep(1, 5),
  cat.fontface = rep("italic", 5),
  cat.fontfamily = rep("serif", 5),
  cat.just = rep(list(c(0.5, 0.5)), 5),
  rotation.degree = 0,
  rotation.centre = c(0.5, 0.5),
  ind = TRUE,
  margin = 0.15
)

dev.off()

singles = df1[grepl("single*", rownames(df1)), ]
uniq_1=sum(singles[, colname1])
uniq_2=sum(singles[, colname2])
uniq_3=sum(singles[, colname3])
uniq_4=sum(singles[, colname4])
uniq_5=sum(singles[, colname5])
orthogroups = df1[grepl("orthogroup*", rownames(df1)), ]
inpara_1 = sum(orthogroups[,colname1] == 1 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 0 & orthogroups[,colname5] == 0)
inpara_2 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 1 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 0 & orthogroups[,colname5] == 0)
inpara_3 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 1 & orthogroups[,colname4] == 0 & orthogroups[,colname5] == 0)
inpara_4 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 1 & orthogroups[,colname5] == 0)
inpara_5 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 0 & orthogroups[,colname5] == 1)
colname1
area1
uniq_1
inpara_1
colname2
area2
uniq_2
inpara_2
colname3
area3
uniq_3
inpara_3
colname4
area4
uniq_4
inpara_4
colname5
area5
uniq_5
inpara_5

warnings()
q()
