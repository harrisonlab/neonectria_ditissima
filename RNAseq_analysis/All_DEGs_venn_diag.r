#!/usr/bin/Rscript

# Plot a 3-way Venn diagram from a tab delimited file containing a matrix showing
# presence of DEGs between 3 timepoints post inoculation.

# This is intended to be used on the output of parse_RNA-Seq.py, after the creation of DEG lists by DeSeq

# The script also requires the colorspace package. This can be downloaded by
# opening R and running the following command:
# options(download.file.method = "wget")
# install.packages("colorspace")

#get config options
library(optparse)
library(colorspace)
library(VennDiagram, lib.loc="/home/armita/R-packages/")
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of DEGs"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

DEGs <-data.frame()
DEGs <- read.table(f, header = TRUE)

TP1=subset(DEGs, DEGs[,"X12hr"] == 1 & DEGs[,"X18hr"] == 0 & DEGs[,"X24hr"] == 0)
TP2=subset(DEGs, DEGs[,"X12hr"] == 0 & DEGs[,"X18hr"] == 1 & DEGs[,"X24hr"] == 0)
TP3=subset(DEGs, DEGs[,"X12hr"] == 0 & DEGs[,"X18hr"] == 0 & DEGs[,"X24hr"] == 1)
# orthologs=subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 1)

# area1=(nrow(nonpath) + nrow(orthologs))
# area2=(nrow(path) + nrow(orthologs))
# area3=(nrow(tomato) + nrow(orthologs))
# area3
# area2
# area1



# Print labels
# label1 <- paste('Path', ' (', area2, ')', sep="" )
# label2 <- paste('NonPath', ' (', area1, ')', sep="" )
# label3 <- paste('FoL', ' (', area3, ')', sep="" )

# Set up labels
label1 <- paste("12hr", sep="" )
label2 <- paste("18hr", sep="" )
label3 <- paste("24hr", sep="" )

n123=nrow(subset(DEGs, DEGs[,"X12hr"] == 1 & DEGs[,"X18hr"] == 1 & DEGs[,"X24hr"] == 1))
n12=n123 + nrow(subset(DEGs, DEGs[,"X12hr"] == 1 & DEGs[,"X18hr"] == 1 & DEGs[,"X24hr"] == 0))
n13=n123 + nrow(subset(DEGs, DEGs[,"X12hr"] == 1 & DEGs[,"X18hr"] == 0 & DEGs[,"X24hr"] == 1))
n23=n123 + nrow(subset(DEGs, DEGs[,"X12hr"] == 0 & DEGs[,"X18hr"] == 1 & DEGs[,"X24hr"] == 1))
summary(n12)
summary(n13)
summary(n23)
summary(n123)

area1=(nrow(TP1) + (n12 - n123) + (n13 - n123) + n123)
area2=(nrow(TP2) + (n12 - n123) + (n23 - n123) + n123)
area3=(nrow(TP3) + (n13 - n123) + (n23 - n123) + n123)
#nrow(nonpath)
nrow(TP1)
nrow(TP2)
nrow(TP3)
n12
n13
n23
n123
area1
#area1 - n12 - n13 + n123
area2
area3

pdf(o)
draw.triple.venn(area1, area2, area3,
    n12, n23, n13,
    n123,
    category = c(label1, label2, label3),
#    rep("", 4),
	    rotation = 1,
    reverse = FALSE,
    lwd = rep(2, 3),
    lty = rep("solid", 3),
    col = rep("black", 3),
    fill = c("purple", "green", "yellow"),
    alpha = rep(0.5, 3),
    label.col = rep("black", 7),
    cex = rep(1.15, 7),
    fontface = rep("plain", 7),
    fontfamily = rep("sans", 7),
    cat.pos = c(-40, 40, 180),
    cat.dist = c(0.05, 0.05, 0.025),
    cat.col = rep("black", 3),
    cat.cex = rep(1.15, 3),
    cat.fontface = rep("plain", 3),
    cat.fontfamily = rep("sans", 3),
    cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
    cat.default.pos = "outer",
    cat.prompts = FALSE,
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE, sep.dist = 0.05, offset = 0,
    )

warnings()
q()
