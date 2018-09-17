#Gene expression of Nd.

This analysis was done repeating the salmon alignment with the option --keepduplicates.

## Episode 1. Gene expression of Nd. GD vs M9. t1 vs t2.

```bash
/home/deakig/R3.4/bin/R
```

```R
setwd("/data/scratch/gomeza/")

#===============================================================================
#       Load libraries (R 3.4 required)
#===============================================================================

library("DESeq2")
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
library(tximport)
library(rjson)
library(readr)
require("pheatmap")
require(data.table)
library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggrepel")

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# First analysis. GD vs M9 and t1 vs t2. Mycelium control removed.

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/trans2gene2.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))


#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

# 1st design
design <- ~ Cultivar + Timepoint
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Library normalisation
dds <- estimateSizeFactors(dds)

# Set reference factor level
#dds$Cultivar<-factor(dds$Cultivar, levels=c("mycelium","GD","M9"))

# Deseq
dds<-DESeq(dds)
resultsNames(dds)
###
[1] "Intercept"          "Cultivar_M9_vs_GD"  "Timepoint_t2_vs_t1"
###

#===============================================================================
#       Results
#===============================================================================

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Cultivar","GD","M9"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
summary(sig.res)
###
out of 19 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 9, 47%
LFC < 0 (down)   : 10, 53%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 13)
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/GD_vs_M9.txt",sep="\t",na="",quote=F)

res <- results(dds, alpha=alpha,contrast=c("Timepoint","t1","t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
summary(sig.res)
###
out of 170 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 12, 7.1%
LFC < 0 (down)   : 158, 93%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/t1_vs_t2.txt",sep="\t",na="",quote=F)

#===============================================================================
#       Exploring and exporting results
#===============================================================================

# Sample Distances

vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Cultivar)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
#heatmap( sampleDistMatrix,
#  trace="none",  # turns off trace lines inside the heat map
#  col=colours, # use on color palette defined earlier
#  margins=c(12,12), # widens margins around plot
#  srtCol=45,
#  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:
rld <- varianceStabilizingTransformation(dds)
rld <- rlog( dds )
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/heatmap_rld.pdf")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Cultivar)
colnames(sampleDistMatrix) <- paste(rld$Timepoint)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# MA-plot
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/plotMA_vst.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Plot counts
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/plotcounts_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="Cultivar")
dev.off()

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/plotcounts2_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup=c("Cultivar","Timepoint"))
dev.off()

# PCA plots
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar","Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup=c("Cultivar","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Cultivar)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()
ggsave("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v4/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```

## Episode 2. Gene expression of Nd. Infection vs Control.

/home/deakig/R3.4/bin/R

```R

setwd("/data/scratch/gomeza/")

#===============================================================================
#       Load libraries (R 3.4 required)
#===============================================================================

library("DESeq2")
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
library(tximport)
library(rjson)
library(readr)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# Second analysis. Infection vs mycelium control.

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/trans2gene2.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

# Group column
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

# design
design <- ~ Condition
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Library normalisation
dds <- estimateSizeFactors(dds)

# Set reference factor level
dds$Condition<-factor(dds$Condition, levels=c("Control","Infected"))

# Deseq
dds<-DESeq(dds)
resultsNames(dds)
###
[1] "Intercept"                      "Condition_Infection_vs_Control"
###

#===============================================================================
#       Results
#===============================================================================

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Condition","Infection","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
summary(sig.res)
###
out of 4288 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 1958, 46%
LFC < 0 (down)   : 2330, 54%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/Infection_vs_Control.txt",sep="\t",na="",quote=F)

#===============================================================================
#       Exploring and exporting results
#===============================================================================

# Sample Distances

vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Condition)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
#heatmap( sampleDistMatrix,
#  trace="none",  # turns off trace lines inside the heat map
#  col=colours, # use on color palette defined earlier
#  margins=c(12,12), # widens margins around plot
#  srtCol=45,
#  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/heatmap_rld.pdf")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Condition)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# MA-plot
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/plotMA_vst.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Plot counts
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/plotcounts_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")
dev.off()

# PCA plotsPl

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Condition","Group"))
dev.off()

# Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Condition", "Group"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup=c("Condition", "Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()
ggsave("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```
