#Salmon tool for quantifying the expression of transcripts using RNA-seq data

Note that it is designed to align to predicted transcripts and not to the whole genome

```bash
for Transcriptome in $(ls gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.cdna.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/rna_seq/*/*); do
  FileNum=$(ls $RNADir/F/*trim.fq.gz | wc -l)
  Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    for num in $(seq 1 $FileNum); do
    printf "\n"
    FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
    FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
    echo $FileF
    echo $FileR
    Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    echo "$Timepoint"
    OutDir=alignment/salmon/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
  done
done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p alignment/salmon/N.ditissima/Hg199/DeSeq2
for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene.txt
done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | grep -v -e 'GD_C1_3' -e 'GD_C2_3' -e 'GD_C3_3' -e 'M9_C1_3' -e 'M9_C2_3' -e 'M9_C3_3' -e 'Hg199_1' -e 'Hg199_2' -e 'Hg199_3'); do
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/N.ditissima/Hg199/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Hg199/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

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
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199/DeSeq2",full.names=F,recursive=F)

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

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
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
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 22 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 14, 64%
LFC < 0 (down)   : 8, 36%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 12)
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199/DeSeq2/GD_vs_M9.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2/GD_vs_M9_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2/GD_vs_M9_down.txt",sep="\t",na="",quote=F)


res <- results(dds, alpha=alpha,contrast=c("Timepoint","t1","t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 176 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 9, 5.1%
LFC < 0 (down)   : 167, 95%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199/DeSeq2/t1_vs_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2/t1_vs_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2/t1_vs_t2_down.txt",sep="\t",na="",quote=F)

#===============================================================================
#       Exploring and exporting results
#===============================================================================

# Sample Distances

vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_vst.pdf", width=12,height=12)
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
rld <- rlog(dds)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_rld.pdf")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Cultivar)
colnames(sampleDistMatrix) <- paste(rld$Timepoint)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# MA-plot
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotMA_vst.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Plot counts
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="Cultivar")
dev.off()

pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts2_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup=c("Cultivar","Timepoint"))
dev.off()

# PCA plots
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar","Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/PCA_rld.pdf")
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
ggsave("alignment/salmon/N.ditissima/Hg199/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

# Gene clustering plots

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new.pdf")
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Timepoint ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new2.pdf")
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Cultivar ]
mat <- assay(rld)[ topVarGenes, ]
#mat <- mat - rowMeans(mat)
#colnames(mat) <- paste0(rld$Cultivar,"-",rld$Timepoint)
heatmap.2(mat, trace="none", col=colors, dendrogram="column",ColSideColors=sidecols,labRow=TRUE, mar=c(10,2), scale="row")
dev.off()
```

### Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Condition)
write.table(raw_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Condition)
write.table(norm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

#robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

tpm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(tpm_counts) <- paste(colData$Group)
write.table(tpm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2/tpm_norm_counts.txt",sep="\t",na="",quote=F)
tpm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(tpm_counts) <- paste(colData$Group)
write.table(tpm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2/tpm_counts.txt",sep="\t",na="",quote=F)
```

## Episode 2. Gene expression of Nd. Infection vs Control.

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC
for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/trans2gene.txt
done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | grep -v -e 'GD_C1_3' -e 'GD_C2_3' -e 'GD_C3_3' -e 'M9_C1_3' -e 'M9_C2_3' -e 'M9_C3_3'); do
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

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
library(genefilter)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# Second analysis. Infection vs mycelium control.

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC",full.names=F,recursive=F)

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

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
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
[1] "Intercept"                     "Condition_Infected_vs_Control"
###

#===============================================================================
#       Results
#===============================================================================

res <- results(dds)
res
summary(res)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Condition","Infected","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4201 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1998, 48%
LFC < 0 (down)   : 2203, 52%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
###
write.table(sig.res,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/Infection_vs_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/Infection_vs_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/Infection_vs_Control_down.txt",sep="\t",na="",quote=F)


#===============================================================================
#       Exploring and exporting results
#===============================================================================

# Sample Distances

vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/heatmap_vst.pdf", width=12,height=12)
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
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/heatmap_rld.pdf")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Condition)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# MA-plot

pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/plotMA_vst.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Plot counts

pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/plotcounts_dds.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")
dev.off()

# PCA plotsPl

pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Condition","Group"))
dev.off()

# Plot using rlog transformation

pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Condition", "Group"))
dev.off()

#Plot using rlog transformation, showing sample names

data <- plotPCA(rld, intgroup=c("Condition", "Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()
ggsave("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

# Gene clustering plots

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/heatmap_new.pdf")
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Condition ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()
```
### Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Condition)
write.table(raw_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Condition)
write.table(norm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/normalised_counts.txt",sep="\t",na="",quote=F)

#robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

tpm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(tpm_counts) <- paste(colData$Group)
write.table(tpm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/tpm_norm_counts.txt",sep="\t",na="",quote=F)
tpm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(tpm_counts) <- paste(colData$Group)
write.table(tpm_counts,"alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/tpm_counts.txt",sep="\t",na="",quote=F)
```


## Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done

for UpFile in $(ls alignment/salmon/N.ditissima/Hg199/DeSeq2/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done
```
```
alignment/salmon/N.ditissima/Hg199/DeSeq2_IvsC/Infection_vs_Control_DEGs.txt
4104
alignment/salmon/N.ditissima/Hg199/DeSeq2/GD_vs_M9_DEGs.txt
20
alignment/salmon/N.ditissima/Hg199/DeSeq2/t1_vs_t2_DEGs.txt
176
```

## Produce a more detailed table of analyses

```bash
for GeneGff in $(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3); do
Strain=Hg199_minion
Organism=N.ditissima
Assembly=$(ls /data/scratch/gomeza/Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/*_contigs_unmasked.fa)
InterPro=$(ls /data/scratch/gomeza/gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv)
SwissProt=$(ls /data/scratch/gomeza/gene_pred/swissprot/N.ditissima/Hg199_minion/swissprot_vJul2016_tophit_parsed.tbl)
OutDir=gene_pred/annotation/v5/$Organism/$Strain
mkdir -p $OutDir
GeneFasta=$(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.pep.fasta)
TFs=$(ls analysis/transcription_factors/N.ditissima/Hg199_minion/Hg199_minion_TF_domains.tsv)
#Antismash=$(ls /data/scratch/gomeza/analysis/secondary_metabolites/antismash/Hg199_minion/fungi-dfc734cf-18aa-414d-b034-0da05c627613/Hg199_minion_antismash_secmet_genes.tsv)
SigP4=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/Hg199_minion_final_sp_no_trans_mem.aa)
effector_total=$(ls analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt)
CAZY_total=$(ls gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt)
TMHMM_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_pos_headers.txt)
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/annotation_tables
Dir1=$(ls -d /data/scratch/gomeza/alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5)
DEG_Files=$(ls \
$Dir1/Infection_vs_Control.txt \
| sed -e "s/$/ /g" | tr -d "\n")
RawCount=$(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/raw_counts.txt)
FPKM=$(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v5/fpkm_counts.txt)
$ProgDir/Nd_annotation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP4 $SigP4 --trans_mem $TMHMM_headers --TFs $TFs --effector_total $effector_total --CAZY_total $CAZY_total --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp.tsv
done
done
```
