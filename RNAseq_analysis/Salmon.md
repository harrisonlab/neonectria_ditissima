#Salmon tool for quantifying the expression of transcripts using RNA-seq data

Note that it is designed to align to predicted transcripts and not to the whole genome

Nd transcripts

```bash
for Transcriptome in $(ls gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta); do
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
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
  done
done
```
Apple transcripts

```bash
for Transcriptome in $(ls apple_genome/transcripts_modified3.fasta); do
#Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
#Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/rna_seq/*/* | grep -v "mycelium"); do
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
    OutDir=alignment/salmon/N.ditissima/Cultivar/$Prefix/$Timepoint/$Sample_Name
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
  done
done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p alignment/N.ditissima/Hg199_minion/salmon/DeSeq2
for File in $(ls alignment/salmon/*/Hg199_minion/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/trans2gene.txt
done
# Put files in a convenient location for DeSeq. Analysis was not performed on
# Apple control samples.
for File in $(ls alignment/salmon/*/Hg199_minion/*/*/*/quant.sf | grep -v -e 'GD_C1_3' -e 'GD_C2_3' -e 'GD_C3_3' -e 'M9_C1_3' -e 'M9_C2_3' -e 'M9_C3_3'); do
  # cat $File | grep 'g6.t1'
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

```bash
mkdir -p alignment/salmon/N.ditissima/Cultivar/DeSeq2
for File in $(ls alignment/salmon/*/Cultivar/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene.txt
done
# Put files in a convenient location for DeSeq. Analysis was not performed on
# mycelium control samples.
for File in $(ls alignment/salmon/*/*/*/*/quant.sf | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_'); do
  # cat $File | grep 'g6.t1'
  Prefix=$(echo $File | cut -f6 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

```bash
export R_LIBS=/home/armita/R/x86_64-pc-linux-gnu-library/3.4
/home/deakig/R3.4/bin/R
```

Method 1: Gene expression of Nd.

```R

setwd("/data/scratch/gomeza/")

#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/RNA-seq/scripts/myfunctions")
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

library(tximport)
library(rjson)
library(readr)

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files	    
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	   

# get the sample names from the folders	    
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))


#==========================================================================================
#       Read sample metadata and annotations
#=========================================================================================

# Read sample metadata

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

colData <- colData[!(colData$Sample=="M9_C1_3"),]      
colData <- colData[!(colData$Sample=="M9_C2_3"),]      
colData <- colData[!(colData$Sample=="M9_C3_3"),]      
colData <- colData[!(colData$Sample=="GD_C1_3"),]      
colData <- colData[!(colData$Sample=="GD_C2_3"),]      
colData <- colData[!(colData$Sample=="GD_C3_3"),]      

dds <- DESeqDataSetFromTximport(txi.genes,colData,design)
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

design <- ~Group
design(dds) <- design

dds <- DESeq(dds,parallel=T)

library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
trace="none",  # turns off trace lines inside the heat map
col=colours, # use on color palette defined earlier
margins=c(12,12), # widens margins around plot
srtCol=45,
srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
dev.off()

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```

#Analysis of gene expression


#set the significance level for BH adjustment	    
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 412 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 382, 93%
LFC < 0 (down)   : 30, 7.3%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 5337 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 4058, 76%
LFC < 0 (down)   : 1279, 24%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 18903 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 9706, 51%
LFC < 0 (down)   : 9197, 49%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 18677 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8620, 46%
LFC < 0 (down)   : 10057, 54%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 15916 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 7590, 48%
LFC < 0 (down)   : 8326, 52%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 11146 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6156, 55%
LFC < 0 (down)   : 4990, 45%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)
tar xvjf filename.tar.bz2

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library("naturalsort",lib.loc="/home/adamst/R/x86_64-pc-linux-gnu-library/3.2/")
library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("/data/scratch/gomeza/apple_genome/transcripts_modified.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)
```





Method 2: Host response gene expression.

```R

setwd("/data/scratch/gomeza/")

#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/RNA-seq/scripts/myfunctions")
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

library(tximport)
library(rjson)
library(readr)

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files	    
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	   

# get the sample names from the folders	    
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))


#==========================================================================================
#       Read sample metadata and annotations
#=========================================================================================

# Read sample metadata

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

colData <- colData[!(colData$Sample=="M9_C1_3"),]      
colData <- colData[!(colData$Sample=="M9_C2_3"),]      
colData <- colData[!(colData$Sample=="M9_C3_3"),]      
colData <- colData[!(colData$Sample=="GD_C1_3"),]      
colData <- colData[!(colData$Sample=="GD_C2_3"),]      
colData <- colData[!(colData$Sample=="GD_C3_3"),]      

dds <- DESeqDataSetFromTximport(txi.genes,colData,design)
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

design <- ~Group
design(dds) <- design

dds <- DESeq(dds,parallel=T)

library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
trace="none",  # turns off trace lines inside the heat map
col=colours, # use on color palette defined earlier
margins=c(12,12), # widens margins around plot
srtCol=45,
srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
dev.off()

pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```
