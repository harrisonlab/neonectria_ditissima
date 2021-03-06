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
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
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
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
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
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done

for File in $(ls alignment/salmon/*/Hg199_minion/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/trans2gene2.txt
done
```

```bash
mkdir -p alignment/salmon/N.ditissima/Cultivar/DeSeq2
for File in $(ls alignment/salmon/*/Cultivar/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene.txt  
done
# Put files in a convenient location for DeSeq. Analysis was not performed on
# mycelium control samples.
for File in $(ls alignment/salmon/N.ditissima/Cultivar/*/*/*/quant.sf | grep -v -e 'Hg199_1' -e 'Hg199_2' -e 'Hg199_3'); do
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/*/Cultivar/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Cultivar/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

```bash
export R_LIBS=/home/armita/R/x86_64-pc-linux-gnu-library/3.4
/home/deakig/R3.4/bin/R
```

#Method 1: Gene expression of Nd.

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


#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata	    
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

Library normalisation

# create dds object from Salmon counts and sample metadata (library size normalisation is taken from the length estimates)

design <- ~Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

dds <- estimateSizeFactors(dds)
Group <- as.factor(dds$Group)

dds <- DESeq(dds,parallel=T)

disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)


Sample Distances

#install and load libraries
require("pheatmap")
require(data.table)
library("RColorBrewer")
#install.packages("gplots")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
#install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/N.ditissima/Hg199_minion/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
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

pdf("alignment/salmon/N.ditissima/Hg199_minion/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Group"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
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

```R
#set the significance level for BH adjustment	    
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 5051 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1961, 39%
LFC < 0 (down)   : 3090, 61%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 2774 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1151, 41%
LFC < 0 (down)   : 1623, 59%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 2128 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1385, 65%
LFC < 0 (down)   : 743, 35%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 5555 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1621, 29%
LFC < 0 (down)   : 3934, 71%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","M9_t1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_M9_t1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_M9_t1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_M9_t1_down.txt",sep="\t",na="",quote=F)

out of 81 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 73, 90%
LFC < 0 (down)   : 8, 9.9%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 10)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 26 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8, 31%
LFC < 0 (down)   : 18, 69%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 35)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 102 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 27, 26%
LFC < 0 (down)   : 75, 74%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 32)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 150 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 5, 3.3%
LFC < 0 (down)   : 145, 97%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

#Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)


 #robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)
```

#Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/*_up.txt); do
  DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
  DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
  # echo $UpFile
  # echo $DownFile
  cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
  echo $DegFile
  cat $DegFile | wc -l
done
```

```
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_GD_t2_DEGs.txt
101
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_M9_t1_DEGs.txt
63
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t1_vs_mycelium_Control_DEGs.txt
4536
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_M9_t2_DEGs.txt
26
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/GD_t2_vs_mycelium_Control_DEGs.txt
2767
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_M9_t2_DEGs.txt
150
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t1_vs_mycelium_Control_DEGs.txt
5029
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/M9_t2_vs_mycelium_Control_DEGs.txt
2128
```

#Inital analysis of tables of DEGs

```bash
  effector_names=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
  CAZY_names=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
  for File in $(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/*_vs_*.txt |)
  do
      Assessment=$(basename $File | sed "s/.txt//g")
      echo $Assessment
      echo "Total number of genes in dataset:"
      cat $File | grep -v 'baseMean' | wc -l
      echo "Total number of effectors in dataset:"
      Effector_File=$(echo $File | sed "s/.txt/_Effector.txt/g")
      cat $File | head -n 1 > $Effector_File
      cat $File | grep -w -f $effector_names >> $Effector_File
      cat $Effector_File | tail -n +2 | wc -l
      echo "Total number of CAZY in dataset:"
      CAZY_File=$(echo $File | sed "s/.txt/_CAZY.txt/g")
      cat $File | head -n 1 > $CAZY_File
      cat $File | grep -w -f $CAZY_names >> $CAZY_File
      cat $CAZY_File | tail -n +2 | wc -l
    done
done
```

#Method 1.2: Gene expression of Nd.

This analysis was done repeating the salmon alignment with the option --keepduplicates.

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
tx2gene <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/trans2gene.txt",header=T,sep="\t")

# import quantification files	    
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	   

# get the sample names from the folders	    
mysamples <- list.dirs("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2",full.names=F,recursive=F)

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

unorderedColData <- read.table("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

Library normalisation

# create dds object from Salmon counts and sample metadata (library size normalisation is taken from the length estimates)

design <- ~Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

dds <- estimateSizeFactors(dds)
Group <- as.factor(dds$Group)

dds <- DESeq(dds,parallel=T)

disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)


Sample Distances

#install and load libraries
require("pheatmap")
require(data.table)
library("RColorBrewer")
#install.packages("gplots")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
#install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/N.ditissima/Hg199_minion/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
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

pdf("alignment/salmon/N.ditissima/Hg199_minion/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Group"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
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

```R
#set the significance level for BH adjustment	    
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 5008 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2024, 40%
LFC < 0 (down)   : 2984, 60%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 2730 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1174, 43%
LFC < 0 (down)   : 1556, 57%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 2198 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1455, 66%
LFC < 0 (down)   : 743, 34%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 5472 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1659, 30%
LFC < 0 (down)   : 3813, 70%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","M9_t1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_M9_t1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_M9_t1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_M9_t1_down.txt",sep="\t",na="",quote=F)

out of 71 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 64, 90%
LFC < 0 (down)   : 7, 9.9%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 13)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 26 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8, 31%
LFC < 0 (down)   : 18, 69%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 34)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 99 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 24, 24%
LFC < 0 (down)   : 75, 76%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 31)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 146 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 5, 3.4%
LFC < 0 (down)   : 141, 97%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
```

#Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/normalised_counts.txt",sep="\t",na="",quote=F)


 #robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/fpkm_counts.txt",sep="\t",na="",quote=F)
```

#Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/*_up.txt); do
  DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
  DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
  # echo $UpFile
  # echo $DownFile
  cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
  echo $DegFile
  cat $DegFile | wc -l
done
```

```
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_GD_t2_DEGs.txt
99
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_M9_t1_DEGs.txt
55
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_DEGs.txt
4553
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_M9_t2_DEGs.txt
26
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_DEGs.txt
2721
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_M9_t2_DEGs.txt
146
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_DEGs.txt
4998
alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_DEGs.txt
2198
```

##Produce a more detailed table of analyses

```bash
for GeneGff in $(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3); do
  Strain=Hg199_minion
  Organism=N.ditissima
  Assembly=$(ls /data/scratch/gomeza/Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/*_contigs_unmasked.fa)
  InterPro=$(ls /data/scratch/gomeza/gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv)
  SwissProt=$(ls /data/scratch/gomeza/gene_pred/swissprot/N.ditissima/Hg199_minion/swissprot_vJul2016_tophit_parsed.tbl)
  OutDir=gene_pred/annotation/salmon2/$Organism/$Strain
  mkdir -p $OutDir
  GeneFasta=$(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.pep.fasta)
  TFs=$(ls analysis/transcription_factors/N.ditissima/Hg199_minion/Hg199_minion_TF_domains.tsv)
  Antismash=$(ls /data/scratch/gomeza/analysis/secondary_metabolites/antismash/Hg199_minion/fungi-dfc734cf-18aa-414d-b034-0da05c627613/Hg199_minion_antismash_secmet_genes.tsv)
  SigP4=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/Hg199_minion_final_sp_no_trans_mem.aa)
  effector_total=$(ls analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt)
  CAZY_total=$(ls gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt)
  TMHMM_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_pos_headers.txt)
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
  Dir1=$(ls -d /data/scratch/gomeza/alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2)
  DEG_Files=$(ls \
  $Dir1/GD_t1_vs_mycelium_Control.txt \
  $Dir1/M9_t1_vs_mycelium_Control.txt \
  $Dir1/GD_t2_vs_mycelium_Control.txt \
  $Dir1/M9_t2_vs_mycelium_Control.txt \
  $Dir1/GD_t1_vs_GD_t2.txt \
  $Dir1/GD_t1_vs_M9_t1.txt \
  $Dir1/GD_t2_vs_M9_t2.txt \
  $Dir1/M9_t1_vs_M9_t2.txt \
  | sed -e "s/$/ /g" | tr -d "\n")
  RawCount=$(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/raw_counts.txt)
  FPKM=$(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/fpkm_counts.txt)
  $ProgDir/Nd_annotation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP4 $SigP4 --trans_mem $TMHMM_headers --TFs $TFs --Antismash $Antismash --effector_total $effector_total --CAZY_total $CAZY_total --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp3.tsv
  done
done
```

##All genes

###All DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Upregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_up.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_up.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_up.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_up.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_up.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_up.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_up.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_up.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_down.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_down.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_down.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_down.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t1_vs_mycelium_Control_down.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t1_vs_mycelium_Control_down.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_t2_vs_mycelium_Control_down.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_t2_vs_mycelium_Control_down.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

#Extract fasta file of all DEGs for BLAST analysis

###All DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_all_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_all_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_all_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_all_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_all_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_all_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_all_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_all_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_all_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_all_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_all_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_all_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

###Upregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

###Downregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_DEGs.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_DEGs_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_DEGs.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

#Investigate enriched functional annotations in DEGs vs all genes

```bash
OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv
```

##GD up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/GD/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_DEGs_names.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_up_DEGs.txt
  Set2Genes=$OutDir/GD_up_genes2.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/GD/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_DEGs_names.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_down_DEGs.txt
  Set2Genes=$OutDir/GD_down_genes2.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##M9 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/M9/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_DEGs_names.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_up_DEGs.txt
  Set2Genes=$OutDir/M9_up_genes2.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/M9/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_DEGs_names.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_down_DEGs.txt
  Set2Genes=$OutDir/M9_down_genes2.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##t1 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/t1/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_DEGs_names.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_up_DEGs.txt
  Set2Genes=$OutDir/t1_up_genes2.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/t1/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_DEGs_names.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_down_DEGs.txt
  Set2Genes=$OutDir/t1_down_genes2.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##t2 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/t2/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_DEGs_names.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_up_DEGs.txt
  Set2Genes=$OutDir/t2_up_genes2.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/t2/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/CultivarExperiment/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_DEGs_names.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_down_DEGs.txt
  Set2Genes=$OutDir/t2_down_genes2.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

#Inital analysis of tables of DEGs

```bash
cat analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt | sed 's/.t.*//g' > analysis/effectorP/N.ditissima/Hg199_minion/EffectorP_gene_headers.txt
cat gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt | sed 's/.t.*//g' > gene_pred/CAZY/N.ditissima/Hg199_minion/CAZY_gene_headers.txt

effector_names=analysis/effectorP/N.ditissima/Hg199_minion/EffectorP_gene_headers.txt
CAZY_names=gene_pred/CAZY/N.ditissima/Hg199_minion/CAZY_gene_headers.txt
TF_names=analysis/transcription_factors/N.ditissima/Hg199_minion/Hg199_minion_TF_geneid_headers.txt

for File in $(ls alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/*vs*.txt | grep -v "genes" | grep -v "countData")
do
Assessment=$(basename $File | sed "s/.txt//g")
echo $Assessment
echo "Total number of genes in dataset:"
cat $File | grep -v 'baseMean' | wc -l

echo "Total number of effectors in dataset:"
Effector_File=$(echo $File | sed "s/.txt/_Effector.txt/g")
cat $File | head -n 1 > $Effector_File
cat $File | grep -w -f $effector_names >> $Effector_File
cat $Effector_File | tail -n +2 | wc -l

echo "Total number of CAZY in dataset:"
CAZY_File=$(echo $File | sed "s/.txt/_CAZY.txt/g")
cat $File | head -n 1 > $CAZY_File
cat $File | grep -w -f $CAZY_names >> $CAZY_File
cat $CAZY_File | tail -n +2 | wc -l

echo "Total number of TF in dataset:"
TF_File=$(echo $File | sed "s/.txt/_TF.txt/g")
cat $File | head -n 1 > $TF_File
cat $File | grep -w -f $TF_names >> $TF_File
cat $TF_File | tail -n +2 | wc -l
done
```

##Effectors

###Upregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t1_vs_mycelium_Control_up_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t2_vs_mycelium_Control_up_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t1_vs_mycelium_Control_up_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t2_vs_mycelium_Control_up_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t1_vs_mycelium_Control_up_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t1_vs_mycelium_Control_up_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t2_vs_mycelium_Control_up_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t2_vs_mycelium_Control_up_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t1_vs_mycelium_Control_down_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t2_vs_mycelium_Control_down_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_down_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t1_vs_mycelium_Control_down_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t2_vs_mycelium_Control_down_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_down_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t1_vs_mycelium_Control_down_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t1_vs_mycelium_Control_down_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_down_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_t2_vs_mycelium_Control_down_Effector.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_t2_vs_mycelium_Control_down_Effector.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_down_Effectors.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

##CAZY

###Upregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t1_vs_mycelium_Control_up_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t2_vs_mycelium_Control_up_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t1_vs_mycelium_Control_up_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t2_vs_mycelium_Control_up_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t1_vs_mycelium_Control_up_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t1_vs_mycelium_Control_up_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t2_vs_mycelium_Control_up_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t2_vs_mycelium_Control_up_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t1_vs_mycelium_Control_down_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t2_vs_mycelium_Control_down_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_down_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t1_vs_mycelium_Control_down_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t2_vs_mycelium_Control_down_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_down_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t1_vs_mycelium_Control_down_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t1_vs_mycelium_Control_down_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_down_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_t2_vs_mycelium_Control_down_CAZY.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_t2_vs_mycelium_Control_down_CAZY.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_down_CAZY.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

##TF

###Upregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t1_vs_mycelium_Control_up_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t2_vs_mycelium_Control_up_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/GD_up_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t1_vs_mycelium_Control_up_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t2_vs_mycelium_Control_up_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/M9_up_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t1_vs_mycelium_Control_up_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t1_vs_mycelium_Control_up_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t1_up_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t2_vs_mycelium_Control_up_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t2_vs_mycelium_Control_up_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/t2_up_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t1_vs_mycelium_Control_down_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t2_vs_mycelium_Control_down_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_down_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t1_vs_mycelium_Control_down_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t2_vs_mycelium_Control_down_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_down_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t1_vs_mycelium_Control_down_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t1_vs_mycelium_Control_down_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_down_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_t2_vs_mycelium_Control_down_TF.txt
inp2=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_t2_vs_mycelium_Control_down_TF.txt
OutDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_down_TF.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

#Extract fasta file of all DEGs for BLAST analysis

##Effectors

###Upregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_up_Effectors.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_up_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_up_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_up_Effectors.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_up_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_up_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_up_Effectors.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_up_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_up_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_up_Effectors.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_up_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_up_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

###Downregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_down_Effector.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_down_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_down_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_down_Effector.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_down_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_down_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_down_Effector.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_down_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_down_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_down_Effector.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_down_Effector_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_down_Effector.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

##CAZY

###Upregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_up_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_up_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_up_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_up_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_up_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_up_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_up_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_up_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_up_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_up_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_up_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_up_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

###Downregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_down_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_down_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_down_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_down_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_down_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_down_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_down_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_down_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_down_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_down_CAZY.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_down_CAZY_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_down_CAZY.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

##TF

###Upregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_up_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_up_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_up_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_up_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_up_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_up_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_up_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_up_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_up_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_up_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_up_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_up_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```

###Downregulated DEGs

```bash
DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_down_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_down_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_down_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_down_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_down_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_down_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_down_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_down_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t1_down_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta

DEGFile=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_down_TF.tsv
DEGNames=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_down_TF_names.txt
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.gene.fasta
DEGFasta=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/t2_down_TF.fa
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
$ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
```


#Investigate enriched functional annotations in DEGs vs all genes

```bash
OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv
```

##GD up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Upregulated_Effectors
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_up_Effector_names.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_up_Effector.txt
  Set2Genes=$OutDir/GD_up_genes2.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Downregulated_Effectors
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/GD_down_Effector_names.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_down_Effector.txt
  Set2Genes=$OutDir/GD_down_genes2.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Upregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_up_CAZY_names.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_up_CAZY.txt
  Set2Genes=$OutDir/GD_up_genes2.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Downregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/GD_down_CAZY_names.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_down_CAZY.txt
  Set2Genes=$OutDir/GD_down_genes2.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Upregulated_TF
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_up_TF_names.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_up_TF.txt
  Set2Genes=$OutDir/GD_up_genes2.txt
  AllGenes=$OutDir/GD_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/GD/Downregulated_TF
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/GD_down_TF_names.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/GD_down_TF.txt
  Set2Genes=$OutDir/GD_down_genes2.txt
  AllGenes=$OutDir/GD_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##M9 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Upregulated_Effectors
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_up_Effector_names.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_up_Effector.txt
  Set2Genes=$OutDir/M9_up_genes2.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Downregulated_Effectors
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/M9_down_Effector_names.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_down_Effector.txt
  Set2Genes=$OutDir/M9_down_genes2.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Upregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_up_CAZY_names.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_up_CAZY.txt
  Set2Genes=$OutDir/M9_up_genes2.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Downregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/M9_down_CAZY_names.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_down_CAZY.txt
  Set2Genes=$OutDir/M9_down_genes2.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Upregulated_TF
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_up_TF_names.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_up_TF.txt
  Set2Genes=$OutDir/M9_up_genes2.txt
  AllGenes=$OutDir/M9_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/M9/Downregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/TF/M9_down_TF_names.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/M9_down_TF.txt
  Set2Genes=$OutDir/M9_down_genes2.txt
  AllGenes=$OutDir/M9_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##t1 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t1/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_up_Effector_names.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_up_Effector.txt
  Set2Genes=$OutDir/t1_up_genes2.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t1/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t1_down_Effector_names.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_down_Effector.txt
  Set2Genes=$OutDir/t1_down_genes2.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t1/Upregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_up_CAZY_names.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_up_CAZY.txt
  Set2Genes=$OutDir/t1_up_genes2.txt
  AllGenes=$OutDir/t1_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t1/Downregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t1_down_CAZY_names.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t1_down_CAZY.txt
  Set2Genes=$OutDir/t1_down_genes2.txt
  AllGenes=$OutDir/t1_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

##t2 up and down

```bash
  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t2/Upregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_up_Effector_names.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_up_Effector.txt
  Set2Genes=$OutDir/t2_up_genes2.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t2/Downregulated
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/Effectors/t2_down_Effector_names.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_down_Effector.txt
  Set2Genes=$OutDir/t2_down_genes2.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t2/Upregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_up_CAZY_names.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_up_CAZY.txt
  Set2Genes=$OutDir/t2_up_genes2.txt
  AllGenes=$OutDir/t2_up_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt

  OutDir=analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/t2/Downregulated_CAZY
  mkdir -p $OutDir
  cp analysis/enrichment/N.ditissima/Hg199_minion/EffectorAnalysis/experiment_all_gene_GO_annots_geneid.tsv $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Gene_enrichment/
  AnnotTable=gene_pred/annotation/salmon2/N.ditissima/Hg199_minion/Hg199_minion_gene_table_incl_exp.tsv
  DEGs=alignment/salmon/N.ditissima/Hg199_minion/DeSeq2_v2/CAZY/t2_down_CAZY_names.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
  Set1Genes=$OutDir/t2_down_CAZY.txt
  Set2Genes=$OutDir/t2_down_genes2.txt
  AllGenes=$OutDir/t2_down_genes.txt
  cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/.t.*//g' | sed -e 's/$/\t1.00/g' > $Set2Genes
  cat $Set1Genes $Set2Genes > $AllGenes
  $ProgDir/GO_enrichment_allGO.r --all_genes $AllGenes --GO_annotations $OutDir/experiment_all_gene_GO_annots_geneid.tsv --out_dir $OutDir > $OutDir/output.txt
```

###Venn diagrams

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
WorkDir=alignment/salmon/N.ditissima/Hg199_minion/DeSeq_v2/Effectors
$ProgDir/All_DEGs_venn_diag_2.r --inp GD_up_Effector.tsv --out GD_up_Effector.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_all_DEGs.tsv --out $WorkDir/t1_all_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_up_DEGs.tsv --out $WorkDir/t1_up_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_down_DEGs.tsv --out $WorkDir/t1_down_DEGs.pdf
```




#Method 2: Host response gene expression.

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
load_all("/home/deakig/pipelines/RNA-seq/scripts/myfunctions")
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
tx2gene <- read.table("alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene.txt",header=T,sep=" ")

# import quantification files	    
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Cultivar/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	   

# get the sample names from the folders	    
mysamples <- list.dirs("alignment/salmon/N.ditissima/Cultivar/DeSeq2",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata

unorderedColData <- read.table("alignment/salmon/N.ditissima/Cultivar/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

design <- ~Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

dds <- estimateSizeFactors(dds)
Group <- as.factor(dds$Group)

dds <- DESeq(dds,parallel=T)

disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
```

#Sample Distances

 #install and load libraries
require("pheatmap")
require(data.table)
library("RColorBrewer")
 #install.packages("gplots")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/N.ditissima/Cultivar/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
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

pdf("alignment/salmon/N.ditissima/Cultivar/DeSeq2/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/salmon/N.ditissima/Cultivar/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Group"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Cultivar/DeSeq2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Group"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()

ggsave("alignment/salmon/N.ditissima/Cultivar/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

```R

#set the significance level for BH adjustment	    
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 16171 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 7677, 47%
LFC < 0 (down)   : 8494, 53%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t2_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t2_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t2_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 13139 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6510, 50%
LFC < 0 (down)   : 6629, 50%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 11954 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6763, 57%
LFC < 0 (down)   : 5191, 43%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 16536 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8721, 53%
LFC < 0 (down)   : 7815, 47%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","GD_t1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_GD_t1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_GD_t1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_GD_t1_down.txt",sep="\t",na="",quote=F)

out of 6182 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 3194, 52%
LFC < 0 (down)   : 2988, 48%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 6661 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 3472, 52%
LFC < 0 (down)   : 3189, 48%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 541 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 469, 87%
LFC < 0 (down)   : 72, 13%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 116 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 111, 96%
LFC < 0 (down)   : 5, 4.3%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

#Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/N.ditissima/Cultivar/DeSeq2/*_up.txt); do
  DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
  DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
  # echo $UpFile
  # echo $DownFile
  cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
  echo $DegFile
  cat $DegFile | wc -l
done
```

```
alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t0_DEGs.txt
11724
alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t1_vs_GD_t2_DEGs.txt
107
alignment/salmon/N.ditissima/Cultivar/DeSeq2/GD_t2_vs_GD_t0_DEGs.txt
9265
alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_GD_t1_DEGs.txt
4766
alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t0_DEGs.txt
11923
alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t1_vs_M9_t2_DEGs.txt
450
alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_GD_t2_DEGs.txt
5113
alignment/salmon/N.ditissima/Cultivar/DeSeq2/M9_t2_vs_M9_t0_DEGs.txt
8043
```
