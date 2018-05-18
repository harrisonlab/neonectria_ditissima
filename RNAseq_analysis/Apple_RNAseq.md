#Method 2: Host response gene expression.

Cultivar folder without duplicates genes

#Create a gene sequence fastafile

Cultivar 2

```bash

cd /data/scratch/gomeza/apple_genome

gffread gene_models_20170612.gff -g GDDH13_1-1_formatted.fasta -w your_transcripts.fasta
cat your_transcripts.fasta | sed 's/gene=.*//g' > transcripts_modified.fasta

gffread gene_models_20170612.gff -g GDDH13_1-1_formatted.fasta -w your_transcripts.fasta
cat transcripts_modified2.fasta | sed 's/ncRNA://g' > transcripts_modified3.fasta
```

#Salmon tool for quantifying the expression of transcripts using RNA-seq data

Note that it is designed to align to predicted transcripts and not to the whole genome

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
mkdir -p alignment/salmon/N.ditissima/Cultivar/DeSeq2
for File in $(ls alignment/salmon/*/Cultivar/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene.txt
done

for File in $(ls alignment/salmon/*/Cultivar/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene2.txt
done

# Put files in a convenient location for DeSeq. Analysis was not performed on
# mycelium control samples.
  for File in $(ls alignment/salmon/N.ditissima/Cultivar/*/*/*/quant.sf | grep -v -e 'Hg199_1' -e 'Hg199_2' -e 'Hg199_3'); do
  Prefix=$(echo $File | cut -f7 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/N.ditissima/Cultivar/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/N.ditissima/Cultivar/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
  done
```

##A file was created with columns referring to experimental treatments

```bash
OutDir=alignment/salmon/N.ditissima/Cultivar2/DeSeq2
mkdir -p $OutDir
printf "Sample.name\tCultivar\tTimepoint\n" > $OutDir/N.dit_Hg199_RNAseq_design.txt
for i in $(seq 1 24)
do
if [ $i == '10' ] || [ $i == '11' ] || [ $i == '12' ]
then
Timepoint='Control'
Cultivar='mycelium'
elif [ $i == '7' ] || [ $i == '8' ] || [ $i == '9' ]
then
Timepoint='t0'
Cultivar='GD'
elif [ $i == '19' ] || [ $i == '20' ] || [ $i == '21' ]
then
Timepoint='t0'
Cultivar='M9'
elif [ $i == '1' ] || [ $i == '4' ] || [ $i == '6' ]
then
Timepoint='t1'
Cultivar='GD'
elif [ $i == '14' ] || [ $i == '16' ] || [ $i == '18' ]
then
Timepoint='t1'
Cultivar='M9'
elif [ $i == '2' ] || [ $i == '3' ] || [ $i == '5' ]
then
Timepoint='t2'
Cultivar='GD'
elif [ $i == '13' ] || [ $i == '15' ] || [ $i == '17' ]
then
Timepoint='t2'
Cultivar='M9'   
fi
#if [ $i == '10' ]
#then
#    printf "Hg199_1\t$Cultivar\t$Timepoint\n"
#elif [ $i == '11' ]
#then
#    printf "Hg199_2\t$Cultivar\t$Timepoint\n"
#elif [ $i == '12' ]
#then
#    printf "Hg199_3\t$Cultivar\t$Timepoint\n"
if [ $i == '7' ]
then
printf "GD_C1_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '8' ]
then
printf "GD_C2_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '9' ]
then
printf "GD_C3_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '19' ]
then
printf "M9_C1_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '20' ]
then
printf "M9_C2_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '21' ]
then
printf "M9_C3_3\t$Cultivar\t$Timepoint\n"
elif [ $i == '1' ]
then
printf "GD_4A4\t$Cultivar\t$Timepoint\n"
elif [ $i == '4' ]
then
printf "GD_5A3\t$Cultivar\t$Timepoint\n"
elif [ $i == '6' ]
then
printf "GD_6A3\t$Cultivar\t$Timepoint\n"
elif [ $i == '14' ]
then
printf "M9_2A3\t$Cultivar\t$Timepoint\n"
elif [ $i == '16' ]
then
printf "M9_5A4\t$Cultivar\t$Timepoint\n"
elif [ $i == '18' ]
then
printf "M9_6A3\t$Cultivar\t$Timepoint\n"
elif [ $i == '2' ]
then
printf "GD_4A5\t$Cultivar\t$Timepoint\n"
elif [ $i == '3' ]
then
printf "GD_5A2\t$Cultivar\t$Timepoint\n"
elif [ $i == '5' ]
then
printf "GD_6A2\t$Cultivar\t$Timepoint\n"
elif [ $i == '13' ]
then
printf "M9_2A2\t$Cultivar\t$Timepoint\n"
elif [ $i == '15' ]
then
printf "M9_5A2\t$Cultivar\t$Timepoint\n"
elif [ $i == '17' ]
then
printf "M9_6A2\t$Cultivar\t$Timepoint\n"
fi
done >> $OutDir/N.dit_Hg199_RNAseq_design.txt
```

```bash
export R_LIBS=/home/armita/R/x86_64-pc-linux-gnu-library/3.4
/home/deakig/R3.4/bin/R
```

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
library(tximport)
library(rjson)
library(readr)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/trans2gene2.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/N.ditissima/Cultivar2/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/N.ditissima/Cultivar2/DeSeq2",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata

unorderedColData <- read.table("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

design <- ~Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

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

pdf("alignment/salmon/N.ditissima/Cultivar2/heatmap_vst.pdf", width=12,height=12)
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

pdf("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/heatmap_rld.pdf")
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
pdf("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Sample.name", "Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Sample.name", "Group"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()

ggsave("alignment/salmon/N.ditissima/Cultivar2/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t1_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t1_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t1_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 16237 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 7710, 47%
LFC < 0 (down)   : 8527, 53%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t2_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t2_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/GD_t2_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 13141 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6510, 50%
LFC < 0 (down)   : 6631, 50%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 11957 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6764, 57%
LFC < 0 (down)   : 5193, 43%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 16539 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8723, 53%
LFC < 0 (down)   : 7816, 47%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_GD_t1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_GD_t1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t1_vs_GD_t1_down.txt",sep="\t",na="",quote=F)

out of 6181 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 3192, 52%
LFC < 0 (down)   : 2989, 48%
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

write.table(sig.res,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/N.ditissima/Cultivar2/DeSeq2/M9_t2_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 6660 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 3471, 52%
LFC < 0 (down)   : 3189, 48%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

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

#Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

 #robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)
```

##Produce a more detailed table of analyses

```bash
for GeneGff in $(ls /data/scratch/gomeza/apple_genome/gene_models_20170612.gff3); do
  Assembly=$(ls /data/scratch/gomeza/apple_genome/GDDH13_1-1_formatted.fasta)
  InterPro=$(ls /data/scratch/gomeza/apple_genome/Interproscan/interpro_1v1.tsv)
  OutDir=gene_pred/annotation/salmon/GDDH13_V_1.1/Apple_vs_Nd
  mkdir -p $OutDir
  GeneFasta=$(ls /data/scratch/gomeza/apple_genome/GDDH13_1-1_prot.fasta)
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
  Dir1=$(ls -d /data/scratch/gomeza/alignment/salmon/N.ditissima/Cultivar/DeSeq2_v2)
  DEG_Files=$(ls \
  $Dir1/GD_t1_vs_GD_t0.txt \
  $Dir1/M9_t1_vs_M9_t0.txt \
  $Dir1/GD_t2_vs_GD_t0.txt \
  $Dir1/M9_t2_vs_M9_t0.txt \
  $Dir1/GD_t1_vs_GD_t2.txt \
  $Dir1/GD_t1_vs_M9_t1.txt \
  $Dir1/GD_t2_vs_M9_t2.txt \
  $Dir1/M9_t1_vs_M9_t2.txt \
  $Dir1/GD_t0_vs_M9_t0.txt \
  | sed -e "s/$/ /g" | tr -d "\n")
  RawCount=$(ls alignment/salmon/N.ditissima/Cultivar/DeSeq2_v2/raw_counts.txt)
  FPKM=$(ls alignment/salmon/N.ditissima/Cultivar/DeSeq2_v2/fpkm_counts.txt)
  $ProgDir/Nd_annotation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP4 $SigP4 --trans_mem $TMHMM_headers --TFs $TFs --Antismash $Antismash --effector_total $effector_total --CAZY_total $CAZY_total --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp2.tsv
  done
done
```







#Method 2: Host response gene expression.

Cultivar

```bash
export R_LIBS=/home/armita/R/x86_64-pc-linux-gnu-library/3.4
/home/deakig/R3.4/bin/R
```

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
library(tximport)
library(rjson)
library(readr)

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/N.ditissima/Cultivar/DeSeq2/trans2gene2.txt",header=T,sep="\t")

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
pdf("alignment/salmon/N.ditissima/Cultivar/PCA_vst.pdf")
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

#Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

#robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'

fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/N.ditissima/Cultivar/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)
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

##Produce a more detailed table of analyses

```bash
for GeneGff in $(ls /data/scratch/gomeza/apple_genome/gene_models_20170612.gff3); do
  Assembly=$(ls /data/scratch/gomeza/apple_genome/GDDH13_1-1_formatted.fasta)
  InterPro=$(ls /data/scratch/gomeza/apple_genome/Interproscan/interpro_1v1.tsv)
  OutDir=gene_pred/annotation/salmon/GDDH13_V_1.1/Apple_vs_Nd
  mkdir -p $OutDir
  GeneFasta=$(ls /data/scratch/gomeza/apple_genome/GDDH13_1-1_prot.fasta)
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
  Dir1=$(ls -d /data/scratch/gomeza/alignment/salmon/N.ditissima/Cultivar/DeSeq2)
  DEG_Files=$(ls \
  $Dir1/GD_t1_vs_GD_t0.txt \
  $Dir1/M9_t1_vs_M9_t0.txt \
  $Dir1/GD_t2_vs_GD_t0.txt \
  $Dir1/M9_t2_vs_M9_t0.txt \
  $Dir1/GD_t1_vs_GD_t2.txt \
  $Dir1/M9_t1_vs_GD_t1.txt \
  $Dir1/M9_t2_vs_GD_t2.txt \
  $Dir1/M9_t1_vs_M9_t2.txt \
  | sed -e "s/$/ /g" | tr -d "\n")
  RawCount=$(ls alignment/salmon/N.ditissima/Cultivar/DeSeq2/raw_counts.txt)
  FPKM=$(ls alignment/salmon/N.ditissima/Cultivar/DeSeq2/fpkm_counts.txt)
  $ProgDir/Apple_annotation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --InterPro $InterPro > $OutDir/Apple_vs_Nd_gene_table_incl_exp2.tsv
done
```
