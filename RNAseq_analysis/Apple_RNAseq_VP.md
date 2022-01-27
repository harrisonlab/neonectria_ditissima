# RNA-Seq analysis

Analysis of the apple reads

## Perform qc on RNA-Seq data

```bash
# Run fastq-mcf
for RNADir in $(ls -d rna_seq/N.ditissima/M9/*); do
FileNum=$(ls $RNADir/F/*_1.fq.gz | wc -l)
    for num in $(seq 1 $FileNum); do
        printf "\n"
        FileF=$(ls $RNADir/F/*.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        IluminaAdapters=../../../home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p himem $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA
    done
done
```

## Decontamination of rRNA reads in RNAseq data

```bash
    for RNADir in $(ls -d qc_rna/N.ditissima/*/*); do
    FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            echo $RNADir
            Strain=$(sed 's/.*\///' <<< $RNADir)
            Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            echo $Sample_Name
            sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned/$Timepoint/$Sample_Name $FileF $FileR $ProgDir $Strain
        done
    done
```

## Salmon 


```bash

cd /projects/neonectria_ditissima/apple_genome/

gffread gene_models_20170612.gff -g GDDH13_1-1_formatted.fasta -w your_transcripts.fasta
cat your_transcripts.fasta | sed 's/gene=.*//g' > transcripts_modified.fasta

gffread gene_models_20170612.gff -g GDDH13_1-1_formatted.fasta -w your_transcripts.fasta
cat transcripts_modified2.fasta | sed 's/ncRNA://g' > transcripts_modified3.fasta

cp transcripts_modified3.fasta GDDH13_1-1_geneID.fasta
```

```bash
    # I have installed the latest version in a new env
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda create -n salmon salmon
```

```bash
    for Transcriptome in $(ls apple_genome/GDDH13_1-1_geneID.fasta); do
        for RNADir in $(ls -d ../../data/scratch/gomeza/qc_rna/N.ditissima/*/t*/cleaned/t*/*); do
        FileNum=$(ls $RNADir/F/*_1_cleaned.fq.gz | wc -l)
            for num in $(seq 1 $FileNum); do
                printf "\n"
                FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
                FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
                echo $FileF
                echo $FileR
                Prefix=$(echo $RNADir | rev | cut -f5 -d '/' | rev)
                Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
                Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
                echo "$Prefix"
                echo "$Timepoint"
                echo "$Sample_Name"
                OutDir=alignment/salmon_VP/Malus/$Prefix/$Timepoint/$Sample_Name
                ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
                sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
            done
        done
    done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p alignment/salmon/Malus/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/Malus/*/*/*/quant.sf | head -n1); do
cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/Malus/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.

for File in $(ls alignment/salmon/Malus/*/*/*/quant.sf); do
Prefix=$(echo $File | cut -f6 -d '/' --output-delimiter '_')
mkdir -p alignment/salmon/Malus/DeSeq2/$Prefix
cp $PWD/$File alignment/salmon/Malus/DeSeq2/$Prefix/quant.sf
# rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/data/scratch/gomeza/neonectria_ditissima")


# Load libraries

library(DESeq2)
library(BiocParallel)
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
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggrepel)

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Malus/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/Malus/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/Malus/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# write.table(txi.genes,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
# write.table(txi.reps,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txireps.txt",sep="\t",na="",quote=F)

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Malus/DeSeq2/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)

# Define the DESeq 'GLM' model
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
###
[1] "Intercept"            "Group_GD_t1_vs_GD_t0" "Group_GD_t2_vs_GD_t0" "Group_M9_t0_vs_GD_t0" "Group_M9_t1_vs_GD_t0"
[6] "Group_M9_t2_vs_GD_t0"
###

# Results

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/Malus/Contrasts/GD_t1_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Malus/Contrasts/GD_t1_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Malus/Contrasts/GD_t1_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/Malus/Contrasts/GD_t2_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Malus/Contrasts/GD_t2_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Malus/Contrasts/GD_t2_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/Malus/Contrasts/M9_t1_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Malus/Contrasts/M9_t1_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Malus/Contrasts/M9_t1_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"alignment/salmon/Malus/Contrasts/M9_t2_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Malus/Contrasts/M9_t2_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Malus/Contrasts/M9_t2_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

#===============================================================================
# Exploring and exporting results
#===============================================================================

# Sample Distances

vst1<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst1), file="alignment/salmon/Malus/DeSeq2/NDmalus_data_vst_T.csv")
vst2<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2), file="alignment/salmon/Malus/DeSeq2/NDmalus_data_vst_F.csv")

pdf("alignment/salmon/Malus/DeSeq2/heatmap_vstT.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Group)
colnames(sampleDistMatrix) <- paste(vst1$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Malus/DeSeq2/heatmap_vstF.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst2$Group)
colnames(sampleDistMatrix) <- paste(vst2$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

# PCA
data <- plotPCA(vst1, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst1))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/Malus/DeSeq2/PCA_vst_true.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst2))) + theme_light()
coord_fixed()
ggsave("alignment/salmon/Malus/DeSeq2/PCA_vst_false.pdf", pca_plot, dpi=300, height=10, width=12)

# Plot using rlog transformation:
rld <- rlog(dds)
pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Media"))
dev.off()
```


## Analysis of DeSeq2 output

```bash
for UpFile in $(ls alignment/salmon/Malus/Contrasts/*_up.txt); do
DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
echo $DegFile
cat $DegFile | wc -l
done
```

```
alignment/salmon/Malus/Contrasts/GD_t1_vs_GD_t0_DEGs.txt
11192
alignment/salmon/Malus/Contrasts/GD_t2_vs_GD_t0_DEGs.txt
8918
alignment/salmon/Malus/Contrasts/M9_t1_vs_M9_t0_DEGs.txt
11327
alignment/salmon/Malus/Contrasts/M9_t2_vs_M9_t0_DEGs.txt
7760
```

# Generating an TSV file with sequencing information

conda activate Python

```bash
for GeneGff in $(ls apple_genome/gene_models_20170612.gff3); do
Assembly=$(ls apple_genome/GDDH13_1-1_formatted.fasta)
InterPro=$(ls apple_genome/Interproscan/interpro_1v1.tsv)
OutDir=gene_pred/annotation/apple_genome/Apple_vs_Nd
mkdir -p $OutDir
GeneFasta=$(ls apple_genome/Genes/GDDH13_1-1_prot.fasta)
Dir1=$(ls -d alignment/salmon/Malus/Contrasts)
DEG_Files=$(ls \
$Dir1/GD_t1_vs_GD_t0.txt \
$Dir1/GD_t2_vs_GD_t0.txt \
$Dir1/M9_t1_vs_M9_t0.txt \
$Dir1/M9_t2_vs_M9_t0.txt \
| sed -e "s/$/ /g" | tr -d "\n")
ProgDir=/home/gomeza/git_repos/scripts/neonectria_ditissima/RNAseq_analysis
python $ProgDir/Apple_annotation_tables.py  --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --InterPro $InterPro > Malus_Nd_gene_table_incl_exp.tsv
done
```

### Add vst to tsv table 

```r
C1<- read.table("Malus_Nd_gene_table_incl_exp.txt",header=T,sep="\t")
C2<- read.table("NDmalus_data_vst_F.txt",header=T,sep="")
C3<- merge(C1,C2,by=gene.id)
```
