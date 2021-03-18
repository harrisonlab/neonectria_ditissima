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
    for RNADir in $(ls -d qc_rna/48DD_experiment2020/WT53/*); do
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
    for Transcriptome in $(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all_parsed.fa); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
        Organism=V.dahliae
        echo "$Organism - $Strain"
        for RNADir in $(ls -d qc_rna/48DD_experiment2020/WT53/*/cleaned/*/*); do
            FileNum=$(ls $RNADir/F/*_1_cleaned.fq.gz | wc -l)
            for num in $(seq 1 $FileNum); do
                printf "\n"
                FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
                FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
                echo $FileF
                echo $FileR
                Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
                Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
                Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
                echo "$Prefix"
                echo "$Timepoint"
                echo "$Sample_Name"
                OutDir=alignment/salmon/48DD_experiment/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
                ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
                sbatch -p himem $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
            done
        done
    done
```

Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/48DD_experiment/V.dahliae/JR2/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/48DD_experiment/V.dahliae/JR2/*/*/*/quant.sf); do
  Prefix=$(echo $File | cut -f8 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/$Prefix
  cp $PWD/$File RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/vertclock2020")

#===============================================================================
#       Load libraries
#===============================================================================

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

library(MetaCycle)
library(docker4seq)



#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files
txi.reps <- tximport(paste(list.dirs("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders
mysamples <- list.dirs("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

write.table(txi.genes,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/txireps.txt",sep="\t",na="",quote=F)




pca(experiment.table="txigenes2.txt", type="TPM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="txigenes_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

#===============================================================================
#       Read sample metadata and annotations
#===============================================================================

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/48DD_experiment/V.dahliae/JR2/DeSeq2/Vd_clocks_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Define the DESeq 'GLM' model
design <- ~ Condition
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
###
[1] "Intercept"                    "Condition_WT53_16_vs_WT53_12"
 [3] "Condition_WT53_20_vs_WT53_12" "Condition_WT53_24_vs_WT53_12"
 [5] "Condition_WT53_28_vs_WT53_12" "Condition_WT53_32_vs_WT53_12"
 [7] "Condition_WT53_36_vs_WT53_12" "Condition_WT53_4_vs_WT53_12" 
 [9] "Condition_WT53_40_vs_WT53_12" "Condition_WT53_44_vs_WT53_12"
[11] "Condition_WT53_48_vs_WT53_12" "Condition_WT53_8_vs_WT53_12" 
###

#Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


write.csv(vst, file="vst.csv")
