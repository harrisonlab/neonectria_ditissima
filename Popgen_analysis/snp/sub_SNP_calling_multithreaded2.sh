#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=4G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (output from pre_SNP_calling_cleanup.sh, filename ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

input=/data/scratch/gomeza/analysis/genome_alignment/bowtie
reference=/data/scratch/gomeza/Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}_temp.vcf"
output2="${filename%.*}.vcf"

gatk=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar \
     -R $reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 6 \
     --allow_potentially_misencoded_quality_scores \
     -I $input/*/Ag02/vs_Hg199_minion/Ag02_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag04/vs_Hg199_minion/Ag04_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag05/vs_Hg199_minion/Ag05_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag06/vs_Hg199_minion/Ag06_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/ND9/vs_Hg199_minion/ND9_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/ND8/vs_Hg199_minion/ND8_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R0905/vs_Hg199_minion/R0905_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R37-15/vs_Hg199_minion/R37-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R39-15/vs_Hg199_minion/R39-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R41-15/vs_Hg199_minion/R41-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R42-15/vs_Hg199_minion/R42-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R45-15/vs_Hg199_minion/R45-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag08/vs_Hg199_minion/Ag08_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/P112/vs_Hg199_minion/P112_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/BGV344/vs_Hg199_minion/BGV344_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/OPC304/vs_Hg199_minion/OPC304_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag09_A/vs_Hg199_minion/Ag09_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag11_A/vs_Hg199_minion/Ag11_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Ag11_C/vs_Hg199_minion/Ag11_C_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/Hg199/vs_Hg199_minion/Hg199_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R6-17-2/vs_Hg199_minion/R6-17-2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/R6-17-3/vs_Hg199_minion/R6-17-3_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -o $output

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

java -jar $gatk/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $reference \
   -V $output \
   -o $output2 \


#####################################
# Notes on GATK parallelisation
#####################################
# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
