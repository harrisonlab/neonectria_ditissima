#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -l h_vmem=4G
#$ -l h=blacklace11.blacklace

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
reference=/data/scratch/gomeza/R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa

filename=$(basename "$reference")
output=analysis/popgen/SNP_calling_R0905/"${filename%.*}_temp.vcf"
output2=analysis/popgen/SNP_calling_R0905/"${filename%.*}.vcf"

gatk=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar \
     -R $reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 24 \
     --allow_potentially_misencoded_quality_scores \
     -I $input/*/*/Ag02/Ag02_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag04/Ag04_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag05/Ag05_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag06/Ag06_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag08/Ag08_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag09_A/Ag09_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_A/Ag11_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_B/Ag11_B_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_C/Ag11_C_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/BGV344/BGV344_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Hg199/Hg199_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/ND8/ND8_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/ND9/ND9_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/OPC304/OPC304_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/P112/P112_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R0905/R0905_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R37-15/R37-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R39-15/R39-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R41-15/R41-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R42-15/R42-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R45-15/R45-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R6-17-2/R6-17-2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R6-17-3/R6-17-3_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R68-17-C2/R68-17-C2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R68-17-C3/R68-17-C3_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/SVK1/SVK1_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/SVK2/SVK2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/NMaj/NMaj_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
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
