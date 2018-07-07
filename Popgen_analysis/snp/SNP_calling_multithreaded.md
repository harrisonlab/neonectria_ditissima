#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

```bash
reference=Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa
input=Hg199_genome/repeat_masked/N.ditissima/Hg199_minion
#input=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs
#reference=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs/N.ditissima_contigs_unmasked.fa
filename=$(basename "$reference")
output="${filename%.*}.dict"
```

##Prepare genome reference indexes required by GATK

```bash
java -jar /home/gomeza/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$reference O=$input/$output
samtools faidx $reference
```

###Copy index file to same folder as BAM alignments

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 ND8 ND9 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R68-17 Ag08 P112 BGV344 OPC304 Ag09_A Ag11_A Ag11_B Ag11_C Hg199 R6-17-2 R6-17-3
do
  Index=Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa.fai
  Directory=analysis/genome_alignment/bowtie/N.ditissima/$Strain/vs_Hg199_minion
  cp $Index $Directory
done
```

##Move to the directory where the output of SNP calling should be placed

```bash
mkdir -p /data/scratch/gomeza/analysis/popgen/SNP_calling
cd /data/scratch/gomeza/analysis/popgen/SNP_calling
#mkdir -p /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling
#cd /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling
```

##Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See inside the submission script below:

```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
qsub $scripts/sub_SNP_calling_multithreaded.sh
```
Removing R68-17 and Ag11_B
```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
qsub $scripts/sub_SNP_calling_multithreaded2.sh
```
