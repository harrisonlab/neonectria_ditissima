#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

```bash
reference=REFERENCE/Hg199_contigs_unmasked.fa
input=REFERENCE/
#input=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs
#reference=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs/N.ditissima_contigs_unmasked.fa
filename=$(basename "$reference")
output="${filename%.*}.dict"
```
```bash
reference=R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa
input=R0905_good/repeat_masked/filtered_contigs/
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
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj   
do
  Index=REFERENCE/Hg199_contigs_unmasked.fa.fai
  Directory=analysis/genome_alignment/bowtie/N.*/$Strain
  cp $Index $Directory
done
```
```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj   
do
  Index=R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa.fai
  Directory=analysis/genome_alignment/bowtie/N.*/R0905/$Strain
  cp $Index $Directory
done
```
##Move to the directory where the output of SNP calling should be placed

```bash
mkdir -p analysis/popgen/SNP_calling3
cd analysis/popgen/SNP_calling3
#mkdir -p /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling
#cd /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling

mkdir -p analysis/popgen/SNP_calling_R0905
cd analysis/popgen/SNP_calling_R0905
```

##Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See inside the submission script below:

All isolates. Complete. Hg199 ref.
```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
qsub $scripts/sub_SNP_calling_multithreaded2.sh
```

All isolates. Complete. R0905 ref.
```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
qsub -R y $scripts/sub_SNP_calling_multithreaded.sh
```
