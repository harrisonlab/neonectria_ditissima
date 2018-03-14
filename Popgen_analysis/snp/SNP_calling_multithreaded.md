#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

```bash
input=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs
reference=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs/N.ditissima_contigs_unmasked.fa

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
for Strain in Ag02 Ag04 Ag05 Ag06 Hg199 ND8 R0905 R37-15 R45-15
do
    Index=repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/filtered_contigs/N.ditissima_contigs_unmasked.fa.fai
    Directory=analysis/genome_alignment/bowtie/*/$Strain/
    cp $Index $Directory
done
``

##Move to the directory where the output of SNP calling should be placed

```bash
mkdir -p /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling
cd /home/groups/harrisonlab/project_files/neonectria_ditissima/SNP_calling
```

##Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See inside the submission script below:

```bash
scripts=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
qsub $scripts/sub_SNP_calling_multithreaded.sh
```
