# SNP calling analysis

This file contains the SNP calling using the traditional method

## Sets up correct formatting for SNP calling analysis

```bash
    cd /data/scratch/gomeza

    # Latest isolates
  for Strain in 118923 118924 226-31 227-31; do
    Reference=$(ls R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa)
    for StrainPath in $(ls -d qc_dna/paired/N.*/$Strain); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=Home/analysis_vAG/genome_alignment/bowtie/vs_R0905/$Organism/R0905/$Strain/
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
    sbatch $ProgDir/bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
    done
  done
```


## Rename input mapping files in each folder by prefixing with the strain ID

```bash
mv analysis/genome_alignment/bowtie/N.ditissima/R0905/R0905_all analysis/genome_alignment/bowtie/N.ditissima/R0905/R0905
```

```bash
for Strain in 118923 118924 226-31 227-31; do
for filename in $(ls -d Home/analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/$Strain); do
Organism=$(echo $filename | rev | cut -f3 -d '/' | rev)
Strain=$(echo $filename | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
mv "$filename/R0905_good_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
mv "$filename/R0905_good_contigs_unmasked.fa_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
mv "$filename/R0905_good_contigs_unmasked.fa_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
mv "$filename/R0905_good_contigs_unmasked.fa_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
done
done
```

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $scripts/sub_pre_snp_calling.sh <SAMPLE_ID> This needs to use samtools 0.1.18 - hash out 1.5 from profile while this is run

```bash
    for Strain in 118923 118924 226-31 227-31; do
        for input in Home/analysis_vAG/genome_alignment/bowtie/vs_*/N.d*/R0905/$Strain/"$Strain"_unmasked.fa_aligned.sam; do
        Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        while [ $Jobs -gt 5 ]; do
        sleep 5m
        printf "."
        Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        done
        printf "\n"
        OutDir=Home/analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/$Strain
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
        sbatch $ProgDir/pre_SNP_calling.sh $input $Strain $OutDir
        done
    done
# Some errors came out related to profile, but it worked
```

## Prepare genome reference indexes required by GATK

<!-- ```bash
reference=R0905_contigs_unmasked.fa
input=./
filename=$(basename "$reference")
output="${filename%.*}.dict"
picard CreateSequenceDictionary R=$reference O=$input/$output
samtools faidx $reference
``` -->

### Copy index file to same folder as BAM alignments

```bash
for Strain in 118923 118924 226-31 227-31; do
Index=R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa.fai
Directory=Home/analysis_vAG/genome_alignment/bowtie/vs_R0905/N.d*/R0905/$Strain
cp $Index $Directory
done
```

## Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See GATK_SNP_calling.sh

```bash
mkdir -p Home/analysis_vAG/SNPs/SNP_calling_R0905_VP
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/GATK_SNP_calling.sh
```
```
less R0905_good_contigs_unmasked_temp.vcf | grep -v '^#' | wc -l
1631836
```
# Post SNP calling

## Filter vcf outputs, only retain biallelic high-quality SNPS with no missing data for genetic analyses.

```bash
    for vcf in $(ls Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/R0905_good_contigs_unmasked_temp.vcf)
    do
    echo $vcf
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
    sbatch $ProgDir/vcf_parser.sh $vcf
    done
```
```
After filtering, kept 32 out of 32 Individuals
Outputting VCF file...
After filtering, kept 885411 out of a possible 1417258 Sites
Run Time = 187.00 seconds
```
```
less R0905_good_contigs_unmasked_temp_filtered.vcf | grep -v '^#' | wc -l
885411
```
## Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls Home/analysis_vAG/SNPs/SNP_calling_R0905_VP//*_filtered.vcf)
do
  scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

# Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/*distance.log)
do
scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
Rscript --vanilla $scripts/distance_matrix.R $log
done
```

# Remove low coverage samples

```bash
#Remove low-coverage isolates. Create a cut-down VCF and filter it
#The isolate Ag11_B has a 22X coverage, so it will be removed and new statistics calculated.
cd Home/analysis_vAG/SNPs/SNP_calling_R0905_VP
mkdir NoAg11_B/
cp *.vcf NoAg11_B/
vcfremovesamples R0905_good_contigs_unmasked_temp.vcf Ag11_B > R0905_good_contigs_unmasked_NoAg11B.vcf
vcfremovesamples R0905_good_contigs_unmasked_temp_filtered.vcf Ag11_B > R0905_good_contigs_unmasked_NoAg11B_filtered.vcf
```
```
less R0905_good_contigs_unmasked_NoAg11B.vcf | grep -v '^#' | wc -l
1631836
less R0905_good_contigs_unmasked_NoAg11B_filtered.vcf | grep -v '^#' | wc -l
885411
```

# Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.vcf)
do
    echo $vcf
    out=$(basename $vcf .vcf)
    echo $out
    vcftools --vcf $vcf --mac 1 --recode --out $out
done

After filtering, kept 31 out of 31 Individuals
Outputting VCF file...
After filtering, kept 727128 out of a possible 885411 Sites
```

# Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.recode.vcf)
do
  scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

# Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls *distance.log)
do
  scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
  Rscript --vanilla $scripts/distance_matrix.R $log
done
```

# Remove outgroup and low coverage samples

```bash
cd Home/analysis_vAG/SNPs/SNP_calling_R0905_VP
mkdir NoNMaj_NoAg11_B/
cp NoAg11_B/R0905_good_contigs_unmasked_NoAg11B_filtered.recode.vcf NoNMaj_NoAg11_B/
cd NoNMaj_NoAg11_B
vcfremovesamples R0905_good_contigs_unmasked_NoAg11B_filtered.recode.vcf NMaj > R0905_good_contigs_unmasked_FINAL.vcf

vcfremovesamples R0905_good_contigs_unmasked_temp.vcf Ag11_B NMaj > R0905_good_contigs_unmasked_FINAL.vcf
vcfremovesamples R0905_good_contigs_unmasked_temp_filtered.vcf Ag11_B NMaj > R0905_good_contigs_unmasked_FINAL_filtered.vcf
```
```
less R0905_good_contigs_unmasked_FINAL.vcf | grep -v '^#' | wc -l
1631836
less R0905_good_contigs_unmasked_FINAL_filtered.vcf | grep -v '^#' | wc -l
885411
```

# Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_FINAL_filtered.vcf)
do
echo $vcf
out=$(basename $vcf .vcf)
echo $out
vcftools --vcf $vcf --mac 1 --recode --out $out
done

After filtering, kept 30 out of 30 Individuals
Outputting VCF file...
After filtering, kept 448001 out of a possible 727128 Sites
```

# Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_FINAL_filtered.recode.vcf)
do
  scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

# Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls *distance.log)
do
scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
Rscript --vanilla $scripts/distance_matrix.R $log
done
```

### Create custom SnpEff genome database


```bash
snpeff=/data/scratch/gomeza/prog/snpEff
nano $snpeff/snpEff.config
```

Add the following lines to the section with databases:
```
Nd_R0905.genome : Nd_R0905
```

```bash
# Collect input files
mkdir -p $snpeff/data/Nd_R0905
cp R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa $snpeff/data/Nd_R0905
cp gene_pred/codingquary/Ref_Genomes_v2/N.ditissima/R0905/final/final_genes_appended_renamed.gff3 $snpeff/data/Nd_R0905

#Rename input files
cd $snpeff/data/Nd_R0905
mv final_genes_appended_renamed.gff3 genes.gff
mv R0905_good_contigs_unmasked.fa sequences.fa

#Build database using GFF3 annotation
java -jar $snpeff/snpEff.jar build -gff3 -v Nd_R0905
```


# Annotate VCF files

```bash
cd Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/NoNMaj
for a in R0905_good_contigs_unmasked_FINAL_filtered.recode.vcf
do
echo $a
filename=$(basename "$a")
java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 Nd_R0905 $a > ${filename%.vcf}_annotated.vcf
mv snpEff_genes.txt snpEff_genes_${filename%.vcf}.txt
mv snpEff_summary.html snpEff_summary__${filename%.vcf}.html
done
```
```
less R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.vcf  | grep -v '^#' | wc -l
448001
```