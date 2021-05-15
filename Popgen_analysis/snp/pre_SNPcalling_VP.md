# SNP calling VP

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
mkdir -p analysis/popgen/SNP_calling_R0905_VP
#cd analysis/popgen/SNP_calling
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/GATK_SNP_calling.sh
sbatch $ProgDir/HaplotypeCaller_GVCF.sh
```

```bash
Reference=R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa 
#for Strain in 118923 226-31 Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C3; do
for Strain in SVK2 118924 227-31 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C2 SVK1 NMaj; do
for Bam in $(ls Home/analysis_vAG/genome_alignment/bowtie/vs_R0905/*/R0905/$Strain/"$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam); do
echo $Strain
Outdir=/projects/neonectria_ditissima/HaplotypeCaller/$Strain
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch -p long $ProgDir/HaplotypeCaller_GVCF.sh $Reference $Strain $Bam $Outdir
done 
done  
```

# SNP calling analysis


### Filter vcf outputs, only retain biallelic high-quality SNPS with no missing data for genetic analyses.

```bash
for vcf in $(ls analysis/popgen/SNP_calling/*_contigs_unmasked.vcf)
do
echo $vcf
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/vcf_parser.sh $vcf
done
```
