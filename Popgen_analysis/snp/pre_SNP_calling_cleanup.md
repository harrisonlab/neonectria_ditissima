## Sets up correct formatting for SNP calling analysis

```bash
input=analysis/genome_alignment/bowtie
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```

## Alignment of raw reads vs the Nd genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Nd genome

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R0905_v2 R0905_all R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj; do
#Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
#New genome version was copied to the REFERENCE folder.
Reference=$(ls REFERENCE/Hg199_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/N.*/$Strain); do
  Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
  Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
  echo "$Organism - $Strain"
  F_Read=$(ls $StrainPath/F/*.fq.gz)
  R_Read=$(ls $StrainPath/R/*.fq.gz)
  echo $F_Read
  echo $R_Read
  OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
  qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
done
```
For R0905, I used R0905_all.

R0905 genome as reference

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905_all R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj; do
Reference=$(ls R0905_good/R0905_contigs_unmasked.fa)
for StrainPath in $(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/qc_dna/paired/N.*/$Strain); do
  Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
  Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
  echo "$Organism - $Strain"
  F_Read=$(ls $StrainPath/F/*.fq.gz)
  R_Read=$(ls $StrainPath/R/*.fq.gz)
  echo $F_Read
  echo $R_Read
  OutDir=analysis/genome_alignment/bowtie/$Organism/R0905/$Strain/
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
  qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
done
```

## Rename input mapping files in each folder by prefixing with the strain ID

```bash
for filename in $(ls -d analysis/genome_alignment/bowtie/N.*/R68-17-C*); do
Organism=$(echo $filename | rev | cut -f2 -d '/' | rev)
Strain=$(echo $filename | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
  mv "$filename/Hg199_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
  mv "$filename/Hg199_contigs_unmasked.fa_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
  mv "$filename/Hg199_contigs_unmasked.fa_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
  mv "$filename/Hg199_contigs_unmasked.fa_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
done
```

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $scripts/sub_pre_snp_calling.sh <SAMPLE_ID> This needs to use samtools 0.1.18 - hash out 1.5 from profile while this is run

```bash
for Strain in R68-17-C2 R68-17-C3; do
#for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17 SVK1 SVK2 NMaj; do
    Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    done
    qsub $scripts/sub_pre_snp_calling_maria.sh $input/N.*/$Strain/"$Strain"_unmasked.fa_aligned.sam $Strain
done
```

## Copy outputs from cleanup to alignment folder

```bash
for Strain in R68-17-C2 R68-17-C3; do
#for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17 SVK1 SVK2 NMaj; do
  #Bam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
  rgBam="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
  Bai="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
  Txt="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
  Directory=analysis/genome_alignment/bowtie/*/$Strain/
  #mv $Bam $Directory
  mv $rgBam $Directory
  mv $Bai $Directory
  mv $Txt $Directory
done
```
