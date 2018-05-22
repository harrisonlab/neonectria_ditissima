## Sets up correct formatting for SNP calling analysis

```bash
input=analysis/genome_alignment/bowtie/N.ditissima
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```

## Alignment of raw reads vs the Nd genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Nd genome

```bash
for Strain in Ag08 Ag11_B R41-15 R6-17-2 R6-17-3; do
#for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
#for Strain in Ag02 Ag05 ND8 R37-15 Ag04 R45-15 R0905 Hg199; do
#Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
Reference=$(ls Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/*_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/N.ditissima/$Strain); do
  Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
  Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
  echo "$Organism - $Strain"
  F_Read=$(ls $StrainPath/F/*.fq.gz)
  R_Read=$(ls $StrainPath/R/*.fq.gz)
  echo $F_Read
  echo $R_Read
  OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Hg199_minion
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
  qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
done
```

## Rename input mapping files in each folder by prefixing with the strain ID

```bash
for filename in $(ls -d analysis/genome_alignment/bowtie/*/*/vs_Hg199_minion); do
Organism=$(echo $filename | rev | cut -f3 -d '/' | rev)
Strain=$(echo $filename | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
  cp "$filename/N.ditissima_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
  cp "$filename/N.ditissima_contigs_unmasked.fa_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
  cp "$filename/N.ditissima_contigs_unmasked.fa_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
  cp "$filename/N.ditissima_contigs_unmasked.fa_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
done
```

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $scripts/sub_pre_snp_calling.sh <SAMPLE_ID> This needs to use samtools 0.1.18 - hash out 1.5 from profile while this is run

```bash
#for Strain in Ag02 Ag04 Ag05 Ag06 Hg199 ND8 R0905 R37-15 R45-15
#for Strain in Ag06 Ag09_A Ag11_A Ag11_C BGV344 ND9 OPC304 P112 R39-15 R42-15 R68-17
for Strain in Ag08 Ag11_B R41-15 R6-17-2 R6-17-3
do
    Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    done
    qsub $scripts/sub_pre_snp_calling_maria.sh $input/$Strain/vs_Hg199_minion/"$Strain"_unmasked.fa_aligned.sam $Strain
done
```

## Copy outputs from cleanup to alignment folder

```bash
for Strain in Ag08 Ag11_B R41-15 R6-17-2 R6-17-3
#for Strain in Ag02 Ag04 Ag05 Ag06 Hg199 ND8 R0905 R37-15 R45-15
#for Strain in Ag06 Ag09_A Ag11_A Ag11_C BGV344 ND9 OPC304 P112 R39-15 R42-15 R68-17
do
  #Bam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
  rgBam="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
  Bai="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
  Txt="$Strain"_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
  Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Hg199_minion/
  mv $Bam $Directory
  mv $rgBam $Directory
  mv $Bai $Directory
  mv $Txt $Directory
done
```
