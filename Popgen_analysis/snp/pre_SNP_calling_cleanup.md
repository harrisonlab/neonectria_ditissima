Sets up correct formatting for SNP calling analysis

Alignment of raw reads vs the Nd genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Nd genome

```bash
  Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/*); do
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
  ```


  ```bash
    Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
    for GenomePath in $(ls -d /data/scratch/gomeza/Nd_genomes/*); do
      Organism=$(echo $GenomePath | rev | cut -f2 -d '/' | rev)
      Strain=$(echo $GenomePath | rev | cut -f1 -d '/' | rev)
      Genome2=$(ls $GenomePath/*nt.gz)
      echo "$Reference - $Genome2"
      OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
      qsub $ProgDir/mauve/run_progressive_mauve.sh $OutDir $Reference $Genome2
    done
    ```
    progressiveMauve --output=analysis/genome_alignment/bowtie/$Organism/$Strain/my_seqs.xmfa repeat_masked/N.ditisima/Ref_Genomes/Hg199_minion/filtered_contigs/N.ditissima_contigs_unmasked.fa /data/scratch/gomeza/Nd_genomes/RS324p/LDPL01.1.fsa_nt

```bash
input=analysis/genome_alignment/bowtie/N.ditissima
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```

## Rename input mapping files in each folder by prefixing with the strain ID

```bash

  for filename in $(ls -d analysis/genome_alignment/bowtie/*/*); do
  Organism=$(echo $filename | rev | cut -f2 -d '/' | rev)
  Strain=$(echo $filename | rev | cut -f1 -d '/' | rev)
  echo "$Organism - $Strain"
      cp "$filename/N.ditissima_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
  done
```

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $scripts/sub_pre_snp_calling.sh <SAMPLE_ID> This needs to use samtools 0.1.18 - hash out 1.5 from profile while this is run

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Hg199 ND8 R0905 R37-15 R45-15
do
    Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    done
    qsub $scripts/sub_pre_snp_calling.sh $input/$Strain/"$Strain"_unmasked.fa_aligned.sam $Strain
done
```


##Copy outputs from cleanup to alignment folder

```bash
for Strain in AgN04 Hg199 R0905 R45-15
do
    #Bam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
    rgBam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
    Bai="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
    Txt="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Hg199_minion/
    mv $Bam $Directory
    mv $rgBam $Directory
    mv $Bai $Directory
    mv $Txt $Directory
done
```
