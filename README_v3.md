Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

The following is a summary of the work presented in this Readme:

Draft Genome assembly Hg199
  * Data qc
  * Genome assembly
  * Repeatmasking

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  for RawData in $(ls raw_dna/paired/*/NG-R0905/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
  Read_F=raw_dna/paired/N.ditissima/NG-R0905/F/*.fastq.gz
  Read_R=raw_dna/paired/N.ditissima/NG-R0905/R/*.fastq.gz
  IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
  qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/*/NG-R0905/*/*.fq.gz); do
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
  qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

```bash
  Reference=$(ls repeat_masked/N.*/R0905_pacbio_canu/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/NG-R0905); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_NG-R0905_pacbio2
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
  ```
