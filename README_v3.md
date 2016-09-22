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
