# Neonectria_Reference_Genome_Assembly
==========

This document details the commands used to assemble and annotate the Hg199 Neonectria genome.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/N.ditissima

# 0. Building of directory structure

screen -a

```bash
  RawDatDir=/home/miseq_data/minion/2017/*_Neonectria_Hg199/fast5/pass
  Organism=N.ditissima
  Strain=Hg199
  Date=17-07-17
  mkdir -p raw_dna/minion/$Organism/$Strain/$Date
  for Fast5Dir in $(ls -d $RawDatDir/*); do
      poretools fastq $Fast5Dir | gzip -cf
    done > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_pass.fastq.gz
  ```

The following is a summary of the work presented in this Readme.

The following processes were applied to Neonectria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/minion/*/*/*.fastq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```
