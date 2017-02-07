# neonectria_ditissima
Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
cd /home/groups/harrisonlab/project_files/neonectria_ditissima
Species=N.ditissima
mkdir -p raw_dna/paired/$Species/ND8/F
mkdir -p raw_dna/paired/$Species/ND8/R
mkdir -p raw_dna/paired/$Species/AgN04/F
mkdir -p raw_dna/paired/$Species/AgN04/R
mkdir -p raw_dna/paired/$Species/R45-15/F
mkdir -p raw_dna/paired/$Species/R45-15/R
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/ND8_S3_L001_R1_001.fastq.gz raw_dna/paired/$Species/ND8/F/
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/ND8_S3_L001_R2_001.fastq.gz raw_dna/paired/$Species/ND8/R/
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/NdAgO4_S1_L001_R1_001.fastq.gz raw_dna/paired/$Species/AgN04/F/
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/NdAgO4_S1_L001_R2_001.fastq.gz raw_dna/paired/$Species/AgN04/R/
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/R45-15_S2_L001_R1_001.fastq.gz raw_dna/paired/$Species/R45-15/F/  
cp /home/miseq_data/2017/RAW/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/R45-15_S2_L001_R2_001.fastq.gz raw_dna/paired/$Species/R45-15/R/
```

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  cd raw_dna/paired/N.ditissima/
  mkdir new_isolates_2017
  mv ND8/ new_isolates_2017/
  mv AgN04/ new_isolates_2017/
  mv R45-15/ new_isolates_2017/
  for RawData in $(ls raw_dna/paired/*/new_isolates_2017/*/*/*.fastq.gz); do
	  echo $RawData;
	  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	  qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
	Read_F=raw_dna/paired/N.ditissima/*/ND8/F/*.fastq.gz
	Read_R=raw_dna/paired/N.ditissima/*/ND8/R/*.fastq.gz
	IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
	qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```
```bash
	Read_F=raw_dna/paired/N.ditissima/*/AgN04/F/*.fastq.gz
	Read_R=raw_dna/paired/N.ditissima/*/AgN04/R/*.fastq.gz
	IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
	qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```
```bash
	Read_F=raw_dna/paired/N.ditissima/*/R45-15/F/*.fastq.gz
	Read_R=raw_dna/paired/N.ditissima/*/R45-15/R/*.fastq.gz
	IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
	qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fastq.gz); do
	echo $RawData;
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
	Trim_F=qc_dna/paired/N.ditissima/NG-R0905/F/NG-R0905_qc_F.fastq.gz
	Trim_R=qc_dna/paired/N.ditissima/NG-R0905/R/NG-R0905_qc_R.fastq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

** Estimated Genome Size is: 48490129

** Esimated Coverage is: 42
