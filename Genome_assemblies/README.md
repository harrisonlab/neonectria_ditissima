# neonectria_ditissima
Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
  cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  Species=N.ditissima
  Strain=NG-R0905
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp /home/groups/harrisonlab/project_files/neonectria/NG-R0905_S4_L001_R1_001.fastq raw_dna/paired/$Species/$Strain/F/.
  cp /home/groups/harrisonlab/project_files/neonectria/NG-R0905_S4_L001_R2_001.fastq raw_dna/paired/$Species/$Strain/R/.
  mkdir -p raw_dna/paired/N.ditissima/Hg199/F
  mkdir -p raw_dna/paired/N.ditissima/Hg199/R
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160725_M04465_0020-AP1N0
  cp $RawDat/Hg199_S1_L001_R1_001.fastq.gz raw_dna/paired/N.ditissima/Hg199/F/.
  cp $RawDat/Hg199_S1_L001_R2_001.fastq.gz raw_dna/paired/N.ditissima/Hg199/R/.
```

RNA data was also copied into the project area:

```bash
  cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  mkdir -p raw_rna/paired/N.ditissima/Hg199/F
  mkdir -p raw_rna/paired/N.ditissima/Hg199/R
  mkdir -p raw_rna/paired/N.ditissima/R0905/F
  mkdir -p raw_rna/paired/N.ditissima/R0905/R
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160718_M04465_0018_ANU02
  cp $RawDat/Hg199_S1_L001_R1_001.fastq.gz raw_rna/paired/N.ditissima/Hg199/F/.
  cp $RawDat/Hg199_S1_L001_R2_001.fastq.gz raw_rna/paired/N.ditissima/Hg199/R/.
  cp $RawDat/R09-05_S2_L001_R1_001.fastq.gz raw_rna/paired/N.ditissima/R0905/F/.
  cp $RawDat/R09-05_S2_L001_R2_001.fastq.gz raw_rna/paired/N.ditissima/R0905/R/.
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq); do
	  echo $RawData;
	  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	  qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
	Read_F=raw_dna/paired/N.ditissima/NG-R0905/F/NG-R0905_S4_L001_R1_001.fastq
	Read_R=raw_dna/paired/N.ditissima/NG-R0905/R/NG-R0905_S4_L001_R2_001.fastq
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

#Assembly
Assembly was performed using: Velvet / Abyss / Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
	F_Read=qc_dna/paired/N.ditissima/NG-R0905/F/NG-R0905_qc_F.fastq.gz
	R_Read=qc_dna/paired/N.ditissima/NG-R0905/R/NG-R0905_qc_R.fastq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/NG-R0905/
	qsub $ProgDir/sumit_SPAdes2.sh $F_Read $R_Read $Outdir correct
```

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/NG-R0905/
	qsub $ProgDir/sumit_SPAdes3.sh $F_Read $R_Read $Outdir correct
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:126848
  * N80:
  * N20:
  * Longest contig:728907
  **

# Filter contigs

The assembled contigs were filtered to remove all contigs shorter than 1kb from
the assembly. This was done using the following commands:

```bash
	InDir=assembly/spades/N.ditissima/NG-R0905
	OutDir=assembly/spades/N.ditissima/NG-R0905_filtered
	mkdir –p $OutDir
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
  	Assembly=$InDir/scaffolds.fasta
  	AssFiltered=$OutDir/scaffolds_filtered_500.fasta
  	$ProgDir/filter_abyss_contigs.py $Assembly 500 > $AssFiltered
  	AssFiltered=$OutDir/scaffolds_filtered_1000.fasta
  	$ProgDir/filter_abyss_contigs.py $Assembly 1000 > $AssFiltered
```

# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss=/assembly/spades/N.ditissima/NG-R0905/scaffolds.fasta
	qsub $ProgDir/rep_modeling.sh $BestAss
	qsub $ProgDir/transposonPSI.sh $BestAss
 ```

** % bases masked by repeatmasker: 12.53 %

** % bases masked by transposon psi: **


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Neonectria genome. This used results from CEGMA as hints for gene models.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
  	Assembly=repeat_masked/spades/N.ditissima/NG-R0905_repmask/N.ditissima_contigs_unmasked.fa
  	qsub $ProgDir/sub_cegma.sh $Assembly dna
```

** Number of cegma genes present and complete: 95.56
** Number of cegma genes present and partial: 97.18

##Gene prediction

Gene prediction was performed for the neonectria genome.
CEGMA genes were used as Hints for the location of CDS.

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/augustus
  	Assembly=repeat_masked/spades/N.ditissima/NG-R0905_repmask/N.ditissima_contigs_unmasked.fa
  	GeneModel=fusarium_graminearum
  	qsub $ProgDir/submit_augustus.sh $GeneModel $Assembly
```

** Number of genes predicted: 13589

#Functional annotation

Interproscan was used to give gene models functional annotations.

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
  	Genes=gene_pred/augustus/spades/N.ditissima/N.ditissima_aug_out.aa
  	$ProgDir/sub_interproscan.sh $Genes
```

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	Genes=gene_pred/augustus/spades/N.ditissima/N.ditissima_aug_out.aa
	InterProRaw=gene_pred/interproscan/spades/N.ditissima/raw
	ProgDir/append_interpro.sh $Genes $InterProRaw
```

#Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/blast
	Query=../../phibase/v3.8/PHI_accessions.fa
	Subject=repeat_masked/spades/N.ditissima/NG-R0905_repmask/N.ditissima_contigs_unmasked.fa
	qsub $ProgDir/blast_pipe.sh $Query protein $Subject
```

Top BLAST hits were used to annotate gene models.

```bash

```

** Blast results of note: **
  * 'Result A'
