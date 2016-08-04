
Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

The following is a summary of the work presented in this Readme:

Draft Genome assembly Hg199
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
  mkdir -p raw_dna/paired/N.ditissima/R0905/F
  mkdir -p raw_dna/paired/N.ditissima/R0905/R
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
    for RawData in $(ls raw_dna/paired/*/Hg199/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
    Read_F=raw_dna/paired/N.ditissima/Hg199/F/*.fastq.gz
    Read_R=raw_dna/paired/N.ditissima/Hg199/R/*.fastq.gz
    IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/*/Hg199/*/*.fq.gz); do
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
  qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
	Trim_F=qc_dna/paired/N.ditissima/Hg199/F/*.fq.gz
	Trim_R=qc_dna/paired/N.ditissima/Hg199/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

mode kmer abundance prior to error correction was reported using the following commands:
```bash
    for File in $(ls qc_dna/kmc/*/Hg199/*_true_kmer_summary.txt); do
        basename $File;
        cat $File | grep -e 'abundance' -e 'size'
    done
```

** Estimated Genome Size is: 42917594

** Esimated Coverage is: 74


#Assembly
Assembly was performed using: Velvet / Abyss / Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
	F_Read=qc_dna/paired/N.ditissima/Hg199/F/*.fq.gz
	R_Read=qc_dna/paired/N.ditissima/Hg199/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/Hg199/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
```

cat contigs.fasta | grep '>' | wc -l
1204

cat contigs_min_500bp.fasta | grep '>' | wc -l
849


```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/Hg199/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:182348
  * N80:
  * N20:
  * Longest contig:614761
  **


# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss=/assembly/spades/N.ditissima/Hg199/filtered_contigs/contigs_min_500bp.fasta
	qsub $ProgDir/rep_modeling.sh $BestAss
	qsub $ProgDir/transposonPSI.sh $BestAss
 ```

** % bases masked by repeatmasker: 10.88% (bases masked:4941726 bp)

** % bases masked by transposon psi: 9.48% (bases masked:4306052 bp)

Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
  for File in $(ls repeat_masked/*/Hg199/filtered_contigs_repmask/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done

    repeat_masked/N.ditissima/Hg199/filtered_contigs_repmask/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
    Number of masked bases:  5062507
```


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
