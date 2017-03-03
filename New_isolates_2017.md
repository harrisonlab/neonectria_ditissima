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
	for RawData in $(ls qc_dna/paired/*/new_isolates_2017/*/*/*.fq.gz); do
	echo $RawData;
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
	Trim_F=qc_dna/paired/N.ditissima/*/ND8/F/*.fq.gz
	Trim_R=qc_dna/paired/N.ditissima/*/ND8/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

** Estimated Genome Size is: 946768573

** Esimated Coverage is: 5

```bash
	Trim_F=qc_dna/paired/N.ditissima/*/AgN04/F/*.fq.gz
	Trim_R=qc_dna/paired/N.ditissima/*/AgN04/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

** Estimated Genome Size is: 48203823

** Esimated Coverage is: 30

```bash
	Trim_F=qc_dna/paired/N.ditissima/*/R45-15/F/*.fq.gz
	Trim_R=qc_dna/paired/N.ditissima/*/R45-15/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

** Estimated Genome Size is: 43409438

** Esimated Coverage is: 86


#Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash
	F_Read=qc_dna/paired/N.ditissima/*/R45-15/F/*.fq.gz
	R_Read=qc_dna/paired/N.ditissima/*/R45-15/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/R45-15/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
```

cat contigs.fasta | grep '>' | wc -l
1427

cat contigs_min_500bp.fasta | grep '>' | wc -l
969

```bash
	F_Read=qc_dna/paired/N.ditissima/*/AgN04/F/*.fq.gz
	R_Read=qc_dna/paired/N.ditissima/*/AgN04/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/AgN04/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
```

cat contigs.fasta | grep '>' | wc -l
1106

cat contigs_min_500bp.fasta | grep '>' | wc -l
728

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/R45-15/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/AgN04/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
  NewSpades=$(ls assembly/spades/*/R45-15/filtered_contigs/contigs_min_500bp.fasta)
  for Assembly in $(ls $NewSpades); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    OutDir=repeat_masked/N.ditissima/R45-15/
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```

** % bases masked by repeatmasker: 10.73% (bases masked:4941726 bp)

** % bases masked by transposon psi: 9.33% (bases masked:4306052 bp)

```bash
  NewSpades2=$(ls assembly/spades/*/AgN04/filtered_contigs/contigs_min_500bp.fasta)
  for Assembly in $(ls $NewSpades2); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    OutDir=repeat_masked/N.ditissima/AgN04/
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```

** % bases masked by repeatmasker: 12.85% (bases masked:4941726 bp)

** % bases masked by transposon psi: 11.49% (bases masked:4306052 bp)

Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
  for File in $(ls repeat_masked/*/*/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done
```
    repeat_masked/N.ditissima/AgN04/AgN04_contigs_softmasked_repeatmasker_TPSI_appended.fa
    Number of masked bases:
    6076978
    repeat_masked/N.ditissima/R0905_merged_2017/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
    Number of masked bases:
    5395566
    repeat_masked/N.ditissima/R45-15/R45-15_contigs_softmasked_repeatmasker_TPSI_appended.fa
    Number of masked bases:
    4967354

# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
	for Genome in $(ls repeat_masked/*/R45-15/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```
** Number of cegma genes present and complete: 95.56%
** Number of cegma genes present and partial: 96.77%

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
	for Genome in $(ls repeat_masked/*/AgN04/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```
** Number of cegma genes present and complete: %
** Number of cegma genes present and partial: %

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/N.*/Hg199/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```
