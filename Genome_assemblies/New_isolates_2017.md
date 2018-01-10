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

mkdir -p raw_dna/paired/$Species/Ag02/F
mkdir -p raw_dna/paired/$Species/Ag02/R
mkdir -p raw_dna/paired/$Species/Ag05/F
mkdir -p raw_dna/paired/$Species/Ag05/R
mkdir -p raw_dna/paired/$Species/R37-15/F
mkdir -p raw_dna/paired/$Species/R37-15/R


cd /data/seq_data/miseq/2017/RAW
Species=N.ditissima

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/AG02_S1_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/AG02_S1_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/AG02_S1_L001_R1_001.fastq.gz > raw_dna/paired/$Species/Ag02/F/AG02_S1_L001_R1_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/AG02_S1_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/AG02_S1_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/AG02_S1_L001_R2_001.fastq.gz > raw_dna/paired/$Species/Ag02/R/AG02_S1_L001_R2_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/AG05_S2_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/AG05_S2_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/AG05_S2_L001_R1_001.fastq.gz > raw_dna/paired/$Species/Ag05/F/AG05_S2_L001_R1_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/AG05_S2_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/AG05_S2_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/AG05_S2_L001_R2_001.fastq.gz > raw_dna/paired/$Species/Ag05/R/AG05_S2_L001_R2_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/ND8_S3_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/ND8_S3_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/ND8_S3_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/ND8_S3_L001_R1_001.fastq.gz > raw_dna/paired/$Species/ND8/F/ND8_S3_L001_R1_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/ND8_S3_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/ND8_S3_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/ND8_S3_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/170203_M04465_0032_000000000-ATMUR/Data/Intensities/BaseCalls/ND8_S3_L001_R2_001.fastq.gz > raw_dna/paired/$Species/ND8/R/ND8_S3_L001_R2_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/R3715_S4_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/R3715_S4_L001_R1_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/R3715_S4_L001_R1_001.fastq.gz > raw_dna/paired/$Species/R37-15/F/R3715_S4_L001_R1_001_allfiles.fastq.gz

cat /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls/R3715_S4_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171205_M04465_0058_000000000-BJ49V/Data/Intensities/BaseCalls/R3715_S4_L001_R2_001.fastq.gz /data/seq_data/miseq/2017/ANALYSIS/171026_M04465_0052_000000000-BFR28/Data/Intensities/BaseCalls/R3715_S4_L001_R2_001.fastq.gz > raw_dna/paired/$Species/R37-15/R/R3715_S4_L001_R2_001_allfiles.fastq.gz

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

  for RawData in $(ls raw_dna/paired/*/R37-15/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
  Read_F=raw_dna/paired/N.ditissima/Ag02/F/*allfiles.fastq.gz
  Read_R=raw_dna/paired/N.ditissima/Ag02/R/*allfiles.fastq.gz
  IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
  qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA

  Read_F=raw_dna/paired/N.ditissima/Ag05/F/*allfiles.fastq.gz
  Read_R=raw_dna/paired/N.ditissima/Ag05/R/*allfiles.fastq.gz
  IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
  qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA

	Read_F=raw_dna/paired/N.ditissima/R37-15/F/*allfiles.fastq.gz
	Read_R=raw_dna/paired/N.ditissima/R37-15/R/*allfiles.fastq.gz
	IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
	qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA

  Read_F=raw_dna/paired/N.ditissima/ND8/F/*allfiles.fastq.gz
  Read_R=raw_dna/paired/N.ditissima/ND8/R/*allfiles.fastq.gz
  IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
  qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/*/*/*/*allfiles_trim.fq.gz); do
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


```bash
  # for Strain in 415 416 A4 SCRP245_v2 Bc23 Nov5 Nov77; do
  for Strain in Ag02; do
    echo $Strain
    Trim_F=$(ls qc_dna/paired/N.*/$Strain/F/*allfiles_trim.fq.gz)
    Trim_R=$(ls qc_dna/paired/N.*/$Strain/R/*allfiles_trim.fq.gz)
    Outdir=qc_dna/kmc/N.ditissima/Ag02/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh 21 $Outdir $Trim_F $Trim_R
  done

  for Strain in Ag02; do
    echo $Strain
    Trim_F=$(ls qc_dna/paired/N.*/$Strain/F/*allfiles_trim.fq.gz)
    Trim_R=$(ls qc_dna/paired/N.*/$Strain/R/*allfiles_trim.fq.gz)
    Outdir=qc_dna/kmc/N.ditissima/paired/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh 41 $Outdir $Trim_F $Trim_R
  done
```

```bash
for Strain in Ag05 ND8 R37-15; do
#for Strain in Ag02; do
  echo $Strain
  Trim_F=$(ls qc_dna/paired/N.*/$Strain/F/*allfiles_trim.fq.gz)
  Trim_R=$(ls qc_dna/paired/N.*/$Strain/R/*allfiles_trim.fq.gz)
  Outdir=qc_dna/kmc/N.ditissima/$Strain
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  qsub $ProgDir/kmc_kmer_counting.sh 21 $Outdir $Trim_F $Trim_R
done
```

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

```bash
  for Strain in Ag02 Ag05 ND8 R37-15; do
  echo $Strain
  F_Read=$(ls qc_dna/paired/N.*/$Strain/F/*allfiles_trim.fq.gz)
  R_Read=$(ls qc_dna/paired/N.*/$Strain/R/*allfiles_trim.fq.gz)
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/$Strain/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
done
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/Ag02/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/Ag05/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/ND8/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/R37-15/filtered_contigs/contigs_min_500bp.fasta); do
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

```bash
  for Strain in Ag02 Ag05 ND8 R37-15; do
  echo $Strain
  Assembly=$(ls assembly/spades/*/$Strain/filtered_contigs/contigs_min_500bp.fasta)
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
  OutDir=repeat_masked/N.ditissima/$Strain/
  qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
  qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

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

  
