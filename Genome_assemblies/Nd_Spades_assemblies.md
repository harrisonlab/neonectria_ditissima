# neonectria_ditissima
Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/neonectria_ditissima

Data is copied from /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls
/data/seq_data/miseq/2018/ANALYSIS/180501_M04465_0078_000000000-BM867/Data/Intensities/BaseCalls

to raw_dna/paired/N.ditissima

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
  #for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
  #for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    RawData=$(ls raw_dna/paired/*/$Strain/F/*.fastq.gz)
	  echo $RawData;
	  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	  qsub $ProgDir/run_fastqc.sh $RawData;
  done

  for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
    #for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
    #for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    RawData=$(ls raw_dna/paired/*/$Strain/R/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
  for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
  #for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
  #for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    Read_F=raw_dna/paired/N.ditissima/$Strain/F/*.fastq.gz
    Read_R=raw_dna/paired/N.ditissima/$Strain/R/*.fastq.gz
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

Data quality was visualised once again following trimming:

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
#for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  RawData=$(ls qc_dna/paired/*/$Strain/F/*.fq.gz)
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
done

for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  RawData=$(ls qc_dna/paired/*/$Strain/R/*.fq.gz)
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
done
```

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
#for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  echo $Strain
  Trim_F=$(ls qc_dna/paired/N.*/$Strain/F/*.fq.gz)
  Trim_R=$(ls qc_dna/paired/N.*/$Strain/R/*.fq.gz)
  Outdir=qc_dna/kmc/N.ditissima/temp/$Strain
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  qsub $ProgDir/kmc_kmer_counting.sh 21 $Outdir $Trim_F $Trim_R
done
```

#Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
  #for Strain in Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  echo $Strain
  F_Read=$(ls qc_dna/paired/N.*/$Strain/F/*.fq.gz)
  R_Read=$(ls qc_dna/paired/N.*/$Strain/R/*.fq.gz)
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/$Strain/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
done
```

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/N.ditissima/$Strain/filtered_contigs/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
done
```

# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
  NewSpades=$(ls assembly/spades/*/$Strain/filtered_contigs/contigs_min_500bp.fasta)
  for Assembly in $(ls $NewSpades); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    OutDir=repeat_masked/N.ditissima/$Strain/
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
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
repeat_masked/N.ditissima/Ag02/Ag02_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4667657
repeat_masked/N.ditissima/Ag04/AgN04_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
6076978
repeat_masked/N.ditissima/Ag05/Ag05_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4812553
repeat_masked/N.ditissima/Hg199/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5062507
repeat_masked/N.ditissima/ND8/ND8_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5608134
repeat_masked/N.ditissima/R0905_canu_2017_v2/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5786417
repeat_masked/N.ditissima/R37-15/R37-15_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4258109
repeat_masked/N.ditissima/R45-15/R45-15_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4967354

# Alignment of raw reads vs the Nd genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Nd genome

```bash
for Strain in Ag02 Ag05 ND8 R37-15 Ag04 R45-15 R0905 Hg199; do
  Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/$Strain); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=/data/scratch/gomeza/analysis/genome_alignment/bowtie/$Organism/$Strain/
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
done
```
Sequence data for isolates with a data from two sequencing runs was aligned against the Fus2 genome

```bash
  Reference=$(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2_canu_new')
  for StrainPath in $(ls -d qc_dna/paired/F.*/* | grep -e 'HB6' -e 'Fus2'); do
    echo $StrainPath
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Fus2_unmasked
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir $Strain
  done
```
