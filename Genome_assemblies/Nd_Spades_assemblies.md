# neonectria_ditissima
Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/neonectria_ditissima

Data is copied from /data/seq_data/miseq/2017/ANALYSIS/171115_M04465_0055_000000000-B87BN/Data/Intensities/BaseCalls
/data/seq_data/miseq/2018/ANALYSIS/180501_M04465_0078_000000000-BM867/Data/Intensities/BaseCalls
/data/seq_data/miseq/2018/ANALYSIS/180823_M04465_0087_000000000-BRM44/Data/Intensities/BaseCalls

to raw_dna/paired/N.ditissima

I cat together the two sequencing runs of the isolate R09/05.
```bash
cat R0905_v2/F/R0905_S1_L001_R1_001.fastq.gz NG-R0905/F/NG-R0905_S4_L001_R1_001.fastq.gz > R0905_all/F/R0905_F_all.fastq.gz
cat R0905_v2/R/R0905_S1_L001_R2_001.fastq.gz NG-R0905/R/NG-R0905_S4_L001_R2_001.fastq.gz > R0905_all/R/R0905_R_all.fastq.gz
```

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
for Strain in R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    RawData=$(ls raw_dna/paired/*/$Strain/F/*.fastq.gz)
	  echo $RawData;
	  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	  qsub $ProgDir/run_fastqc.sh $RawData;
  done

for Strain in R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    RawData=$(ls raw_dna/paired/*/$Strain/R/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
for Strain in R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
    Read_F=raw_dna/paired/*/$Strain/F/*.fastq.gz
    Read_R=raw_dna/paired/*/$Strain/R/*.fastq.gz
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

Data quality was visualised once again following trimming:

```bash
for Strain in R0905_all R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  RawData=$(ls qc_dna/paired/*/$Strain/F/*.fq.gz)
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
done

for Strain in R0905_all R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  RawData=$(ls qc_dna/paired/*/R0905_all/R/*.fq.gz)
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
	qsub $ProgDir/run_fastqc.sh $RawData;
done
```

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
for Strain in R0905_all R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
  #for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do  echo $Strain
  Trim_F=$(ls qc_dna/paired/N.*/$Strain/F/*.fq.gz)
  Trim_R=$(ls qc_dna/paired/N.*/$Strain/R/*.fq.gz)
  Outdir=qc_dna/kmc/N.ditissima/$Strain
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  qsub $ProgDir/kmc_kmer_counting.sh $Outdir $Trim_F $Trim_R
done
```
```bash
for DataDir in $(ls -d qc_dna/paired/N.ditissima/R45-15); do
   F_Read=$(ls $DataDir/F/*.gz)
   R_Read=$(ls $DataDir/R/*.gz)
   Strain=$(echo $DataDir | rev | cut -f1 -d ‘/’ | rev)
   Organism=$(echo $DataDir | rev | cut -f2 -d ‘/’ | rev)
   WorkDir=tmp_dir/R45-15
   mkdir -p $WorkDir
   cp -r $F_Read $WorkDir
   cp -r $R_Read $WorkDir
   cd $WorkDir
   Read1=*R1*
   Read2=*R2*
   gunzip $Read1
   gunzip $Read2
   Sub1=*R1*.fq
   Sub2=*R2*.fq
   echo “$Organism - $Strain”
   count_nucl.pl -i $Sub1 -i $Sub2 -g 45
   cd /home/groups/harrisonlab/project_files/neonectria_ditissima
done
```
```
Ag02 144.71
Ag04 46.19
Ag05 125.05
Ag06 73.52
Ag08 45.42
Ag09_A 67.41
Ag11_A 45.39
Ag11_B 22.29
Ag11_C 87.03
ND8 242.97
ND9 63.88
Hg199 89.92
R6-17-2 52.60
R6-17-3 65.28
R37-15 117.17
R39-15 51.68
R41-15 51.77
R42-15 66.00
R45-15 118.89
R68-17-C3 65.59
P112 54.03
BGV344 55.05
OPC304 37.73
R0905 67.06
R09/05_v2 48.69
R68-17-C2 52.09
NMaj 57.09
SVK1 57.17
SVK2 65.08
```

#Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash
for Strain in NMaj; do
#for Strain in R0905_all R0905_v2 R68-17-C2 SVK1 SVK2 Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag08 Ag11-B R6-17-2 R6-17-3 R41-15; do
  echo $Strain
  F_Read=$(ls qc_dna/paired/N.*/$Strain/F/*.fq.gz)
  R_Read=$(ls qc_dna/paired/N.*/$Strain/R/*.fq.gz)
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.major/$Strain/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
done
```

```bash
for Strain in NMaj; do
#for Strain in R0905_all R0905_v2 SVK1 SVK2 R68-17-C2 Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17; do
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/N.*/$Strain/filtered_contigs/contigs_min_500bp.fasta); do
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
for Strain in R0905_all R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
#for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag11_C BGV344 ND9 OPC304 P112 Ag08 Ag11_B R41-15 R6-17-2 R6-17-3; do  
  for Assembly in $(ls assembly/spades/*/$Strain/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    OutDir=repeat_masked/$Organism/$Strain/
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
done
```

Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
for Strain in R0905_all R0905_v2 R68-17-C2 NMaj SVK1 SVK2; do
#for Strain in Ag11_C BGV344 ND9 OPC304 P112; do
#for Strain in Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag11_C BGV344 ND9 OPC304 P112 Ag08 Ag11_B R41-15 R6-17-2 R6-17-3; do
  for File in $(ls repeat_masked/*/$Strain/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done
  done
```

# Alignment of raw reads vs the Nd genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Nd genome

This was done again using the newest version of the Hg199 genome.

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 OPC304 P112 R0905 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj; do
#Reference=$(ls repeat_masked/N.*/*/Hg199_minion/*/*_contigs_unmasked.fa)
#New genome version was copied to the REFERENCE folder.
Reference=$(ls REFERENCE/Hg199_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/N.*/$Strain); do
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
done
```
Sequence data for isolates with a data from two sequencing runs

```bash
  Reference=$(ls repeat_masked/*/*/*/*_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.*/*); do
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
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Nd_unmasked
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir $Strain
  done
```
