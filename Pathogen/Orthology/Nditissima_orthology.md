#!/bin/bash
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

# Orthology analysis between Neonectria ditissima isolates using Orthomcl

```bash
  ProjDir=/data/scratch/gomeza
  cd $ProjDir
  IsolateAbrv=Nd_all_isolates
  WorkDir=analysis/orthology/Orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 1.a) Format fasta files

```bash

  #Working directory used for OrthoFinder
  WorkDir=analysis_vAG/orthology/OrthoFinder
  mkdir -p $WorkDir

  Taxon_code=199R
  Fasta_file=$(ls gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R09R
  Fasta_file=$(ls gene_pred/codingquary/Ref_Genomes_v2/N.ditissima/R0905/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag02
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag02/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag04
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag04/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag05
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag05/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag06
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag06/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag08
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag08/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A09A
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag09_A/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11A
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_A/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11B
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_B/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11C
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_C/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=BGV3
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/BGV344/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  #Taxon_code=H199
  #Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Hg199/final/final_genes_combined.pep.fasta)
  #Id_field=1
  #orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  #mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=ND8
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/ND8/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=ND9
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/ND9/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=OPC3
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/OPC304/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=P112
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/P112/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R37
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R37-15/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R39
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R39-15/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R41
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R41-15/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R42
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R42-15/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R45
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R45-15/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R602
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R6-17-2/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R603
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R6-17-3/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=68C2
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R68-17-C2/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=68C3
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R68-17-C3/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R305
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/RS305p/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R324
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/RS324p/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=SKV1
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/SVK1/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=SVK2
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/SVK2/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=NMAJ
  Fasta_file=$(ls gene_pred/codingquary/N.major/NMaj/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 2) Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 3.1) Perform an all-vs-all blast of the proteins

```bash

IsolateAbrv=Nd_all_isolates
WorkDir=analysis/orthology/Orthomcl/$IsolateAbrv
BlastDB=$WorkDir/blastall/$IsolateAbrv.db

makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
BlastOut=$WorkDir/all-vs-all_results.tsv
mkdir -p $WorkDir/splitfiles

SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Orthology
for File in $(find $WorkDir/splitfiles)
do
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
BlastOut=$(echo $File | sed 's/.fa/.tab/g')
qsub -h $ProgDir/blast_500.sh $BlastDB $File $BlastOut

JobID=$(qstat | grep 'blast_500' | tail -n 1 | cut -d ' ' -f1)
Queue_Status=$(qstat | grep 'blast_500' | grep 'hqw' | wc -l)
while (($Queue_Status > 0))
do
    Queue_Status=$(qstat | grep 'blast_500' | grep 'hqw' | wc -l)
    load05=$(qstat -u "*" | grep 'blacklace05'| grep 'blast_500' | wc -l)
    load06=$(qstat -u "*" | grep 'blacklace06'| grep 'blast_500' | wc -l)
    load07=$(qstat -u "*" | grep 'blacklace07'| grep 'blast_500' | wc -l)
    load08=$(qstat -u "*" | grep 'blacklace08'| grep 'blast_500' | wc -l)
    load09=$(qstat -u "*" | grep 'blacklace09'| grep 'blast_500' | wc -l)
    load10=$(qstat -u "*" | grep 'blacklace10'| grep 'blast_500' | wc -l)
  if (($load05 < 3))
    then
    qalter $JobID -l h=blacklace05.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 5"
  elif (($load06 < 3))
    then
    qalter $JobID -l h=blacklace06.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 6"
  elif (($load07 < 3))
    then
    qalter $JobID -l h=blacklace07.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 7"
  elif (($load08 < 3))
    then
    qalter $JobID -l h=blacklace08.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 8"
  elif (($load09 < 3))
    then
    qalter $JobID -l h=blacklace09.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 9"
    elif (($load10 < 3))
    then
    qalter $JobID -l h=blacklace10.blacklace
    sleep 5s
    qalter $JobID -h U
    sleep 5s
    echo "Submitted to node 10"
    else
    echo "all nodes full, waiting ten minutes"
    sleep 10m
    fi
    done    
done
done
done
```

## 3.2) Merge the all-vs-all blast results  
```bash
IsolateAbrv=Nd_all_isolates
WorkDir=analysis/orthology/Orthomcl/$IsolateAbrv

  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.fa | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4) Perform ortholog identification

```bash
# Need to create a mysql database before run orthomcl.
# It seems better results are shown in OrthoFinder, therefore I am running this now.

  ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

# Identification of orthologs amongst selected genomes using the OrthoFinder.

## Get whole-genome protein sequences and rename the files

scripts=/home/sobczm/bin/popgen/clock/motif_discovery

```bash
#Create a directory for OrthoFinder run, copy input FASTA files and run orthofinder
#This is not needed
#for Strain in Ag02 Ag08 Ag11_C ND8 R0905_canu_2017_v2 R42-15 R68-17 Ag04 Ag09_A BGV344 ND9 #R37-15 R45-15 RS305p Ag05 Ag11_A Hg199 OPC304 R39-15 R6-17-2 RS324p Ag06 Ag11_B Hg199_minion #P112 R41-15 R6-17-3; do
#cp gene_pred/codingquary/N.ditissima/$Strain/final/final_genes_combined.pep.fasta analysis/orthology/OrthoFinder/"$Strain".pep.fa
#done
#cd analysis/orthology/OrthoFinder
#qsub $scripts/sub_orthofinder.sh OrthoFinder
#It turns out that the pipeline terminates for no reason at different random
#points when executing through SGE. Therefore, will run it on a worker node in a screen.
```
```bash
#IMPORTANT: Input files have to be the pep files of each sample, renamed with maximum 4 letters name + .fasta. This step is done in Orthomcl section.

cp -r ../orthology/Orthomcl/Nd_all_isolates/formatted/ OrthoFinder2/
screen -a

#login in a node with 16 threads and give 1G of memory to each. This is going to take a while, so try to avoid using blacklace11 and leave it for other users.
qlogin -pe smp 24 -l virtual_free=4G
cd analysis_vAG/orthology/OrthoFinder2/formatted
#16 threads used
orthofinder -f ./ -t 24
```

## Orthofinder results:

```bash
#IsolateAbrv=Nd_all_isolates
WorkDir=analysis_vAG/orthology/OrthoFinder
cd $WorkDir/formatted/Results_Mar08
```
```
Number of genes 430235
Number of genes in orthogroups  424576
Number of unassigned genes      5659
Percentage of genes in orthogroups      98.7
Percentage of unassigned genes  1.3
Number of orthogroups   17360
Number of species-specific orthogroups  8
Number of genes in species-specific orthogroups 22
Percentage of genes in species-specific orthogroups     0.0
Mean orthogroup size    24.5
Median orthogroup size  30.0
```

## 5) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.

The genes unique to N.ditissima were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis_vAG/orthology/OrthoFinder
  Orthogroups=$(ls $WorkDir/formatted/Results_Mar08/Orthogroups.txt)
```
sed -i 's/SKV1/SVK1/g' Orthogroups.txt > Orthogroups_edited.txt
```bash
for num in 1
do
echo "The total number of orthogroups is:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | wc -l
echo "The total number of genes in orthogroups is:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -o '|' | wc -l
echo "The number of orthogroups common to N. ditissima is:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e '199R|' | grep -e '68C2|' | grep -e '68C3|' | grep -e '199R|' | grep -e 'Ag02|' | grep -e 'Ag04|' | grep -e 'Ag05|' | grep -e 'Ag06|' | grep -e 'Ag08|' | grep -e 'A09A|' | grep -e 'A11A|' | grep -e 'A11B|' | grep -e 'A11C|' | grep -e 'BGV3|' | grep -e 'ND8|' | grep -e 'ND9|' | grep -e 'OPC3|' | grep -e 'P112|' | grep -e 'R09R|' | grep -e 'R305|' | grep -e 'R324|' | grep -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -e 'R602|' | grep -e 'R603|' | grep -e 'SVK1|' | grep -e 'SVK2|' | wc -l
echo "This represents the following number of genes:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e '68C2|' | grep -e '68C3|' | grep -e '199R|' | grep -e 'Ag02|' | grep -e 'Ag04|' | grep -e 'Ag05|' | grep -e 'Ag06|' | grep -e 'Ag08|' | grep -e 'A09A|' | grep -e 'A11A|' | grep -e 'A11B|' | grep -e 'A11C|' | grep -e 'BGV3|' | grep -e 'ND8|' | grep -e 'ND9|' | grep -e 'OPC3|' | grep -e 'P112|' | grep -e 'R09R|' | grep -e 'R305|' | grep -e 'R324|' | grep -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -e 'R602|' | grep -e 'R603|' | grep -e 'SVK1|' | grep -e 'SVK2|' | grep -o '|' | wc -l

echo "The number of orthogroups common to Neonectria is:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e '199R|' | grep -e '68C2|' | grep -e '68C3|' | grep -e 'Ag02|' | grep -e 'Ag04|' | grep -e 'Ag05|' | grep -e 'Ag06|' | grep -e 'Ag08|' | grep -e 'A09A|' | grep -e 'A11A|' | grep -e 'A11B|' | grep -e 'A11C|' | grep -e 'BGV3|' | grep -e 'ND8|' | grep -e 'ND9|' | grep -e 'OPC3|' | grep -e 'P112|' | grep -e 'R09R|' | grep -e 'R305|' | grep -e 'R324|' | grep -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -e 'R602|' | grep -e 'R603|' | grep -e 'SVK1|' | grep -e 'SVK2|' | grep -e 'NMAJ|' | wc -l
echo "This represents the following number of genes:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e '199R|' | grep -e '68C2|' | grep -e '68C3|' | grep -e 'Ag02|' | grep -e 'Ag04|' | grep -e 'Ag05|' | grep -e 'Ag06|' | grep -e 'Ag08|' | grep -e 'A09A|' | grep -e 'A11A|' | grep -e 'A11B|' | grep -e 'A11C|' | grep -e 'BGV3|' | grep -e 'ND8|' | grep -e 'ND9|' | grep -e 'OPC3|' | grep -e 'P112|' | grep -e 'R09R|' | grep -e 'R305|' | grep -e 'R324|' | grep -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -e 'R602|' | grep -e 'R603|' | grep -e 'SVK1|' | grep -e 'SVK2|' | grep -e 'NMAJ|' | grep -o '|' | wc -l

echo "The number of orthogroups common to high pathogenic isolates tested and absent in low pathogenic is:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e 'Ag02|' | grep -e 'Ag06|' | grep -e '199R|' | grep -v -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -v -e 'R602|' | grep -e 'R603|' | wc -l
echo "This represents the following number of genes:"
cat $WorkDir/formatted/Results_Mar08/Orthogroups.txt | grep -e 'Ag02|' | grep -e 'Ag06|' | grep -e '199R|' | grep -v -e 'R37|' | grep -e 'R39|' | grep -e 'R41|' | grep -e 'R42|' | grep -e 'R45|' | grep -v -e 'R602|' | grep -e 'R603|' | grep -o '|' | wc -l
done
```
```
The total number of orthogroups is:
23019
The total number of genes in orthogroups is:
430235
The number of orthogroups common to N. ditissima is:
11607
This represents the following number of genes:
356487

echo "The number of orthogroups common to Neonectria is:"
10830
This represents the following number of genes
332917

The number of orthogroups common to high pathogenic isolates tested and absent in low pathogenic is
5
This represents the following number of genes
119
```

## 6) Plot venn diagrams:

Orthofinder output:

```bash
IsolateAbrv=Nd_all_isolates
WorkDir=analysis/orthology/OrthoFinder2
GoodProts=analysis/orthology/Orthomcl/Nd_all_isolates/goodProteins/goodProteins.fasta
ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
$ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/Results_Jul25/Orthogroups.txt > $WorkDir/formatted/Results_Jul25/"$IsolateAbrv"_orthogroups.tab

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Orthology
$ProgDir/Ndit_Ortholog_venn.r --inp $WorkDir/formatted/Results_Jul25/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/formatted/Results_Jul25/"$IsolateAbrv"_orthogroups.pdf
```
```
Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name
total number of orthogroups
number of unique singleton genes
number of unique groups of inparalogs
```

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


## Extracting fasta files for orthogroups:

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupsTxt=$(ls $WorkDir/Results_Jul05/Orthogroups.txt)
  GoodProts=analysis/orthology/Orthomcl/Nd_all_isolates/goodProteins/goodProteins.fasta
  OutDir=$WorkDir/orthogroups_fasta
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
```


# Identification of orthologs amongst reference genomes using the OrthoFinder.

## Get whole-genome protein sequences and rename the files

```bash
Taxon_code=199R
Fasta_file=$(ls gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.pep.fasta )
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta analysis/orthology/OrthoFinderRef/"$Taxon_code".fasta

Taxon_code=R09R
Fasta_file=$(ls gene_pred/codingquary/Ref_Genomes/N.ditissima/R0905/final/final_genes_appended_renamed.pep.fasta )
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta analysis/orthology/OrthoFinderRef/"$Taxon_code".fasta
```

## Run orthofinder using qlogin

```bash
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
#login in a node with 16 threads and give 1G of memory to each. This is going to take a while, so try to avoid using blacklace11 and leave it for other users.
qlogin -pe smp 16 -l virtual_free=1G
cd /data/scratch/gomeza/analysis/orthology/OrthoFinderRef/
#8 threads used
orthofinder -f ./ -t 3 -a 3
```


#!/bin/bash
scripts=/home/sobczm/bin/popgen/clock/motif_discovery

# Orthology analysis between Neonectria ditissima isolates using Orthomcl

```bash
  ProjDir=/data/scratch/gomeza
  cd $ProjDir
  IsolateAbrv=Nd_Ref_isolates
  WorkDir=analysis/orthology/OrthomclRef/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  

  cp analysis/orthology/OrthoFinderRef/*.fasta analysis/orthology/OrthomclRef/Nd_Ref_isolates/formatted/
```

## 2) Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 3.1) Perform an Hg-vs-R9 blast of the proteins

```bash

IsolateAbrv=Nd_Ref_isolates
WorkDir=analysis/orthology/OrthomclRef/$IsolateAbrv
BlastDB=$WorkDir/blastall/$IsolateAbrv.db

makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
BlastOut=$WorkDir/Hg-vs-R9_results.tsv
mkdir -p $WorkDir/splitfiles

SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Orthology
for File in $(find $WorkDir/splitfiles)
do
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
BlastOut=$(echo $File | sed 's/.fa/.tab/g')
qsub -h $ProgDir/blast_500.sh $BlastDB $File $BlastOut
done
```

## 3.2) Merge the all-vs-all blast results  
```bash
IsolateAbrv=Nd_all_isolates
WorkDir=analysis/orthology/Orthomcl/$IsolateAbrv

  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.fa | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4) Perform ortholog identification

```bash
# Need to create a mysql database before run orthomcl.
# It seems better results are shown in OrthoFinder, therefore I am running this now.

  ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```
