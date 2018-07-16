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
  Taxon_code=Ag02
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag02/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag04
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag04/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag05
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag05/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag06
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag06/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=Ag08
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag08/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A09A
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag09_A/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11A
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_A/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11B
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_B/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=A11C
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Ag11_C/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=BGV3
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/BGV344/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=H199
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/Hg199/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=ND8
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/ND8/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=ND9
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/ND9/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=OPC3
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/OPC304/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=P112
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/P112/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R09
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R0905_canu_2017_v2/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R37
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R37-15/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R39
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R39-15/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R41
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R41-15/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R42
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R42-15/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R45
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R45-15/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R602
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R6-17-2/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R603
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R6-17-3/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

  Taxon_code=R68
  Fasta_file=$(ls gene_pred/codingquary/N.ditissima/R68-17/final/final_genes_combined.pep.fasta)
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

for Strain in Ag02 Ag08 Ag11_C ND8 R0905_canu_2017_v2 R42-15 R68-17 Ag04 Ag09_A BGV344 ND9 R37-15 R45-15 RS305p Ag05 Ag11_A Hg199 OPC304 R39-15 R6-17-2 RS324p Ag06 Ag11_B Hg199_minion P112 R41-15 R6-17-3; do
cp gene_pred/codingquary/N.ditissima/$Strain/final/final_genes_combined.pep.fasta analysis/orthology/OrthoFinder/"$Strain".pep.fa
done

cd analysis/orthology/OrthoFinder
qsub $scripts/sub_orthofinder.sh OrthoFinder

#It turns out that the pipeline terminates for no reason at different random
#points when executing through SGE. Therefore, will run it on a worker node in a screen.
#Input files have to be the pep files of each sample, renamed with maximum 4 letters name + .fasta

cp -r ../orthology/Orthomcl/Nd_all_isolates/formatted/ OrthoFinder2/
screen -a

#login in a node with 16 threads and give 1G of memory to each. This is going to take a while, so try to avoid using blacklace11 and leave it for other users.
qlogin -pe smp 16 -l virtual_free=1G
cd analysis/orthology/OrthoFinder2/formatted
#16 threads used
orthofinder -f ./ -t 3 -a 3
```
## Orthofinder results:

WorkDir=analysis/orthology/OrthoFinder  
IsolateAbrv=Nd_all_isolates

output files are in:
```bash
ls $WorkDir/Results_Jul05
```
```
Number of genes 389858
Number of genes in orthogroups  385064
Number of unassigned genes      4794
Percentage of genes in orthogroups      98.8
Percentage of unassigned genes  1.2
Number of orthogroups   17163
Number of species-specific orthogroups  8
Number of genes in species-specific orthogroups 22
Percentage of genes in species-specific orthogroups     0.0
Mean orthogroup size    22.4
Median orthogroup size  27.0
```










## 5) Plot venn diagrams:

Orthofinder output:

```bash
  GoodProts=analysis/orthology/Orthomcl/Nd_all_isolates/goodProteins/goodProteins.fasta
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/Results_Jul05/Orthogroups.txt > $WorkDir/Results_Jul05/"$IsolateAbrv"_orthogroups.tab
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.pdf
```



Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name
total number of orthogroups
number of unique singleton genes
number of unique groups of inparalogs

```
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

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) Clade unique gene families

The genes unique to A. tenuissima apple pathotypes were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/At_Aa_Ag_all_isolates
  Orthogroups=$(ls $WorkDir/formatted/Results_May04/Orthogroups.txt)
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  # Braker_genes=
```ENG <- c("Ag02", "Ag04", "Ag05", "Ag06", "Hg199", "R0905", "R6-17-2","R6-17-3")
NI <- c("Ag11_A", "Ag11_B", "Ag11_C")
IRL <- c("Ag08", "Ag09_A")
BRA <- c("ND8", "ND9")
ESP <- c("BGV344", "P112","OPC304")
BEL <- c("R37-15", "R39-15")
NL <- c("R41-15", "R42-15","R45-15")

#### 6.1.a ) Orthologroups only containing A. tenuissima genes were extracted:
```bash
  for num in 1
  do
      echo "The total number of orthogroups is:"
      cat $WorkDir/Results_Jul05/Orthogroups.txt | wc -l
      echo "The total number of genes in orthogroups is:"
      cat $WorkDir/Results_Jul05/Orthogroups.txt | grep -o '|' | wc -l
      echo "The number of orthogroups common to N. ditissima is:"
      cat $WorkDir/Results_Jul05/Orthogroups.txt | grep -e 'Ag02|' | grep -e 'Ag04|' | grep -e 'Ag05|' | grep -e 'Ag06|' grep -e 'Hg199|' grep -e 'R0905|' grep -e 'R6-17-2|' grep -e 'R6-17-3|' grep -e 'Ag11_A|' grep -e 'Ag11_B|' grep -e 'Ag11_C|' grep -e 'Ag08|' grep -e 'Ag09_A|' grep -e 'ND8|' grep -e 'ND9|' grep -e 'BGV344|' grep -e 'P112|' grep -e 'OPC304|' grep -e 'R37-15|' grep -e 'R39-15|' grep -e 'R41-15|' grep -e 'R42-15|' grep -e 'R45-15|' | wc -l

    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -e 'SCRP249|' | grep -e 'SCRP324|' | grep -e 'SCRP333|' | grep -v -e 'A4|' -e 'Bc1|' -e 'Bc16|' -e 'Bc23|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'Nov77|' -e 'Nov9|' -e 'ONT3|' -e 'SCRP245_v2|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to P. fragariae is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'A4|' | grep -e 'Bc1|' | grep -e 'Bc16|' | grep -e 'Bc23|' | grep -e 'Nov5|' | grep -e 'Nov27|' | grep -e 'Nov71|' | grep -e 'Nov77|' | grep -e 'Nov9|' | grep -e 'ONT3|' | grep -e 'SCRP245_v2|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'A4|' | grep -e 'Bc1|' | grep -e 'Bc16|' | grep -e 'Bc23|' | grep -e 'Nov5|' | grep -e 'Nov27|' | grep -e 'Nov71|' | grep -e 'Nov77|' | grep -e 'Nov9|' | grep -e 'ONT3|' | grep -e 'SCRP245_v2|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to UK1 isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov27|' -e 'Nov71|' -e 'Nov77|' -e 'Nov9|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc1|' | grep -e 'Nov5|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov27|' -e 'Nov71|' -e 'Nov77|' -e 'Nov9|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc1|' | grep -e 'Nov5|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to UK2 isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'Bc1|' -e 'Bc23|' -e 'Nov27|' -e 'Nov71|' -e 'Nov77|' -e 'Nov9|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'Nov5|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc16|' | grep -e 'A4|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'Bc1|' -e 'Bc23|' -e 'Nov27|' -e 'Nov71|' -e 'Nov77|' -e 'Nov9|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'Nov5|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc16|' | grep -e 'A4|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to UK3 isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Nov27|' | grep -e 'Nov71|' | grep -e 'Nov9|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'ONT3|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Nov27|' | grep -e 'Nov71|' | grep -e 'Nov9|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to CA4 isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'Nov9|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'ONT3|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'Nov9|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'ONT3|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to CA5 isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'ONT3|' -e 'Nov9|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc23|' | grep -e 'Nov77|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'A4|' -e 'Bc16|' -e 'ONT3|' -e 'Nov9|' -e 'SCRP245_v2|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'Bc23|' | grep -e 'Nov77|' | grep -o '|' | wc -l
    echo "The number of orthogroups common to Unknown race isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'ONT3|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'Nov9|' -e 'A4|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'SCRP245_v2|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'ONT3|' -e 'Bc16|' -e 'Bc23|' -e 'Nov77|' -e 'Nov9|' -e 'A4|' -e 'Bc1|' -e 'Nov5|' -e 'Nov27|' -e 'Nov71|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' | grep -e 'SCRP245_v2|' | grep -o '|' | wc -l
    echo "The number of orthogroups with only seven highly conserved target isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'ONT3|' -e 'Bc23|' -e 'Nov77|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' -e 'SCRP245_v2|' | grep -e 'Nov5|' | grep -e 'Nov27' | grep -e 'Nov71' | grep -e 'Bc16' | grep -e 'Nov9' | grep -e 'Bc1' | grep -e 'A4' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -v -e 'ONT3|' -e 'Bc23|' -e 'Nov77|' -e 'SCRP249|' -e 'SCRP324|' -e 'SCRP333|' -e 'SCRP245_v2|' | grep -e 'Nov5|' | grep -e 'Nov27' | grep -e 'Nov71' | grep -e 'Bc16' | grep -e 'Nov9' | grep -e 'Bc1' | grep -e 'A4' | grep -o '|' | wc -l
    echo "The number of orthogroups containing all seven highly conserved target isolates is:"
    cat $WorkDir/Orthogroups.txt | grep -e 'Nov5|' | grep -e 'Nov27' | grep -e 'Nov71' | grep -e 'Bc16' | grep -e 'Nov9' | grep -e 'Bc1' | grep -e 'A4' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/Orthogroups.txt | grep -e 'Nov5|' | grep -e 'Nov27' | grep -e 'Nov71' | grep -e 'Bc16' | grep -e 'Nov9' | grep -e 'Bc1' | grep -e 'A4' | grep -o '|' | wc -l
done


```bash
for num in 1; do
  # AtenUniq
AtenUniqDir=$WorkDir/A.tenuissima_unique
Uniq_Aten_groups=$AtenUniqDir/A.tenuissima_unique.txt
mkdir -p $AtenUniqDir
cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'Aa' -e 'Ag' > $Uniq_Aten_groups
echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
cat $Uniq_Aten_groups | wc -l
echo "The following number genes from isolate 648 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_1' | wc -l
echo "The following number genes from isolate 1082 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_2' | wc -l
echo "The following number genes from isolate 1164 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_3' | wc -l
echo "The following number genes from isolate 24350 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_4' | wc -l
echo "The following number genes from isolate 648 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_5' | wc -l
echo "The following number genes from isolate 743 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_6' | wc -l
echo "The following number genes from isolate 1166 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_7' | wc -l
echo "The following number genes from isolate 1177 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
51
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 1082 are contained in these orthogorups:
51
The following number genes from isolate 1164 are contained in these orthogorups:
51
The following number genes from isolate 24350 are contained in these orthogorups:
52
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 743 are contained in these orthogorups:
51
The following number genes from isolate 1166 are contained in these orthogorups:
51
The following number genes from isolate 1177 are contained in these orthogorups:
51
```

#### 6.1.b ) Orthologroups only containing A. arborescens genes were extracted:

```bash
for num in 1; do
  # AarbUniq
  AarbUniqDir=$WorkDir/A.arborescens_unique
  Uniq_Aarb_groups=$AarbUniqDir/A.arborescens_unique.txt
  mkdir -p $AarbUniqDir
  cat $Orthogroups | grep -e 'Aa_1' | grep -e 'Aa_2' | grep -e 'Aa_3' | grep -v -e 'At' -e 'Ag' > $Uniq_Aarb_groups
  echo "The number of orthogroups unique to A.arborescens are:"
  cat $Uniq_Aarb_groups | wc -l
  echo "The following number genes from isolate 675 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_1' | wc -l
  echo "The following number genes from isolate 97.0013 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_2' | wc -l
  echo "The following number genes from isolate 97.0016 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_3' | wc -l
done
```

```
The number of orthogroups unique to A.arborescens are:
150
The following number genes from isolate 675 are contained in these orthogorups:
155
The following number genes from isolate 97.0013 are contained in these orthogorups:
154
The following number genes from isolate 97.0016 are contained in these orthogorups:
156
```

#### 6.1.c ) Orthologroups only containing A. gaisen genes were extracted:

```bash
for num in 1; do
  # AgaiPathUniq
  AgaiUniqDir=$WorkDir/A.gaisen_unique
  Uniq_Agai_groups=$AgaiUniqDir/A.gaisen_unique.txt
  mkdir -p $AgaiUniqDir
  cat $Orthogroups | grep -e 'Ag_1' | grep -v  -e 'At' -e 'Aa' > $Uniq_Agai_groups
  echo "The number of orthogroups unique to A.gaisen pear pathotype pathotype are:"
  cat $Uniq_Agai_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Agai_groups | grep -o -e 'Ag_1' | wc -l
done
```

```
The number of orthogroups unique to A.gaisen pear pathotype pathotype are:
327
The following number genes from isolate 650 are contained in these orthogorups:
329
```

#### 6.1.d ) Orthologroups only containing A. tenuissima non pathotype genes were extracted:

```bash
for num in 1; do
  # AtenNonPathUniq
  AtenNonPathUniqDir=$WorkDir/A.tenuissima_non_pathotype_unique
  Uniq_AtenNonPath_groups=$AtenNonPathUniqDir/A.tenuissima_non_pathotype_unique.txt
  mkdir -p $AtenNonPathUniqDir
  cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -v -e 'At_5' -e 'At_6' -e 'At_7' -e 'At_8' -e 'Aa' -e 'Ag' > $Uniq_AtenNonPath_groups
  echo "The number of orthogroups unique to A.tenuissima non-apple pathotype are:"
  cat $Uniq_AtenNonPath_groups | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_1' | wc -l
  echo "The following number genes from isolate 1082 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_2' | wc -l
  echo "The following number genes from isolate 1164 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_3' | wc -l
  echo "The following number genes from isolate 24350 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_4' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima non-apple pathotype are:
0
The following number genes from isolate 648 are contained in these orthogorups:
0
The following number genes from isolate 1082 are contained in these orthogorups:
0
The following number genes from isolate 1164 are contained in these orthogorups:
0
The following number genes from isolate 24350 are contained in these orthogorups:
0
```

#### 6.1.d ) Orthologroups only containing A. tenuissima apple pathotype genes were extracted:

```bash
for num in 1; do
  # AtenPathUniq
  AtenPathUniqDir=$WorkDir/A.tenuissima_apple_pathotype_unique
  Uniq_AtenPath_groups=$AtenPathUniqDir/A.tenuissima_apple_pathotype_unique.txt
  mkdir -p $AtenPathUniqDir
  cat $Orthogroups | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' -e 'Ag' > $Uniq_AtenPath_groups
  echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
  cat $Uniq_AtenPath_groups | wc -l
  echo "The following number genes from isolate 635 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
49
The following number genes from isolate 635 are contained in these orthogorups:
56
The following number genes from isolate 743 are contained in these orthogorups:
57
The following number genes from isolate 1166 are contained in these orthogorups:
72
The following number genes from isolate 1177 are contained in these orthogorups:
55
```

```bash
cat $Uniq_AtenPath_groups | grep -o -P "At_7\|\S+" | cut -f2 -d '|' > tmp.txt
Isolate='1166'
AnnotTab=$(ls gene_pred/annotation/A.*/$Isolate/${Isolate}_annotation_ncbi.tsv)
OutFile=$(echo $AnnotTab | sed 's/_annotation_ncbi.tsv/_apple_pathotype_unique_orthogroups.tsv/g')
for Gene in $(cat tmp.txt); do
  cat $AnnotTab | grep -w "^$Gene"
done > $OutFile
```


#### 6.1.e ) Orthologroups only containing genes from apple and pear pathotypes were extracted:

```bash
for num in 1; do
  # PathUniq
  PathUniqDir=$WorkDir/Pathotype_unique
  Uniq_Path_groups=$PathUniqDir/Path_unique.txt
  mkdir -p $PathUniqDir
  cat $Orthogroups | grep -e 'Ag' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' > $Uniq_Path_groups
  echo "The number of orthogroups unique to apple and pear pathotypes are:"
  cat $Uniq_Path_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'Ag_1' | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to apple and pear pathotypes are:
48
The following number genes from isolate 650 are contained in these orthogorups:
62
The following number genes from isolate 648 are contained in these orthogorups:
58
The following number genes from isolate 743 are contained in these orthogorups:
58
The following number genes from isolate 1166 are contained in these orthogorups:
54
The following number genes from isolate 1177 are contained in these orthogorups:
60
```
