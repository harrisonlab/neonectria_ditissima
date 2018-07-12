# Second stage of summary_stats analysis

The R scripts called for this analysis have to be changed for every population.

## Nd analysis - No R68

### Set initial variables

```bash
input=/data/scratch/gomeza/analysis/summary_stats
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
```

In order to calculate different statistics in Popgenome, the input has to be
arranged in a particular way.
The input directory should contain two folders.
Folder No. 1: named "gff", contains GFF files for all the contigs output from
the split_gff_contig.sh script
Folder No. 2: named "contigs", contains subfolders, each subfolder named with
exact contig name and containing one individual contig FASTA file, also named
with exact contig name, as output from vcf_to_fasta.py

### Create this directory structure

```bash
cd $input/NoR68/all
#This folder contaings only contig FASTA files
#So create a new "contigs" directory to hold those files:
mkdir contigs
mv *.fasta ./contigs
```

### Split the master GFF file into one contig --> one GFF file. Required for analyses in Popgenome below.

```bash

cd $input
mkdir -p NoR68/gff
sh $scripts/summary_stats/split_gff_contig.sh final_genes_appended.gff3
mv *.gff NoR68/gff
```


In the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
cd NoR68/all/contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
```

Navigate to the input folder holding the two folders: "contigs" and "gff" to proceed with Popgenome run.

Lastly, test if all contigs have a matching GFF file. In some cases, no genes are predicted on a given contig, and GFF file for it will be missing

Check for orphan contigs with no matching gff file, which need to be removed prior to the run.

```bash
cd $input/all
for a in $PWD/contigs/*/*.fasta
do
filename=$(basename "$a")
expected_gff="$PWD/gff/${filename%.fa*}.gff"
if [ ! -f "$expected_gff" ];
then
   rm -rf $(dirname $a)
fi
done
```

The R script used below is custom-made for each run (see first few lines of it)
It requires custom definition of populations, and individual assignment to them.
The example below calculates nucleotide diversity within (Pi) and between (Dxy) populations/
Other scripts (sub_calculate_neutrality_stats.sh) are used in analogous manner.

```bash
qsub $scripts/summary_stats/sub_calculate_nucleotide_diversity.sh
qsub $scripts/summary_stats/sub_calculate_neutrality_stats.sh
qsub $scripts/summary_stats/sub_calculate_fst.sh
#qsub $scripts/sub_calculate_haplotype_based_stats.sh
```





#This calculation was done over all sites. Now going to proceed for site subsets:
#synonymous, non-synonymous and four-fold degenerate (silent), in the respective folders

#four_fold_degenerate (analogous to above, for all sites)
cd $input/silent
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
cd $input/silent

qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh

#For synonymous and non-synonymous have to create FASTA input first, as done
#for silent and all sites in fus_variant_annotation.sh
##synonymous
cd $input
ref_genome=/home/sobczm/popgen/summary_stats/Fus2_canu_contigs_unmasked.fa
python $scripts/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_syn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir syn
mv *.fasta ./syn

##non-synonymous
cd $input
ref_genome=/home/sobczm/popgen/summary_stats/Fus2_canu_contigs_unmasked.fa
python $scripts/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_nonsyn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir nonsyn
mv *.fasta ./nonsyn

## And now back to creating dir structure and carrying out Popgenome analysis
cd $input/syn
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/syn
qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh

cd $input/nonsyn
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/nonsyn
qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh














```bash
scripts2=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis/popgenome_scripts
qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```

This calculation was done over all sites. Now going to proceed for site subsets:
synonymous, non-synonymous and four-fold degenerate (silent)

### four_fold_degenerate (analogous to above, for all sites), called silent

```bash
cd $input/silent
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done
cd $input/silent

# Check all contigs have a matching gff and remove any that do not

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```

## For synonymous and non-synonymous have to create FASTA input first

### Synonymous sites

```bash
cd $input
ref_genome=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa
ProgDir=/home/adamst/git_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py \
polished_contigs_unmasked_UK123_filtered.recode_syn.vcf $ref_genome 2
mkdir syn
mv *.fasta ./syn
```

### Non-synonymous sites

```bash
cd $input
ref_genome=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa
ProgDir=/home/adamst/git_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py \
polished_contigs_unmasked_UK123_filtered.recode_nonsyn.vcf $ref_genome 2
mkdir nonsyn
mv *.fasta ./nonsyn
```

## Create directory structure and carry out Popgenome analysis

```bash
cd $input/syn
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done

cd $input/syn

# Check all contigs have a matching gff and remove any that do not

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh

cd $input/nonsyn
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done

# Check all contigs have a matching gff and remove any that do not

cd $input/nonsyn

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```

## Pf analysis

### Set initial variables for Pf

```bash
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats
scripts=/home/adamst/git_repos/scripts/popgen_analysis
```

```
In order to calculate different statistics in Popgenome, the input has to be
arranged in a particular way.
The input directory should contain two folders.
Folder No. 1: named "gff", contains GFF files for all the contigs output from
the split_gff_contig.sh script
Folder No. 2: named "contigs", contains subfolders, each subfolder named with
exact contig name and containing one individual contig FASTA file, also named
with exact contig name, as output from vcf_to_fasta.py
```

### Create directory structure

```bash
mkdir -p $input/all_Pf
cd $input/all_Pf
```

### This folder contains only contig FASTA files for Pf

So create a new "contigs" directory to hold those files

```bash
mkdir contigs
mv *.fasta ./contigs
```

### Copy the "gff" folder containing gff files for Pf

```bash
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
```

### In the folder "contigs" create subfolders, each to hold one FASTA file

```bash
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done
```

## Navigate to the input folder and proceed with Popgenome run

```bash
cd $input/all_Pf
```

### Test if all contigs have a matching gff and remove any which do not

```bash
for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done
```

```
The R script used below is custom-made for each run (see first few lines of it).
It requires custom definition of populations, and individual assignment to them.
The example below calculates nucleotide diversity within (Pi) and between (Dxy) populations.
Other scripts (sub_calculate_neutrality_stats.sh) are used in analogous manner.
Vcf of all Pf strains, bar NOV-77 has been phased, run for haplotype-based stats.
```

```bash
scripts2=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis/popgenome_scripts
qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```

## This calculation was done over all sites. Proceed for site subsets

### synonymous, non-synonymous and four-fold degenerate (silent)

#### four_fold_degenerate (analogous to above, for all sites) for Pf

```bash
cd $input/silent_Pf
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./

cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done

cd $input/silent_Pf

# Check all contigs have a matching gff and remove any that do not

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```

## For synonymous and non-synonymous have to create FASTA input first for Pf

### Synonymous sites for Pf

```bash
cd $input
ref_genome=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa
ProgDir=/home/adamst/git_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py \
polished_contigs_unmasked_Pf_filtered.recode_syn.vcf $ref_genome 2
mkdir syn_Pf
mv *.fasta ./syn_Pf
```

### Non-synonymous

```bash
cd $input
ref_genome=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa
ProgDir=/home/adamst/git_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py \
polished_contigs_unmasked_Pf_filtered.recode_nonsyn.vcf $ref_genome 2
mkdir nonsyn_Pf
mv *.fasta ./nonsyn_Pf
```

## Create directory structure and carry out Popgenome analysis for Pf

```bash
cd $input/syn_Pf
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done

cd $input/syn_Pf

# Check all contigs have a matching gff and remove any that do not

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh

cd $input/nonsyn_Pf
mkdir contigs
mv *.fasta ./contigs
cp -r \
/home/groups/harrisonlab/project_files/phytophthora_fragariae/summary_stats/gff ./
cd contigs
for f in *.fasta
do
    folder=${f%.fasta}
    mkdir $folder
    mv $f $folder
done

cd $input/nonsyn_Pf

# Check all contigs have a matching gff and remove any that do not

for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

qsub $scripts2/sub_calculate_nucleotide_diversity.sh
qsub $scripts2/sub_calculate_neutrality_stats.sh
qsub $scripts2/sub_calculate_fst.sh
qsub $scripts2/sub_calculate_4_gamete_test.sh
```
