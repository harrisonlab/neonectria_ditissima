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
qsub $scripts/summary_stats/sub_calculate_haplotype_based_stats.sh
qsub $scripts/summary_stats/sub_calculate_4gt.sh
```

This calculation was done over all sites. Now going to proceed for site subsets: synonymous, non-synonymous and four-fold degenerate (silent), in the respective folders

###four_fold_degenerate (analogous to above, for all sites), called silent

```bash
cd $input/NoR68/silent
mkdir contigs
mv *.fasta ./contigs
cp -r ./all/gff/ ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
cd $input/NoR68/silent

qsub $scripts/summary_stats/sub_calculate_nucleotide_diversity.sh
failed
qsub $scripts/summary_stats/sub_calculate_neutrality_stats.sh
qsub $scripts/summary_stats/sub_calculate_fst.sh
qsub $scripts/summary_stats/sub_calculate_haplotype_based_stats.sh
qsub $scripts/summary_stats/sub_calculate_4gt.sh
```

###synonymous
For synonymous and non-synonymous have to create FASTA input first, as done
for silent and all sites in fus_variant_annotation.sh

```bash
cd $input/NoR68
ref_genome=/data/scratch/gomeza/analysis/summary_stats/N.ditissima_contigs_unmasked.fa
python $scripts/summary_stats/vcf_to_fasta.py N.ditissima_contigs_unmasked_noR68_filtered.recode_syn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir syn
mv *.fasta ./syn
```

##non-synonymous
```bash
ref_genome=/data/scratch/gomeza/analysis/summary_stats/N.ditissima_contigs_unmasked.fa
python $scripts/summary_stats/vcf_to_fasta.py N.ditissima_contigs_unmasked_noR68_filtered.recode_nonsyn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir nonsyn
mv *.fasta ./nonsyn
```

##And now back to creating dir structure and carrying out Popgenome analysis

```bash
cd $input/NoR68/syn
mkdir contigs
mv *.fasta ./contigs
cp -r ../all/gff/ ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/NoR68/syn
for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

cd $input/NoR68/syn
qsub $scripts/summary_stats/sub_calculate_nucleotide_diversity.sh
failed
qsub $scripts/summary_stats/sub_calculate_neutrality_stats.sh
qsub $scripts/summary_stats/sub_calculate_fst.sh
qsub $scripts/summary_stats/sub_calculate_haplotype_based_stats.sh
qsub $scripts/summary_stats/sub_calculate_4gt.sh
```
```bash
cd $input/NoR68/nonsyn
mkdir contigs
mv *.fasta ./contigs
cp -r ../all/gff/ ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/NoR68/nonsyn
for a in $PWD/contigs/*/*.fasta
do
    filename=$(basename "$a")
    expected_gff="$PWD/gff/${filename%.fa*}.gff"
    if [ ! -f "$expected_gff" ]
    then
       rm -rf $(dirname $a)
    fi
done

Need to run
cd $input/NoR68/nonsyn
qsub $scripts/summary_stats/sub_calculate_nucleotide_diversity.sh
qsub $scripts/summary_stats/sub_calculate_neutrality_stats.sh
qsub $scripts/summary_stats/sub_calculate_fst.sh
qsub $scripts/summary_stats/sub_calculate_haplotype_based_stats.sh
qsub $scripts/summary_stats/sub_calculate_4gt.sh
```
