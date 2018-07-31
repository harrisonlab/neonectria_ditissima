#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 20
#$ -l virtual_free=1G
#$ -l h=blacklace11.blacklace

input=$1
cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $input $temp_dir
input_f=$(basename "$input")
cd $temp_dir

bayescan=/home/gomeza/bin/bayescan2.1/binaries/bayescan_2.1
$bayescan -threads 16 -od ./ $input_f

rm $input_f
cp -r * $cpath
rm -rf $temp_dir

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
cd $cpath
Rscript --vanilla $scripts/plot_bayescan.R ${input%vcf}_fst.txt
