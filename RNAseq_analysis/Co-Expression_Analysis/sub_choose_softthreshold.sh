#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=8G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis

/home/adamst/prog/R/R-3.2.5/bin/Rscript --vanilla $scripts/choose_softthreshold.R --out_dir $1 --max_SFT $2
